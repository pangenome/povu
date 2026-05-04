#include "ita/variation/rov.hpp"

#include <cctype>   // for isdigit
#include <cstddef>  // for size_t
#include <cstdlib>  // for exit, EXIT_FAILURE
#include <optional> // for optional, operator==
#include <string>   // for basic_string, string
#include <vector>   // for vector

#include <liteseq/refs.h>     // for ref_walk
#include <log.h>	      // for log_error, log_fatal
#include <oza/graph/pvst.hpp> // for Tree, VertexBase
#include <quilt/shim.hpp>     // for format
#include <quilt/types.hpp>    // for qt

#include "ita/graph/graph.hpp"	   // for RoV, find_walks, pgt
#include "ita/variation/color.hpp" // for color, color_e

namespace ita::rov
{
namespace pgg = povu::genomics::graph;
namespace pvst = oza::pvst;
namespace lq = liteseq;

std::ostream &operator<<(std::ostream &os, var_type_e vt)
{
	return os << to_string_view(vt);
}

var_type_e covariant(var_type_e a) noexcept
{
	switch (a) {
	case var_type_e::ins:
		return var_type_e::del;
	case var_type_e::del:
		return var_type_e::ins;
	case var_type_e::sub:
		return var_type_e::sub;
	case var_type_e::subr:
		return var_type_e::subr;
	case var_type_e::inv:
		return var_type_e::inv;
	}

	log_fatal("Unknown variant type");
	std::exit(EXIT_FAILURE);
}

/**
 * Parse a genomic region string in the format: ref:start-end
 * Supports optional reference prefix (e.g., GRCh38#chr17:start-end)
 */
std::optional<genomic_region>
parse_genomic_region(const std::string &region_str)
{
	// Find the colon separator
	std::size_t colon_pos = region_str.find(':');
	if (colon_pos == std::string::npos) {
		log_error("Invalid region format '%s': missing colon "
			  "separator",
			  region_str.c_str());
		return std::nullopt;
	}

	std::string ref_name = region_str.substr(0, colon_pos);
	if (ref_name.empty()) {
		log_error("Invalid region format '%s': empty reference name",
			  region_str.c_str());
		return std::nullopt;
	}

	std::string range_str = region_str.substr(colon_pos + 1);
	std::size_t dash_pos = range_str.find('-');
	if (dash_pos == std::string::npos) {
		log_error("Invalid region format '%s': missing dash "
			  "separator in range",
			  region_str.c_str());
		return std::nullopt;
	}

	std::string start_str = range_str.substr(0, dash_pos);
	std::string end_str = range_str.substr(dash_pos + 1);

	// Validate that start and end are numeric
	for (char c : start_str) {
		if (!std::isdigit(c)) {
			log_error("Invalid region format '%s': non-numeric "
				  "start position",
				  region_str.c_str());
			return std::nullopt;
		}
	}
	for (char c : end_str) {
		if (!std::isdigit(c)) {
			log_error("Invalid region format '%s': non-numeric "
				  "end position",
				  region_str.c_str());
			return std::nullopt;
		}
	}

	try {
		qt::idx_t start = std::stoull(start_str);
		qt::idx_t end = std::stoull(end_str);

		if (start >= end) {
			log_error("Invalid region '%s': start position (%ul) "
				  "must be less than end position (%ul)",
				  region_str.c_str(), start, end);
			return std::nullopt;
		}

		return genomic_region(ref_name, start, end);
	}
	catch (const std::exception &e) {
		log_error("Failed to parse region '%s': %s", region_str.c_str(),
			  e.what());
		return std::nullopt;
	}
}

/**
 * Check if a PVST vertex overlaps with a genomic region
 * A vertex overlaps if any of its constituent graph vertices have
 * positions within the region on the specified reference path
 */
bool pvst_vertex_overlaps_region(const bd::VG &g,
				 const pvst::VertexBase *pvst_v,
				 const genomic_region &region, qt::id_t ref_id)
{
	// Get the route parameters which contain source and sink vertex
	// IDs
	auto opt_route = pvst_v->get_route_params();
	if (!opt_route.has_value()) {
		return false; // vertices without route params are
			      // skipped
	}

	const pvst::route_params_t &route = opt_route.value();
	qt::id_t start_v_id = route.start.v_id;
	qt::id_t end_v_id = route.end.v_id;

	// Convert vertex IDs to indices
	qt::idx_t start_v_idx = g.v_id_to_idx(start_v_id);
	qt::idx_t end_v_idx = g.v_id_to_idx(end_v_id);

	// Get step indices for these vertices on the reference
	const std::vector<qt::idx_t> &start_steps =
		g.get_vertex_ref_idxs(start_v_idx, ref_id);
	const std::vector<qt::idx_t> &end_steps =
		g.get_vertex_ref_idxs(end_v_idx, ref_id);

	// If either vertex is not on this reference, skip
	if (start_steps.empty() || end_steps.empty())
		return false;

	// Get the reference walk to access genomic positions (loci)
	const lq::ref *ref_ptr = g.get_ref_vec(ref_id);
	const lq::ref_walk *ref_w = ref_ptr->walk;
	qt::idx_t ref_step_count = lq::get_step_count(ref_ptr);

	// Check if any genomic position falls within the region
	// Use actual genomic coordinates (loci) instead of step indices
	for (qt::idx_t step : start_steps) {
		if (step < ref_step_count) {
			qt::idx_t genomic_pos = ref_w->loci[step];
			if (genomic_pos >= region.start &&
			    genomic_pos < region.end) {
				return true;
			}
		}
	}

	for (qt::idx_t step : end_steps) {
		if (step < ref_step_count) {
			qt::idx_t genomic_pos = ref_w->loci[step];
			if (genomic_pos >= region.start &&
			    genomic_pos < region.end) {
				return true;
			}
		}
	}

	return false;
}

void find_pvst_rovs(const bd::VG &g, const pvst::Tree &pvst,
		    const std::set<qt::u32> &colored_vtxs, std::vector<RoV> &rs,
		    const std::optional<genomic_region> &region = std::nullopt,
		    const std::optional<qt::id_t> &region_ref_id = std::nullopt)
{

	for (qt::u32 i : colored_vtxs) { // i is pvst_v_idx
		const pvst::VertexBase *v = pvst.get_vertex_const_ptr(i);

		// Apply region filtering if specified
		if (region.has_value() && region_ref_id.has_value()) {
			if (!pvst_vertex_overlaps_region(
				    g, v, region.value(),
				    region_ref_id.value())) {
				continue; // Skip vertices outside the
					  // region
			}
		}

		RoV r{v};

		qt::status_t s = povu::genomics::graph::find_walks(g, r);

		if (r.size() < 3) {
			log_warn("RoV too small (size=%ul): %s. Skipping.",
				 r.size(), r.as_str().c_str());
			continue;
		}

		if (s != 0) {
			std::cerr << "Skipping RoV: " << r.as_str() << "\n";
			continue;
		}

		rs.emplace_back(std::move(r));
	}
}

/**
 * find walks in the graph based on the leaves of the pvst
 * initialize RoVs from flubbles
 * @param region Optional genomic region to filter RoVs
 */
std::vector<RoV> gen_rov(const std::vector<pvst::Tree> &pvsts, const bd::VG &g,
			 const std::set<qt::id_t> &to_call_ref_ids,
			 const std::optional<genomic_region> &region)
{
	// the set of RoVs to return
	std::vector<RoV> rs;
	rs.reserve(pvsts.size());

	// Resolve region reference ID if region is specified
	std::optional<qt::id_t> region_ref_id = std::nullopt;
	if (region.has_value()) {
		region_ref_id = g.get_ref_id(region.value().ref_name);
		if (!region_ref_id.has_value()) {
			std::string err = qs::format(
				"Reference path '{}' not found in graph",
				region.value().ref_name);
			log_error("%s", err.c_str());
			log_warn("No RoVs will be generated due to invalid "
				 "region");
			return rs; // Return empty result
		}
		std::string inf =
			qs::format("Filtering RoVs to region {}:{}-{}",
				   region.value().ref_name,
				   region.value().start, region.value().end);
		log_info("%s", inf.c_str());
	}

	for (qt::idx_t i{}; i < pvsts.size(); i++) { // for each pvst
		const pvst::Tree &pvst = pvsts.at(i);
		std::set<qt::u32> colored_vtxs =
			ic::color_pvst(g, pvst, to_call_ref_ids);

		find_pvst_rovs(g, pvst, colored_vtxs, rs, region,
			       region_ref_id);
	}

	if (region.has_value()) {
		std::string inf =
			qs::format("Generated {} RoVs in region {}:{}-{}",
				   rs.size(), region.value().ref_name,
				   region.value().start, region.value().end);
		log_info("%s", inf.c_str());
	}

	return rs;
}

} // namespace ita::rov
