#include "mto/to_structure_export.hpp"

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string_view>

#include "ita/variation/rov.hpp"
#include "liteseq/refs.h"
#include "povu/algorithms/flubbles.hpp"
#include "povu/common/constants.hpp"
#include "povu/refs/refs.hpp"

namespace mto::to_structure_export
{
namespace fs = std::filesystem;
namespace ir = ita::rov;
namespace pc = povu::constants;
namespace pgt = povu::types::graph;
namespace pr = povu::refs;

namespace
{
void write_json_string(std::ostream &os, std::string_view value)
{
	os << '"';
	for (char ch : value) {
		switch (ch) {
		case '"':
			os << "\\\"";
			break;
		case '\\':
			os << "\\\\";
			break;
		case '\b':
			os << "\\b";
			break;
		case '\f':
			os << "\\f";
			break;
		case '\n':
			os << "\\n";
			break;
		case '\r':
			os << "\\r";
			break;
		case '\t':
			os << "\\t";
			break;
		default:
			os << ch;
			break;
		}
	}
	os << '"';
}

void write_key(std::ostream &os, std::string_view key)
{
	write_json_string(os, key);
	os << ':';
}

void write_key_string(std::ostream &os, std::string_view key,
		      std::string_view value)
{
	write_key(os, key);
	write_json_string(os, value);
}

void write_nullable_string(std::ostream &os,
			   const std::optional<std::string> &value)
{
	if (value.has_value())
		write_json_string(os, *value);
	else
		os << "null";
}

std::string route_name(pvst::rt_e route)
{
	return route == pvst::rt_e::s2e ? "s2e" : "e2s";
}

std::string oriented_step(const pgt::id_or_t &step)
{
	return step.as_str();
}

std::string endpoint_label(const pvst::route_params_t &route)
{
	return oriented_step(route.start) + oriented_step(route.end);
}

std::string gfa_from_orientation(pgt::v_end_e end)
{
	return end == pgt::v_end_e::r ? "+" : "-";
}

std::string gfa_to_orientation(pgt::v_end_e end)
{
	return end == pgt::v_end_e::l ? "+" : "-";
}

std::string join_genotype_column(const std::vector<std::string> &column)
{
	bool untraversed = std::all_of(
		column.begin(), column.end(),
		[](const std::string &value) { return value == "."; });
	if (untraversed)
		return ".";

	std::ostringstream os;
	for (std::size_t idx{}; idx < column.size(); ++idx) {
		if (idx > 0)
			os << '|';
		os << column[idx];
	}
	return os.str();
}

std::string allele_frequency_ratio(pt::idx_t alternate_count,
				   pt::idx_t total_alleles)
{
	return std::to_string(alternate_count) + "/" +
	       std::to_string(total_alleles);
}

std::optional<std::string> maybe_node_id(
	const std::map<const pvst::VertexBase *, std::string> &node_ids,
	const pvst::VertexBase *vertex)
{
	if (vertex == nullptr)
		return std::nullopt;
	auto it = node_ids.find(vertex);
	if (it == node_ids.end())
		return std::nullopt;
	return it->second;
}

void write_flubble_debug(std::ostream &out, const core::config &app_config)
{
	write_key(out, "flubble_debug");
	out << '{';
	write_key_string(out, "schema", "povu.flubble-debug.v1");
	out << ',';
	write_key(out, "frames");
	out << '[';

	bool first = true;
	if (app_config.has_structure_export_path()) {
		const fs::path sidecar = povu::flubbles::debug_sidecar_path(
			*app_config.get_structure_export_path());
		std::ifstream in(sidecar);
		std::string line;
		while (std::getline(in, line)) {
			if (line.empty())
				continue;
			if (!first)
				out << ',';
			first = false;
			out << line;
		}
	}

	out << ']';
	out << '}';
}
} // namespace

Writer::Writer(const core::config &app_config, const bd::VG &graph,
	       const std::vector<pvst::Tree> &pvsts)
    : graph_(graph), pvsts_(pvsts)
{
	const auto &path = app_config.get_structure_export_path();
	if (!path.has_value()) {
		throw std::logic_error(
			"structure export writer requires an output path");
	}

	out_.open(*path);
	if (!out_.is_open()) {
		throw std::runtime_error("could not open structure export path: " +
					 path->string());
	}

	index_pvst_nodes();
	write_prefix(app_config);
}

void Writer::index_pvst_nodes()
{
	for (std::size_t component_idx{}; component_idx < pvsts_.size();
	     ++component_idx) {
		const pvst::Tree &tree = pvsts_[component_idx];
		for (std::size_t node_idx{}; node_idx < tree.vtx_count();
		     ++node_idx) {
			pvst_node_ids_[tree.get_vertex_const_ptr(node_idx)] =
				node_id(component_idx, node_idx);
		}
	}
}

std::string Writer::node_id(std::size_t component_idx,
			    std::size_t node_idx) const
{
	return std::to_string(component_idx + 1) + ":" +
	       std::to_string(node_idx);
}

std::string Writer::node_id(const pvst::VertexBase *vertex) const
{
	auto id = maybe_node_id(pvst_node_ids_, vertex);
	return id.value_or("");
}

void Writer::write_prefix(const core::config &app_config)
{
	out_ << '{';
	write_key_string(out_, "schema", "povu.lean4.structure.v1");
	out_ << ',';
	write_key_string(out_, "producer", "current-povu-cli");
	out_ << ',';
	write_accepted_gfa(app_config);
	out_ << ',';

	write_key(out_, "reference_prefixes");
	out_ << '[';
	const std::vector<std::string> &prefixes =
		app_config.get_ref_name_prefixes();
	for (std::size_t idx{}; idx < prefixes.size(); ++idx) {
		if (idx > 0)
			out_ << ',';
		write_json_string(out_, prefixes[idx]);
	}
	out_ << ']';
	out_ << ',';

	write_pvst_boundary_candidates();
	out_ << ',';
	write_pvst_nodes();
	out_ << ',';
	write_flubble_debug(out_, app_config);
	out_ << ',';
	write_key(out_, "variant_calls");
	out_ << '[';
}

void Writer::write_accepted_gfa(const core::config &app_config)
{
	write_key(out_, "accepted_gfa");
	out_ << '{';
	write_key_string(out_, "input_name",
			 fs::path(app_config.get_input_gfa()).filename().string());
	out_ << ',';

	write_key(out_, "segments");
	out_ << '[';
	for (std::size_t v_idx{}; v_idx < graph_.vtx_count(); ++v_idx) {
		if (v_idx > 0)
			out_ << ',';
		const bd::Vertex &vertex = graph_.get_vertex_by_idx(v_idx);
		out_ << '{';
		write_key(out_, "order");
		out_ << v_idx << ',';
		write_key(out_, "id");
		out_ << vertex.id() << ',';
		write_key_string(out_, "sequence", vertex.get_label());
		out_ << '}';
	}
	out_ << ']';
	out_ << ',';

	write_key(out_, "links");
	out_ << '[';
	for (std::size_t edge_idx{}; edge_idx < graph_.edge_count();
	     ++edge_idx) {
		if (edge_idx > 0)
			out_ << ',';
		const bd::Edge &edge = graph_.get_edge(edge_idx);
		out_ << '{';
		write_key(out_, "order");
		out_ << edge_idx << ',';
		write_key(out_, "from");
		out_ << graph_.v_idx_to_id(edge.get_v1_idx()) << ',';
		write_key_string(out_, "from_orient",
				 gfa_from_orientation(edge.get_v1_end()));
		out_ << ',';
		write_key(out_, "to");
		out_ << graph_.v_idx_to_id(edge.get_v2_idx()) << ',';
		write_key_string(out_, "to_orient",
				 gfa_to_orientation(edge.get_v2_end()));
		out_ << ',';
		write_key_string(out_, "overlap", "0M");
		out_ << '}';
	}
	out_ << ']';
	out_ << ',';

	write_key(out_, "paths");
	out_ << '[';
	for (pt::idx_t ref_id{}; ref_id < graph_.get_hap_count(); ++ref_id) {
		if (ref_id > 0)
			out_ << ',';
		const pr::Ref &ref = graph_.get_ref_by_id(ref_id);
		const liteseq::ref_walk *walk = graph_.get_ref_vec(ref_id)->walk;
		out_ << '{';
		write_key(out_, "order");
		out_ << ref_id << ',';
		write_key_string(out_, "name", ref.tag());
		out_ << ',';
		write_key_string(out_, "sample", ref.get_sample_name());
		out_ << ',';
		write_key(out_, "steps");
		out_ << '[';
		for (pt::idx_t step_idx{}; step_idx < ref.get_length();
		     ++step_idx) {
			if (step_idx > 0)
				out_ << ',';
			std::string step;
			step.push_back(pr::lq_strand_to_char(
				walk->strands[step_idx]));
			step += std::to_string(walk->v_ids[step_idx]);
			write_json_string(out_, step);
		}
		out_ << ']';
		out_ << '}';
	}
	out_ << ']';
	out_ << '}';
}

void Writer::write_pvst_boundary_candidates()
{
	write_key(out_, "boundary_candidates");
	out_ << '[';
	std::size_t order{};
	bool first = true;
	for (std::size_t component_idx{}; component_idx < pvsts_.size();
	     ++component_idx) {
		const pvst::Tree &tree = pvsts_[component_idx];
		for (std::size_t node_idx{}; node_idx < tree.vtx_count();
		     ++node_idx) {
			const pvst::VertexBase &node = tree.get_vertex(node_idx);
			std::optional<pvst::route_params_t> route =
				node.get_route_params();
			if (!route.has_value())
				continue;
			if (!first)
				out_ << ',';
			first = false;
			out_ << '{';
			write_key(out_, "order");
			out_ << order++ << ',';
			write_key_string(out_, "kind",
					 "pvst-route-boundary");
			out_ << ',';
			write_key_string(out_, "node_id",
					 node_id(component_idx, node_idx));
			out_ << ',';
			write_key_string(out_, "family",
					 pvst::to_str(node.get_fam()));
			out_ << ',';
			write_key_string(out_, "start",
					 oriented_step(route->start));
			out_ << ',';
			write_key_string(out_, "end", oriented_step(route->end));
			out_ << ',';
			write_key_string(out_, "route", route_name(route->route));
			out_ << ',';
			write_key_string(out_, "boundary",
					 endpoint_label(*route));
			out_ << '}';
		}
	}
	out_ << ']';
}

void Writer::write_pvst_nodes()
{
	write_key(out_, "pvst_nodes");
	out_ << '[';
	std::size_t order{};
	bool first = true;
	for (std::size_t component_idx{}; component_idx < pvsts_.size();
	     ++component_idx) {
		const pvst::Tree &tree = pvsts_[component_idx];
		for (std::size_t node_idx{}; node_idx < tree.vtx_count();
		     ++node_idx) {
			if (!first)
				out_ << ',';
			first = false;
			const pvst::VertexBase &node = tree.get_vertex(node_idx);
			std::optional<pvst::route_params_t> route =
				node.get_route_params();
			out_ << '{';
			write_key(out_, "order");
			out_ << order++ << ',';
			write_key(out_, "component");
			out_ << component_idx + 1 << ',';
			write_key(out_, "local_index");
			out_ << node_idx << ',';
			write_key_string(out_, "node_id",
					 node_id(component_idx, node_idx));
			out_ << ',';
			write_key_string(out_, "family",
					 pvst::to_str(node.get_fam()));
			out_ << ',';
			write_key_string(out_, "type", pvst::to_str(node.get_fam()));
			out_ << ',';
			write_key_string(out_, "label", node.as_str());
			out_ << ',';
			write_key(out_, "start");
			write_nullable_string(
				out_, route.has_value()
					      ? std::optional<std::string>(
							oriented_step(route->start))
					      : std::nullopt);
			out_ << ',';
			write_key(out_, "end");
			write_nullable_string(
				out_, route.has_value()
					      ? std::optional<std::string>(
							oriented_step(route->end))
					      : std::nullopt);
			out_ << ',';
			write_key(out_, "route");
			write_nullable_string(
				out_, route.has_value()
					      ? std::optional<std::string>(
							route_name(route->route))
					      : std::nullopt);
			out_ << ',';
			write_key(out_, "boundary");
			write_nullable_string(
				out_, route.has_value()
					      ? std::optional<std::string>(
							endpoint_label(*route))
					      : std::nullopt);
			out_ << ',';
			write_key(out_, "parent");
			if (node_idx == tree.root_idx())
				out_ << "null";
			else
				write_json_string(
					out_, node_id(component_idx,
						      tree.get_parent_idx(
							      node_idx)));
			out_ << ',';
			write_key(out_, "children");
			out_ << '[';
			const std::vector<pt::idx_t> &children =
				tree.get_children(node_idx);
			for (std::size_t child_idx{}; child_idx < children.size();
			     ++child_idx) {
				if (child_idx > 0)
					out_ << ',';
				write_json_string(
					out_, node_id(component_idx,
						      children[child_idx]));
			}
			out_ << ']';
			out_ << ',';
			write_key(out_, "tree_depth");
			out_ << node.get_height();
			out_ << '}';
		}
	}
	out_ << ']';
}

void Writer::write_variant_calls(iv::VcfRecIdx &vcf_recs)
{
	for (auto &[ref_id, recs] : vcf_recs.get_recs_mut()) {
		std::string chrom = graph_.get_ref_by_id(ref_id).tag();
		for (const iv::VcfRec &record : recs)
			write_variant_call(chrom, record);
	}
}

void Writer::write_variant_call(const std::string &chrom,
				const iv::VcfRec &record)
{
	if (wrote_variant_)
		out_ << ',';
	wrote_variant_ = true;

	const pvst::VertexBase *source = record.get_source_pvst_vtx();
	std::optional<std::string> source_node_id =
		maybe_node_id(pvst_node_ids_, source);
	const bool is_subr = record.get_var_type() == ir::var_type_e::subr;
	const pt::idx_t an = record.get_an();

	out_ << '{';
	write_key(out_, "order");
	out_ << next_variant_order_++ << ',';
	write_key(out_, "source");
	out_ << '{';
	write_key_string(out_, "kind", is_subr ? "hairpin" : "flubble");
	out_ << ',';
	write_key(out_, "node_id");
	write_nullable_string(out_, source_node_id);
	out_ << ',';
	write_key(out_, "family");
	write_nullable_string(
		out_, source == nullptr
			      ? std::nullopt
			      : std::optional<std::string>(
					pvst::to_str(source->get_fam())));
	out_ << ',';
	write_key_string(out_, "endpoint_id",
			 is_subr ? record.get_id() : record.get_enc_flubble());
	out_ << ',';
	write_key(out_, "enclosing_site");
	if (is_subr)
		out_ << "null";
	else
		write_json_string(out_, record.get_enc_flubble());
	out_ << ',';
	write_key(out_, "level");
	if (is_subr)
		out_ << "null";
	else
		out_ << (record.get_height() > 0 ? record.get_height() - 1 : 0);
	out_ << '}';
	out_ << ',';

	write_key_string(out_, "chrom", chrom);
	out_ << ',';
	write_key(out_, "contig_order");
	out_ << record.get_ref_id() << ',';
	write_key(out_, "pos");
	out_ << record.get_pos() << ',';
	write_key_string(out_, "id", record.get_id());
	out_ << ',';
	write_key_string(out_, "ref", record.get_ref_as_dna_str(graph_));
	out_ << ',';
	write_key_string(out_, "ref_traversal",
			 record.get_ref_slice().as_str(record.get_var_type()));
	out_ << ',';

	write_key(out_, "alternates");
	out_ << '[';
	for (pt::idx_t alt_idx{}; alt_idx < record.get_alt_allele_group_count();
	     ++alt_idx) {
		if (alt_idx > 0)
			out_ << ',';
		const pt::idx_t count = record.get_alt_allele_count(alt_idx);
		out_ << '{';
		write_key(out_, "index");
		out_ << alt_idx + 1 << ',';
		write_key_string(out_, "dna",
				 record.get_alt_as_dna_str(graph_, alt_idx));
		out_ << ',';
		write_key_string(out_, "traversal",
				 record.get_alt_as_str(alt_idx));
		out_ << ',';
		write_key(out_, "count");
		out_ << count;
		out_ << '}';
	}
	out_ << ']';
	out_ << ',';

	write_key_string(out_, "variant_type",
			 ir::to_string_view(record.get_var_type()));
	out_ << ',';
	write_key(out_, "tangled");
	out_ << (record.is_tangled() ? "true" : "false") << ',';
	write_key_string(out_, "qual", record.get_qual());
	out_ << ',';
	write_key_string(out_, "filter", record.get_filter());
	out_ << ',';
	write_key_string(out_, "format", record.get_format());
	out_ << ',';
	write_key(out_, "reference_allele_count");
	out_ << record.get_reference_allele_count() << ',';

	write_key(out_, "ac");
	out_ << '[';
	for (pt::idx_t alt_idx{}; alt_idx < record.get_alt_allele_group_count();
	     ++alt_idx) {
		if (alt_idx > 0)
			out_ << ',';
		out_ << record.get_alt_allele_count(alt_idx);
	}
	out_ << ']';
	out_ << ',';

	write_key(out_, "af");
	out_ << '[';
	for (pt::idx_t alt_idx{}; alt_idx < record.get_alt_allele_group_count();
	     ++alt_idx) {
		if (alt_idx > 0)
			out_ << ',';
		write_json_string(
			out_, allele_frequency_ratio(
				      record.get_alt_allele_count(alt_idx), an));
	}
	out_ << ']';
	out_ << ',';
	write_key(out_, "an");
	out_ << an << ',';
	write_key(out_, "ns");
	out_ << record.get_ns() << ',';

	write_key(out_, "genotypes");
	out_ << '[';
	const std::vector<std::string> &sample_names =
		graph_.get_genotype_col_names();
	const std::vector<std::vector<std::string>> &genotype_cols =
		record.get_genotype_cols();
	for (std::size_t sample_idx{}; sample_idx < sample_names.size();
	     ++sample_idx) {
		if (sample_idx > 0)
			out_ << ',';
		out_ << '{';
		write_key_string(out_, "sample", sample_names[sample_idx]);
		out_ << ',';
		write_key_string(out_, "value",
				 join_genotype_column(genotype_cols[sample_idx]));
		out_ << '}';
	}
	out_ << ']';
	out_ << '}';
}

void Writer::finish()
{
	if (finished_)
		return;
	out_ << "]}";
	out_.close();
	finished_ = true;
}

} // namespace mto::to_structure_export
