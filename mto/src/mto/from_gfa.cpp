#include <cstddef> // for size_t
#include <set>	   // for set
#include <string>  // for basic_string
#include <thread>  // for thread, sleep_for
#include <vector>  // for vector

#include <liteseq/gfa.h>	 // for gfa_config, gfa...
#include <liteseq/refs.h>	 // for get_step_count
#include <log/log.h>		 // for log_warn, log_info
#include <quilt/graph_types.hpp> // for v_end_e, side_n_id_t, side_n_idx_t
#include <quilt/types.hpp>	 // for qt

#include "mto/from_gfa.hpp"

namespace mto::from_gfa
{
namespace lq = liteseq;
namespace pgt = quilt::types::graph;

inline lq::gfa_config gen_lq_conf(const core::config &app_config,
				  std::string &gfa_fp)
{
	gfa_fp = app_config.get_input_gfa();

	lq::gfa_config_cpp lq_conf(
		gfa_fp.c_str(),		     // file path
		app_config.inc_vtx_labels(), // include vertex labels
		app_config.inc_refs()	     // include references
	);

	return lq_conf;
}

/**
 * Read GFA into a variation graph represented as a bidirected graph
 *
 * @param [in] filename The GFA file to read
 * @param [in] app_config The application configuration
 * @return A VariationGraph object from the GFA file
 */
bd::VG *to_bd(const core::config &app_config)
{
	if (app_config.verbosity() > 1)
		log_info("Processing GFA");

	/* initialize a liteseq gfa */
	std::vector<const char *> refs;
	std::string gfa_fp;
	lq::gfa_config conf = gen_lq_conf(app_config, gfa_fp);
	lq::gfa_props *gfa = nullptr;

	std::thread get_gfa_async([&]() { gfa = lq::gfa_new(&conf); });

	get_gfa_async.join();

	qt::idx_t vtx_count = gfa->vtx_arr_size;
	qt::idx_t edge_count = gfa->l_line_count;
	qt::idx_t ref_count = gfa->ref_count;

	/* initialize a povu bidirected graph */
	auto vg = new bd::VG(gfa); // vg is bd::VG *

	/* add vertices */
	if (app_config.verbosity() > 0)
		log_info("Adding Vertices");

	for (qt::idx_t i{}; i < vtx_count; ++i) {
		lq::vtx *v = lq::get_vtx(gfa, i);
		if (v == nullptr) // skip uninitialized vertices
			continue;

		std::size_t v_id = v->id;
		std::string label = app_config.inc_vtx_labels()
					    ? std::string(v->seq)
					    : std::string();
		vg->add_vertex(v_id, label);
	}

	/* add edges */
	if (app_config.verbosity() > 1)
		log_info("Adding edges");

	for (std::size_t i{}; i < edge_count; ++i) {
		std::size_t v1 = gfa->e[i].v1_id;
		pgt::v_end_e v1_end = gfa->e[i].v1_side == lq::vtx_side_e::LEFT
					      ? pgt::v_end_e::l
					      : pgt::v_end_e::r;

		std::size_t v2 = gfa->e[i].v2_id;
		pgt::v_end_e v2_end = gfa->e[i].v2_side == lq::vtx_side_e::LEFT
					      ? pgt::v_end_e::l
					      : pgt::v_end_e::r;

		vg->add_edge(v1, v1_end, v2, v2_end);
	}

	/* refs */
	if (app_config.inc_refs()) {
		if (app_config.verbosity() > 1)
			log_info("Adding refs");

		vg->set_refs_meta(gfa->refs, ref_count);

		for (qt::idx_t ref_idx{}; ref_idx < ref_count; ref_idx++) {
			lq::ref *ref = lq::get_ref(gfa, ref_idx);
			qt::idx_t N = lq::get_step_count(ref);
			for (qt::idx_t step_idx{}; step_idx < N; step_idx++) {
				qt::id_t v_id = ref->walk->v_ids[step_idx];
				vg->set_vtx_ref_idx(v_id, ref_idx, step_idx);
			}
		}
	}

	/* populate tips */
	if (app_config.verbosity() > 1)
		log_info("Processing tips");

	for (std::size_t v_idx{}; v_idx < vg->vtx_count(); ++v_idx) {
		const bd::Vertex &v = vg->get_vertex_by_idx(v_idx);

		// TODO: [c] improve logic on handling isolated vertices
		if (v.get_edges_l().empty() && v.get_edges_r().empty()) {
			if (app_config.verbosity() > 2)
				log_warn("isolated vertex %ul", v.id());
			vg->add_tip(v.id(), pgt::v_end_e::l);
		}
		else if (v.get_edges_l().empty()) {
			vg->add_tip(v.id(), pgt::v_end_e::l);
		}
		else if (v.get_edges_r().empty()) {
			vg->add_tip(v.id(), pgt::v_end_e::r);
		}
	}

	return vg;
}
}; // namespace mto::from_gfa
