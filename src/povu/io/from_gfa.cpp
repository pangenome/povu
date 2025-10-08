#include <chrono>	 // for milliseconds
#include <cstddef>	 // for size_t
#include <liteseq/gfa.h> // for gfa_config, gfa...
#include <optional>	 // for optional
#include <set>		 // for set
#include <string>	 // for basic_string
#include <thread>	 // for thread, sleep_for
#include <vector>	 // for vector

#include "fmt/core.h"				     // for format
#include "indicators/indeterminate_progress_bar.hpp" // for IndeterminatePr...
#include "indicators/multi_progress.hpp"	     // for MultiProgress
#include "indicators/progress_bar.hpp"		     // for ProgressBar
#include "indicators/setting.hpp"		     // for PostfixText
#include "liteseq/refs.h"			     // for get_step_count
#include "povu/common/core.hpp"			     // for pt, idx_t, id_t
#include "povu/common/log.hpp"			     // for WARN
#include "povu/common/progress.hpp"		     // for set_progress_ba...
#include "povu/graph/types.hpp"			     // for v_end_e
#include "povu/io/from_gfa.hpp"

namespace povu::io::from_gfa
{
namespace lq = liteseq;
namespace bd = povu::bidirected;
namespace pgt = povu::types::graph;
using namespace povu::progress;

inline lq::gfa_config gen_lq_conf(const core::config &app_config,
				  std::string &gfa_fp)
{
	gfa_fp = app_config.get_input_gfa();
	pt::idx_t ref_count = 0;

	bool read_all_refs =
		app_config.inc_refs() && app_config.inc_vtx_labels();

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
	bool show_prog = app_config.show_progress();

	IndeterminateProgressBar prep_bar;
	prep_bar.set_option(
		indicators::option::PostfixText{"Preparing GFA..."});
	set_progress_bar_ind(&prep_bar);

	/* initialize a liteseq gfa */
	std::vector<const char *> refs;
	std::string gfa_fp;
	lq::gfa_config conf = gen_lq_conf(app_config, gfa_fp);
	lq::gfa_props *gfa = nullptr;

	std::thread get_gfa_async(
		[&]()
		{
			gfa = lq::gfa_new(&conf);
			if (show_prog) {
				prep_bar.set_option(
					indicators::option::PostfixText{
						"GFA prepared."});
				prep_bar.mark_as_completed();
			}
		});

	if (show_prog) {
		while (!prep_bar.is_completed()) {
			prep_bar.tick();
			std::this_thread::sleep_for(
				std::chrono::milliseconds(100));
		}
	}

	get_gfa_async.join();

	pt::idx_t vtx_count = gfa->vtx_arr_size;
	pt::idx_t edge_count = gfa->l_line_count;
	pt::idx_t ref_count = gfa->ref_count;

	/* initialize a povu bidirected graph */
	bd::VG *vg = new bd::VG(gfa);

	/* set up progress bars */
	const std::size_t PROG_BAR_COUNT{3};
	ProgressBar vtx_bar, edge_bar, ref_bar;
	{
		set_progress_bar_common_opts(&vtx_bar, vtx_count);
		set_progress_bar_common_opts(&edge_bar, edge_count);
		set_progress_bar_common_opts(&ref_bar, ref_count);
	}
	MultiProgress<ProgressBar, PROG_BAR_COUNT> bars(vtx_bar, edge_bar,
							ref_bar);
	const size_t VTX_BAR_IDX{0};
	const size_t EDGE_BAR_IDX{1};
	const size_t REF_BAR_IDX{2};

	/* add vertices */
	for (pt::idx_t i{}; i < vtx_count; ++i) {
		if (show_prog) { // update progress bar
			std::string prog_msg = fmt::format(
				"Loading vertices ({}/{})", i + 1, vtx_count);
			vtx_bar.set_option(
				indicators::option::PostfixText{prog_msg});
			bars.set_progress<VTX_BAR_IDX>(
				static_cast<size_t>(i + 1));
		}

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
	for (std::size_t i{}; i < edge_count; ++i) {
		if (show_prog) { // update progress bar
			std::string prog_msg = fmt::format(
				"Loading edges ({}/{})", i + 1, edge_count);
			edge_bar.set_option(
				indicators::option::PostfixText{prog_msg});
			bars.set_progress<EDGE_BAR_IDX>(
				static_cast<size_t>(i + 1));
		}

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
		vg->add_all_refs(gfa->refs, ref_count);

		for (pt::idx_t ref_idx{}; ref_idx < ref_count; ref_idx++) {
			lq::ref *ref = lq::get_ref(gfa, ref_idx);
			pt::idx_t N = lq::get_step_count(ref);
			for (pt::idx_t step_idx{}; step_idx < N; step_idx++) {
				pt::id_t v_id = ref->walk->v_ids[step_idx];
				vg->set_vtx_ref_idx(v_id, ref_idx, step_idx);
			}
		}

		vg->gen_genotype_metadata();
	}

	/* populate tips */
	for (std::size_t v_idx{}; v_idx < vg->vtx_count(); ++v_idx) {
		const bd::Vertex &v = vg->get_vertex_by_idx(v_idx);

		// TODO: [c] improve logic on handling isolated vertices
		if (v.get_edges_l().empty() && v.get_edges_r().empty()) {
			if (app_config.verbosity() > 2)
				WARN("isolated vertex {}", v.id());
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
}; // namespace povu::io::from_gfa
