#include "./io.hpp"
#include <cstddef>

namespace povu::io::from_gfa {

/**
 * Read GFA into a variation graph represented as a bidirected graph
 *
 * @param [in] filename The GFA file to read
 * @param [in] app_config The application configuration
 * @return A VariationGraph object from the GFA file
 */
bd::VG *to_bd(const char* filename, const core::config& app_config) {
  std::string fn_name { std::format("[povu::io::{}]", __func__) };

  /* initialize a liteseq graph */
  lq::vg *ls_g = lq::vg_new();
  lq::vg_props ls_cfg = {
    .inc_vtx_labels = app_config.inc_vtx_labels(),
    .inc_refs = app_config.inc_refs()
  };

  lq::gfa_to_vg(filename, ls_g, &ls_cfg); // read the GFA file into a liteseq graph

  std::size_t vtx_count = ls_g->vtx_count;
  std::size_t edge_count = ls_g->edge_count;
  std::size_t ref_count = ls_g->ref_count;

  bd::VG *vg = new bd::VG(vtx_count, edge_count); // initialize a bidirected graph

  //bd::VG vg(vtx_count, edge_count); // initialize a bidirected graph

  /* add vertices */
  for (size_t i {}; i < ls_g->vtx_count; ++i) {
    std::size_t v_id = ls_g->v[i].id;
    std::string label = app_config.inc_vtx_labels() ? ls_g->v[i].seq : std::string();
    vg->add_vertex(v_id, label);
  }

  /* add edges */
  for (std::size_t i {}; i < ls_g->edge_count; ++i) {
    std::size_t v1 = ls_g->e[i].v1_id;
    pgt::v_end_e v1_end = ls_g->e[i].v1_side == lq::vtx_side::LEFT ? pgt::v_end_e::l : pgt::v_end_e::r;

    std::size_t v2 = ls_g->e[i].v2_id;
    pgt::v_end_e v2_end = ls_g->e[i].v2_side == lq::vtx_side::LEFT ? pgt::v_end_e::l : pgt::v_end_e::r;

    vg->add_edge(v1, v1_end, v2, v2_end);
  }

  /* add references if necessary */
  // TODO: to parallise run in parallel for each vertex
  if (app_config.inc_refs()) {
    std::size_t path_pos {}; // the position of a base in a reference path
    for (std::size_t ref_idx {}; ref_idx < ref_count; ++ref_idx) {

      vg->add_ref(ls_g->rs->names[ref_idx]);
      path_pos = 1; // this is 1 indexed

      // color each vertex in the path
      for (std::size_t lq_v_idx {}; lq_v_idx < vtx_count; ++lq_v_idx) {
        bool h = lq::vec_has_ref(ls_g->rs->x, lq_v_idx, ref_idx);
        if (!h) { continue; }

        lq::strand s = lq::get_ref_strand(ls_g->rs->s, lq_v_idx, ref_idx);
        pgt::or_t o = (s == lq::strand::FORWARD) ? pgt::or_t::forward : pgt::or_t::reverse;

        bd::Vertex& v = vg->get_vertex_mut_by_id(ls_g->v[lq_v_idx].id);

        v.add_ref(ref_idx, o, path_pos);
        path_pos += v.get_label().length();
      }
    }
  }

  lq::vg_free(ls_g, &ls_cfg); // very important: free the liteseq graph

  /* populate tips */
  for (std::size_t v_idx{}; v_idx < vg->vtx_count(); ++v_idx) {
    const bd::Vertex &v = vg->get_vertex_by_idx(v_idx);

    if (v.get_edges_l().empty() && v.get_edges_r().empty()) {
      if (app_config.verbosity() > 2 ){
        std::cerr << std::format(" {} WARN isolated node {} \n", fn_name, v.id());
      }
      vg->add_tip(v.id(), pgt::v_end_e::l);
    }
    else if (v.get_edges_l().empty()) { vg->add_tip(v.id(), pgt::v_end_e::l); }
    else if (v.get_edges_r().empty()) { vg->add_tip(v.id(), pgt::v_end_e::r); }
  }

  return vg;
}

};
