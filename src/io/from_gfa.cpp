#include "./from_gfa.hpp"


namespace povu::io::from_gfa {

lq::gfa_config gen_lq_conf(const core::config &app_config,
                           std::vector<const char *> &refs,
                           std::string &gfa_fp) {
  gfa_fp = app_config.get_input_gfa();
  pt::idx_t ref_count = 0;

  bool read_all_refs = app_config.inc_refs() && app_config.inc_vtx_labels();

  // if (app_config.inc_refs()) {
  //   ref_count = app_config.get_reference_paths().size();
  //   refs.reserve(ref_count);
  //   for (const std::string &r : app_config.get_reference_paths()) {
  //     refs.push_back(r.c_str());
  //   }
  // }


  lq::gfa_config_cpp lq_conf(
      gfa_fp.c_str(), // file path
      app_config.inc_vtx_labels(), // include vertex labels
      app_config.inc_refs(), // include references
      read_all_refs, // read all references
      ref_count, // reference count
      NULL // reference names
  );

  // lq::gfa_config conf = {.fp = gfa_fp.c_str(),
  //                        .inc_vtx_labels = app_config.inc_vtx_labels(),
  //                        .inc_refs = app_config.inc_refs(),
  //                        .read_all_refs = true,
  //                        .ref_count = 0,
  //                        .ref_names = NULL,
  // };

  return lq_conf;
}



/**
 * Read GFA into a variation graph represented as a bidirected graph
 *
 * @param [in] filename The GFA file to read
 * @param [in] app_config The application configuration
 * @return A VariationGraph object from the GFA file
 */
bd::VG *to_bd(const core::config& app_config) {
  std::string fn_name { pv_cmp::format("{}::{}]", MODULE, __func__) };

  /* initialize a liteseq gfa */
  std::vector<const char *> refs;
  std::string gfa_fp;
  lq::gfa_config conf = gen_lq_conf(app_config, refs, gfa_fp);
  lq::gfa_props *gfa = lq::gfa_new(&conf);

  pt::idx_t vtx_count = gfa->s_line_count;
  pt::idx_t edge_count = gfa->l_line_count;
  pt::idx_t ref_count = gfa->p_line_count;

  /* initialize a povu bidirected graph */
  bd::VG *vg = new bd::VG(vtx_count, edge_count, app_config.inc_refs());

  /* add vertices */
  for (size_t i {}; i < vtx_count; ++i) {
    std::size_t v_id = gfa->v[i].id;
    std::string label = app_config.inc_vtx_labels() ? std::string(gfa->v[i].seq) : std::string();
    vg->add_vertex(v_id, label);
  }

  /* add edges */
  for (std::size_t i {}; i < edge_count; ++i) {
    std::size_t v1 = gfa->e[i].v1_id;
    pgt::v_end_e v1_end = gfa->e[i].v1_side == lq::vtx_side_e::LEFT ? pgt::v_end_e::l : pgt::v_end_e::r;

    std::size_t v2 = gfa->e[i].v2_id;
    pgt::v_end_e v2_end = gfa->e[i].v2_side == lq::vtx_side_e::LEFT ? pgt::v_end_e::l : pgt::v_end_e::r;

    vg->add_edge(v1, v1_end, v2, v2_end);
  }

  /* add references if necessary */
  // TODO: to parallise run in parallel for each vertex
  if (app_config.inc_refs()) {

    //std::cerr << "inc ref\n";

    std::size_t path_pos {}; // the position of a base in a reference path
    pt::id_t curr_ref_id {pc::INVALID_ID};
    for (pt::idx_t ref_idx {}; ref_idx < ref_count; ++ref_idx) {
      const std::string &label = gfa->refs[ref_idx].name;

      //std::cerr << "reading ref " << label << "\n";

      char delim = '#';
      curr_ref_id = vg->add_ref(label, delim);
      //pt::id_t vg_ref_id = vg->add_ref(label, delim);
      path_pos = 1; // this is 1 indexed

      // color each vertex in the path
      for (pt::idx_t step_idx{}; step_idx < gfa->refs[ref_idx].step_count; ++step_idx) {
        pt::id_t v_id = gfa->refs[ref_idx].steps[step_idx].v_id;
        lq::strand_e s = gfa->refs[ref_idx].steps[step_idx].s;
        pgt::or_e o = (s == lq::strand_e::FORWARD) ? pgt::or_e::forward : pgt::or_e::reverse;

        bd::Vertex& v = vg->get_vertex_mut_by_id(v_id);
        v.add_ref(curr_ref_id, o, path_pos);
        path_pos += v.get_label().length();
      }

      // set the length of the reference
      pgt::Ref &ref = vg->get_ref_by_id_mut(curr_ref_id);
      ref.set_length(path_pos - 1);
    }
  }

  gfa_free(gfa); // very important free the gfa props

  /* populate tips */
  for (std::size_t v_idx{}; v_idx < vg->vtx_count(); ++v_idx) {
    const bd::Vertex &v = vg->get_vertex_by_idx(v_idx);

    if (v.get_edges_l().empty() && v.get_edges_r().empty()) {
      if (app_config.verbosity() > 2 ){
        std::cerr << pv_cmp::format(" {} WARN isolated node {} \n", fn_name, v.id());
      }
      vg->add_tip(v.id(), pgt::v_end_e::l);
    }
    else if (v.get_edges_l().empty()) { vg->add_tip(v.id(), pgt::v_end_e::l); }
    else if (v.get_edges_r().empty()) { vg->add_tip(v.id(), pgt::v_end_e::r); }
  }

  return vg;
}

};
