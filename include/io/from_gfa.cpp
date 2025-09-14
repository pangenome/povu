#include "./from_gfa.hpp"


namespace povu::io::from_gfa {

inline lq::gfa_config gen_lq_conf(const core::config &app_config, std::string &gfa_fp) {
  gfa_fp = app_config.get_input_gfa();
  pt::idx_t ref_count = 0;

  bool read_all_refs = app_config.inc_refs() && app_config.inc_vtx_labels();

  lq::gfa_config_cpp lq_conf(
      gfa_fp.c_str(), // file path
      app_config.inc_vtx_labels(), // include vertex labels
      app_config.inc_refs(), // include references
      read_all_refs, // read all references
      ref_count, // reference count
      NULL // reference names
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
bd::VG *to_bd(const core::config &app_config) {

  bool show_prog = app_config.show_progress();

  pv_prog::IndeterminateProgressBar prep_bar;
  prep_bar.set_option(indicators::option::PostfixText{"Preparing GFA..."});
  set_progress_bar_ind(&prep_bar);


  /* initialize a liteseq gfa */
  std::vector<const char *> refs;
  std::string gfa_fp;
  lq::gfa_config conf = gen_lq_conf(app_config, gfa_fp);
  lq::gfa_props *gfa = nullptr;
  //lq::gfa_props *gfa = lq::gfa_new(&conf);
  std::thread get_gfa_async([&]() {
    gfa = lq::gfa_new(&conf);
    if (show_prog) {
      prep_bar.set_option(indicators::option::PostfixText{"GFA prepared."});
      prep_bar.mark_as_completed();
    }
  });

  if (show_prog) {
    while (!prep_bar.is_completed()) {
      prep_bar.tick();
      std::this_thread::sleep_for(std::chrono::milliseconds(100));
    }
  }

  get_gfa_async.join();


  pt::idx_t vtx_count = gfa->s_line_count;
  pt::idx_t edge_count = gfa->l_line_count;
  pt::idx_t ref_count = gfa->p_line_count;

  /* initialize a povu bidirected graph */
  bd::VG *vg = new bd::VG(vtx_count, edge_count, ref_count);

  /* set up progress bars */
  const std::size_t PROG_BAR_COUNT {3};
  ProgressBar vtx_bar, edge_bar, ref_bar;
  {
    set_progress_bar_common_opts(&vtx_bar, vtx_count);
    set_progress_bar_common_opts(&edge_bar, edge_count);
    set_progress_bar_common_opts(&ref_bar, ref_count);
  }
  MultiProgress<ProgressBar, PROG_BAR_COUNT> bars(vtx_bar, edge_bar, ref_bar);
  const size_t VTX_BAR_IDX{0};
  const size_t EDGE_BAR_IDX{1};
  const size_t REF_BAR_IDX{2};

  /* add vertices */
  for (size_t i{}; i < vtx_count; ++i) {

    if (show_prog) { // update progress bar
      std::string prog_msg = fmt::format("Loading vertices ({}/{})", i + 1, vtx_count);
      vtx_bar.set_option(indicators::option::PostfixText{prog_msg});
      bars.set_progress<VTX_BAR_IDX>(static_cast<size_t>(i + 1));
    }

    std::size_t v_id = gfa->v[i].id;
    std::string label = app_config.inc_vtx_labels() ? std::string(gfa->v[i].seq) : std::string();
    vg->add_vertex(v_id, label);
  }

  /* add edges */
  for (std::size_t i{}; i < edge_count; ++i) {
    if (show_prog) { // update progress bar
      std::string prog_msg = fmt::format("Loading edges ({}/{})", i + 1, edge_count);
      edge_bar.set_option(indicators::option::PostfixText{prog_msg});
      bars.set_progress<EDGE_BAR_IDX>(static_cast<size_t>(i + 1));
    }

    std::size_t v1 = gfa->e[i].v1_id;
    pgt::v_end_e v1_end = gfa->e[i].v1_side == lq::vtx_side_e::LEFT ? pgt::v_end_e::l : pgt::v_end_e::r;

    std::size_t v2 = gfa->e[i].v2_id;
    pgt::v_end_e v2_end = gfa->e[i].v2_side == lq::vtx_side_e::LEFT ? pgt::v_end_e::l : pgt::v_end_e::r;

    vg->add_edge(v1, v1_end, v2, v2_end);
  }

  /* add references if necessary */
  // TODO: to parallise run in parallel for each vertex
  if (app_config.inc_refs()) {

    pt::idx_t path_pos {}; // the position of a base in a reference path
    pt::id_t curr_ref_id{pc::INVALID_ID};

    // TODO: [C] set up ref names independently?

    for (pt::idx_t ref_idx{}; ref_idx < ref_count; ++ref_idx) {

      if (show_prog) { // update progress bar
        std::string prog_msg =fmt::format("Loading references ({}/{})", ref_idx + 1, ref_count);
        ref_bar.set_option(indicators::option::PostfixText{prog_msg});
        bars.set_progress<REF_BAR_IDX>(static_cast<size_t>(ref_idx + 1));
      }

      const std::string &label = gfa->refs[ref_idx].name;

      char delim = '#';
      curr_ref_id = vg->add_ref(label, delim);
      path_pos = 1; // this is 1 indexed

      pgt::ref_walk_t &ref_vector = vg->get_ref_vec_mut(curr_ref_id);
      ref_vector.reserve(gfa->refs[ref_idx].step_count);

      // color each vertex in the path
      for (pt::idx_t step_idx{}; step_idx < gfa->refs[ref_idx].step_count; ++step_idx) {
        pt::id_t v_id = gfa->refs[ref_idx].steps[step_idx].v_id;
        lq::strand_e s = gfa->refs[ref_idx].steps[step_idx].s;
        pgt::or_e o = (s == lq::strand_e::FORWARD) ? pgt::or_e::forward : pgt::or_e::reverse;

        // add step to the reference walk
        pt::idx_t idx_in_ref_vec = static_cast<pt::idx_t>(ref_vector.size());
        ref_vector.emplace_back(pgt::ref_step_t{v_id, o, path_pos});
        vg->set_vtx_ref_idx(v_id, curr_ref_id, idx_in_ref_vec);

        bd::Vertex& v = vg->get_vertex_mut_by_id(v_id);
        path_pos += v.get_label().length();
      }

      // set the length of the reference
      pr::Ref &ref = vg->get_ref_by_id_mut(curr_ref_id);
      ref.set_length(path_pos - 1);
    }

    vg->gen_genotype_metadata();
  }

  gfa_free(gfa); // very important free the gfa props

  /* populate tips */
  for (std::size_t v_idx{}; v_idx < vg->vtx_count(); ++v_idx) {
    const bd::Vertex &v = vg->get_vertex_by_idx(v_idx);

    if (v.get_edges_l().empty() && v.get_edges_r().empty()) {
      if (app_config.verbosity() > 2) {
        WARN("isolated vertex {}", v.id());
      }
      vg->add_tip(v.id(), pgt::v_end_e::l);
    }
    else if (v.get_edges_l().empty()) { vg->add_tip(v.id(), pgt::v_end_e::l); }
    else if (v.get_edges_r().empty()) { vg->add_tip(v.id(), pgt::v_end_e::r); }
  }

  return vg;
}
};
