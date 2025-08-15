#include "./variants.hpp"


namespace povu::variants {

/**
 * Remove walks that are prefixes of other walks in the same Itn
 * This is done to avoid redundant walks in the RoV
 * A walk is a prefix of another walk if it has fewer steps and starts with the same step
 */
void remove_prefix_walks(pvt::Itn &itn) {
  std::string fn_name{pv_cmp::format("[{}::{}]", MODULE, __func__)};

  // we add walk idx which are prefixes is_prefix
  std::set<pt::idx_t> to_remove;

  for (pt::idx_t qry_w_idx {}; qry_w_idx < itn.at_count(); ++qry_w_idx)    {
    const pvt::AW &qry_aw = itn.get_at(qry_w_idx);
    for (pt::idx_t txt_w_idx = {}; txt_w_idx < itn.at_count(); ++txt_w_idx) {
      const pvt::AW &txt_aw = itn.get_at(txt_w_idx);

      if(qry_aw.step_count() < txt_aw.step_count() && qry_aw.get_steps().front() == txt_aw.get_steps().front()) {
        to_remove.insert(qry_w_idx);
      }
    }
  }

  // remove the walks that are prefixes
  itn.remove_aws(to_remove);
}

/**
 * Associate walks in an RoV with references
 */
void gen_rov_ref_walks(const bd::VG &g, const pvt::RoV &rov, std::vector<pvt::Exp> &ref_walks_vec) {

  const std::vector<pvt::walk_t> &walks = rov.get_walks();
  const pvst::VertexBase *pvst_v_ptr = rov.get_pvst_vtx();

  // create an expedition object for the RoV
  pvt::Exp ref_walks {pvst_v_ptr};

  // a walk is a single traversal bounded by start to the end of an RoV
  for (pt::idx_t w_idx{}; w_idx < rov.walk_count(); w_idx++) {
    const pvt::walk_t &w = walks[w_idx];
    pgu::variants::comp_itineraries(g, w, w_idx, ref_walks);
  }

  for (pt::idx_t ref_id : ref_walks.get_ref_ids()) {
    remove_prefix_walks(ref_walks.get_itn_mut(ref_id));
  }

  for (const pt::idx_t ref_id : ref_walks.get_ref_ids()) {
    if (ref_walks.get_itn(ref_id).at_count() > 1) {
      // if any ref has more than one walk, the RoV is tangled
      ref_walks.set_tangled(true);
      // no need to check other refs, we know the RoV is tangled
      break;
    }
  }

  for (pt::id_t ref_id : ref_walks.get_ref_ids()) {
    pvt::Itn &itn = ref_walks.get_itn_mut(ref_id);
    for (pvt::AW &aw : itn.get_ats_mut()) {
      aw.add_ref_id(ref_id); // add the ref id to each walk
    }
  }

  if (ref_walks.is_tangled()) {
    put::untangle_ref_walks(ref_walks);
  }


  //{ 
  //   std::cerr << fn_name
  //             << " " << pvst_v_ptr->as_str()
  //             << " ref count " << ref_walks.get_ref_ids().size()
  //             << "\n";

  //   for (pt::idx_t ref_id : ref_walks.get_ref_ids()) {
  //     std::cerr << fn_name << " " << g.get_ref_label(ref_id) << "\n";

  //     auto itn = ref_walks.get_itn(ref_id);
  //     for (pt::idx_t i{}; i < itn.at_count(); ++i) {
  //       const pvt::AW &aw = itn.get_at(i);
  //       std::cerr << " aw: " << aw.as_str() << "\n";
  //       std::cerr << " ref count " << aw.get_ref_ids().size() << "\n";
  //     }
  //   }
  // }

  //  if (pvst_v_ptr->as_str() == ">11>13") {
  //    exit(1);
  //  }

  ref_walks_vec.push_back(std::move(ref_walks));

  return;
}


/**
 * Check if a vertex in the pvst is a flubble leaf
 * A flubble leaf is a vertex that has no children that are also flubbles
 */
bool is_fl_leaf(const pvtr::Tree &pvst, pt::idx_t pvst_v_idx) noexcept {
  const pvst::VertexBase *pvst_v_ptr = pvst.get_vertex_const_ptr(pvst_v_idx);

  if (!pvst::is_fl_like(pvst_v_ptr->get_type())) {
    return false; // not a flubble
  }

  for (pt::idx_t v_idx : pvst.get_children(pvst_v_idx)) {
    if (pvst::is_fl_like(pvst.get_vertex_const_ptr(v_idx)->get_type())) {
      return false;
    }
   }

  return true;
}

/**
 * find walks in the graph based on the leaves of the pvst
 * initialize RoVs from flubbles
 */
std::vector<pvt::RoV> gen_rov(const std::vector<pvtr::Tree> &pvsts, const bd::VG &g) {
  // the set of RoVs to return
  std::vector<pvt::RoV> rs;
  rs.reserve(pvsts.size());

  for (const pvtr::Tree &pvst : pvsts) { // for each pvst
    // loop through each tree

    for (pt::idx_t pvst_v_idx{}; pvst_v_idx < pvst.vtx_count(); pvst_v_idx++) {
      // call variants on the leaves only
      const pvst::VertexBase *pvst_v_ptr = pvst.get_vertex_const_ptr(pvst_v_idx);

      pvst::traversal_params_t p = pvst_v_ptr->get_traversal_params();

      if (p.traversable && (is_fl_leaf(pvst, pvst_v_idx) )) {
        pvt::RoV r { pvst_v_ptr };

        // get the set of walks for the RoV
        pgu::graph::find_walks(g, r);

        rs.push_back(std::move(r));
      }
    }
  }

  return rs;
}

pvt::VcfRecIdx gen_vcf_rec_map(const std::vector<pvtr::Tree> &pvsts, const bd::VG &g) {

  std::vector<pvt::RoV> all_rovs = gen_rov(pvsts, g);

  std::vector<pvt::Exp> all_ref_walks;
  all_ref_walks.reserve(all_rovs.size());
  for (const pvt::RoV &r : all_rovs) {
    gen_rov_ref_walks(g, r, all_ref_walks);
  }

  pvt::VcfRecIdx rs = pgv::gen_vcf_records(g, all_ref_walks);

  return rs;
}
} // namespace povu::genomics
