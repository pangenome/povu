#include "./genomics.hpp"



namespace povu::genomics {

/**
 * Remove walks that are prefixes of other walks in the same Itn
 * This is done to avoid redundant walks in the RoV
 * A walk is a prefix of another walk if it has fewer steps and starts with the same step
 */
void remove_prefix_walks(pga::Itn &itn) {
  std::string fn_name{pv_cmp::format("[{}::{}]", MODULE, __func__)};

  // we add walk idx which are prefixes is_prefix
  std::set<pt::idx_t> to_remove;

  for (pt::idx_t qry_w_idx {}; qry_w_idx < itn.at_count(); ++qry_w_idx)    {
    const pga::AW &qry_aw = itn.get_at(qry_w_idx);
    for (pt::idx_t txt_w_idx = {}; txt_w_idx < itn.at_count(); ++txt_w_idx) {
      const pga::AW &txt_aw = itn.get_at(txt_w_idx);

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
pga::Exp exp_frm_rov(const bd::VG &g, const pgg::RoV &rov) {

  // create an expedition object for the RoV
  const std::vector<pgt::walk_t> &walks = rov.get_walks();
  const pvst::VertexBase *pvst_v_ptr = rov.get_pvst_vtx();
  pga::Exp ref_walks{pvst_v_ptr};

  // compute the itineraries for each reference in the RoV
  std::map<pt::id_t, pga::Itn> &ref_map = ref_walks.get_ref_itns_mut();
  pga::comp_itineraries(g, walks, ref_map);

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
    pga::Itn &itn = ref_walks.get_itn_mut(ref_id);
    for (pga::AW &aw : itn.get_ats_mut()) {
      aw.add_ref_id(ref_id); // add the ref id to each walk
    }
  }

  if (ref_walks.is_tangled()) {
    put::untangle_ref_walks(ref_walks);
  }

  return ref_walks;
}

std::vector<pga::Exp> comp_expeditions(const bd::VG &g,
                                       const std::vector<pgg::RoV> &all_rovs,
                                       std::size_t thread_count) {
  std::vector<pga::Exp> all_exp(all_rovs.size());

  //std::size_t thread_count = 16; // default
  auto [num_threads, chunk_size] = pu::compute_thread_allocation(thread_count, all_rovs.size());

  std::vector<std::thread> threads(num_threads);
  std::size_t start, end;
  for (unsigned int thread_idx{}; thread_idx < num_threads; ++thread_idx) {
    start = thread_idx * chunk_size;
    end = (thread_idx == num_threads - 1) ? all_rovs.size() : (thread_idx + 1) * chunk_size;

    threads[thread_idx] = std::thread([&, start, end]() {
          for (std::size_t i{start}; i < end; i++) {
            const pgg::RoV &r = all_rovs[i];
            pga::Exp rov_rws = exp_frm_rov(g, r);
            all_exp[i] = std::move(rov_rws);
          }});
  }

  // Wait for all threads to finish
  for (auto &thread : threads) {
    thread.join();
  }

  return all_exp;
}

/**
 * Check if a vertex in the pvst is a flubble leaf
 * A flubble leaf is a vertex that has no children that are also flubbles
 */
bool is_fl_leaf(const pvtr::Tree &pvst, pt::idx_t pvst_v_idx) noexcept {
  const pvst::VertexBase *pvst_v_ptr = pvst.get_vertex_const_ptr(pvst_v_idx);

  // we assume that the vertex has a clan
  pvst::vf_e prt_fam = pvst_v_ptr->get_fam();
  if (pvst::to_clan(prt_fam).value() != pvst::vc_e::fl_like) {
    return false; // not a flubble
  }

  for (pt::idx_t v_idx : pvst.get_children(pvst_v_idx)) {
    pvst::vf_e c_fam = pvst.get_vertex_const_ptr(v_idx)->get_fam();
    if (pvst::to_clan(c_fam) == pvst::vc_e::fl_like) {
      return false;
    }
   }

  return true;
}

/**
 * find walks in the graph based on the leaves of the pvst
 * initialize RoVs from flubbles
 */
std::vector<pgg::RoV> gen_rov(const std::vector<pvtr::Tree> &pvsts, const bd::VG &g) {
  // the set of RoVs to return
  std::vector<pgg::RoV> rs;
  rs.reserve(pvsts.size());

  // true when the vertex is a flubble leaf or a leaf in the pvst
  auto should_call = [&](const pvtr::Tree &pvst, const pvst::VertexBase *pvst_v_ptr,
                     pt::idx_t pvst_v_idx) -> bool {
    if (std::optional<pvst::route_params_t> opt_rp = pvst_v_ptr->get_route_params()) {
      return is_fl_leaf(pvst, pvst_v_idx) || pvst.is_leaf(pvst_v_idx);
    }
    return false;
  };

  for (const pvtr::Tree &pvst : pvsts) { // for each pvst
    // loop through each tree

    for (pt::idx_t pvst_v_idx{}; pvst_v_idx < pvst.vtx_count(); pvst_v_idx++) {
      const pvst::VertexBase *pvst_v_ptr = pvst.get_vertex_const_ptr(pvst_v_idx);
      if (should_call(pvst, pvst_v_ptr, pvst_v_idx)) {
        pgg::RoV r{pvst_v_ptr};

        // get the set of walks for the RoV
        povu::genomics::graph::find_walks(g, r);

        if (r.get_walks().size() == 0) {
          // no walks found, skip this RoV
          continue;
        }

        rs.push_back(std::move(r));
        }
    }
  }

  return rs;
}

void gen_ref_idxs(bd::VG &g, const std::vector<pgg::RoV> &all_rovs) {
  for (const pgg::RoV &r : all_rovs) {
    const std::vector<pgt::walk_t> &walks = r.get_walks();
    for (const pgt::walk_t& w : walks) {
      for (const auto &[v_id, _] : w) {
        g.gen_vtx_ref_idx(v_id);
      }
    }
  }
}

pgv::VcfRecIdx gen_vcf_rec_map(const std::vector<pvtr::Tree> &pvsts, bd::VG &g,
                               std::size_t thread_count) {

  std::vector<pgg::RoV> all_rovs = gen_rov(pvsts, g);
  gen_ref_idxs(g, all_rovs);
  std::vector<pga::Exp> exps = comp_expeditions(g, all_rovs, thread_count);
  pgv::VcfRecIdx rs = pgv::gen_vcf_records(g, exps);

  return rs;
  }
} // namespace povu::genomics
