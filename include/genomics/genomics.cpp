#include "./genomics.hpp"
#include <vector>


namespace povu::genomics {

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
pvt::Exp comp_expeditions(const bd::VG &g, const pvt::RoV &rov) {
#ifdef DEBUG
  std::chrono::duration<double> timeRefRead;
  auto t0 = pt::Time::now();
#endif

  const std::vector<pvt::walk_t> &walks = rov.get_walks();
  const pvst::VertexBase *pvst_v_ptr = rov.get_pvst_vtx();

  // create an expedition object for the RoV
  pvt::Exp ref_walks {pvst_v_ptr};

  // a walk is a single traversal bounded by start to the end of an RoV
  for (pt::idx_t w_idx{}; w_idx < rov.walk_count(); w_idx++) {
    const pvt::walk_t &w = walks[w_idx];
    povu::genomics::allele::comp_itineraries(g, w, w_idx, ref_walks);
  }


#ifdef DEBUG
  if (true) {
    timeRefRead = pt::Time::now() - t0;
    INFO("ROV {}", rov.as_str());
    INFO("ref count in exp {} walk count {}", ref_walks.ref_count(), rov.walk_count());
    INFO("Time spent finding walks in itns {:.2f} sec", timeRefRead.count());
    t0 = pt::Time::now();
  }
#endif

  for (pt::idx_t ref_id : ref_walks.get_ref_ids()) {
    remove_prefix_walks(ref_walks.get_itn_mut(ref_id));
  }

#ifdef DEBUG
  if (true) {
    timeRefRead = pt::Time::now() - t0;
    INFO("Time spent removing prefixes {:.2f} sec", timeRefRead.count());
    t0 = pt::Time::now();
  }
#endif

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

#ifdef DEBUG
  if (true) {
    timeRefRead = pt::Time::now() - t0;
    INFO("Time spent bla {:.2f} sec", timeRefRead.count());
    t0 = pt::Time::now();
  }
#endif

  if (ref_walks.is_tangled()) {
    put::untangle_ref_walks(ref_walks);
  }

#ifdef DEBUG
  if (true) {
    timeRefRead = pt::Time::now() - t0;
    INFO("Time spent untangling {:.2f} sec", timeRefRead.count());
    t0 = pt::Time::now();
  }
#endif

  return ref_walks;
}

// std::vector<pvt::Exp> do_ref_walks(const bd::VG &g, const std::vector<pvt::RoV> &all_rovs) {
//   std::vector<pvt::Exp> all_exp(all_rovs.size());
//   //all_exp.reserve(all_rovs.size());

//   std::size_t thread_count = 16; // default
//   auto [num_threads, chunk_size] = pu::compute_thread_allocation(thread_count, all_rovs.size());


//   INFO("Using {} threads to compute reference walks for {} RoVs", num_threads, all_rovs.size());

//   std::vector<std::thread> threads(num_threads);
//   std::size_t start, end;
//   for (unsigned int thread_idx{}; thread_idx < num_threads; ++thread_idx) {
//     start = thread_idx * chunk_size;
//     end = (thread_idx == num_threads - 1) ? all_rovs.size() : (thread_idx + 1) * chunk_size;

//     INFO("Thread {} processing RoVs from {} to {}", thread_idx, start, end);

//     threads[thread_idx] = std::thread([&, start, end]() {
//           for (std::size_t i{start}; i < end; i++) {
//             const pvt::RoV &r = all_rovs[i];
//             pvt::Exp rov_rws = comp_expeditions(g, r);
//             all_exp[i] = std::move(rov_rws);
//             //all_ref_walks.push_back(std::move(rov_rws));
//           }
//         });
//   }

//   // Wait for all threads to finish
//   for (auto &thread : threads) {
//     thread.join();
//   }

//   return all_exp;
// }

// Assuming bd::VG is thread-safe for concurrent read-only access,
// and gen_rov_ref_walks(g, r) returns pvt::Exp.
std::vector<pvt::Exp> do_ref_walks_pool(const bd::VG &g,
                                        const std::vector<pvt::RoV> &all_rovs,
                                        pu::ThreadPool &pool) {
  const std::size_t N = all_rovs.size();
  std::vector<pvt::Exp> out(N);
  if (N == 0) {
    return out;
  }

  parallel_for(pool, N, [&](std::size_t i) {
    INFO("Processing RoV {}/{}", i + 1, N);
    out[i] = comp_expeditions(g, all_rovs[i]); // disjoint writes → no locks
  });

  return out;
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
std::vector<pvt::RoV> gen_rov(const std::vector<pvtr::Tree> &pvsts, const bd::VG &g) {
  // the set of RoVs to return
  std::vector<pvt::RoV> rs;
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
        pvt::RoV r{pvst_v_ptr};

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

void gen_ref_idxs(bd::VG &g, const std::vector<pvt::RoV> &all_rovs) {
  for (const pvt::RoV &r : all_rovs) {
    const std::vector<pvt::walk_t> &walks = r.get_walks();
    for (const pvt::walk_t& w : walks) {
      for (const auto &[v_id, _] : w) {
        g.gen_vtx_ref_idx(v_id);
      }
    }
  }
}
pvt::VcfRecIdx gen_vcf_rec_map(const std::vector<pvtr::Tree> &pvsts, bd::VG &g) {
#ifdef DEBUG
  std::chrono::duration<double> timeRefRead;
  auto t0 = pt::Time::now();
#endif
  std::vector<pvt::RoV> all_rovs = gen_rov(pvsts, g);

#ifdef DEBUG
  if (true) {
    timeRefRead = pt::Time::now() - t0;
    INFO("Time spent finding walks in RoVs {:.2f} sec", timeRefRead.count());
    t0 = pt::Time::now();
  }
#endif
  gen_ref_idxs(g, all_rovs);



#ifdef DEBUG
  if (true) {
    timeRefRead = pt::Time::now() - t0;
    INFO("Time spent gen idxs {:.2f} sec", timeRefRead.count());
    t0 = pt::Time::now();
  }
#endif

  const std::size_t thread_count = 16; // default
  pu::ThreadPool pool(thread_count); // default thread count
  std::vector<pvt::Exp> all_ref_walks_ = do_ref_walks_pool(g, all_rovs, pool);
  //std::vector<pvt::Exp> all_ref_walks;
  // for (const pvt::RoV &r : all_rovs) {
  //   pvt::Exp rov_rws = gen_rov_ref_walks(g, r);
  //   all_ref_walks.push_back(std::move(rov_rws));
  // }
  //
  //   gen_rov_ref_walks(g, r, all_ref_walks);
  // }

#ifdef DEBUG
      if (true) {
    timeRefRead = pt::Time::now() - t0;
    INFO("Time spent computing walks {:.2f} sec", timeRefRead.count());
    t0 = pt::Time::now();
  }
#endif

  pvt::VcfRecIdx rs = pgv::gen_vcf_records(g, all_ref_walks_);

#ifdef DEBUG
  if (true) {
    timeRefRead = pt::Time::now() - t0;
    INFO("Time gen VCFs {:.2f} sec", timeRefRead.count());
    t0 = pt::Time::now();
  }
#endif

  return rs;
  }
} // namespace povu::genomics
