#include "./genomics.hpp"
#include "indicators/dynamic_progress.hpp"
#include "vcf.hpp"
#include <cstddef>
#include <unordered_set>

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
pga::Exp exp_frm_rov(const bd::VG &g, const pgg::RoV &rov, povu::thread::thread_pool &pool) {

  // create an expedition object for the RoV
  const std::vector<pgt::walk_t> &walks = rov.get_walks();
  const pvst::VertexBase *pvst_v_ptr = rov.get_pvst_vtx();
  pga::Exp ref_walks{pvst_v_ptr};

  // compute the itineraries for each reference in the RoV
  std::map<pt::id_t, pga::Itn> &ref_map = ref_walks.get_ref_itns_mut();
  pga::comp_itineraries(g, walks, ref_map, pool);

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

std::vector<pga::Exp> comp_expeditions_work_steal(
    const bd::VG &g, const std::vector<pgg::RoV> &all_rovs, pt::idx_t start,
    pt::idx_t count, povu::thread::thread_pool &pool,
    std::size_t outer_concurrency,
    std::size_t reserve_for_inner
    ) {
  const std::size_t base = static_cast<std::size_t>(start);
  const std::size_t want = static_cast<std::size_t>(count);

  const std::size_t N = count;
  std::vector<pga::Exp> all_exp(N);

  // choose K workers; leave some for inner tasks
  const std::size_t pool_sz = pool.size();
  if (reserve_for_inner >= pool_sz)
    reserve_for_inner = 0;
  const std::size_t K = std::max<std::size_t>(
      1, std::min({outer_concurrency, pool_sz - reserve_for_inner, N}));

  std::atomic<std::size_t> next{0}; // counts *offsets* within [0, N)
  const std::size_t block = 64;

  povu::thread::task_group tg(pool);

  for (std::size_t k = 0; k < K; ++k) {
    tg.run([&, base] {
      for (;;) {
        // offset within the chunk
        const std::size_t off =
            next.fetch_add(block, std::memory_order_relaxed);
        if (off >= N)
          break;

        // global index range [g_begin, g_end)
        const std::size_t g_begin = base + off;
        const std::size_t g_end = std::min(g_begin + block, base + N);

        // fill results; local index = global - base
        for (std::size_t gi = g_begin; gi < g_end; ++gi) {
          const std::size_t li = gi - base;
          all_exp[li] = exp_frm_rov(g, all_rovs[gi], pool);
        }
      }
    });
  }

  tg.wait(); // rethrows first exception if any

  return all_exp;
}

/**
 * Check if a vertex in the pvst is a flubble leaf
 * A flubble leaf is a vertex that has no children that are also flubbles
 */
bool is_fl_leaf(const pvst::Tree &pvst, pt::idx_t pvst_v_idx) noexcept {
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
std::vector<pgg::RoV> gen_rov(const std::vector<pvst::Tree> &pvsts,
                              const bd::VG &g, const core::config &app_config) {

  // the set of RoVs to return
  std::vector<pgg::RoV> rs;
  rs.reserve(pvsts.size());

  // true when the vertex is a flubble leaf or a leaf in the pvst
  auto should_call = [&](const pvst::Tree &pvst, const pvst::VertexBase *pvst_v_ptr,
                     pt::idx_t pvst_v_idx) -> bool {
    if (std::optional<pvst::route_params_t> opt_rp = pvst_v_ptr->get_route_params()) {
      return is_fl_leaf(pvst, pvst_v_idx) || pvst.is_leaf(pvst_v_idx);
    }
    return false;
  };

  // reuse msg buffer for progress bar
  DynamicProgress<ProgressBar> bars;
  std::string prog_msg;
  prog_msg.reserve(128);

  for (pt::idx_t i{}; i < pvsts.size(); i++) { // for each pvst
    const pvst::Tree &pvst = pvsts[i];
    // loop through each tree
    const pt::idx_t total = pvst.vtx_count();

    // reset progress bar
    ProgressBar bar;
    set_progress_bar_common_opts(&bar);
    std::size_t bar_idx = bars.push_back(bar);
    set_progress_bar_common_opts(&bar, pvst.vtx_count());

    for (pt::idx_t pvst_v_idx{}; pvst_v_idx < pvst.vtx_count(); pvst_v_idx++) {

      if (app_config.show_progress()) { // update progress bar
        prog_msg.clear();
        fmt::format_to(std::back_inserter(prog_msg), "Generating RoVs for PVST {} ({}/{})", i+1, pvst_v_idx + 1, total);
        bars[bar_idx].set_option(option::PostfixText{prog_msg});
        bars[bar_idx].set_progress(pvst_v_idx + 1);
      }

      const pvst::VertexBase *pvst_v_ptr = pvst.get_vertex_const_ptr(pvst_v_idx);
      if (should_call(pvst, pvst_v_ptr, pvst_v_idx)) {
        pgg::RoV r{pvst_v_ptr};

        // get the set of walks for the RoV
        povu::genomics::graph::find_walks(g, r);

        if (r.get_walks().size() == 0) { // no walks found, skip this RoV
          continue;
        }

        rs.push_back(std::move(r));
      }

      if (app_config.show_progress()) {
        bars[bar_idx].mark_as_completed();
      }
    }
  }

  return rs;
}

struct ThreadSplit {
  std::size_t outer;
  std::size_t inner;
};

// Prefer inner >> outer; cap outer to a small number (default 2).
inline ThreadSplit
split_threads(std::size_t total, std::size_t outer_cap = 2,
              double outer_ratio = 0.25) // use ≤25% of pool for outer
{
  total = std::max<std::size_t>(1, total);
  if (total == 1)
    return {1, 0}; // outer runs; inner is serial

  // Suggested outer from ratio, but no more than outer_cap
  std::size_t suggested_outer = std::max<std::size_t>(
      1, static_cast<std::size_t>(std::ceil(total * outer_ratio)));
  std::size_t outer = std::min(outer_cap, suggested_outer);

  // Always leave at least one thread for inner
  if (outer >= total)
    outer = total - 1;

  std::size_t inner = total - outer; // guaranteed ≥ 1
  return {outer, inner};
}

void gen_vcf_rec_map(const std::vector<pvst::Tree> &pvsts, bd::VG &g,
                     pbq::bounded_queue<pgv::VcfRecIdx> &q,
                     DynamicProgress<ProgressBar> &prog, std::size_t bar_idx,
                     std::size_t thread_count, const core::config &app_config) {

  povu::thread::thread_pool pool(thread_count);
  std::vector<pgg::RoV> all_rovs = gen_rov(pvsts, g, app_config);

  const std::size_t CHUNK_SIZE = app_config.get_chunk_size();
  const std::size_t N = all_rovs.size();
  const std::size_t CHUNK_COUNT = (N + CHUNK_SIZE - 1) / CHUNK_SIZE;
  std::string prog_msg;
  prog_msg.reserve(128);

  std::vector<pga::Exp> exps;
  exps.reserve(CHUNK_SIZE);

  // decide on thread split
  auto [outer, inner] = split_threads(pool.size());

  try {

  for (pt::idx_t base{}; base < all_rovs.size(); base += CHUNK_SIZE) {
    const pt::idx_t end = std::min(base + CHUNK_SIZE, N);
    const pt::idx_t count = end - base;

    if (app_config.show_progress()) { // update progress bar
      prog_msg.clear();
      pt::idx_t chunk_num = (base / CHUNK_SIZE) + 1;
      fmt::format_to(std::back_inserter(prog_msg), "Processing Chunk ({}/{})", chunk_num, CHUNK_COUNT);
      prog[bar_idx].set_option(option::PostfixText{prog_msg});
      prog[bar_idx].set_progress(pu::comp_prog(chunk_num+1, CHUNK_COUNT));
    }

    exps = comp_expeditions_work_steal(g, all_rovs, base, count, pool, outer, inner);

    pgv::VcfRecIdx rs = pgv::gen_vcf_records(g, exps);
    if (!q.push(std::move(rs))) {
      break; // queue was closed early
    }
    exps.clear();
  }
  }
  catch (...) {
    q.close(); // make sure consumers wake up on errors
    throw;
  }

  q.close(); // we're done
}
} // namespace povu::genomics
