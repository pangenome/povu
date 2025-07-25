#include "./variants.hpp"

namespace povu::variants {


// void populate_walks(const bd::VG &g, std::vector<pvt::RoV> &rovs) {
//   std::string fn_name{std::format("[{}::{}]", MODULE, __func__)};

//   for (pt::idx_t i {}; i < rovs.size(); i++) {
//     bd::populate_walks(g, rovs[i], MAX_FLUBBLE_STEPS);
//   }

//   return;
// }

// std::vector<pvt::RoV> par_populate_walks(const bd::VG &g, std::vector<pvt::RoV> &rs) {
//   const std::string fn_name{std::format("[{}::{}]", MODULE, __func__)};

//   std::size_t cfl_count = rs.size();

//   uint8_t num_threads = std::thread::hardware_concurrency();
//   if (num_threads == 0) {
//     num_threads = 1;
//   }
//   // avoid creating more threads than flubbles
//   if (static_cast<std::size_t>(num_threads) > cfl_count) {
//     num_threads = cfl_count;
//   }

//   std::vector<std::thread> threads;
//   threads.reserve(num_threads);

//   //std::size_t fl_count = can_fls.size();
//   std::size_t chunk_size = cfl_count / num_threads;
//   std::size_t remainder = cfl_count % num_threads;
//   std::size_t start { 0 };

//   auto worker = [&](std::size_t start, std::size_t end) {
//     for (std::size_t i = start; i < end; ++i) {
//       bd::populate_walks(g, rs[i], MAX_FLUBBLE_STEPS);
//     }
//   };

//   for (uint8_t i = 0; i < num_threads; ++i) {
//     std::size_t end = start + chunk_size + (i < remainder ? 1 : 0);
//     threads.push_back(std::thread(worker, start, end));
//     start = end;
//   }

//   for (auto &t : threads) {
//     t.join();
//   }

//   return rs;
// }

 /* filter out RoVs whose walk count is less than 2 */
// inline void filter_invalid(std::vector<pvt::RoV> &r) {
//   const std::string fn_name{std::format("[{}::{}]", MODULE, __func__)};
//   int y = 0;
//   int z = 0;

//   for (auto x: r) {
//     if( x.walk_count() < 2) {
//       z += 1;
//     }
//     else {
//       y += 1;
//     }
//   }

//   std::cerr << "removed " << z << " keppt " << y << "\n";

//   return;
// };

/* initialize RoVs from flubbles */
inline std::vector<pvt::RoV> init_rovs(const bd::VG &g,
                                       const std::vector<pvtr::Tree> &pvsts) {
  std::vector<pvt::RoV> rs;
  rs.reserve(pvsts.size());

  for (const pvtr::Tree &pvst : pvsts) {
    for (pt::idx_t ft_v_idx{}; ft_v_idx < pvst.vtx_count(); ft_v_idx++) {
      if (pvst.is_leaf(ft_v_idx)) {
        const pvst::VertexBase *pvst_v_ptr = pvst.get_vertex_const_ptr(ft_v_idx);
        pvt::RoV r { pvst_v_ptr };
        std::vector<pgt::Walk> walks = graph_utils::get_walks(g, pvst_v_ptr);
        r.set_walks(std::move(walks));
        rs.push_back(std::move(r));
      }
    }
  }

  return rs;
}

std::vector<pvt::RoV> gen_rov(const std::vector<pvtr::Tree> &pvsts,
                              const bd::VG &g,
                              const core::config &app_config) {
  std::string fn_name {std::format("[{}::{}]", MODULE, __func__)};

  std::vector<pvt::RoV> rs = init_rovs(g, pvsts);
  //populate_walks(g, rs);
  //filter_invalid(rs);

  return rs;
}

} // namespace povu::variants
