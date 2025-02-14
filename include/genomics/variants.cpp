#include "./variants.hpp"

namespace povu::variants {
#define MODULE "povu::variants"

std::vector<std::vector<pgt::walk>>
find_flubble_paths(const std::vector<pgt::flubble> &canonical_flubbles,
                  const bd::VG &bd_vg) {
  std::string fn_name{std::format("[{}::{}]", MODULE, __func__)};

  std::vector<std::vector<pgt::walk>> all_paths;
  const pt::idx_t MAX_STEPS { 20 };

  // TODO: perform a parallel for loop
  for (std::size_t i{}; i < canonical_flubbles.size(); ++i) {

    const auto &[entry, exit] = canonical_flubbles[i];
    // a walk is a vector of ID and orientation
    std::vector<pgt::walk> paths = bd::get_walks(bd_vg, entry, exit, MAX_STEPS);
    if (paths.size() < 2) {
      std::cerr << std::format("{} WARN: Bubble {} {} has {} paths\n", fn_name,
                               entry.as_str(), exit.as_str(), paths.size());
    }

    all_paths.push_back(paths);
  }
  return all_paths;
}

void call_variants(const std::vector<pgt::flubble>& canonical_flubbles,
                   const bd::VG& bd_vg,
                   const core::config& app_config) {
  std::string fn_name{std::format("[{}::{}]", MODULE, __func__)};
  
  std::vector<std::vector<pgt::walk>> flubble_paths =
    find_flubble_paths(canonical_flubbles, bd_vg);

  
}


} // namespace povu::variants
