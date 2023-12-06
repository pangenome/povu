#include <vector>
#include <tuple>

#include "../core/core.hpp"
#include "../graph/tree.hpp"


namespace pvst {
/**
 */
//tree::Tree compute_pvst(std::vector<std::tuple<std::size_t , std::size_t, std::size_t>> const& v);

  tree::Tree compute_pvst(std::vector<std::pair<std::size_t, std::size_t>> v, const core::config& app_config);
} // namespace pvst
