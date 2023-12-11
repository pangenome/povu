#include <cstddef>
#include <vector>
#include <tuple>

#include "../core/core.hpp"
#include "../graph/tree.hpp"


namespace pvst {
/**
 */
//tree::Tree compute_pvst(std::vector<std::tuple<std::size_t , std::size_t, std::size_t>> const& v);

  
// Eq class and node id
struct eq_n_id_t {
  std::size_t v_id;
  std::size_t eq_class;
};

tree::Tree compute_pvst(std::vector<std::pair<std::size_t, std::size_t>> v,
						const core::config& app_config);

tree::Tree compute_pvst(std::vector<eq_n_id_t> v,
						const core::config& app_config);

} // namespace pvst
