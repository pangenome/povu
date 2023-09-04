#include <vector>
#include <tuple>



#include "../graph/tree.hpp"


namespace pvst {
/**
 */
tree::Tree compute_pvst(std::vector<std::tuple<std::size_t , std::size_t, std::size_t>> const& v);
  
/**
 * @brief compute_pvst
 */
tree::Tree compute_pst(std::vector<std::size_t> classes);
} // namespace pvst
