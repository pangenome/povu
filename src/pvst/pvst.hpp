#include <vector>
#include <tuple>



#include "../graph/tree.hpp"


namespace pvst {
/**
 */
tree::Tree compute_pvst(std::vector<std::tuple<std::size_t , std::size_t, std::size_t>> const& v);
} // namespace pvst
