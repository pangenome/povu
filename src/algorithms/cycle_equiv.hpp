#ifndef VST_HPP
#define VST_HPP

#include <cstddef>
#include <vector>

#include "../graph/spanning_tree.hpp"


#include "../common/common.hpp"

namespace algorithms {
struct id_n_cls {
    std::size_t id;
    std::size_t cls;
};


void cycle_equiv(spanning_tree::Tree &t);

// TODO: rename to something related to reflect return type
std::vector<graph_types::canonical_sese> find_seses(spanning_tree::Tree& t);

} // namespace algorithms

#endif
