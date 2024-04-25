#ifndef VST_HPP
#define VST_HPP

#include "../graph/spanning_tree.hpp"
#include "../common/constants.hpp"

namespace algorithms {
struct id_n_cls {
    std::size_t id;
    std::size_t cls;
};


void cycle_equiv(spanning_tree::Tree &t);

void find_seses(spanning_tree::Tree& t);

} // namespace algorithms

#endif
