#ifndef VST_HPP
#define VST_HPP

#include <cstddef>
#include <vector>

#include "../graph/spanning_tree.hpp"


#include "../common/typedefs.hpp"

namespace algorithms {
struct id_n_cls {
    std::size_t id;
    std::size_t cls;
};


void cycle_equiv(spanning_tree::Tree &t);

std::vector<common::typedefs::size_t_pair> find_seses(spanning_tree::Tree& t);

} // namespace algorithms

#endif
