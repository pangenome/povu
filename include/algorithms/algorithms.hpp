#ifndef PV_ALGOS_HPP
#define PV_ALGOS_HPP

#include "../graph/spanning_tree.hpp"

namespace povu::algorithms {
namespace pst = povu::spanning_tree;
using namespace povu::graph_types;
namespace pc = povu::constants;
namespace pt = povu::types;

/**
 * @brief
 */
void eulerian_cycle_equiv(pst::Tree &t);

void simple_cycle_equiv(pst::Tree &t);

} // namespace povu::algorithms
#endif
