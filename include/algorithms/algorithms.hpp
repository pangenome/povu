#ifndef PV_ALGOS_HPP
#define PV_ALGOS_HPP

#include "../../include/common/types.hpp"
#include "../../src/cli/app.hpp"
#include "../graph/spanning_tree.hpp"
#include "../graph/bidirected.hpp"
#include "../graph/tree.hpp"


namespace povu::algorithms {
#define MODULE "povu::algorithms"

namespace pst = povu::spanning_tree;

using namespace povu::graph_types;
namespace pgt = povu::graph_types;
namespace pvtr = povu::tree;

namespace pc = povu::constants;
namespace pt = povu::types;
namespace bd = povu::bidirected;

/**
 * @brief
 */
void eulerian_cycle_equiv(pst::Tree &t);

void simple_cycle_equiv(pst::Tree &t, const core::config &app_config);

} // namespace povu::algorithms
#endif
