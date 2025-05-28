#include <iostream>
#include <optional>
#include <stack>
#include <unordered_set>
#include <vector>
#include <map>

#include "../../include/common/types.hpp"
#include "../graph/spanning_tree.hpp"
#include "../graph/tree.hpp"

namespace povu::hubbles {

#define MODULE "povu::hubbles"

namespace pc = povu::constants;
namespace pgt = povu::graph_types;
namespace pvtr = povu::tree;
namespace pt = povu::types;
namespace pst = povu::spanning_tree;

void find_hubbles(pst::Tree &st, pvtr::Tree<pgt::flubble> &ft);

} // namespace povu::hubbles

