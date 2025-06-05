#ifndef PV_SLUBBLES_HPP
#define PV_SLUBBLES_HPP

#include <any>
#include <cassert>
#include <cstddef>
#include <cstdio>
#include <format>
#include <iostream>
#include <map>
#include <optional>
#include <queue>
#include <stack>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "../graph/spanning_tree.hpp"
#include "../graph/tree.hpp"
#include "../common/types.hpp"
#include "../common/utils.hpp"


namespace povu::slubbles {

#define MODULE "povu::slubbles"

using namespace povu::graph_types;
namespace pgt = povu::graph_types;

namespace pc = povu::constants;
namespace pvtr = povu::tree;
namespace pt = povu::types;
namespace pst = povu::spanning_tree;
namespace pvst= povu::types::pvst;
namespace pc = povu::constants;
namespace pt = povu::types;
namespace pvtr = povu::tree;

void find_hubbles(pst::Tree &st, pvtr::Tree<pvst::Vertex> &ft);

} // namespace povu::slubbles
#endif // PV_SLUBBLES_HPP
