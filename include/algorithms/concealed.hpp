#ifndef PV_CONCEALED_HPP
#define PV_CONCEALED_HPP

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
#include "../common/tree_utils.hpp"

namespace povu::concealed {

#define MODULE "povu::concealed"

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
namespace ptu = povu::tree_utils;

void find_concealed(const pst::Tree &st, pvtr::Tree &ft, const ptu::tree_meta &tm);
} // namespace povu::concealed
#endif // PV_CONCEALED_HPP
