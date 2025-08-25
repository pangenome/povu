#ifndef PV_TINY_HPP
#define PV_TINY_HPP

#include <any>
#include <cassert>
#include <cstddef>
#include <cstdio>
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

#include "../common/tree_utils.hpp"
#include "../common/compat.hpp"
#include "../common/utils.hpp"
#include "../graph/spanning_tree.hpp"
#include "../graph/pvst.hpp"


namespace povu::tiny {
inline constexpr std::string_view MODULE = "povu::tiny";

namespace pc = povu::constants;
namespace ptu = povu::tree_utils;
namespace pst = povu::spanning_tree;
namespace pvst = povu::pvst;
namespace pgt = povu::types::graph;


void find_tiny(const pst::Tree &st, pvst::Tree &pvst, const ptu::tree_meta &tm);

} // namespace povu::parallel
#endif // PV_TINY_HPP
