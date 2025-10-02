#ifndef PV_MISC_HPP
#define PV_MISC_HPP

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

#include "../common/compat.hpp"
#include "../common/utils.hpp"
#include "../graph/pvst.hpp"
#include "../graph/spanning_tree.hpp"
#include "../graph/tree_utils.hpp"

namespace povu::midi
{
inline constexpr std::string_view MODULE = "povu::misc";

namespace pgt = povu::types::graph;
namespace pc = povu::constants;
namespace pst = povu::spanning_tree;
namespace pvst = povu::pvst;
namespace pc = povu::constants;
namespace ptu = povu::tree_utils;

void find_midi(const pst::Tree &st, pvst::Tree &pvst, const ptu::tree_meta &tm);
} // namespace povu::midi

#endif // PV_MISC_HPP
