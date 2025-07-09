#ifndef PV_SMOTHERED_HPP
#define PV_SMOTHERED_HPP

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

#include "../common/tree_utils.hpp"
#include "../common/types/types.hpp"
#include "../common/utils.hpp"
#include "../graph/spanning_tree.hpp"
#include "../graph/tree.hpp"

namespace povu::smothered {
  //#define MODULE "povu::smothered"
inline constexpr std::string_view MODULE = "povu::smothered";

namespace pgt = povu::types::graph;

namespace pc = povu::constants;
namespace pvtr = povu::tree;
namespace pt = povu::types;
namespace pst = povu::spanning_tree;
namespace pvst = povu::types::pvst;
namespace pc = povu::constants;
namespace pt = povu::types;
namespace pvtr = povu::tree;
namespace ptu = povu::tree_utils;

void find_smothered(const pst::Tree &st, pvtr::Tree &ft, const ptu::tree_meta &tm);
} // namespace povu::smothered

#endif // PV_SMOTHERED_HPP
