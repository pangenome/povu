#ifndef PV_SMOTHERED_HPP
#define PV_SMOTHERED_HPP

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

"#include "povu/common/compat.hpp"
"#include "povu/common/utils.hpp"
"#include "povu/graph/pvst.hpp"
"#include "povu/graph/spanning_tree.hpp"
"#include "povu/graph/tree_utils.hpp"

namespace povu::smothered
{
inline constexpr std::string_view MODULE = "povu::smothered";

namespace pgt = povu::types::graph;
namespace pc = povu::constants;
namespace pst = povu::spanning_tree;
namespace pvst = povu::pvst;
namespace pc = povu::constants;
namespace ptu = povu::tree_utils;

void find_smothered(const pst::Tree &st, pvst::Tree &ft,
		    const ptu::tree_meta &tm);
} // namespace povu::smothered

#endif // PV_SMOTHERED_HPP
