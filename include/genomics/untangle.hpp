#ifndef POVU_UNTANGLE_HPP
#define POVU_UNTANGLE_HPP

#include <cstddef>
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

#include "../../app/cli/app.hpp"
#include "../common/types/genomics.hpp"
#include "../common/compat.hpp"
#include "../graph/bidirected.hpp"
#include "../align/align.hpp"
#include "../common/types/pvst.hpp"


namespace povu::genomics::untangle {
namespace pgt = povu::types::graph;
namespace pvt = povu::types::genomics;
namespace bd = povu::bidirected;
namespace pc = povu::constants;
namespace pa = povu::align;
namespace pvst = povu::types::pvst;

inline constexpr std::string_view MODULE = "povu::untangle";

void untangle_ref_walks(pvt::Exp &rt);
} // namespace povu::untangle

#endif // POVU_UNTANGLE_HPP
