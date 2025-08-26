#ifndef POVU_UNTANGLE_HPP
#define POVU_UNTANGLE_HPP

#include <cstddef>
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

#include "../../app/cli/app.hpp"
#include "../common/compat.hpp"
#include "../graph/bidirected.hpp"
#include "../align/align.hpp"
#include "../graph/pvst.hpp"
#include "./allele.hpp"

namespace povu::genomics::untangle {
namespace pgt = povu::types::graph;
namespace bd = povu::bidirected;
namespace pc = povu::constants;
namespace pa = povu::align;
namespace pvst = povu::pvst;
namespace pga = povu::genomics::allele;


inline constexpr std::string_view MODULE = "povu::genomics::untangle";

void untangle_ref_walks(pga::Exp &rt);
} // namespace povu::untangle

#endif // POVU_UNTANGLE_HPP
