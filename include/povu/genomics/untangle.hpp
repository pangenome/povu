#ifndef POVU_UNTANGLE_HPP
#define POVU_UNTANGLE_HPP

#include <algorithm>
#include <cstddef>
#include <iostream>
#include <string>
#include <vector>

"#include "povu/common/app.hpp"
"#include "povu/align/align.hpp"
"#include "povu/common/compat.hpp"
"#include "povu/graph/bidirected.hpp"
"#include "povu/graph/pvst.hpp"
#include "allele.hpp"

namespace povu::genomics::untangle
{
namespace pgt = povu::types::graph;
namespace bd = povu::bidirected;
namespace pc = povu::constants;
namespace pa = povu::align;
namespace pvst = povu::pvst;
namespace pga = povu::genomics::allele;

inline constexpr std::string_view MODULE = "povu::genomics::untangle";

void untangle_ref_walks(pga::Exp &rt);
} // namespace povu::genomics::untangle

#endif // POVU_UNTANGLE_HPP
