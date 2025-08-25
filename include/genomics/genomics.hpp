#ifndef POVU_GENOMICS_HPP
#define POVU_GENOMICS_HPP

#include "../graph/bidirected.hpp"
#include "../common/compat.hpp"
#include "./allele.hpp"
#include "./graph.hpp"
#include "./untangle.hpp"
#include "./vcf.hpp"
#include "../common/utils.hpp"
#include "../common/log.hpp"

namespace povu::genomics {
inline constexpr std::string_view MODULE = "povu::genomics";

namespace put = povu::genomics::untangle;
namespace pgt = povu::types::graph;
namespace pgv = povu::genomics::vcf;
namespace pvst = povu::pvst;
namespace pga = povu::genomics::allele;
namespace pgg = povu::genomics::graph;

pgv::VcfRecIdx gen_vcf_rec_map(const std::vector<pvst::Tree> &pvsts, bd::VG &g,
                               std::size_t thread_count);
} // namespace povu::genomics

#endif
