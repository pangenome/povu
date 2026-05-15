#ifndef ZIEN_COMMON_HPP
#define ZIEN_COMMON_HPP

#include <set>
#include <string>

#include "mto/from_vcf.hpp"	     // for VCFile
#include "povu/common/core.hpp"	     // for pt
#include "povu/graph/bidirected.hpp" // for bidirected

namespace zien::common
{
constexpr povu::refs::ref_format_e PN = povu::refs::ref_format_e::PANSN;

std::set<pt::id_t> get_ref_ids(const bd::VG &g,
			       const mto::from_vcf::VCFile &vcf_file,
			       pt::u32 sample_idx, pt::u32 phase_idx);

char to_char(liteseq::strand s);
} // namespace zien::common

#endif // ZIEN_COMMON_HPP
