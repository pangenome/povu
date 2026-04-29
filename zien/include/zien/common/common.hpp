#ifndef ZIEN_COMMON_HPP
#define ZIEN_COMMON_HPP

#include <set>

#include <quilt/types.hpp> // for qt

#include "mto/from_vcf.hpp" // for VCFile
// #include "povu/common/core.hpp"	     // for pt
#include "povu/graph/bidirected.hpp" // for bidirected

#include <povu/refs/refs.hpp> // for ref_format_e

namespace zien::common
{
constexpr oza::refs::ref_format_e PN = oza::refs::ref_format_e::PANSN;

std::set<qt::id_t> get_ref_ids(const bd::VG &g,
			       const mto::from_vcf::VCFile &vcf_file,
			       qt::u32 sample_idx, qt::u32 phase_idx);

char to_char(liteseq::strand s);
} // namespace zien::common

#endif // ZIEN_COMMON_HPP
