#ifndef IT_GENOMICS_HPP
#define IT_GENOMICS_HPP

#include <set>	       // for set
#include <string_view> // for string_view
#include <vector>      // for vector

#include "ita/genomics/vcf.hpp" // for VcfRecIdx

#include "povu/common/app.hpp"		 // for config
#include "povu/common/bounded_queue.hpp" // for pbq, bounded_queue
#include "povu/common/core.hpp"		 // for id_t, pt
#include "povu/graph/bidirected.hpp"	 // for VG, bd
#include "povu/graph/pvst.hpp"		 // for Tree

namespace ita::genomics
{
inline constexpr std::string_view MODULE = "povu::genomics";

void gen_vcf_rec_map(const std::vector<pvst::Tree> &pvsts, bd::VG &g,
		     const std::set<pt::id_t> &to_call_ref_ids,
		     pbq::bounded_queue<iv::VcfRecIdx> &q,
		     const core::config &app_config);
} // namespace ita::genomics

namespace ig = ita::genomics;

#endif // IT_GENOMICS_HPP
