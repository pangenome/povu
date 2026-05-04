#ifndef IT_GENOMICS_HPP
#define IT_GENOMICS_HPP

#include <set>	       // for set
#include <string_view> // for string_view
#include <vector>      // for vector

#include <oza/common/app.hpp>	    // for config
#include <oza/graph/bidirected.hpp> // for VG, bd
#include <oza/graph/pvst.hpp>	    // for Tree
#include <quilt/types.hpp>	    // for qt

#include "ita/genomics/vcf.hpp"	       // for VcfRecIdx
#include "ita/queue/bounded_queue.hpp" // for bounded_queue

namespace ita::genomics
{
inline constexpr std::string_view MODULE = "povu::genomics";

void gen_vcf_rec_map(const std::vector<pvst::Tree> &pvsts, bd::VG &g,
		     const std::set<qt::id_t> &to_call_ref_ids,
		     bq::bounded_queue<iv::VcfRecIdx> &q,
		     const core::config &app_config);
} // namespace ita::genomics

namespace ig = ita::genomics;

#endif // IT_GENOMICS_HPP
