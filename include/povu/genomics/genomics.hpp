#ifndef POVU_GENOMICS_HPP
#define POVU_GENOMICS_HPP

// #include <cstddef>     // for size_t
#include <set>	       // for set
#include <string_view> // for string_view
#include <vector>      // for vector

// #include "indicators/dynamic_progress.hpp" // for DynamicProgress
// #include "indicators/progress_bar.hpp"	   // for ProgressBar
#include "povu/common/app.hpp"		 // for config
#include "povu/common/bounded_queue.hpp" // for pbq, bounded_queue
#include "povu/common/core.hpp"		 // for id_t, pt
// #include "povu/common/progress.hpp"	 // for progress
#include "povu/genomics/vcf.hpp"     // for VcfRecIdx
#include "povu/graph/bidirected.hpp" // for VG, bd
#include "povu/graph/pvst.hpp"	     // for Tree

namespace povu::genomics
{
inline constexpr std::string_view MODULE = "povu::genomics";
// using namespace povu::progress;
namespace pgv = povu::genomics::vcf;

void gen_vcf_rec_map(const std::vector<pvst::Tree> &pvsts, bd::VG &g,
		     const std::set<pt::id_t> &to_call_ref_ids,
		     pbq::bounded_queue<pgv::VcfRecIdx> &q,
		     const core::config &app_config);
} // namespace povu::genomics

#endif
