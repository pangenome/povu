#ifndef POVU_GENOMICS_HPP
#define POVU_GENOMICS_HPP

#include <algorithm>
#include <chrono>
#include <unordered_set>
#include <vector>

#include "../../app/cli/app.hpp"
#include "../common/bounded_queue.hpp"
#include "../common/compat.hpp"
#include "../common/log.hpp"
#include "../common/progress.hpp"
#include "../common/thread.hpp"
#include "../common/utils.hpp"
#include "../graph/bidirected.hpp"
#include "./allele.hpp"
#include "./graph.hpp"
#include "./untangle.hpp"
#include "./vcf.hpp"
namespace povu::genomics {
inline constexpr std::string_view MODULE = "povu::genomics";

using namespace povu::progress;

namespace put = povu::genomics::untangle;
namespace pgt = povu::types::graph;
namespace pgv = povu::genomics::vcf;
namespace pvst = povu::pvst;
namespace pga = povu::genomics::allele;
namespace pgg = povu::genomics::graph;

void gen_vcf_rec_map(const std::vector<pvst::Tree> &pvsts, bd::VG &g,
                     const std::set<pt::id_t> &to_call_ref_ids,
                     pbq::bounded_queue<pgv::VcfRecIdx> &q,
                     DynamicProgress<ProgressBar> &prog, std::size_t prog_idx,
                     const core::config &app_config);
} // namespace povu::genomics

#endif
