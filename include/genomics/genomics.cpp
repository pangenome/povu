#include "./genomics.hpp"
#include "vcf.hpp"

namespace povu::genomics {
void to_vcf(const std::vector<pgt::flubble> &canonical_flubbles,
            const bd::VG &g,
            const core::config &app_config) {
  std::vector<pvt::RoV> rovs = pv::gen_rov(canonical_flubbles, g, app_config);
  std::vector<pvt::RefWalks> rts = put::untangle_flb_rovs(g, rovs);
  pgv::gen_vcfs(g, rts);
}
} // namespace povu::genomics
