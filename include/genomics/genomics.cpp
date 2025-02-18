#include "./genomics.hpp"

namespace povu::genomics {
void to_vcf(const std::vector<pgt::flubble> &canonical_flubbles,
            const bd::VG &g, const core::config &app_config) {
  std::vector<pvt::RoV> rovs = pv::gen_rov(canonical_flubbles, g, app_config);
  put::untangle_flb_rovs(g, rovs);
}
} // namespace povu::genomics
