#include "./genomics.hpp"
#include "vcf.hpp"

namespace povu::genomics {
pvt::VcfRecIdx gen_vcf_rec_map(const std::vector<pgt::flubble> &canonical_flubbles,
                               const bd::VG &g,
                               const core::config &app_config) {
  std::vector<pvt::RoV> rovs = pv::gen_rov(canonical_flubbles, g, app_config);
  std::vector<pvt::RefWalks> rts = put::untangle_flb_rovs(g, rovs);
  pvt::VcfRecIdx rs = pgv::gen_vcf_records(rts);

  return rs;
}
} // namespace povu::genomics
