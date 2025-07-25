#include "./genomics.hpp"

namespace povu::genomics {
pvt::VcfRecIdx
gen_vcf_rec_map(const std::vector<pvtr::Tree> &pvsts, const bd::VG &g, const core::config &app_config) {
  std::vector<pvt::RoV> rovs = pv::gen_rov(pvsts, g, app_config);
  std::vector<pvt::RefWalks> rts = put::untangle_flb_rovs(g, rovs);
  pvt::VcfRecIdx rs = pgv::gen_vcf_records(g, rts);
  return rs;
}
} // namespace povu::genomics
