#include "./genomics.hpp"
#include "vcf.hpp"

namespace povu::genomics {
pvt::VcfRecIdx
gen_vcf_rec_map(const std::vector<pgt::flubble> &canonical_flubbles,
                const bd::VG &g,
                const std::set<pt::id_t> &ref_ids,
                const core::config &app_config) {
  std::vector<pvt::RoV> rovs = pv::gen_rov(canonical_flubbles, g, app_config);

  std::cerr << "Extracted " << rovs.size() << " regions of variation\n";

  //std::cerr << "Found " << rovs.size() << "rovs\n";

  std::vector<pvt::RefWalks> rts = put::untangle_flb_rovs(g, rovs, ref_ids);



  std::cerr << "untangled " << rts.size() << " ref walks \n";
  pvt::VcfRecIdx rs = pgv::gen_vcf_records(g, rts);

  std::cerr << "generated vcf recs " << rs.get_recs().size() << " recs \n";

  //std::cerr << "left with " << rs.get_recs().size() << " recs \n";
  return rs;
}
} // namespace povu::genomics
