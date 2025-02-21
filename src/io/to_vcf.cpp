#include "./to_vcf.hpp"

namespace povu::io::to_vcf {

inline void write_header(const std::string &chrom , std::ostream &os) {

  os << "##fileformat=VCFv4.2\n";
  os << "##fileDate=" << pu::today() << std::endl;
  os << "##source=povu\n";
  os << "##reference=" << chrom << "\n";
  os << "##INFO=<ID=AT,Number=R,Type=String,Description=\"Allele "
              "Traversal as path in graph\">\n";
  os << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\n";
}


void write_vcf(const bd::VG &g, const std::string &chrom,
               std::vector<pvt::VcfRec> recs, std::ostream &os) {

  auto get_label = [&](const pvt::Step &s) -> std::string {
    if (s.get_o() == pgt::or_e::forward) {
      return g.get_vertex_by_id(s.get_v_id()).get_label();
    } else {
      return g.get_vertex_by_id(s.get_v_id()).get_rc_label();
    }
  };

  auto at_to_str = [&](const std::vector<pvt::AT> &ats) -> std::string {
    std::string at_str = "";
    for (std::size_t step_idx {}; step_idx< ats.size(); ++step_idx) {
      const auto &at = ats[step_idx];
      for (auto &s: at.get_steps()) {
        at_str += get_label(s);
      }

      if (step_idx < ats.size() - 1) {
        at_str += ",";
      }
    }
    return at_str;
  };

  const std::string qual = "60";

  write_header(chrom, os);

  for (const auto &r: recs) {
    os << chrom
       << "\t" << r.get_pos()
       << "\t" << "."
       << "\t" << at_to_str( std::vector<pvt::AT>{r.get_ref_at()})
       << "\t" << at_to_str( r.get_alt_ats())
       << "\t" << qual
       << "\t" << "."
       << "\t" << "."
       << "\t" << "GT"
       << "\n";
  }

  return;
}

void write_vcfs(const pvt::VcfRecIdx &vcf_recs, const bd::VG &g, const core::config &app_config) {

  std::string out_dir = std::string(app_config.get_output_dir());


  for (const auto &[ref_id, recs] : vcf_recs.get_recs()) {

    std::string ref_name = g.get_ref_name(ref_id);
    std::string vcf_fp = std::format("{}/{}.vcf", out_dir, ref_name);
    std::ofstream os(vcf_fp);
    write_vcf(g, ref_name, recs, os);
  }

  return;
}

} // namespacepovu::io::vcf
