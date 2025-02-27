#include "./to_vcf.hpp"
#include <vector>

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

  auto at_to_dna_str = [&](const std::vector<pvt::AT> &ats) -> std::string {
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


  auto at_as_str = [](const pvt::AT &at) -> std::string {
    std::string str;

    for (auto &s: at.get_steps()) {
      str += std::format("{}{}",
                         s.get_o() == pgt::or_e::forward ? ">" : "<",
                         s.get_v_id() );
    }
    return str;
  };

  auto ats_as_str = [&](const std::vector<pvt::AT> & ats) -> std::string {
    std::string str;
    // add a comma between ats
    for (std::size_t i {}; i < ats.size(); ++i) {
      str += at_as_str(ats[i]);
      if (i < ats.size() - 1) {
        str += ",";
      }
    }
    return str;
  };

  auto fmt_field = [&](const pvt::VcfRec &r) -> std::string {
    std::string s;
    s += ats_as_str(std::vector<pvt::AT>{r.get_ref_at()});
    s += ",";
    s += ats_as_str(r.get_alt_ats());

    return s;
  };

  auto ref_dna_str = [&](const pvt::VcfRec &r) -> std::string {
    return at_to_dna_str(std::vector<pvt::AT>{r.get_ref_at()});
  };

  auto alt_dna_strs = [&](const pvt::VcfRec &r) {
    return at_to_dna_str(r.get_alt_ats());
  };

  write_header(chrom, os);

  for (const pvt::VcfRec &r: recs) {
    os << chrom
       << "\t" << r.get_pos()
       << "\t" << r.get_id()
       << "\t" << ref_dna_str(r)
       << "\t" << alt_dna_strs(r)
       << "\t" << qual
       << "\t" << "."
       << "\t" << "AT=" << fmt_field(r)
       << "\t" << "GT"
       << "\n";
  }

  return;
}

void write_vcfs(const pvt::VcfRecIdx &vcf_recs, const bd::VG &g, const core::config &app_config) {

  std::string out_dir = std::string(app_config.get_output_dir());

  std::cerr << "called write_vcfs\n";

  for (const auto &[ref_id, recs] : vcf_recs.get_recs()) {

    std::string ref_name = g.get_ref_name(ref_id);
    std::string vcf_fp = std::format("{}/{}.vcf", out_dir, ref_name);

    std::cerr << "writing " << vcf_fp << "\n";

    std::ofstream os(vcf_fp);
    write_vcf(g, ref_name, recs, os);
  }

  return;
}

} // namespacepovu::io::vcf
