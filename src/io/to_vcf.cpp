#include "./to_vcf.hpp"
#include <string>
#include <utility>
#include <vector>

namespace povu::io::to_vcf {

inline void write_header(const std::string &chrom , std::ostream &os) {

  os << "##fileformat=VCFv4.2\n";
  os << "##fileDate=" << pu::today() << std::endl;
  os << "##source=povu\n";
  os << "##reference=" << chrom << "\n";
  os << "##INFO=<ID=AT,Number=R,Type=String,Description=\"Allele Traversal as path in graph\">\n";
  os << "##INFO=<ID=VARTYPE,Number=R,Type=String,Description=\"type of variation ins for insertion, del for deletion and sub for substitution\">\n";
  os << "##INFO=<ID=TANGLED,Number=R,Type=String,Description=\"if the variant is in a tangled part of the graph T if tangled and F if not tangled\n";
  os << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\n";
}

void write_vcf_rec(const bd::VG &g, const pvt::VcfRec &r, const std::string &chrom, std::ostream &os) {
  std::string fn_name{std::format("[{}::{}]", MODULE, __func__)};

  const std::string qual = "60";
  pvt::var_type_e var_typ = r.get_var_type();

  auto get_label = [&](const pvt::AS &s) -> std::string {
    if (s.get_o() == pgt::or_e::forward) {
      return g.get_vertex_by_id(s.get_v_id()).get_label();
    }
    else {
      return g.get_vertex_by_id(s.get_v_id()).get_rc_label();
    }
  };

  auto at_to_dna_str = [&](const pvt::AW &at) -> std::string {
    std::string at_str = "";

    // 1) Anchor base for deletions & insertions
    switch (var_typ) {
    case pvt::var_type_e::del:
    case pvt::var_type_e::ins: {
      // grab the first step’s label, take its last character
      auto lbl = get_label(at.get_step(0));
      at_str.push_back(lbl.back());
      break;
    }
    default:
      break;
    }

    // 2) Middle steps (for all types) does nothing for deletions
    //    steps [1 .. step_count()-2], inclusive

    for (pt::idx_t step_idx{1}; step_idx < at.step_count() - 1; ++step_idx) {
      const pvt::AS &s = at.get_step(step_idx);
      at_str += get_label(s);
    }

    // 3) If nothing got added, use “.”
    return at_str.empty() ? std::string{"."} : at_str;

  };

  // TODO: [A] look into the prefix here
  auto ats_to_dna_str = [&](std::string &&prefix, const std::vector<pvt::AW> &allele_walks) -> std::string {
    std::string allele_walks_str = "";
    for (std::size_t w_idx {}; w_idx < allele_walks.size(); ++w_idx) {
      const auto &aw = allele_walks[w_idx];
      allele_walks_str += prefix;
      allele_walks_str += at_to_dna_str(aw);

      if (w_idx < allele_walks.size() - 1) {
        allele_walks_str += ",";
      }
    }
    return allele_walks_str;
  };

  auto at_as_str = [](const pvt::AW &at) -> std::string {
    std::string str;

    for (auto &s : at.get_steps()) {
      str += std::format("{}{}", s.get_o() == pgt::or_e::forward ? ">" : "<", s.get_v_id());
    }

    return str;
  };

  auto ats_as_str = [&](const std::vector<pvt::AW> & ats) -> std::string {
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
    s += at_as_str(r.get_ref_at());
    s += ",";
    s += ats_as_str(r.get_alt_ats());

    return s;
  };

  pt::idx_t pos = var_typ == pvt::var_type_e::del || var_typ == pvt::var_type_e::ins ? r.get_pos() - 1 : r.get_pos();
  std::string ref_dna = at_to_dna_str(r.get_ref_at());

  // When you delete one or more bases, you must include at least one “anchor”
  // base so that neither REF nor ALT becomes empty
  // the prefix is the anchor base for deletions and insertions
  std::string anchor_base = var_typ == pvt::var_type_e::del ? ref_dna : "";
  std::string alt_dna = ats_to_dna_str(std::move(anchor_base), r.get_alt_ats());
  std::string info_field = std::format("AT={},VARTYPE={},TANGLED={}", fmt_field(r), pvt::to_string_view(var_typ), r.is_tangled() ? "T" : "F");

  os << chrom
     << "\t" << pos
     << "\t" << r.get_id()
     << "\t" << ref_dna
     << "\t" << alt_dna
     << "\t" << qual
     << "\t" << "."
     << "\t" << info_field
     << "\t" << "GT"
     << "\n";
}

void write_vcf(const bd::VG &g, const std::string &chrom,
               std::vector<pvt::VcfRec> recs, std::ostream &os) {
  std::string fn_name{std::format("[{}::{}]", MODULE, __func__)};

  write_header(chrom, os);
  for (const pvt::VcfRec &r: recs) {
    write_vcf_rec(g, r, chrom, os);
  }

  return;
}



void write_vcfs(const pvt::VcfRecIdx &vcf_recs, const bd::VG &g, const core::config &app_config) {
  std::string fn_name{std::format("[{}::{}]", MODULE, __func__)};

  std::string out_dir = std::string(app_config.get_output_dir());

  for (const auto &[ref_id, recs] : vcf_recs.get_recs()) {
    std::string ref_name = g.get_ref_name(ref_id);
    std::string vcf_fp = std::format("{}/{}.vcf", out_dir, ref_name);
    std::ofstream os(vcf_fp);
    write_vcf(g, ref_name, recs, os);
    std::cerr << "wrote " << vcf_fp << "\n";
  }

  return;
}

} // namespacepovu::io::vcf
