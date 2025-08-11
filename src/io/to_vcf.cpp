#include "./to_vcf.hpp"
#include <algorithm>
#include <cmath>
#include <set>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace povu::io::to_vcf {

inline void write_header(const std::vector<std::pair<std::string, pt::idx_t>> &contigs, std::ostream &os) {
  os << "##fileformat=VCFv4.2\n";
  os << "##fileDate=" << pu::today() << std::endl;
  os << "##source=povu\n";
  os << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
  os << "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Total number of alternate alleles in called genotypes\">\n";
  os << "##INFO=<ID=AT,Number=R,Type=String,Description=\"Allele traversal path through the graph\">\n";
  os << "##INFO=<ID=AN,Number=1,Type=String,Description=\"Total number of alleles in called genotypes\">\n";
  os << "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele frequency in the population\">\n";
  os << "##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of samples with data\">\n";
  os << "##INFO=<ID=VARTYPE,Number=1,Type=String,Description=\"Type of variation: INS (insertion), DEL (deletion), SUB (substitution)\">\n";
  os << "##INFO=<ID=TANGLED,Number=1,Type=String,Description=\"Variant lies in a tangled region of the graph: T or F\">\n";
  os << "##INFO=<ID=LV,Number=1,Type=Integer,Description=\"Level in the snarl tree (0=top level)\">\n";
  os << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
  
  for (const auto &[chrom, len] : contigs) {
    os << std::format("##contig=<ID={},length={}>\n", chrom, len);
  }

  return;
}

inline void write_single_header(const std::string &chrom, pt::idx_t len, std::ostream &os) {
  write_header({{chrom, len}}, os);
}


inline void write_combined_header(const std::vector<std::pair<std::string, pt::idx_t>> &contigs, std::ostream &os) {
  write_header(contigs, os);
}

inline void write_col_header(pvt::genotype_data_t gtd, std::ostream &os) {
  os << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t";
  for (std::size_t i = 0; i < gtd.genotype_cols.size(); ++i) {
    os << gtd.genotype_cols[i];
    if (i < gtd.genotype_cols.size() - 1) {
      os << "\t";
    }
  }
  os << "\n";
}


// ns and the genotype columns are generated from the genotype data
std::pair<pt::idx_t, std::string> gen_genotype_cols(const bd::VG &g, pt::id_t ref_count,
                                                    const pvt::genotype_data_t &gtd,
                                                    pt::id_t ref_hap_id,
                                                    const pvt::AW &ref_aw,
                                                    const std::vector<pvt::AW> &alt_aws) {
  std::string fn_name{std::format("[{}::{}]", MODULE, __func__)};

  typedef std::vector<std::string> gt_col;

  // initialise all with empty vectors
  std::vector<gt_col> cols(gtd.genotype_cols.size());

  pt::idx_t ns {}; // no of samples with data

  // matching ref
  for (pt::id_t ref_allele_hap_id : ref_aw.get_ref_ids()) {

    pt::idx_t hap_col_idx = gtd.ref_id_to_col_idx.at(ref_allele_hap_id);
    gt_col &col = cols[hap_col_idx];

    const std::set<pt::id_t> &sample_refs = g.get_shared_samples(ref_allele_hap_id);

    if (col.empty()) {
      ns++; // increment the number of samples with data
      col.resize(sample_refs.size(), ".");
    }

    pt::idx_t col_col_idx{};
    for ( auto it = sample_refs.begin(); it != sample_refs.end(); it++, col_col_idx++) {
      if (*it == ref_allele_hap_id) {
        col[col_col_idx] = "0";
      }
    }
  }

  for (pt::idx_t alt_w_idx{}; alt_w_idx < alt_aws.size(); alt_w_idx++) {

    const pvt::AW &alt_w = alt_aws[alt_w_idx];
    for (pt::id_t alt_hap_id : alt_w.get_ref_ids()){

      pt::idx_t alt_hap_col_idx = gtd.ref_id_to_col_idx.at(alt_hap_id);
      gt_col &alt_col = cols[alt_hap_col_idx];

      const std::set<pt::id_t> &alt_col_sample_refs = g.get_shared_samples(alt_hap_id);

      if (alt_col.empty()) {
        ns++; // increment the number of samples with data
        alt_col.resize(alt_col_sample_refs.size(), ".");
      }

      pt::idx_t col_col_idx{};
      for (auto it = alt_col_sample_refs.begin(); it != alt_col_sample_refs.end(); it++, col_col_idx++) {
        if (*it == alt_hap_id) {
          alt_col[col_col_idx] = std::to_string(alt_w_idx + 1);
        }
      }
    }
  }

  // gt_col to string
  std::string s;
  for (auto col : cols) {
    if (col.empty()) {
      // append . and a tab
      s.append(".\t");
    }
    else {
      // join the col with a comma
      std::string joined_col = pu::concat_with(col, '|');
      s.append(joined_col);
      s.append("\t");
    }
  }

  return {ns, s};
}


void write_vcf_rec(const bd::VG &g, const pvt::genotype_data_t &gtd,
                     const pvt::VcfRec &r, const std::string &chrom,
                     std::ostream &os) {
    std::string fn_name{std::format("[{}::{}]", MODULE, __func__)};


    const std::string qual = "60";
    const std::string filter = "PASS";
    pvt::var_type_e var_typ = r.get_var_type();

    auto get_label = [&](const pvt::AS &s) -> std::string {
      if (s.get_o() == pgt::or_e::forward) {
        return g.get_vertex_by_id(s.get_v_id()).get_label();
      } else {
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

    auto ats_to_dna_str = [&](const std::vector<pvt::AW> &allele_walks) -> std::string {
      std::string allele_walks_str = "";
      for (std::size_t w_idx{}; w_idx < allele_walks.size(); ++w_idx) {

        const auto &aw = allele_walks[w_idx];
        // convert to DNA string and takes care of indels
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
        str += std::format("{}{}", s.get_o() == pgt::or_e::forward ? ">" : "<",
                           s.get_v_id());
      }

      return str;
    };

    auto ats_as_str = [&](const std::vector<pvt::AW> &ats) -> std::string {
      std::string str;
      // add a comma between ats
      for (std::size_t i{}; i < ats.size(); ++i) {
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

    // Total number of called alleles (ref + alt), across all samples
    auto comp_AN = [](const pvt::AW &ref_aw, const std::vector<pvt::AW> &alt_aws) -> pt::idx_t {
      pt::idx_t sum {};
      sum = ref_aw.get_ref_ids().size();
      for (const auto &aw : alt_aws) {
        sum += aw.get_ref_ids().size();
      }
      return sum;
    };

    auto comp_AC = [](const std::vector<pvt::AW> &alt_aws) -> std::vector<pt::idx_t> {
      std::vector<pt::idx_t> res;
      for (const auto &aw : alt_aws) {
        res.push_back(aw.get_ref_ids().size());
      }
      return res;
    };

    auto comp_AF = [](pt::idx_t an,const std::vector<pt::idx_t> &ac) -> std::vector<double> {
      std::vector<double> res;
      res.reserve(ac.size()); // avoid reallocations
      for (pt::idx_t x : ac) {
        // cast to double so we get a real fraction
        res.push_back(static_cast<double>(x) / static_cast<double>(an));
      }
      return res;
    };

    pt::idx_t pos = var_typ == pvt::var_type_e::del || var_typ == pvt::var_type_e::ins
      ? r.get_pos() - 1 : r.get_pos();

    // Both ref_dna and alt_dna are plain std::string values over the DNA letters {A, C, G, T}.
    //   ref_dna  is a single contiguous sequence.
    //   alt_dna  may contain one or more such sequences, separated by commas.
    std::string ref_dna = at_to_dna_str(r.get_ref_at());
    std::string alt_dna = ats_to_dna_str(r.get_alt_ats());

    // INFO fields & genotypes
    pt::idx_t an = comp_AN(r.get_ref_at(), r.get_alt_ats());
    std::string an_str = std::to_string(an);
    std::vector<pt::idx_t> ac = comp_AC(r.get_alt_ats());
    std::string ac_str = pu::concat_with(ac, ',');
    std::vector<double> af = comp_AF(an, ac);
    std::string af_str = pu::concat_with(af, ',');

    // gen genotype cols
    auto [ns_str, gt_cols] = gen_genotype_cols(g, gtd.ref_id_to_col_idx.size(),
                                               gtd, r.get_ref_id(),
                                               r.get_ref_at(), r.get_alt_ats());

    std::ostringstream info_field;
    info_field << "AC=" << ac_str
               << ";AF=" << af_str
               << ";AN=" << an_str
               << ";NS=" << ns_str
               << ";AT=" << fmt_field(r)
               << ";VARTYPE=" << pvt::to_string_view(var_typ)
               << ";TANGLED=" << (r.is_tangled() ? "T" : "F")
               << ";LV=" << (r.get_height() - 1);

    // std::string info_field_ = std::format(
    //     "AC={},AF={},AN={},NS={},AT={},VARTYPE={},TANGLED={}", ac_str, af_str, an_str, ns_str,
    //     fmt_field(r), pvt::to_string_view(var_typ), r.is_tangled() ? "T" : "F");

    os << chrom
       << "\t" << pos
       << "\t" << r.get_id()
       << "\t" << ref_dna
       << "\t" << alt_dna
       << "\t" << qual
       << "\t" << filter
       << "\t" << info_field.str()
       << "\t" << "GT"
       << "\t" << gt_cols
       << "\n";
  }

void write_vcf(const bd::VG &g, pt::id_t ref_id, const std::string &chrom,
               const pvt::genotype_data_t &gtd, std::vector<pvt::VcfRec> recs,
               std::ostream &os) {
    std::string fn_name{std::format("[{}::{}]", MODULE, __func__)};

    write_single_header(chrom, g.get_ref_by_id(ref_id).get_length(), os);
    write_col_header(gtd, os);
    for (const pvt::VcfRec &r : recs) {
      write_vcf_rec(g, gtd, r, chrom, os);
    }

    return;
  }

void write_vcfs(const pvt::VcfRecIdx &vcf_recs, const bd::VG &g, const core::config &app_config) {
  std::string fn_name{std::format("[{}::{}]", MODULE, __func__)};

  const pvt::genotype_data_t &gtd = vcf_recs.get_genotype_data();
  const std::vector<std::string> &ref_paths = app_config.get_reference_paths();

  auto is_in_ref_paths = [&](const std::string &ref_name) -> bool {
    return std::find(ref_paths.begin(), ref_paths.end(), ref_name) != ref_paths.end();
  };

  if (app_config.get_stdout_vcf()) {
    // Write to stdout
    write_combined_vcf_to_stdout(vcf_recs, g, app_config);
  } else {
    // Write to separate files
    std::string out_dir = std::string(app_config.get_output_dir());
    
    for (const auto &[ref_id, recs] : vcf_recs.get_recs()) {
      std::string ref_name = g.get_ref_label(ref_id);

      if (!is_in_ref_paths(ref_name)) {
        continue;
      }

      std::string vcf_fp = std::format("{}/{}.vcf", out_dir, ref_name);
      std::ofstream os(vcf_fp);
      write_vcf(g, ref_id, ref_name, gtd, recs, os);
      std::cerr << "wrote " << vcf_fp << "\n";
    }
  }

  return;
}

void write_combined_vcf_to_stdout(const pvt::VcfRecIdx &vcf_recs, const bd::VG &g, const core::config &app_config) {
  std::string fn_name{std::format("[{}::{}]", MODULE, __func__)};

  const pvt::genotype_data_t &gtd = vcf_recs.get_genotype_data();
  const std::vector<std::string> &ref_paths = app_config.get_reference_paths();

  auto is_in_ref_paths = [&](const std::string &ref_name) -> bool {
    return std::find(ref_paths.begin(), ref_paths.end(), ref_name) != ref_paths.end();
  };

  // Collect all contigs for combined header
  std::vector<std::pair<std::string, pt::idx_t>> contigs;
  for (const auto &[ref_id, recs] : vcf_recs.get_recs()) {
    std::string ref_name = g.get_ref_label(ref_id);
    if (is_in_ref_paths(ref_name)) {
      contigs.emplace_back(ref_name, g.get_ref_by_id(ref_id).get_length());
    }
  }
  
  // Write combined header to stdout
  write_combined_header(contigs, std::cout);
  write_col_header(gtd, std::cout);
  
  // Write all records to stdout
  for (const auto &[ref_id, recs] : vcf_recs.get_recs()) {
    std::string ref_name = g.get_ref_label(ref_id);
    if (!is_in_ref_paths(ref_name)) {
      continue;
    }
    for (const pvt::VcfRec &r : recs) {
      write_vcf_rec(g, gtd, r, ref_name, std::cout);
    }
  }

  return;
}

} // namespacepovu::io::vcf
