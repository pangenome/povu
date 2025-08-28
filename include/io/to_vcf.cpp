#include "./to_vcf.hpp"
#include <algorithm>

#include <fstream>
#include <ostream>
#include <set>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace povu::io::to_vcf {

void write_header_common(std::ostream &os) {
  os << "##fileformat=VCFv4.2\n";
  os << pv_cmp::format("##fileDate={}\n", pu::today());
  os << "##source=povu\n";
  os << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
  os << "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Total number of alternate alleles in called genotypes\">\n";
  os << "##INFO=<ID=AT,Number=R,Type=String,Description=\"Allele traversal path through the graph\">\n";
  os << "##INFO=<ID=AN,Number=1,Type=String,Description=\"Total number of alleles in called genotypes\">\n";
  os << "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele frequency in the population\">\n";
  os << "##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of samples with data\">\n";
  os << "##INFO=<ID=VARTYPE,Number=1,Type=String,Description=\"Type of variation: INS (insertion), DEL (deletion), SUB (substitution)\">\n";
  os << "##INFO=<ID=TANGLED,Number=1,Type=String,Description=\"Variant lies in a tangled region of the graph: T or F\">\n";
  os << "##INFO=<ID=LV,Number=1,Type=Integer,Description=\"Level in the PVST (0=top level)\">\n";
  os << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";

  return;
}

void write_header_contig_line(const pr::Ref &r, std::ostream &os) {
  os << pv_cmp::format("##contig=<ID={},length={}>\n", r.tag(), r.get_length());
  return;
}


void write_col_header(pgv::genotype_data_t gtd, std::ostream &os) {
  os << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t";
  for (std::size_t i = 0; i < gtd.genotype_cols.size(); ++i) {
    os << gtd.genotype_cols[i];
    if (i < gtd.genotype_cols.size() - 1) {
      os << "\t";
    }
  }
  os << "\n";
}

// TODO: should this be the same across all VCFs?
povu::genomics::vcf::genotype_data_t comp_gt(const bd::VG &g) {
  povu::genomics::vcf::genotype_data_t gd;

  std::set<pt::id_t> handled;

  for (pt::id_t ref_id = 0; ref_id < g.ref_count(); ++ref_id) {
    if (pv_cmp::contains(handled, ref_id)) {
      continue;
    }

    handled.insert(ref_id);

    const std::set<pt::id_t> &sample_refs =  g.get_shared_samples(ref_id);

    std::string col_name = g.get_ref_by_id(ref_id).get_sample_name();

    if (sample_refs.empty()) {
      // throw an exeption
      ERR("No sample names found for ref_id {}", ref_id);
      std::exit(EXIT_FAILURE);
    }
    else if (sample_refs.size() == 1) {
      // throw an exeption
      gd.ref_id_to_col_idx[ref_id] = gd.genotype_cols.size();
      gd.genotype_cols.push_back(col_name);
    }
    else if (sample_refs.size() > 1) {
      for (pt::id_t ref_id_ : sample_refs) {
        gd.ref_id_to_col_idx[ref_id_] = gd.genotype_cols.size();
        handled.insert(ref_id_);
      }
      gd.genotype_cols.push_back(col_name);
    }
  }

  return gd;
}


void init_vcfs(bd::VG &g, const std::vector<std::string> &sample_names, VcfOutput &vout) {

  vout.for_each_stream([&](std::ostream &os) {
    write_header_common(os); // write common header lines
  });

  // add contig lines
  for (const auto &sample_name : sample_names) {
    std::ostream &os = vout.stream_for(sample_name);
    std::set<pt::id_t> ref_ids = g.get_refs_in_sample(sample_name);
    for (pt::id_t ref_id : ref_ids) {
      const pr::Ref &ref = g.get_ref_by_id(ref_id);
      write_header_contig_line(ref, os);
    }
  }

  vout.for_each_stream([&](std::ostream &os) {
    write_col_header(comp_gt(g), os); // write column header
  });

  vout.flush_all();

  return;
}
// ns and the genotype columns are generated from the genotype data
std::pair<pt::idx_t, std::string> gen_genotype_cols(const bd::VG &g,
                                                    const pgv::genotype_data_t &gtd,
                                                    const pga::AW &ref_aw,
                                                    const std::vector<pga::AW> &alt_aws) {
  std::string fn_name{pv_cmp::format("[{}::{}]", MODULE, __func__)};

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

    const pga::AW &alt_w = alt_aws[alt_w_idx];
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


void write_vcf_rec(const bd::VG &g, const pgv::genotype_data_t &gtd,
                     const pgv::VcfRec &r, const std::string &chrom,
                     std::ostream &os) {
    std::string fn_name{pv_cmp::format("[{}::{}]", MODULE, __func__)};


    const std::string qual = "60";
    const std::string filter = "PASS";
    pgv::var_type_e var_typ = r.get_var_type();

    auto get_label = [&](const pga::AS &s) -> std::string {
      if (s.get_o() == pgt::or_e::forward) {
        return g.get_vertex_by_id(s.get_v_id()).get_label();
      } else {
        return g.get_vertex_by_id(s.get_v_id()).get_rc_label();
      }
    };

    auto at_to_dna_str = [&](const pga::AW &at) -> std::string {
      std::string at_str = "";

      // 1) Anchor base for deletions & insertions
      switch (var_typ) {
      case pgv::var_type_e::del:
      case pgv::var_type_e::ins: {
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
        const pga::AS &s = at.get_step(step_idx);
        at_str += get_label(s);
      }

      // 3) If nothing got added, use “.”
      return at_str.empty() ? std::string{"."} : at_str;
    };

    auto ats_to_dna_str = [&](const std::vector<pga::AW> &allele_walks) -> std::string {
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

    auto at_as_str = [](const pga::AW &at) -> std::string {
      std::string str;

      for (auto &s : at.get_steps()) {
        str += pv_cmp::format("{}{}", s.get_o() == pgt::or_e::forward ? ">" : "<",
                           s.get_v_id());
      }

      return str;
    };

    auto ats_as_str = [&](const std::vector<pga::AW> &ats) -> std::string {
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

    // takes r by reference from outer scope
    auto fmt_field = [&]() -> std::string {
      std::string s;
      s += at_as_str(r.get_ref_at());
      s += ",";
      s += ats_as_str(r.get_alt_ats());

      return s;
    };

    // Total number of called alleles (ref + alt), across all samples
    auto comp_AN = [](const pga::AW &ref_aw, const std::vector<pga::AW> &alt_aws) -> pt::idx_t {
      pt::idx_t sum {};
      sum = ref_aw.get_ref_ids().size();
      for (const auto &aw : alt_aws) {
        sum += aw.get_ref_ids().size();
      }
      return sum;
    };

    auto comp_AC = [](const std::vector<pga::AW> &alt_aws) -> std::vector<pt::idx_t> {
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

    pt::idx_t pos = var_typ == pgv::var_type_e::del || var_typ == pgv::var_type_e::ins
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
    auto [ns_str, gt_cols] = gen_genotype_cols(g, gtd, r.get_ref_at(), r.get_alt_ats());

    std::string at_str = fmt_field();

    std::ostringstream info_field;
    info_field << "AC=" << ac_str
               << ";AF=" << af_str
               << ";AN=" << an_str
               << ";NS=" << ns_str
               << ";AT=" << at_str
               << ";VARTYPE=" << pgv::to_string_view(var_typ)
               << ";TANGLED=" << (r.is_tangled() ? "T" : "F")
               << ";LV=" << (r.get_height() - 1);

    // std::string info_field_ = pv_cmp::format(
    //     "AC={},AF={},AN={},NS={},AT={},VARTYPE={},TANGLED={}", ac_str, af_str, an_str, ns_str,
    //     fmt_field(r), pgv::to_string_view(var_typ), r.is_tangled() ? "T" : "F");

    os << chrom
       << "\t" << pos
       << "\t" << r.get_id()
       << "\t" << ref_dna
       << "\t" << alt_dna
       << "\t" << qual
       << "\t" << filter
       << "\t" << info_field.str()
       << "\t" << "GT" // TODO: [c] make this a const
       << "\t" << gt_cols
       << "\n";
    return;
}



void write_vcf(const bd::VG &g, pt::id_t ref_id, const std::string &chrom,
               const pgv::genotype_data_t &gtd, std::vector<pgv::VcfRec> recs,
               std::ostream &os) {
  for (const pgv::VcfRec &r : recs) {
    write_vcf_rec(g, gtd, r, chrom, os);
  }

  return;
}



void write_vcfs(const pgv::VcfRecIdx &vcf_recs, const bd::VG &g,
                std::set<pt::id_t> vcf_ref_ids,
                VcfOutput &vout, const core::config &app_config) {

  const pgv::genotype_data_t &gtd = vcf_recs.get_genotype_data();
  const std::vector<std::string> &ref_paths = app_config.get_reference_paths();

  auto is_in_ref_paths = [&](const std::string &ref_name) -> bool {
    return std::find(ref_paths.begin(), ref_paths.end(), ref_name) != ref_paths.end();
  };

  if (app_config.get_stdout_vcf()) {
    // Write to stdout
    std::ostream &os = vout.stream_for("");
    const pgv::genotype_data_t &gtd = vcf_recs.get_genotype_data();
    for (const auto &[ref_id, recs] : vcf_recs.get_recs()) {
      if (!pv_cmp::contains(vcf_ref_ids, ref_id)) {
        continue;
      }
      const std::string &ref_name = g.get_sample_name(ref_id);
      write_vcf(g, ref_id, ref_name, gtd, recs, os);
    }
  }
  else {

    for (const auto &[ref_id, recs] : vcf_recs.get_recs()) {
      if (!pv_cmp::contains(vcf_ref_ids, ref_id)) {
        continue;
      }
      std::string sample_name = g.get_sample_name(ref_id);
      std::ostream &os = vout.stream_for(sample_name);
      write_vcf(g, ref_id, sample_name, gtd, recs, os);
    }
  }

  return;
}
} // namespacepovu::io::vcf
