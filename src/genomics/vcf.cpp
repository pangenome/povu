#include <algorithm>
#include <cstddef>
#include <cstdlib>
#include <map>
#include <tuple>
#include <utility>
#include <vector>
#include <set>
#include <format>
#include <string>

#include "./genomics.hpp"
#include "../graph/bidirected.hpp"


namespace povu::genomics::vcf {
using povu::graph_types::walk;
namespace bd = povu::bidirected;
namespace pt = povu::types;
namespace put = povu::untangle; // TODO: rename

typedef std::tuple < Bubble,
                     std::map<pt::id_t, std::vector<put::exp_cnv>>,
                     std::map<pt::id_t, std::set<std::size_t>>> meta_bub;

vcf::vcf_record gen_vcf_rec(const bd::VG& bd_vg,
                            const genomics::Bubble& c_bub,
                            const std::vector<genomics::variant_type>& variant_cats,
                            const std::string& record_id,
                            const put::exp_cnv& exp_cnv) {
  std::string fn_name{std::format("[povu::genomics::vcf::{}]", __func__)};

  auto [pos, w_idx, sp] = exp_cnv;
  const genomics::Path& p { c_bub.get_path(w_idx) };

  vcf::vcf_record vcf_rec;

  // ---
  // pos
  // ---
  bool has_sub { std::find(variant_cats.begin(), variant_cats.end(), genomics::variant_type::SUB) != variant_cats.end() };
  vcf_rec.pos = (has_sub ? pos : pos - 1 );

  // ---
  // id
  // ---
  vcf_rec.id = record_id;

  // ---
  // ref
  // ---
  vcf_rec.ref = p.as_DNA_str(bd_vg, sp, variant_cats, record_id);

  // ---
  // alt
  // ---
  std::vector<std::string> alts;
  for (std::size_t i{}; i < c_bub.size(); i++) {
    if (i == w_idx) { continue; }
    const genomics::Path& alt_p { c_bub.get_path(i) };
    std::string alt_str { alt_p.as_DNA_str(bd_vg) };

    if (variant_cats[i] == genomics::variant_type::DEL) {
      alt_str = alt_p.as_DNA_str(bd_vg, 0).back();
      vcf_rec.ref = alt_str + vcf_rec.ref;
    }
    else if (variant_cats[i] == genomics::variant_type::INS) {
      alt_str = vcf_rec.ref + alt_str;
    }
    else if (vcf_rec.ref.empty()) {
      vcf_rec.ref = alt_p.as_DNA_str(bd_vg, 0).back();
      alt_str = vcf_rec.ref + alt_str;
    }

    alts.push_back(alt_str);
  }
  vcf_rec.alt = alts;

  // ------
  // format
  // ------

  // .....................
  //
  // AT (allele traversal)
  // .....................
  std::string allele_traversal {"AT="};
  allele_traversal += p.as_str(bd_vg, sp);
  for (std::size_t i {}; i < c_bub.size(); i++ ) {
    if (i == w_idx) { continue; }
    allele_traversal += ",";
    allele_traversal += c_bub.get_path(i).as_str(bd_vg);
  }
  vcf_rec.format += allele_traversal;

  return vcf_rec;
}

std::vector<genomics::variant_type> categorize_variants(const genomics::Bubble& c_bub, pt::id_t w_idx) {
  std::string fn_name{std::format("[povu::genomics::vcf::{}]", __func__)};

  const genomics::Path& p { c_bub.get_path(w_idx) };
  walk const& w { p.get_walk() };
  std::vector<genomics::variant_type> variant_types;
  variant_types.reserve(c_bub.size() );

  for (std::size_t i {}; i < c_bub.size(); i++ ) {
    if (i == w_idx) {
      variant_types.push_back(genomics::variant_type::INVALID);
      continue;
    }

    genomics::Path const& alt_p { c_bub.get_path(i) };

    if (w.size() == 2 ) {
      variant_types.push_back(genomics::variant_type::INS);
    }
    else if (alt_p.get_walk().size() == 2) {
      variant_types.push_back(genomics::variant_type::DEL);
    }
    else if (alt_p.get_walk().size() > 2 && w.size() > 2) {
      variant_types.push_back(genomics::variant_type::SUB);
    }
    else {
      variant_types.push_back(genomics::variant_type::UNDEFINED);
    }
  }

  return variant_types;
}


/**
 * @brief Generate VCF records for a bubble
 */
void populate_vcf_rec(const bd::VG& bd_vg,
                      std::map<std::size_t, std::vector<vcf::vcf_record>>& vcf_records,
                      const meta_bub& meta_c_bub,
                      const std::string& record_id,
                      std::size_t hap_id) {
  std::string fn_name { std::format("[povu::genomics::vcf::{}]",  __func__) };

  std::vector<vcf::vcf_record> vcf_recs;
  std::vector<genomics::variant_type> variant_cats;

  auto [c_bub, unrolled, keeping ] = meta_c_bub;
  const std::vector<put::exp_cnv>& un = unrolled[hap_id];
  const std::set<std::size_t>& k = keeping[hap_id];

  vcf::vcf_record vcf_rec;
  for (std::size_t i{}; i < un.size(); i++) {
    if (k.find(i) == k.end()) { continue; }
    auto [_, w_idx, __] = un[i];
    variant_cats = categorize_variants(c_bub, w_idx);
    vcf_rec = gen_vcf_rec(bd_vg, c_bub, variant_cats, record_id, un[i]);
    vcf_records[hap_id].push_back(vcf_rec);
  }
}


/**
 * @brief Generate VCF records for a bubble as a map of haplotype IDs to VCF records
 */
std::map<std::size_t, std::vector<vcf::vcf_record>> gen_bub_vcf_records(const bd::VG &bd_vg,
                                                                        const meta_bub& meta_c_bub,
                                                                        const core::config& app_config) {
  std::string fn_name { std::format("[povu::genomics::vcf::{}]",  __func__) };

  auto [c_bub, _, __ ] = meta_c_bub;

  // all vcf records in this bubble will have the same record ID
  // ID col in the VCF file
  std::string record_id { std::format("{}{}", c_bub.start().as_str(), c_bub.end().as_str()) };

  // all our VCF files for this bubble will be in this map
  std::map<std::size_t, std::vector<vcf::vcf_record>> vcf_records;
  // for each ref that goes through the bubble
  for (std::size_t ref_id : c_bub.get_hap_ids()) {
    populate_vcf_rec(bd_vg, vcf_records, meta_c_bub, record_id, ref_id);
  }

  return vcf_records;
}

std::map<std::size_t, std::vector<vcf::vcf_record>> gen_vcf_records(const bd::VG &bd_vg,
                                                                    const std::vector<meta_bub>& c_bubs,
                                                                    const core::config& app_config) {
  std::string fn_name { std::format("[povu::genomics::vcf::{}]",  __func__) };

  std::map<std::size_t, std::vector<vcf::vcf_record>> all_vcf_records;
  // map of haplotype IDs to VCF records for this bubble
  std::map<std::size_t, std::vector<vcf::vcf_record>> b_vcf_records;
  for (const auto& m : c_bubs) { // for each bubble
    b_vcf_records = gen_bub_vcf_records(bd_vg, m, app_config);

    for (auto& [hap_id, vcf_recs] : b_vcf_records) {
      all_vcf_records[hap_id].insert(all_vcf_records[hap_id].end(), vcf_recs.begin(), vcf_recs.end());
    }
  }

  return all_vcf_records;
}

}; // namespace genomics::vcf
