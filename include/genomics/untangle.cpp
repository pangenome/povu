#include <algorithm>
#include <cmath>  // For std::round
//#include <concepts>
#include <cstddef>
#include <cstdlib>
#include <format>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "WFAligner.hpp"

#include "../algorithms/algorithms.hpp"
#include "./genomics.hpp"
#include "../graph/bidirected.hpp"


namespace povu::untangle {
namespace pgt = povu::graph_types;


/**
 * @brief Get the path metadata for a haplotype/reference in a bubble
 *
 * @param c_bub The bubble
 * @param ref_id The haplotype/reference ID
 */
std::vector<ref_meta> find_ref_meta(const pg::Bubble& c_bub, std::size_t ref_id) {
  std::string fn_name{std::format("[povu::untangle::vcf::{}]", __func__)};

  std::vector<ref_meta> m;
  const std::set<std::size_t>& w_idxs { c_bub.get_walk_idxs(ref_id) };

  // is a span a prefix of another span?
  auto is_pre = [](const pgt::walk& qry_w, const pt::span& qry_span, const pgt::walk& w, const pt::span& s) -> bool {
    auto [qs, ql] = qry_span;
    auto [ts, tl] = s;

    std::size_t len = std::min(ql, tl);

    for (std::size_t i {qs}, j {ts}; i < len && j < len; j++, i++) {
      if (qry_w[i] != w[j]) { return false; }
    }
    if (ql > tl) { return false; }

    return true;
  };

  std::set<std::size_t> prefix_spans;
  for (std::size_t qry_w_idx : w_idxs) {
    const pgt::walk& qry_wlk = c_bub.get_path(qry_w_idx).get_walk() ;
    const pt::span& qry_span = c_bub.get_path(qry_w_idx).get_hap_span(ref_id);

    for (std::size_t idx : w_idxs) {
      if (idx == qry_w_idx) { continue; }
      const pgt::walk& w = c_bub.get_path(idx).get_walk();
      const pt::span& s = c_bub.get_path(idx).get_hap_span(ref_id);

      if (is_pre(qry_wlk, qry_span, w, s)) { prefix_spans.insert(qry_w_idx); }
    }
  }

  // filter out spans that are prefixes of existing spans
  std::set<std::size_t> diff;
  // Perform the set difference (A - B)
  std::set_difference(w_idxs.begin(), w_idxs.end(), prefix_spans.begin(), prefix_spans.end(),
                      std::inserter(diff, diff.end()));

  // ------------------------------
  // extract metadata for each walk
  // ------------------------------
  for (std::size_t w_idx : diff) {
    const pg::Path& p { c_bub.get_path(w_idx) };
    const pt::span& sp { p.get_hap_span(ref_id) };
    m.push_back({p, w_idx, sp});
  }

  return m;
}


/**
 * @brief untangling i.e. get the position(s) of a variant in the given reference
 *
 * invalidate some of the positions in a walk if they cannot be extended
 * The position as pertains the VCF format
 *
 * @param bd_vg The bidirected graph
 * @param meta The reference metadata
 * @param hap_id The haplotype/reference ID
 * @return The position(s) of the variant in the reference
 */
std::vector<pt::idx_t> get_positions(const bd::VG &bd_vg, const ref_meta& meta, pt::id_t hap_id) {
  std::string fn_name{std::format("[povu::untangle::vcf::{}]", __func__)};

  auto &[p, _, __] = meta;
  pt::span sp = meta.w_span;

  const std::vector<pgt::id_n_orientation_t>& w = p.get_walk();
  auto [strt_v_idx, ___] = w[sp.start];

  // invalidate some of the positions in a walk if they cannot be extended

  //bd::Vertex w_start_vtx = bd_vg.get_vertex(strt_v_idx);

  bd::Vertex v = bd_vg.get_vertex(strt_v_idx);
  std::vector<bd::PathInfo> ps { v.get_paths() };

  // true if the position is valid, false otherwise.
  auto pred = [&](const bd::PathInfo& pi, std::size_t len) -> bool {

    auto [ref_id, si, _] = pi;

    if (ref_id != hap_id) { return false; }

    std::size_t curr_si = si;

    // for each vertex in the span of the walk after the first
    for (std::size_t i { sp.start + 1 }; i < sp.start + sp.length; i++) {

      auto [curr_v_idx, __] = w[i];
      const bd::Vertex& curr_v = bd_vg.get_vertex(curr_v_idx);
      const std::vector<bd::PathInfo>& curr_ps {  curr_v.get_paths() };

      bool found { false };

      for (const bd::PathInfo& cpi : curr_ps) {

        auto [c_ref_id, cv_si, _] = cpi;
        if (c_ref_id != ref_id) { continue; }

        if (cv_si == curr_si + len) {
          found = true;
          break;
        }
      }

      if (!found) { return false; }

      curr_si += len;
      len = curr_v.get_label().length();
    }

    return true;
  };


  // check that the path id is valid and that
  // the position is valid for the walk i.e the position can be extended in the
  // next vertex in the walk
  std::size_t l = v.get_label().length();
  ps.erase(std::remove_if(ps.begin(), ps.end(), [&](bd::PathInfo& pi) { return !pred(pi, l); } ), ps.end());

  for (bd::PathInfo& pi : ps) { pi.step_index += l; }

  std::vector<pt::idx_t> positions;
  for (const bd::PathInfo& pi : ps) { positions.push_back(pi.step_index); }

  return positions;
}


/**
 * @brief Convert a step to an ASCII character
 *
 * @param bd_vg The bidirected graph
 * @param step The step
 * @param bub_min The minimum vertex ID in the bubble
 * @param bub_max The maximum vertex ID in the bubble
 * @return The ASCII representation of the step
 */
char step_to_char(const bd::VG& bd_vg, pgt::id_or step, pt::id_t bub_min, pt::id_t bub_max) {
  std::string fn_name{std::format("[povu::untangle::{}]", __func__)};

  // Define the ASCII printable range
  const int ascii_min = 65; // 32;
  const int ascii_max = 90; // 126;

  auto [v_idx, o] = step;
  std::size_t num = bd_vg.idx_to_id(v_idx);

  if (num < bub_min || num > bub_max) {
    throw std::runtime_error(std::format("{} vertex {} out of range min {} max {}", fn_name, num, bub_min, bub_max));
  }

  // Perform linear mapping
  double mapped_value = ascii_min + (static_cast<double>(num - bub_min) * (ascii_max - ascii_min)) / (bub_max - bub_min);
  // Round the result to the nearest integer
  int val = static_cast<int>(std::round(mapped_value));

  return static_cast<char>(val);
}


/**
 * @brief expects the expanded cnvs are already sorted by position
 *
 * @param bd_vg The bidirected graph
 * @param c_bub The bubble
 * @param unrolled The expanded CNVs
 * @return a map of reference IDs to the ASCII representation of the expanded CNVs
 */
std::map<pt::id_t, std::string> exp_cnv_to_ascii(const bd::VG &bd_vg,
                                                 const pg::Bubble &c_bub,
                                                 const std::map<pt::id_t, std::vector<exp_cnv>>& unrolled) {
  std::string fn_name{std::format("[povu::untangle::{}]", __func__)};

  auto walk_to_ascii = [&](pt::idx_t w_idx, pt::span sp) -> std::string {
    const pgt::walk w = c_bub.get_path(w_idx).get_walk();

    std::string at_str;
    for (std::size_t i{sp.start}; i < sp.length; ++i) {
      at_str += step_to_char(bd_vg, w[i], c_bub.start().v_idx, c_bub.end().v_idx);
    }

    return at_str;
  };

  std::map<pt::id_t, std::string> ascii;

  for (pt::id_t ref_id : c_bub.get_hap_ids()) {

    const std::vector<std::tuple<pt::idx_t, pt::idx_t, pt::span>>& expanded = unrolled.at(ref_id);

    std::string s;

    for (auto &[_, w_idx, sp] : expanded) {
      s += walk_to_ascii(w_idx, sp);
    }

    ascii[ref_id] = s;
  }

  if (c_bub.start().v_idx == 6 && c_bub.end().v_idx == 8) {
      for (pt::id_t ref_id : c_bub.get_hap_ids()) {
        std::cerr << "name: " << bd_vg.get_path(ref_id).name << "  ";
        // std::vector<cnv>& cnv_ = bub_cnvs[ref_id];
        std::cerr << ascii[ref_id];
        std::cerr << std::endl;
      }
    }

  return ascii;
}


std::map<pt::id_t, std::vector<exp_cnv>> linearise(const bd::VG &bd_vg,
                                                   const pg::Bubble &c_bub,
                                                   std::map<pt::id_t, std::vector<cnv>> bub_cnvs) {

  std::map<pt::id_t, std::vector<exp_cnv>> unrolled;

  for (pt::id_t ref_id : c_bub.get_hap_ids()) {
    std::vector<cnv>& cnv_ = bub_cnvs[ref_id];

    std::vector<exp_cnv> expanded;

    for (auto& [positions, w_idx, sp] : cnv_) {
      for (auto pos : positions) {
        expanded.push_back(std::make_tuple(pos, w_idx, sp));
      }
    }

    // sort by position
    std::sort(expanded.begin(), expanded.end(), [](const auto& a, const auto& b) {
      return std::get<0>(a) < std::get<0>(b);
    });

    unrolled[ref_id] = expanded;
  }

  if (c_bub.start().v_idx == 6 && c_bub.end().v_idx == 8) {
    for (auto [ref_id, v] : unrolled) {
      //std::cerr << "ref: " << bd_vg.get_ref(ref_id).name << " count: " << v.size() << std::endl;
      std::cerr << "ref: " << bd_vg.get_ref(ref_id).name << " pos: ";
      for (auto [pos, _, __] : v) {
        std::cerr << pos << " ";
      }
      std::cerr << std::endl;
    }
  }

  return unrolled;
}


/**
 * @brief Generate VCF records for a bubble
 */
std::map<pt::id_t, std::vector<cnv>> untangle_bub(const bd::VG& bd_vg, const pg::Bubble& c_bub) {
  std::string fn_name { std::format("[povu::untangle::{}]",  __func__) };

  std::map<pt::id_t, std::vector<cnv>> bub_cnvs;
  for (pt::id_t ref_id : c_bub.get_hap_ids()) {
    std::vector<ref_meta> m { find_ref_meta(c_bub, ref_id) };
    std::vector<pt::idx_t> positions; // TODO: rename to coordinates?
    for (const auto& i : m) {
      positions = get_positions(bd_vg, i, ref_id);
      bub_cnvs[ref_id].push_back({positions, i.w_idx, i.w_span});
    }
  }

  return bub_cnvs;
}


/**
 * @brief given an alignment, flip the query and text
 *
 * @param es edit script or edit transcript
 * @return the "flipped" edit transcript
 */
inline std::string flip_aln(const std::string& es) {
  std::string fes;
  for (char c : es) {
    switch (c) {
    case 'M':
      fes.push_back(c);
      break;
    case 'X':
      fes.push_back(c);
      break;
    case 'D':
      fes.push_back('I');
      break;
    case 'I':
      fes.push_back('D');
      break;
    default:
      break;
    }
  }

  return fes;
}


std::map<pt::Pair<pt::id_t>, std::string> run_align(const bd::VG& bd_vg,
                                                    std::map<pt::id_t, std::string>& cnvs,
                                                    wfa::WFAlignerGapAffine& aligner,
                                                    const core::config& app_config,
                                                    const pg::Bubble& c_bub) {
  std::string fn_name{std::format("[povu::untangle::{}]", __func__)};

  // the key is q,t pair and the value is the cigar string of the alignment
  std::map<pt::Pair<pt::id_t>, std::string> aln_res;

  for (const std::string& s : app_config.get_reference_paths()) {
    pt::id_t ref_id = bd_vg.get_ref(s).id;
    // make sure that the ref path exists in the bubble
    if (cnvs.find(ref_id) == cnvs.end()) { continue; }
    const std::string& t = cnvs.at(ref_id);

    for (const auto& [cnv_ref_id, q] : cnvs) {
      if (ref_id == cnv_ref_id) { continue; }

      // update the map with the edit script of the alignment
      if (aln_res.find(pt::Pair<pt::id_t>{ref_id, cnv_ref_id}) != aln_res.end()) {
        aln_res[pt::Pair<pt::id_t>{cnv_ref_id, ref_id}]
          = flip_aln(aln_res[pt::Pair<pt::id_t>{ref_id, cnv_ref_id}]);
      }
      else {
        aln_res[pt::Pair<pt::id_t>{cnv_ref_id, ref_id}]
          = povu::align::wfa2(aligner, std::move(q), std::move(t));
      }
    }
  }

  if (c_bub.start().v_idx == 6 && c_bub.end().v_idx == 8) {
    for (auto [ref_id, aln] : aln_res) {
      std::cerr << std::format("ref: {} {} \n {}\n",
                               bd_vg.get_ref(ref_id.second).name, bd_vg.get_ref(ref_id.first).name, aln);
    }
  }

  return aln_res;
}


/**
 * @brief Filter the expanded CNVs based on the alignment results
 *
 * @param exp The expanded CNVs
 * @param aln The alignment results
 * return A set of indices of the expanded CNVs to keep
*/
inline std::set<std::size_t> is_keep(std::vector<exp_cnv> exp, const std::string &aln) {
  std::set<char> targetChars { 'D', 'X'};

  std::set<std::size_t> s;
  pt::idx_t i{}; pt::idx_t j{};
  while (i < aln.length() && j < exp.size()) {
    auto [_, __, sp] = exp[j];

    bool keep =  std::any_of(aln.begin() + i, aln.begin() + i + sp.length, [&targetChars](char c) {
      return targetChars.contains(c);
    });
    if (keep) { s.insert(j); }

    j++;
    i += sp.length;
  }

  return s;
}


/**
 * @brief Filter the expanded CNVs based on the alignment results
 *
 * @param unrolled The expanded CNVs
 * @param aln_res The alignment results
 * @param c_bub The bubble
 * @param bd_vg The bidirected graph
 * return A map of reference IDs to a set of indices of the expanded CNVs to keep
*/
std::map<std::size_t, std::set<std::size_t>> filter(std::map < pt::id_t,
                                                    std::vector<exp_cnv>> unrolled,
                                                    std::map<pt::Pair<pt::id_t>, std::string> aln_res,
                                                    const pg::Bubble& c_bub, const bd::VG& bd_vg) {
  std::string fn_name{std::format("[povu::untangle::{}]", __func__)};

  std::map<std::size_t, std::set<std::size_t>> to_keep;
  for (auto [ref_id, e_cnv] : unrolled) {
    for (auto [p, aln] : aln_res) {
      if (p.first == ref_id) {
        // append the results of to keep to bar
        auto s = is_keep(e_cnv, aln);
        to_keep[ref_id].insert(s.begin(), s.end());
      }
    }
  }

  if (false && c_bub.start().v_idx == 2701 && c_bub.end().v_idx == 2703) {
    // print to_keep
    std::cerr << fn_name << "\n";
    for (auto [ref_id, s] : to_keep) {
      std::cerr << "ref: " << bd_vg.get_ref(ref_id).name << " keeping: ";
      for (auto i : s) {
        std::cerr << i << " ";
      }
      std::cerr << std::endl;
    }
  }

  return to_keep;
}


/**
 * @brief does the bubble contain a CNV?
 *
 * @param unrolled The expanded CNVs
 * @return true if the bubble contains a CNV, false otherwise
 */
bool has_cnv(const std::map<pt::id_t, std::vector<exp_cnv>> &unrolled) {
  for (const auto& [_, v] : unrolled) {
    if (v.size() > 1) { return true; }
  }

  return false;
}

/**
 * @brief untangle ...
 *
 * @param bd_vg The bidirected graph
 */
std::vector<std::tuple <pg::Bubble, std::map<pt::id_t, std::vector<exp_cnv>>, std::map<std::size_t, std::set<std::size_t>>>>
untangle(const bd::VG& bd_vg, const std::vector<pg::Bubble>& c_bubs, const core::config& app_config) {
  std::string fn_name { std::format("[povu::untangle::vcf::{}]",  __func__) };

  std::vector<std::tuple <pg::Bubble, std::map<pt::id_t, std::vector<exp_cnv>>, std::map<std::size_t, std::set<std::size_t>>>> res;

  // Create a WFAligner object
  wfa::WFAlignerGapAffine aligner(4,6,2,wfa::WFAligner::Alignment,wfa::WFAligner::MemoryHigh);

  for (pg::Bubble const& c_bub : c_bubs) { // for each bubble
    if (c_bub.size() < 2) { continue; } // TODO: avoid these cases from getting here
    if (c_bub.get_hap_ids().empty()) { continue; } // if the bubble has no reference paths

    std::map<pt::id_t, std::vector<cnv>> bub_cnvs = untangle_bub(bd_vg, c_bub);
    std::map<pt::id_t, std::vector<exp_cnv>> unrolled = linearise(bd_vg, c_bub, bub_cnvs);

    if (has_cnv(unrolled)) {
      std::map<pt::id_t, std::string> cnv_as_strs = exp_cnv_to_ascii(bd_vg, c_bub, unrolled);
      std::map<pt::Pair<pt::id_t>, std::string> aln_res = run_align(bd_vg, cnv_as_strs, aligner, app_config, c_bub);
      std::map<std::size_t, std::set<std::size_t>> keeping = filter(unrolled, aln_res, c_bub, bd_vg);
      res.push_back(std::make_tuple(c_bub, unrolled, keeping));
    }
    else {
      std::map<std::size_t, std::set<std::size_t>> keeping;
      // keep the first CNV, that is, everything.
      for (const auto& [ref_id, _] : unrolled) { keeping[ref_id].insert(0); }
      res.push_back(std::make_tuple(c_bub, unrolled, keeping));
    }
  }

  return res;
}

} // namespace untangle
