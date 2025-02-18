#include <algorithm>
#include <cstddef>
#include <ctime>
#include <format>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <sys/types.h>
#include <utility>
#include <vector>

#include "../common/types.hpp"
#include "./genomics.hpp"
#include "../graph/bidirected.hpp"
#include "../io/io.hpp"


namespace povu::genomics {
namespace bd = povu::bidirected;
namespace pt = povu::types;
namespace pc = povu::constants;
namespace pgt = povu::graph_types;
namespace pu = povu::untangle;

// TODO: replace with stride
typedef std::pair<std::size_t, std::size_t> range; // start and length covered by the haplotype

// Overload the << operator
std::ostream &operator<<(std::ostream &os, variant_type vt) {
  switch (vt) {
  case povu::genomics::variant_type::SNP:     os << "SNP"; break;
  case povu::genomics::variant_type::NP:      os << "NP"; break;
  case povu::genomics::variant_type::INDEL:   os << "INDEL"; break;
  case povu::genomics::variant_type::DEL:     os << "DEL"; break;
  case povu::genomics::variant_type::INS:     os << "INS"; break;
  case povu::genomics::variant_type::SUB:     os << "SUB"; break;
  case povu::genomics::variant_type::INV:     os << "INV"; break;
  case povu::genomics::variant_type::DUP:     os << "DUP"; break;
  case povu::genomics::variant_type::CNV:     os << "CNV"; break;
  case povu::genomics::variant_type::BND:     os << "BND"; break;
  case povu::genomics::variant_type::INVALID: os << "INVALID"; break;
  default:                              os << "UNDEFINED"; break;
  }

  return os;
}


/**
 * @brief Find the haplotypes that exist within a walk
 *
 * a reference exists once, the longest one, within a walk
 * a haplotype exists within a walk, W, iff it exists within a continuous range of vertices in W
 */
// TODO: remove
std::map<pt::id_t, pt::Stride> find_walk_refs(const bd::VG& bd_vg,
                                              const pgt::walk& w,
                                              const std::set<pt::id_t>& ref_ids) {
  std::string fn_name{std::format("[povu::genomics::{}]", __func__)};

  // key is the haplotype id and value are the ranges of the haplotype
  std::map<pt::id_t, pt::Stride> m;

  auto get_v_refs = [&](std::size_t v_idx) -> std::set<std::size_t> {
    std::set<std::size_t> refs;
    for (const auto& p : bd_vg.get_vertex_by_idx(v_idx).get_refs()) {
      refs.insert(p.path_id);
    }

    return refs;
  };

  // find the stride of the reference paths in the walk
  for (pt::id_t r : ref_ids) {

    std::vector<pt::Stride> strides;

    const povu::types::Stride INVALID_STRIDE {pc::INVALID_IDX, 1};
    povu::types::Stride stride { INVALID_STRIDE };

    // for each step in the walk
    for (std::size_t s{1}; s < w.size(); ++s) {
      auto v1 = w[s-1];
      std::set<std::size_t> v1_refs { get_v_refs(v1.v_idx) };

      auto v2 = w[s];
      std::set<std::size_t> v2_refs { get_v_refs(v2.v_idx) };

      const bd::Edge& e = bd_vg.get_edge(v1, v2);
      const std::set<std::size_t> &e_refs = e.get_refs();

      if (v1_refs.count(r) && v2_refs.count(r) && e_refs.count(r)) {
        if (stride.start == pc::INVALID_IDX) {
          stride.start = s-1;
        }
        ++stride.length;
      }
      else if (v1_refs.count(r)) {
        if (stride.start != pc::INVALID_IDX) {
          strides.push_back(stride);
          stride = INVALID_STRIDE;
        }
        else {
          stride = {s-1, 1};
          strides.push_back(stride);
          stride = INVALID_STRIDE;
        }
      }
      else {
        if (stride.start != pc::INVALID_IDX) {
          strides.push_back(stride);
          stride = INVALID_STRIDE;
        }
      }

    }

    // get the longest stride
    if (strides.size() > 0) {
      auto it = std::max_element(strides.begin(), strides.end(), [](const pt::Stride& a, const pt::Stride& b) {
        return a.length < b.length;
      });

      stride = *it;
    }

    //std::cerr << std::format("stride: {} {}\n", stride.start, stride.length);

    // check that the stride is valid
    // i.e. it covers entire bubble or more than just the bubble ends
    if ((stride.start == 0 && stride.length > 1) ||
        (stride.start > 0 && stride.start + stride.length == w.size())
        ) {
      m[r] = stride;
    }
  }

  return m;
}



/**
 * @brief Find the reference paths that are in each path of the bubble
 *
 * creates a vector of genomics::Path
 *
 * @param bd_vg The graph
 * @param bub_walks The walks in the bubble
 * @param relevant_refs The reference paths that are in the graph & in the app config
 * @param app_config The application configuration
 * @return The reference paths that are in each path of the bubble
 */
std::vector<Path> find_bubble_refs(const bd::VG& bd_vg,
                                   const std::vector<pgt::walk>& bub_walks,
                                   const std::set<std::size_t>& ref_ids) {
  std::vector<Path> bub_paths_vec;

  // for each walk in the bubble
  for (std::size_t i{}; i < bub_walks.size(); ++i) {
    const pgt::walk& w = bub_walks[i];
    std::map<std::size_t, pt::Stride> w_refs = find_walk_refs(bd_vg, w, ref_ids);
    bub_paths_vec.push_back( { w, w_refs });
  }

  return bub_paths_vec;
}

/**
 * Associate each bubble with a set of haplotypes
 *
 */
std::vector<Bubble> find_haplotypes(
    const bd::VG &bd_vg, const std::vector<std::vector<pgt::walk>> &all_paths,
    const std::vector<pgt::flubble>& canonical_flubbles,
    const std::set<pt::id_t>& ref_ids) {
  std::string fn_name{std::format("[povu::genomics::{}]", __func__)};

  std::vector<Bubble> bubs;
  for (std::size_t b_idx{}; b_idx < all_paths.size(); ++b_idx) { // for each bubble
    auto [entry, exit] = canonical_flubbles[b_idx];
    std::vector<Path> v { find_bubble_refs(bd_vg, all_paths[b_idx], ref_ids) };
    bubs.push_back({entry, exit, v});
  }

  return bubs;
}

// TODO: remove
/**
 * @brief Get the side of a start or end vertex that is incident to the SESE
 *
 * the alternative side faces the SESE
 * if the incident side is the left side then the vertex is in forward orientation and vice versa
 * Being an SESE we expect that all vertices connected to one side are in the SESE
 * and all vertices connected to the other side are not in the SESE
 *
 * @param in_sese The set of vertices in the SESE
 * @param start_id The id of the vertex to get the orientation for
 */
//bd::VertexEnd
pgt::v_end_e get_boundary_incidence(const bd::VG& g,
                                    const std::set<std::size_t>& in_sese,
                                    std::size_t v_id, std::size_t alt_id) {
  std::string fn_name = std::format("[povu::genomics::]", __func__);

  std::size_t v_idx = g.v_id_to_idx(v_id);
  const bd::Vertex& v = g.get_vertex_by_idx(v_idx);

  auto foo = [&](std::size_t e_idx) ->bool {
    auto [side, alt_v_idx] = g.get_edge(e_idx).get_other_vtx(v_idx);
    return in_sese.count(g.v_idx_to_id(alt_v_idx)) || g.v_idx_to_id(alt_v_idx) == alt_id ;
  };

  bool allLeftInSet = std::any_of(v.get_edges_l().begin(), v.get_edges_l().end(), foo);

  bool allRightNotInSet = allLeftInSet ? !std::none_of(v.get_edges_r().begin(), v.get_edges_r().end(), foo)
                                       : std::any_of(v.get_edges_r().begin(), v.get_edges_r().end(), foo);
  if (!(allLeftInSet ^ allRightNotInSet)) {
    throw std::runtime_error(std::format("{} {} {}", v_idx, allLeftInSet, allRightNotInSet));
  }

  pgt::v_end_e o = allLeftInSet ? pgt::v_end_e::l : pgt::v_end_e::r;

  return o;
}


std::pair<bd::id_n_orientation_t, bd::id_n_orientation_t>
foo(const bd::VG& g, const pgt::canonical_sese& sese) {
  std::string fn_name = std::format("[povu::genomics::{}]", __func__);
  if (false) { std::cerr << fn_name << "\n"; }

  auto [start_id, stop_id, in_sese] = sese;

  /*
    assume e is the incoming vertex end
   */
  auto v_end_to_orientation = [&](pgt::v_end_e e) -> bd::orientation_t {
    if (e == pgt::v_end_e::l ) {
      return bd::orientation_t::forward;
    }

    return bd::orientation_t::reverse;
  };

  /*
    Determine start side and orientation
   */
  pgt::v_end_e start_end;
  try {
    start_end = get_boundary_incidence(g, in_sese, start_id, stop_id);
  }
  catch (std::exception& e) {
    std::cerr << std::format("{} WARN: SESE ({} {}) start boundary: {}\n", fn_name, start_id, stop_id, e.what());
    return{};
  }

  pgt::v_end_e stop_end;
  try {
    stop_end = get_boundary_incidence(g, in_sese, stop_id, start_id);
  }
  catch (std::exception& e) {
    std::cerr << std::format("{} WARN: SESE ({} {}) stop boundary: {}\n", fn_name, start_id, stop_id, e.what());
    return{};
  }


  //std::cerr << start_orientation << " " << stop_orientation << "\n";

    // TODO: remove when sort is implemented
  // if SESE is flipped then we need to reverse the start and stop
  if (start_id > stop_id && start_end == pgt::v_end_e::l && stop_end == pgt::v_end_e::r) {
    std::swap(start_id, stop_id);
    std::swap(start_end, stop_end);
  }

  bd::orientation_t start_o =
    start_end == pgt::v_end_e::r ? bd::orientation_t::forward : bd::orientation_t::reverse;

  bd::orientation_t stop_o =
    stop_end == pgt::v_end_e::l ? bd::orientation_t::forward : bd::orientation_t::reverse;

  bd::id_n_orientation_t entry { start_id, start_o };
  bd::id_n_orientation_t exit { stop_id, stop_o };

  return { entry, exit };
}

/**
 * @brief loop over bubbles and the paths
 */
void find_bubble_paths(const std::vector<pgt::flubble>& canonical_flubbles,
                       const bd::VG& bd_vg,
                       std::vector<std::vector<pgt::walk>>& all_paths) {
  std::string fn_name { std::format("[povu::genomics::{}]" , __func__) };

  for (std::size_t i{}; i < canonical_flubbles.size(); ++i) {

    const auto& [entry, exit] = canonical_flubbles[i];
    std::vector<pgt::walk> paths = bd_vg.get_paths(entry, exit);
    if (paths.size() < 2) {
      std::cerr << std::format("{} WARN: Bubble {} {} has {} paths\n", fn_name, entry.as_str(), exit.as_str(), paths.size());
    }

    all_paths.push_back(paths);
  }
}
/**
  * @brief reference paths that are in the graph & in the app config
  *
  * basically an intersection of the app config (user requested refs) and
component
  *
  * @param bd_vg The graph
  * @param app_config The application configuration
  * @return The relevant reference paths
  */
std::set<pt::id_t> find_relevant_refs(const bd::VG& bd_vg, const core::config& app_config) {
  std::set<pt::id_t> relevant_refs {};
  for (pgt::path_t& hap : bd_vg.get_haplotypes()) {
    auto [name, hap_id, _] = hap;
    std::vector<std::string> const& ref_paths = app_config.get_reference_paths();
    if (std::find(ref_paths.begin(), ref_paths.end(), name) != ref_paths.end()) {
      relevant_refs.insert(hap.id);
    }
  }

  return relevant_refs;
}

/**
 *
 *
 */
void call_variants(const std::vector<pgt::flubble>& canonical_flubbles,
                   const bd::VG& bd_vg,
                   const core::config& app_config) {
  std::string fn_name { std::format("[povu::genomics::{}]" , __func__) };

  std::vector<std::vector<pgt::walk>> all_paths;
  //std::vector<std::pair<pgt::id_n_orientation_t, pgt::id_n_orientation_t>> bubble_boundaries;
  find_bubble_paths(canonical_flubbles, bd_vg, all_paths);

  std::set<pt::id_t> ref_ids { find_relevant_refs(bd_vg, app_config) };

  std::vector<Bubble> c_bubs { // canonical bubbles
    find_haplotypes(bd_vg, all_paths, canonical_flubbles, ref_ids)
  };

  std::vector<std::tuple<Bubble, std::map<pt::id_t, std::vector<pu::exp_cnv>>, std::map<std::size_t, std::set<std::size_t>>>>
      res { povu::untangle::untangle(bd_vg, c_bubs, app_config) };


  std::map<std::size_t, std::vector<vcf::vcf_record>> vcf_records {
    genomics::vcf::gen_vcf_records(bd_vg, res, app_config)
  };


  //return;
  io::vcf::write_vcfs(vcf_records, bd_vg, app_config);
}

} // namespace genomics
