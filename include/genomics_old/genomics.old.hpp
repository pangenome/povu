#ifndef GENOMICS_HPP
#define GENOMICS_HPP

#include <cstddef>
#include <set>
#include <string>
#include <vector>

#include "../../src/cli/app.hpp"
#include "../common/types.hpp"
#include "../common/types.hpp"
#include "../graph/bidirected.hpp"


namespace povu::variants {


}

namespace povu::genomics {
namespace bd = povu::bidirected;
namespace pt = povu::types;
namespace pgt = povu::graph_types;

// TODO: make use of this or delete
enum class variant_type {
  SNP,
  NP,    // nucleotide polymorphism (single or multiple nucleotide)
  INDEL, // insertion or deletion, specifics unknown at this point
  DEL,   // deletion
  INS,   // insertion
  SUB,   // substitution
  INV,
  DUP,
  CNV,
  BND,
  INVALID,
  UNDEFINED,
};

// Overload the << operator
std::ostream& operator<<(std::ostream& os, variant_type vt);

// TODO which version of VCF is best?
enum class output_format {
    VCF, //  currently outputs v4.2
    PAF, // not yet supported
};

class Path {
  pgt::walk walk_; // a path is a sequence of vertices and their orientation

  // haplotypes associated with the path
  // key is the haplotype id and value is the range of the haplotype in the path
  std::map<std::size_t, pt::Stride> haps_; // key is the haplotype id and value are the ranges of the haplotype in walk_ above

  // the refs that are contained in the walk and are valid for variant calling as ref paths
  std::vector<std::size_t> refs_;

private:
  std::string to_g_id (const bd::VG& bd_vg, pgt::id_n_orientation_t const& x) const {
    return std::format("{}{}", or_to_str(x.orientation), bd_vg.v_idx_to_id(x.v_idx));
  }

public:
  // --------------
  // constructor(s)
  // --------------
  Path() = default;
  Path(std::vector<pgt::id_or> path, std::map<std::size_t, pt::Stride> haps)
    : walk_{path}, haps_{haps} {}

  // ---------
  // getter(s)
  // ---------
  const std::vector<pgt::id_n_orientation_t>& get_walk() const {
    return this->walk_;
  }

  const std::map<std::size_t, pt::Stride>& get_haps() const {
    return this->haps_;
  }

  const pt::Stride& get_hap_span(std::size_t hap_id) const {
    return this->haps_.at(hap_id);
  }

  std::string as_DNA_str(const bd::VG& bd_vg, std::size_t w_idx) const {
    // skip the first and last vertices in the path
    auto [idx, o] = this->walk_[w_idx];
    bd::Vertex const& v = bd_vg.get_vertex_by_idx(idx);
    return o == pgt::orientation_t::forward ? v.get_label() : v.get_rc_label();
  }

  std::string as_DNA_str(const bd::VG& bd_vg) const {
    std::string dna_str;
    // skip the first and last vertices in the path
    for (std::size_t i {1}; i < this->walk_.size() - 1; ++i) {
      dna_str += this->as_DNA_str(bd_vg, i);
    }
    return dna_str;
  }

  std::string as_DNA_str(const bd::VG& bd_vg,
                         pt::Stride sp,
                         const std::vector<variant_type> &variant_cats, const std::string &record_id) const {
    std::string dna_str;

    bool has_ins{std::find(variant_cats.begin(), variant_cats.end(), genomics::variant_type::INS) != variant_cats.end()};



    if (has_ins) {
      std::string s_str { this->as_DNA_str(bd_vg, 0) };
      dna_str += s_str.back();
      for (std::size_t i { 1 }; i < sp.length - 1; ++i) {
        dna_str += this->as_DNA_str(bd_vg, i);
      }


      return dna_str;
    }

    std::size_t i { sp.start == 0 ? 1 : sp.start }; // start idx
    std::size_t last_idx { sp.start + sp.length == this->walk_.size() ? sp.length - 1 : sp.length };
    for (; i < last_idx; ++i) {
      dna_str += this->as_DNA_str(bd_vg, i);
    }




    return dna_str;
  }

  std::string as_str(const bd::VG& bd_vg) const {
    std::string at_str;
    for ( const pgt::id_n_orientation_t& x : this->walk_) {
      at_str += to_g_id(bd_vg, x);
    }
    return at_str;
  }

  std::string as_str(const bd::VG& bd_vg, pt::Stride sp) const {
    std::string at_str;
    for (std::size_t i { sp.start }; i < sp.length; ++i) {
      at_str += to_g_id(bd_vg, this->walk_[i]);
    }
    return at_str;
  }
};

/**
  * @brief A class to hold the paths within a single bubble
  *
  * A bubble has multiple paths
  */
class Bubble {
  pgt::id_n_orientation_t start_;
  pgt::id_n_orientation_t end_;
  std::vector<Path> bubble_paths_;
  // set of haplotype ids in the bubble
  std::set<std::size_t> haps_;
  // a map of haplotype id to the set of Paths in bubble_paths_ that contain the haplotype
  std::map<std::size_t, std::set<std::size_t>> hap_paths_;
  // haplotype id to the path index in bubble_paths_ that contains subpath
  std::map<std::size_t, std::size_t> longest_hap_walk_;
  // the haplotypes/references contained in a walk
  // key is walk id and value is a vector of ref ids
  //std::map<std::size_t, <std::size_t>> walk_refs_;

public:
  // --------------
  // constructor(s)
  // --------------
  Bubble() = default;
  Bubble(pgt::id_n_orientation_t start, pgt::id_n_orientation_t end)
    : start_{start}, end_{end} {}

  Bubble(const std::vector<Path> &bubble_paths) : bubble_paths_{bubble_paths} {}

  Bubble(pgt::id_n_orientation_t start,
         pgt::id_n_orientation_t end,
         const std::vector<Path>& bubble_paths)
    : start_{start}, end_{end}, bubble_paths_{bubble_paths} {

    for (std::size_t i {}; i < bubble_paths.size(); ++i) {
      for (auto const& [hap_id, _] : bubble_paths[i].get_haps()) {
        this->haps_.insert(hap_id);
        this->hap_paths_[hap_id].insert(i);
      }
    }

    for (auto const& hap_id : this->haps_) {
      std::size_t longest {};
      for (std::size_t walk_idx : this->hap_paths_[hap_id]) {
        if (this->bubble_paths_[walk_idx].get_hap_span(hap_id).length > longest) {
          longest_hap_walk_[hap_id] = walk_idx;
          longest = this->bubble_paths_[walk_idx].get_hap_span(hap_id).length;
        }
      }
    }
  }


  // ---------
  // getter(s)
  // ---------
  std::size_t size() const { return this->bubble_paths_.size(); }
  pgt::id_n_orientation_t start() const { return this->start_; }
  pgt::id_n_orientation_t end() const { return this->end_; }
  // get the Path at the given index in bubble_paths_
  const Path& get_path(std::size_t walk_idx) const { return this->bubble_paths_[walk_idx]; }
  const std::vector<Path>& get_paths() const { return this->bubble_paths_; }
  const std::set<std::size_t>& get_hap_ids() const { return this->haps_; }
  const std::set<std::size_t> &get_walk_idxs(std::size_t hap_id) const {
    if (this->hap_paths_.find(hap_id) == this->hap_paths_.end()) {
      throw std::runtime_error("Haplotype id not found in bubble");
    }

    return this->hap_paths_.at(hap_id);
  }
  std::size_t get_longest_walk_idx(std::size_t hap_id) const { return this->longest_hap_walk_.at(hap_id); }
};

void call_variants(const std::vector<pgt::flubble>& canonical_flubbles,
                   const bd::VG& bd_vg,
                   const core::config& app_config);
} // namespace genomics


namespace povu::untangle {
namespace bd = bidirected;
namespace pg = genomics;
namespace pt = povu::types;
namespace pgt = povu::graph_types;
using povu::graph_types::walk;


struct ref_meta {
  const pg::Path& path;
  pt::idx_t w_idx; // walk index in the bubble
  const pt::span& w_span; // walk span in path
};

// (nested) copy number variant
struct cnv {
  std::vector<pt::idx_t> positions;
  pt::idx_t w_idx; // walk index in the bubble
  const pt::span& w_span; // walk span in path // TODO: should this be ref instead?
};


// expanded cnv
// position, walk index, span
typedef std::tuple<pt::idx_t, pt::idx_t, pt::span> exp_cnv;

std::vector<std::tuple<pg::Bubble, std::map<pt::id_t, std::vector<exp_cnv>>, std::map<std::size_t, std::set<std::size_t>>>>
untangle(const bd::VG& bd_vg, const std::vector<genomics::Bubble>& c_bubs, const core::config& app_config);
} // namespace povu::untangle


namespace povu::genomics::vcf {
struct vcf_record {
  std::string chrom;
  std::size_t pos;
  std::string id;
  std::string ref;
  std::vector<std::string> alt;
  //std::string qual;
  //std::string filter;
  //std::string info;
  std::string format;
};

std::map<std::size_t, std::vector<vcf::vcf_record>> gen_vcf_records(
  const bd::VG& bd_vg,
  const std::vector<std::tuple<Bubble, std::map<pt::id_t, std::vector<povu::untangle::exp_cnv>>, std::map<std::size_t, std::set<std::size_t>>>>& c_bubs,
  const core::config& app_config);
} // namespace genomics::vcf
#endif
