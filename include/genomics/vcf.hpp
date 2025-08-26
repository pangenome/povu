#ifndef POVU_GENOMICS_VCF_HPP
#define POVU_GENOMICS_VCF_HPP

#include <algorithm>
#include <map>
#include <string>
#include <utility>
#include <vector>


#include "../common/compat.hpp"
#include "../common/types/graph.hpp"
#include "../graph/bidirected.hpp"
#include "./allele.hpp"

namespace povu::genomics::vcf {
inline constexpr std::string_view MODULE = "povu::genomics::vcf";

namespace bd = povu::bidirected;
namespace pgt = povu::types::graph;
namespace pvst = povu::pvst;
namespace pga = povu::genomics::allele;

enum class var_type_e {
  del, // deletion
  ins, // insertion
  sub, // substitution
  und  // undetermined
};

// add operator << for var_type_e
inline std::ostream &operator<<(std::ostream &os, var_type_e vt) {
  switch (vt) {
  case var_type_e::del:
    os << "DEL";
    break;
  case var_type_e::ins:
    os << "INS";
    break;
  case var_type_e::sub:
    os << "SUB";
    break;
  case var_type_e::und:
    os << "UND";
    break;
  }
  return os;
}

// 1) A helper to turn enum → string_view:
constexpr std::string_view to_string_view(var_type_e vt) noexcept {
  switch (vt) {
  case var_type_e::del:
    return "DEL";
  case var_type_e::ins:
    return "INS";
  case var_type_e::sub:
    return "SUB";
  case var_type_e::und:
    return "UND";
  }
  // optional: handle out-of-range
  return "??";
}

class VcfRec {
  pt::id_t ref_id_; // chrom
  pt::idx_t pos_; // 1-based step idx
  std::string id_; // start and end of a RoV e.g >1>4
  pga::AW ref_at_; //
  std::vector<pga::AW> alt_ats_;
  // std::string qual; // fixed at 60
  // std::string filter; // fixed as pass
  // std::string info; // computed later
  // std::string format; // GT for now
  pt::idx_t height_; // height of the pvst node in the tree
  var_type_e var_type_; // type of the variant, e.g. del, ins, sub, und
  bool is_tangled_ = false; // is true when tangling exists, i.e. when a walk traverses an RoV more than once

public:

  // --------------
  // constructor(s)
  // --------------

  VcfRec(pt::id_t ref_id, pt::idx_t pos, std::string id, pga::AW ref_at,
         std::vector<pga::AW> alt_ats, pt::idx_t height, var_type_e variant_type, bool is_tangled)
    : ref_id_(ref_id), pos_(pos), id_(id), ref_at_(ref_at), alt_ats_(alt_ats),
      height_(height), var_type_(variant_type), is_tangled_(is_tangled) {}

  // ---------
  // getter(s)
  // ---------
  pt::id_t get_ref_id() const { return this->ref_id_; }
  pt::idx_t get_pos() const { return this->pos_; }
  std::string get_id() const { return this->id_; }
  const pga::AW &get_ref_at() const { return this->ref_at_; }
  const std::vector<pga::AW> &get_alt_ats() const { return this->alt_ats_; }
  std::vector<pga::AW> &get_alt_ats_mut() { return this->alt_ats_; }
  pt::idx_t get_height() const { return this->height_; }
  var_type_e get_var_type() const { return this->var_type_; }
  bool is_tangled() const { return this->is_tangled_; }

  std::set<pt::id_t> get_ref_ids() const {
    std::set<pt::id_t> ref_ids;
    ref_ids.insert(this->ref_id_);
    for (const auto &aw : this->alt_ats_) {
      for (const auto &ref_id : aw.get_ref_ids()) {
        ref_ids.insert(ref_id);
      }
    }
    return ref_ids;
  }


  // ---------
  // setter(s)
  // ---------
  pt::idx_t append_alt_at(pga::AW &&alt_at) {
    this->alt_ats_.emplace_back(std::move(alt_at));
    return this->alt_ats_.size() - 1;
  }
};


struct genotype_data_t {
  std::map<pt::id_t, pt::idx_t> ref_id_to_col_idx;

  // the size is the number of columns in the genotype data
  // string contains the sample name or label
  std::vector<std::string> genotype_cols;
};

// TODO: [c] find a better name
class VcfRecIdx {
  std::map<pt::idx_t, std::vector<VcfRec>> vcf_recs_;
  genotype_data_t gd_;

public:

  // --------------
  // constructor(s)
  // --------------

  VcfRecIdx() : vcf_recs_() {}

  // ---------
  // getter(s)
  // ---------
  const std::map<pt::idx_t, std::vector<VcfRec>> &get_recs() const {
    return this->vcf_recs_;
  }

  const genotype_data_t &get_genotype_data() const {
    return this->gd_;
  }

  // ---------
  // setter(s)
  // ---------

  void add_rec(pt::idx_t ref_idx, VcfRec &&vcf_rec) {
    if (this->vcf_recs_.find(ref_idx) == this->vcf_recs_.end()) {
      this->vcf_recs_[ref_idx] = std::vector<VcfRec>{vcf_rec};
    } else {
      this->vcf_recs_[ref_idx].emplace_back(vcf_rec);
    }
  }

  void set_genotype_data(genotype_data_t &&gd) {
    this->gd_ = std::move(gd);
  }

};


VcfRecIdx gen_vcf_records(const bd::VG &g, const std::vector<pga::Exp> &ref_walks);

} // namespace povu::genomics::vcf

#endif // POVU_GENOMICS_VCF_HPP
