#ifndef POVU_TYPES_GRAPH_HPP
#define POVU_TYPES_GRAPH_HPP

#include <map>
#include <optional>
#include <set>
#include <string>
#include <vector>

#include "./compat.hpp"
#include "./core.hpp"
#include "../utils.hpp"

namespace povu::types::graph {
namespace pt = povu::types;
namespace pu = povu::utils;

// should this be renamed to clr_e or color_e?
enum class color {
  gray,
  black
};
typedef color color_e;
std::ostream& operator<<(std::ostream& os, const color_e& c);

// Eq class and node id
struct eq_n_id_t {
  std::size_t eq_class;
  std::size_t v_id;
};

struct id_n_cls {
  std::size_t id;
  std::size_t cls;
};


/**
  * l (left) 5' or +
  * r (right) 3' or -
  */
// TODO: replace with struct or class to allow methods like complement
enum class v_end_e {
  l,
  r
};
std::ostream& operator<<(std::ostream& os, const v_end_e& vt);
constexpr v_end_e complement(v_end_e s) {
  return s == v_end_e::l ? v_end_e::r : v_end_e::l;
};

/**
 * l (left) 5' or +
 * r (right) 3' or -
 */
enum class v_type_e {
    l,
    r,
    dummy
};
std::ostream& operator<<(std::ostream& os, const v_type_e& vt);

// TODO: remove
// Merge path_t and biedged PathInfo into one namespace
struct path_t {
  std::string name; // name as pertains the GFA file
  std::size_t id; // numerical id to be associated with handle
  bool is_circular; // is the path circular?
};

struct side_n_id_t {
  v_end_e v_end;
  pt::id_t v_idx; // TODO [B] rename v_idx to v_id

  // -------
  // methods
  // -------
  friend bool operator<(const side_n_id_t& lhs, const side_n_id_t& rhs);
  side_n_id_t complement() const;
};
std::ostream& operator<<(std::ostream& os, const side_n_id_t& x);
typedef side_n_id_t side_n_idx_t; // to use when the id is an index

struct canonical_sese {
  std::size_t start;
  std::size_t end;
  std::set<std::size_t> in_sese; // set of sese ids that contain this sese excluding the start and end
};

enum class or_e {
  forward,
  reverse
};
std::ostream& operator<<(std::ostream& os, const or_e& o);
std::string or_to_str (or_e o);

// TODO: move to povu::types::variation
struct id_or_t {
  pt::id_t v_id; // TODO change type and name to id to pt::id_t
  or_e orientation;


  std::string as_str() const {
    return pv_cmp::format("{}{}", or_to_str(this->orientation) , this->v_id);
  }
};

std::ostream& operator<<(std::ostream& os, const id_or_t& x);
bool operator!=(const id_or_t & lhs, const id_or_t& rhs);
bool operator==(const id_or_t & lhs, const id_or_t& rhs);
bool operator<(const id_or_t& lhs, const id_or_t& rhs);

// a walk is a sequence of vertices also a path
typedef std::vector<id_or_t> walk;


struct flubble {
  id_or_t start_;
  id_or_t end_;

/* constructor(s) */
flubble(id_or_t start, id_or_t end) : start_(start), end_(end) {}

flubble(const std::string& s) {

  // split s based on > and < signs and using s.substr
  // s is in the form of >1<2 or >1>2 or <1<2 or <1>2

  // find the first > or <
  auto first = s.find_first_of("><");
  auto last = s.find_last_of("><");

  // substring based on first and last occurences and store them as size_t

  this->start_.v_id = std::stoull(s.substr(first + 1, last - first - 1));
  this->start_.orientation = s[first] == '>' ? or_e::forward : or_e::reverse;


  this->end_.v_id = std::stoull(s.substr(last + 1, s.size() - last - 1));
  this->end_.orientation = s[last] == '>' ? or_e::forward : or_e::reverse;
}

/* other(s) */
std::string as_str() const {
  std::string s;
  s += this->start_.as_str();
  s += this->end_.as_str();
  return s;
  //return pv_cmp::format("{}{}", this->start_.as_str(), this->end_.as_str());
}
};

typedef  flubble flubble_t ;

bool operator<(const flubble_t &lhs, const flubble_t &rhs);

struct pan_sn {
  std::string sample_name_; // not unique
  pt::id_t haplotype_id_;   // is unique
  std::string contig_name_; // is unique contig or scaffold name
};

std::optional<pan_sn> label_to_pan_sn(const std::string &label, char delim);

// applies to a single line in the P line of a GFA file

class Ref {
  pt::id_t ref_id_;
  std::string label_; // the entire label in the P line

  pt::idx_t len_; // length of the reference/contig/scaffold

  // PanSN spec
  bool has_pansn_data_ = false; // true if the label is in the PanSN format
  pan_sn pansn_data_; // holds the PanSN data if has_pansn_data_ is true
  // std::string sample_name_; // not unique
  // pt::id_t haplotype_id_; // is unique
  // std::string contig_name_; // is unique contig or scaffold name

public:

  // --------------
  // constructor(s)
  // --------------

  static Ref parse_pansn_label(pt::id_t ref_id, const std::string &label, char delim) {
    Ref info;
    info.ref_id_ = ref_id;
    info.label_ = label;

    std::optional<pan_sn> maybe_pansn = label_to_pan_sn(label, delim);
    if (maybe_pansn.has_value()) {
      info.has_pansn_data_ = true;
      info.pansn_data_ = maybe_pansn.value();
    }
    else {
      info.has_pansn_data_ = false;
    }

    return info;
  }

  // ---------
  // getter(s)
  // ---------

  pt::id_t id() const {
    return this->ref_id_;
  }

  bool has_pansn_data() const {
    return this->has_pansn_data_;
  }

  const std::string &get_label() const {
    return this->label_;
  }

  const std::string &get_sample_name() const {
    return this->pansn_data_.sample_name_;
  }

  std::string get_col_name() const {
    if (this->has_pansn_data_) {
      return this->pansn_data_.sample_name_;
    }
    else {
      return this->label_;
    }
  }

  pt::id_t get_haplotype_id() const { return this->pansn_data_.haplotype_id_; }

  const std::string &get_contig_name() const {
    return this->pansn_data_.contig_name_;
  }

  pt::idx_t get_length() const {
    return this->len_;
  }

  // ---------
  // setter(s)
  // ---------

  void set_length(pt::idx_t len) {
    this->len_ = len;
  }
};

class Refs {
  std::vector<Ref> refs_;

  // the ref id is the index in the vector and the value is the set of ref ids
  // which share the sample name with the ref id at that index
  std::vector<std::set<pt::id_t>> ref_ids_;

  // the string is the sample name and the values are the ref ids which share
  // a sample name
  std::map<std::string, std::set<pt::id_t>> sample_to_ref_ids_;

  pu::TwoWayMap<std::string, pt::id_t> label_and_ref_id_;

public:

  // --------------
  // constructor(s)
  // --------------
  Refs() = default;

  // ---------
  // getter(s)
  // ---------

  const Ref &get_ref(pt::id_t ref_id) const {
    return this->refs_.at(ref_id);
  }

  Ref &get_ref_mut(pt::id_t ref_id) {
    return this->refs_.at(ref_id);
  }

  const std::string &get_ref_label(pt::id_t ref_id) const {
    return this->refs_.at(ref_id).get_label();
  }

  pt::id_t get_ref_id(const std::string &ref_label) const {
    return this->label_and_ref_id_.get_value(ref_label);
  }

  pt::id_t ref_id_count() const {
    return this->refs_.size();
  }

  /**
   * Get the set of ref ids that share the same sample name as the ref id
   */
  const std::set<pt::id_t> &get_shared_samples(pt::id_t ref_id) const {
    return this->ref_ids_.at(ref_id);
  }

  // ---------
  // setter(s)
  // ---------

  pt::id_t add_ref(const std::string &label, char delim) {
    pt::id_t ref_id = this->refs_.size();
    Ref ref = Ref::parse_pansn_label(ref_id, label, delim);
    this->refs_.push_back(ref);
    this->label_and_ref_id_.insert(label, ref_id);

    const std::string &sample_name = ref.get_col_name();

    this->sample_to_ref_ids_.try_emplace(sample_name, std::set<pt::id_t>{});
    this->sample_to_ref_ids_[sample_name].insert(ref_id);

    // TODO [B] use already pre computed ref count
    this->ref_ids_.resize(this->refs_.size());

    // update for all refs
    for (pt::id_t r_id : this->sample_to_ref_ids_[sample_name]) {
      this->ref_ids_[r_id] = this->sample_to_ref_ids_[sample_name];
    }

    return ref_id;
  }
};

} // namespace povu::graph_types

#endif
