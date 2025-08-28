#ifndef POVU_REFS_HPP
#define POVU_REFS_HPP

#include <optional>
#include <string>
#include <string_view>
#include <variant>
#include <set>
#include <charconv>

#include "../common/core.hpp"
#include "../common/constants.hpp"

namespace povu::refs {
inline constexpr std::string_view MODULE = "povu::refs";

enum class ref_format_e {
  PANSN = 1,
  UNDEFINED = 0,
};

const uint8_t PANSN_ELEMENT_COUNT = 3;
const char    EXPECTED_PANSN_DELIM = '#';

struct pan_sn {
private:
  std::string sample_name_;
  pt::id_t    haplotype_id_;
  std::string contig_name_;

public:
  pan_sn(std::string_view sample_name, pt::id_t hap_id,
         std::string_view contig_name)
      : sample_name_(sample_name), haplotype_id_(hap_id),
        contig_name_(contig_name) {}

  const std::string &sample() const { return this->sample_name_; }
  pt::id_t hap() const { return this->haplotype_id_; }
  const std::string &contig() const { return this->contig_name_; }
};

// The sum type: either a parsed format, or the original string when unknown.
using ref_variant = std::variant <
  std::string, // unknown / passthrough
  pan_sn
  // add future types here
>;

class Ref {
  ref_format_e format_ = ref_format_e::UNDEFINED;
  // TODO find a better way and name to store the data_
  ref_variant  data_; // holds either the parsed data or the original string
  pt::id_t ref_id_ = pc::INVALID_ID; // the ref id assigned in Refs
  pt::idx_t len_;                    // length of the contig/scaffold

  // ---------------
  // private methods
  // ---------------
  /**
   * Split a string_view into exactly 3 parts based on a delimiter.
   * Returns false if there are not exactly 3 parts.
   *
   * @param s The input string_view to split.
   * @param delim The delimiter character.
   * @param out An array of 3 string_views to hold the split parts.
   */
  [[nodiscard]] static bool
  split3_tag(std::string_view s, char delim,
             std::string_view (&out)[PANSN_ELEMENT_COUNT]) {
    size_t pos = 0, idx = 0;
    while (idx < PANSN_ELEMENT_COUNT - 1) {
      size_t next = s.find(delim, pos);
      if (next == std::string_view::npos) {
        return false; // not enough parts
      }
      out[idx++] = s.substr(pos, next - pos);
      pos = next + 1;
    }
    // last part is the remainder; disallow extra delimiters
    if (s.find(delim, pos) != std::string_view::npos) {
      return false; // too many parts
    }
    out[idx] = s.substr(pos);

    return true;
  }

  [[nodiscard]] static std::optional<pan_sn>
  parse_pan_sn(std::string_view tag, char delim) {
    std::string_view p[PANSN_ELEMENT_COUNT];
    if (!split3_tag(tag, delim, p)) {
      return std::nullopt;
    }
    if (p[0].empty() || p[1].empty() || p[2].empty()) {
      return std::nullopt;
    }

    // parse haplotype with from_chars (no throw)
    unsigned long long tmp = 0;
    auto res = std::from_chars(p[1].data(), p[1].data() + p[1].size(), tmp, 10);
    if (res.ec != std::errc{} || res.ptr != p[1].data() + p[1].size())
      return std::nullopt;

    // range-check for pt::id_t
    using id_lim = std::numeric_limits<std::make_unsigned_t<pt::id_t>>;
    if (tmp > static_cast<unsigned long long>(id_lim::max()))
      return std::nullopt;

    pan_sn out(p[0], static_cast<pt::id_t>(tmp), p[2]);
    // pan_sn out();
    // out.sample_name_.assign(p[0]);
    // out.haplotype_id_ = static_cast<pt::id_t>(tmp);
    // out.contig_name_.assign(p[2]);
    return out;
  }

  [[nodiscard]] static std::pair<ref_format_e, ref_variant>
  determine_format(std::string_view tag, char delim) {
  if (auto psn = parse_pan_sn(tag, delim)) {
    return {ref_format_e::PANSN, std::move(*psn)};
  }

  // unknown → keep original
  return {ref_format_e::UNDEFINED, std::string(tag)};
  }

  // make default constructor private
  Ref() = default;

public:
  // --------------
  // constructor(s)
  // --------------
  static Ref parse(pt::id_t ref_id, std::string_view tag, char delim) {
    Ref r;
    r.ref_id_ = ref_id;
    auto [f, d] = determine_format(tag, delim);
    r.data_ = d;
    r.format_ = f;
    return r;
  }

  // ---------
  // getter(s)
  // ---------
  std::string tag() const {
    if (std::holds_alternative<std::string>(this->data_)) {
      return std::get<std::string>(this->data_);
    }
    else if (std::holds_alternative<pan_sn>(this->data_)) {
      const pan_sn &psn = std::get<pan_sn>(this->data_);
      return psn.sample() + EXPECTED_PANSN_DELIM + std::to_string(psn.hap()) +
             EXPECTED_PANSN_DELIM + psn.contig();
    }

    // should not get here
    throw std::runtime_error("Could not generate tag for Ref. Invalid ref_variant state.");
  }

  ref_format_e get_format() const {
    return this->format_;
  }

  const std::string &get_sample_name() const {
    switch (this->format_) {
    case ref_format_e::PANSN:
      return std::get<pan_sn>(this->data_).sample();
    case ref_format_e::UNDEFINED:
      return std::get<std::string>(this->data_);
    default:
      throw std::runtime_error("Unknown ref format");
    }
  }

  pt::id_t get_length() const {
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

  Refs() = default;

public:
  // --------------
  // constructor(s)
  // --------------
  Refs (pt::idx_t initial_capacity) {
    this->refs_.reserve(initial_capacity);
  }

  // ---------
  // getter(s)
  // ---------

  const Ref &get_ref(pt::id_t ref_id) const {
    return this->refs_.at(ref_id);
  }

  Ref &get_ref_mut(pt::id_t ref_id) {
    return this->refs_[ref_id];
  }

  const std::string &get_sample_name(pt::id_t ref_id) const {
    return this->refs_[ref_id].get_sample_name();
  }

  std::optional<pt::id_t> get_ref_id(std::string_view tag) const {
    for (pt::id_t ref_id{}; ref_id < this->refs_.size(); ref_id++) {
      const Ref &r = this->refs_[ref_id];
      if (r.tag() == tag) {
        return ref_id;
      }
    }

    return std::nullopt;
  }



  // the sample name could also be referred to as a prefix
  std::set<pt::id_t> get_refs_in_sample(std::string_view sample_name) const {

    auto is_prefix = [](std::string_view pre, std::string_view str) -> bool {
      return str.compare(0, pre.size(), pre) == 0;
    };

    std::set<pt::id_t> in_sample;
    for (pt::id_t ref_id{}; ref_id < this->refs_.size(); ref_id++) {
      const Ref &r = this->refs_[ref_id];

      switch (r.get_format()) {
      case ref_format_e::PANSN:
        if (r.get_sample_name() == sample_name) { in_sample.insert(ref_id);}
        break;
      case ref_format_e::UNDEFINED:
        if (is_prefix(sample_name, r.get_sample_name())) { in_sample.insert(ref_id); }
        break;
      default:
        ;
      }
    }

    return in_sample;
  }

  std::set<pt::id_t> get_shared_samples(pt::id_t ref_id) const {
    const Ref &r = this->refs_[ref_id];
    return this->get_refs_in_sample(r.get_sample_name());
  }

  // ---------
  // setter(s)
  // ---------

  pt::id_t add_ref(const std::string &label, char delim) {
    pt::id_t ref_id = static_cast<pt::id_t>(this->refs_.size());
    Ref ref = Ref::parse(ref_id, label, delim);
    this->refs_.emplace_back(std::move(ref));
    return ref_id;
  }

  pt::id_t ref_count() const {
    return static_cast<pt::id_t>(this->refs_.size());
  }
};

}; // namespace povu::refs::pansn

// NOLINTNEXTLINE(misc-unused-alias-decls)
namespace pr = povu::refs;

#endif // POVU_REFS_HPP
