#ifndef POVU_TYPES_PVST_HPP
#define POVU_TYPES_PVST_HPP

#include "./core.hpp"
#include "./constants.hpp"
#include "./graph.hpp"
#include <algorithm>
#include <format>

/* === PVST pangenome variation structure tree === */

namespace povu::types::pvst {
namespace pc = povu::constants;
namespace pt = povu::types;
namespace pgt = povu::types::graph;

// short for VertexType
enum class vt_e {
  dummy,
  /* types of flubbles */
  flubble, // TODO: aka generic flubble, rename to generic_flubble
  tiny,        // for SNPs
  parallel,     //
  /* types of bubbles */
  slubble, // rename to concealed
  smothered,   // rename to smothered
  midi,
};

enum class sl_type_e {
  ai_trunk,
  ai_branch,
  zi_trunk,
  zi_branch
};

// type of the concealed vertex
enum class cn_type_e {
  g,
  s
};

struct bounds_t {
  pt::idx_t upper;
  pt::idx_t lower; // when invalid all leaves are the lower boundaries
};

/* an abstract class for vertices  */
class VertexBase {
  povu::types::id_t idx_; // idx of the vertex in the vst
  vt_e type_;

public:
  VertexBase(povu::types::id_t idx, vt_e type)
      : idx_(idx), type_(type) {}

  povu::types::id_t get_idx() const { return this->idx_; }
  vt_e get_type() const { return this->type_; }

  void set_idx(povu::types::id_t idx) { this->idx_ = idx; }
  void set_type(vt_e type) { this->type_ = type; }
  virtual std::string as_str() const = 0;

  // TODO: add get_bounds
};


class Dummy : public VertexBase {
public:
  Dummy() : VertexBase(pc::INVALID_IDX, vt_e::dummy) {}

  std::string as_str() const override { return "."; }
};


class Flubble : public VertexBase {
  pgt::id_or_t a_; // start
  pgt::id_or_t z_; // end

  pt::idx_t ai_; // a_i
  pt::idx_t zi_; // z_i

  pt::idx_t m_; // m
  pt::idx_t n_; // n

public:
  /*
   == constructor ==
   */

  Flubble(pgt::id_or_t a, pgt::id_or_t z, pt::idx_t ai, pt::idx_t zi)
    : VertexBase(pc::INVALID_IDX, vt_e::flubble), a_(a), z_(z), ai_(ai), zi_(zi) {}

  /*
   == getters ==
   */
  pgt::id_or_t get_a() const { return this->a_; }
  pgt::id_or_t get_z() const { return this->z_; }
  pt::idx_t get_ai() const { return this->ai_; }
  pt::idx_t get_zi() const { return this->zi_; }
  pt::idx_t get_m() const { return this->m_; }
  pt::idx_t get_n() const { return this->n_; }
  bounds_t get_bounds() const { return {this->get_ai(), this->get_zi()}; }

  /*
   == setters ==
   */
  void set_m(pt::idx_t m) { this->m_ = m; }
  void set_n(pt::idx_t n) { this->n_ = n; }

  /*
   == other(s) ==
   */
  std::string as_str() const override {
    return std::format("{}{}", this->a_.as_str(), this->z_.as_str());
  }
};


class Concealed : public VertexBase {
  pt::idx_t fl_idx;
  sl_type_e sl_type_; // type of the slubble (trunk or branch)
  pt::idx_t sl_st_idx_; // idx in the spanning tree for slubble

  // b for boundary
  pgt::id_or_t fl_b_; // a or z
  pgt::id_or_t cn_b_; // g or s

  bounds_t bounds_;

private:
  bool with_ai() const {
    return this->sl_type_ == sl_type_e::ai_trunk || this->sl_type_ == sl_type_e::ai_branch;
  }

  bool with_zi() const {
    return (this->sl_type_ == sl_type_e::zi_trunk || this->sl_type_ == sl_type_e::zi_branch);
  }

public:
  /*
   == constructor ==
   */
  Concealed(pgt::id_or_t fl_b, pgt::id_or_t cn_b, bounds_t bounds, pt::idx_t fl_idx,
            sl_type_e sl_type, pt::idx_t sl_st_idx)
    : VertexBase(pc::INVALID_IDX, vt_e::slubble), fl_idx(fl_idx),
      sl_type_(sl_type), sl_st_idx_(sl_st_idx), fl_b_(fl_b), cn_b_(cn_b), bounds_(bounds) {}

  /*
   == getters ==
   */
  pt::idx_t get_fl_idx() const { return this->fl_idx; }
  sl_type_e get_sl_type() const { return this->sl_type_; }
  pt::idx_t get_sl_st_idx() const { return this->sl_st_idx_; }
  pgt::id_or_t get_fl_b() const { return this->fl_b_; }
  pgt::id_or_t get_cn_b() const { return this->cn_b_; }
  bounds_t get_bounds() const { return this->bounds_; }

  /*
   == others ==
   */
  std::string as_str() const override {
    if (with_ai()) { // formed with a
      return std::format("{}{}", this->fl_b_.as_str(), this->cn_b_.as_str());
    }
    else  { // formed with z
      return std::format("{}{}", this->cn_b_.as_str(), this->fl_b_.as_str());
    }
  }
};


class Smothered : public VertexBase {
  pt::idx_t cn_idx; // idx of the concealed vertex
  pt::idx_t sm_st_idx; // idx in the spanning tree for smothered vertex

  // b for boundary
  pgt::id_or_t cn_b_; // g or s
  pgt::id_or_t sm_b_; // e or w

  // is true when cn_b_ is an ancestor of sm_b_
  bool cn_b_is_ans_;

  cn_type_e cn_type_; // type of the concealed vertex (g or s)

  bounds_t bounds_;

public:
  Smothered(pgt::id_or_t cn_b, pgt::id_or_t sm_b, pt::idx_t cn_idx,
            bool cn_b_is_ans, pt::idx_t sm_st_idx, cn_type_e sm_type, bounds_t bounds)
      : VertexBase(pc::INVALID_IDX, vt_e::smothered),
        cn_idx(cn_idx), sm_st_idx(sm_st_idx), cn_b_(cn_b), sm_b_(sm_b),
        cn_b_is_ans_(cn_b_is_ans), cn_type_(sm_type), bounds_(bounds) {}

  pt::idx_t get_cn_idx() const { return this->cn_idx; }
  pt::idx_t get_sm_st_idx() const { return this->sm_st_idx; }
  bounds_t get_bounds() const { return this->bounds_; }
  cn_type_e get_cn_type() const { return this->cn_type_; }
  bool is_cn_b_ancestor() const { return this->cn_b_is_ans_; }

  std::string as_str() const override {
    if (this->cn_type_ == cn_type_e::g) { // g
      if (this->cn_b_is_ans_) { // cn_b is ancestor of sm_b
        return std::format("{}{}", this->cn_b_.as_str(), this->sm_b_.as_str());
      }
      else { // sm_b is ancestor of cn_b
        return std::format("{}{}", this->sm_b_.as_str(), this->cn_b_.as_str());
      }
      //return std::format("{}{}", this->cn_b_.as_str(), this->sm_b_.as_str());
    }
    else { // s
      return std::format("{}{}", this->sm_b_.as_str(), this->cn_b_.as_str());
    }
  }
};


class MidiBubble : public VertexBase {
  pt::idx_t g_cn_idx_;
  pt::idx_t s_cn_idx_;
  pgt::id_or_t g_;
  pgt::id_or_t s_;

public:
  MidiBubble(pt::idx_t g_cn_idx, pgt::id_or_t g, pt::idx_t s_cn_idx,
             pgt::id_or_t s)
      : VertexBase(pc::INVALID_IDX, vt_e::midi), g_cn_idx_(g_cn_idx),
        s_cn_idx_(s_cn_idx), g_(g), s_(s) {}

  bounds_t get_bounds() const {
    return bounds_t { std::min(g_cn_idx_, s_cn_idx_), std::max(g_cn_idx_, s_cn_idx_) };
  }

  std::string as_str() const {
    return std::format("{}{}", this->g_.as_str(), this->s_.as_str());
  }
};


class Vertex {
  // base
  povu::types::id_t idx_; // idx of the vertex in the vst
  vt_e type_;

  // flubble
  pgt::id_or_t a_;
  pgt::id_or_t z_;
  pt::idx_t ai_;
  pt::idx_t zi_;

  // only applies when is slubble
  //fl_vtx_type_e fl_vtx_type_; // type of the flubble vertex (ai or zi)
  sl_type_e fl_type_; // type of the flubble (trunk or branch)
  //pt::idx_t fl_st_idx_; // idx in the spanning tree for flubble
  pt::idx_t sl_st_idx_; // idx in the spanning tree for slubble

  pt::idx_t sm_st_idx_; // idx in the spanning tree for smothered vertex

private:
  // constructor for a concealed
  Vertex(pgt::id_or_t start, pgt::id_or_t end, pt::idx_t sl_st_idx, sl_type_e t, Vertex fl)
      : idx_(pc::INVALID_IDX), type_(vt_e::slubble), a_(start), z_(end),
        ai_(fl.get_ai()), zi_(fl.get_zi()), fl_type_(t), sl_st_idx_(sl_st_idx) {}

  // constructor for a smothered vertex
  Vertex(pgt::id_or_t start, pgt::id_or_t end, pt::idx_t sm_st_idx, Vertex fl)
    : type_(vt_e::smothered), a_(start), z_(end), sl_st_idx_(fl.get_sl_st_idx()),
      sm_st_idx_(sm_st_idx)  {}

  // constructor for a flubble
  Vertex(pgt::id_or_t a, pgt::id_or_t z, pt::idx_t ai, pt::idx_t zi)
    : idx_(pc::INVALID_IDX), type_(vt_e::flubble), a_(a), z_(z), ai_(ai), zi_(zi) {}

  // constructor for a dummy vertex
  Vertex()
      : idx_(pc::INVALID_IDX), type_(vt_e::dummy), a_(pgt::id_or_t{pc::INVALID_ID, pgt::or_e::forward}),
        z_(pgt::id_or_t{pc::INVALID_ID, pgt::or_e::forward}),
        ai_(pc::INVALID_IDX), zi_(pc::INVALID_IDX),
        sl_st_idx_(pc::INVALID_IDX) {}

public:
  // --------------
  // constructor(s)
  // --------------
  static Vertex make_flubble(pgt::id_or_t a, pgt::id_or_t z, pt::idx_t ai, pt::idx_t zi) {
    return Vertex(a, z, ai, zi);
  }

  static Vertex make_slubble(pgt::id_or_t start, pgt::id_or_t end,
                             pt::idx_t sl_st_idx, sl_type_e sl_t, Vertex fl) {
    return Vertex(start, end, sl_st_idx, sl_t, fl);
  }

  static Vertex make_smothered(pgt::id_or_t start, pgt::id_or_t end,
                               pt::idx_t sm_st_idx, Vertex fl) {
    return Vertex(start, end, sm_st_idx, fl);
  }

  static Vertex make_dummy() {
    return Vertex();
  }

  // ---------
  // getter(s)
  // ---------
  povu::types::id_t get_idx() const { return this->idx_; }
  pgt::id_or_t get_start() const { return this->a_; }
  pgt::id_or_t get_end() const { return this->z_; }
  vt_e get_type() const { return this->type_; }
  pt::idx_t get_sl_st_idx() const { return this->sl_st_idx_; }

  sl_type_e get_sl_type() const { return this->fl_type_; }

  // get the idx of ui in the spanning tree
  pt::idx_t get_ai() const { return this->ai_; }
  pt::idx_t get_zi() const { return this->zi_; }

  // ---------
  // setter(s)
  // ---------
  void set_type(vt_e t) { this->type_ = t; }
  void set_v_idx(pt::idx_t v_idx) { this->idx_ = v_idx; }

  // ---------
  // other(s)
  // ---------
  std::string as_str() const {

    if (this->type_ == vt_e::dummy) {
      return ".";
    }

    std::string s;
    s += this->a_.as_str();
    s += this->z_.as_str();
    return s;
  }
};
} // namespace povu::types::pvst

#endif
