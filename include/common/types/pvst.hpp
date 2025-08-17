#ifndef POVU_TYPES_PVST_HPP
#define POVU_TYPES_PVST_HPP

#include <algorithm>

#include "./core.hpp"
#include "./constants.hpp"
#include "./graph.hpp"
#include "../compat.hpp"

/* === PVST pangenome variation structure tree === */

namespace povu::types::pvst {
namespace pc = povu::constants;
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

std::ostream& operator<<(std::ostream& os, vt_e t);



constexpr bool is_fl_like(pvst::vt_e t) noexcept {
  switch (t) {
  case pvst::vt_e::flubble:
  case pvst::vt_e::tiny:
  case pvst::vt_e::parallel:
    return true;
  default:
    return false;
  }
}

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

// bounds in the spanning tree
struct bounds_t {
  pt::idx_t upper;
  pt::idx_t lower; // when invalid all leaves are the lower boundaries
};

// bounds in the bidirected graph
struct traversal_params_t {
  pgt::id_or_t start;
  pgt::id_or_t end;

  bool traversable;
  // if true traverse start to end else end to start
  bool s2e; // start to end or end to start
};

traversal_params_t null_tp();


/* an abstract class for vertices  */
class VertexBase {
  pt::id_t idx_; // idx of the vertex in the vst
  pt::idx_t height_; // height of the vertex in the tree
  vt_e type_;

public:
  // ——— constructors ———
  VertexBase(pt::id_t idx, vt_e type)
    : idx_(idx), height_(pc::INVALID_IDX), type_(type) {}

  // ——— getters ———
  pt::id_t get_idx() const { return this->idx_; }
  vt_e get_type() const { return this->type_; }
  pt::idx_t get_height() const { return this->height_; }


  // ——— pure virutal functions ———
  virtual std::string as_str() const = 0;
  virtual traversal_params_t get_traversal_params() const = 0;
  virtual ~VertexBase() = default;

  // ——— setters ———
  void set_idx(pt::id_t idx) { this->idx_ = idx; }
  void set_type(vt_e type) { this->type_ = type; }
  void set_height(pt::idx_t height) { this->height_ = height; }
};


class Dummy : public VertexBase {
public:
  // ——— constructors ———
  Dummy() : VertexBase(pc::INVALID_IDX, vt_e::dummy) {}

  // ——— getters ———
  std::string as_str() const override { return "."; }
  traversal_params_t get_traversal_params() const override {
    // dummy vertex does not have any walks
    return null_tp();
  }
};


class Flubble : public VertexBase {
  pgt::id_or_t a_; // start
  pgt::id_or_t z_; // end

  pt::idx_t ai_; // a_i
  pt::idx_t zi_; // z_i

  pt::idx_t m_; // m
  pt::idx_t n_; // n

public:

  // ——— constructors ———
  Flubble(vt_e typ, pgt::id_or_t a, pgt::id_or_t z)
      : VertexBase(pc::INVALID_IDX, typ), a_(a), z_(z),
        ai_(pc::INVALID_IDX), zi_(pc::INVALID_IDX),
        m_(pc::INVALID_IDX), n_(pc::INVALID_IDX) {}

  Flubble(pgt::id_or_t a, pgt::id_or_t z, pt::idx_t ai, pt::idx_t zi)
    : VertexBase(pc::INVALID_IDX, vt_e::flubble), a_(a), z_(z), ai_(ai), zi_(zi),
      m_(pc::INVALID_IDX), n_(pc::INVALID_IDX) {}

  // ——— getters ———
  pgt::id_or_t get_a() const { return this->a_; }
  pgt::id_or_t get_z() const { return this->z_; }
  pt::idx_t get_ai() const { return this->ai_; }
  pt::idx_t get_zi() const { return this->zi_; }
  pt::idx_t get_m() const { return this->m_; }
  pt::idx_t get_n() const { return this->n_; }
  bounds_t get_bounds() const { return {this->get_ai(), this->get_zi()}; }
  //bdg_bounds_t get_bdg_bounds() const { return {this->get_a(), this->get_z()}; }
  traversal_params_t get_traversal_params() const override {
    // dummy vertex does not have any walks
    return traversal_params_t{
      this->get_a(), this->get_z(), true, true
    };
  }

  // ——— setters ———
  void set_m(pt::idx_t m) { this->m_ = m; }
  void set_n(pt::idx_t n) { this->n_ = n; }

  // ——— others ———
  std::string as_str() const override {
    return pv_cmp::format("{}{}", this->a_.as_str(), this->z_.as_str());
  }
};


class Concealed : public VertexBase {
  pt::idx_t fl_idx_; // v idx of the parent flubble in the PVST
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
  // ——— constructors ———
  //Concealed(pgt::id_or_t fl_b, pgt::id_or_t cn_b)
  //  : VertexBase(pc::INVALID_IDX, vt_e::slubble), fl_b_(fl_b), cn_b_(cn_b) {}

  Concealed(pgt::id_or_t fl_b, pgt::id_or_t cn_b, bounds_t bounds,
            pt::idx_t fl_idx, sl_type_e sl_type, pt::idx_t sl_st_idx)
    : VertexBase(pc::INVALID_IDX, vt_e::slubble), fl_idx_(fl_idx),
      sl_type_(sl_type), sl_st_idx_(sl_st_idx), fl_b_(fl_b), cn_b_(cn_b),
      bounds_(bounds) {}

  // ——— getters ———
  pt::idx_t get_fl_idx() const { return this->fl_idx_; }
  sl_type_e get_sl_type() const { return this->sl_type_; }
  pt::idx_t get_sl_st_idx() const { return this->sl_st_idx_; }
  pgt::id_or_t get_fl_b() const { return this->fl_b_; }
  pgt::id_or_t get_cn_b() const { return this->cn_b_; }
  bounds_t get_bounds() const { return this->bounds_; }
  traversal_params_t get_traversal_params() const override {
    // dummy vertex does not have any walks
    return traversal_params_t{
      this->get_fl_b(), this->get_cn_b(), true, this->with_ai()
    };
  }

  // ——— others ———
  std::string as_str() const override {
    if (with_ai()) { // formed with a
      return pv_cmp::format("{}{}", this->fl_b_.as_str(), this->cn_b_.as_str());
    }
    else  { // formed with z
      return pv_cmp::format("{}{}", this->cn_b_.as_str(), this->fl_b_.as_str());
    }
  }
};


class Smothered : public VertexBase {
  pt::idx_t cn_idx_; // idx of the concealed vertex
  pt::idx_t sm_st_idx_; // idx in the spanning tree for smothered vertex

  // b for boundary
  pgt::id_or_t cn_b_; // g or s
  pgt::id_or_t sm_b_; // e or w

  // is true when cn_b_ is an ancestor of sm_b_
  bool cn_b_is_ans_;

  cn_type_e cn_type_; // type of the concealed vertex (g or s)

  bounds_t bounds_;

public:

  // ——— constructors ———
  // Smothered(pgt::id_or_t cn_b, pgt::id_or_t sm_b)
  //     : VertexBase(pc::INVALID_IDX, vt_e::smothered), cn_b_(cn_b), sm_b_(sm_b) {}

  Smothered(pgt::id_or_t cn_b, pgt::id_or_t sm_b, pt::idx_t cn_idx,
            bool cn_b_is_ans, pt::idx_t sm_st_idx, cn_type_e sm_type, bounds_t bounds)
      : VertexBase(pc::INVALID_IDX, vt_e::smothered),
        cn_idx_(cn_idx), sm_st_idx_(sm_st_idx), cn_b_(cn_b), sm_b_(sm_b),
        cn_b_is_ans_(cn_b_is_ans), cn_type_(sm_type), bounds_(bounds) {}

  // ——— getters ———
  pt::idx_t get_cn_idx() const { return this->cn_idx_; }
  pt::idx_t get_sm_st_idx() const { return this->sm_st_idx_; }
  bounds_t get_bounds() const { return this->bounds_; }
  cn_type_e get_cn_type() const { return this->cn_type_; }
  bool is_cn_b_ancestor() const { return this->cn_b_is_ans_; }
  traversal_params_t get_traversal_params() const override {
    // dummy vertex does not have any walks
    return traversal_params_t{
      this->sm_b_, this->cn_b_, true, this->cn_b_is_ans_
    };
  }

  std::string as_str() const override {
    if (this->cn_type_ == cn_type_e::g) { // g
      if (this->cn_b_is_ans_) { // cn_b is ancestor of sm_b
        return pv_cmp::format("{}{}", this->cn_b_.as_str(), this->sm_b_.as_str());
      }
      else { // sm_b is ancestor of cn_b
        return pv_cmp::format("{}{}", this->sm_b_.as_str(), this->cn_b_.as_str());
      }
    }
    else { // s
      return pv_cmp::format("{}{}", this->sm_b_.as_str(), this->cn_b_.as_str());
    }
  }
};


class MidiBubble : public VertexBase {
  pt::idx_t g_cn_idx_;
  pt::idx_t s_cn_idx_;
  pgt::id_or_t g_;
  pgt::id_or_t s_;

public:

  // ——— constructors ———
  // MidiBubble(pgt::id_or_t g, pgt::id_or_t s)
  //     : VertexBase(pc::INVALID_IDX, vt_e::midi), g_(g), s_(s) {}

  MidiBubble(pt::idx_t g_cn_idx, pgt::id_or_t g, pt::idx_t s_cn_idx, pgt::id_or_t s)
      : VertexBase(pc::INVALID_IDX, vt_e::midi), g_cn_idx_(g_cn_idx),
        s_cn_idx_(s_cn_idx), g_(g), s_(s) {}

  // ——— getters ———
  bounds_t get_bounds() const {
    return bounds_t { std::min(g_cn_idx_, s_cn_idx_), std::max(g_cn_idx_, s_cn_idx_) };
  }
  pgt::id_or_t get_g() const { return this->g_; }
  pgt::id_or_t get_s() const { return this->s_; }

  traversal_params_t get_traversal_params() const override {
    // dummy vertex does not have any walks
    return traversal_params_t{this->get_g(), this->get_s(), true, true};
  }

    std::string as_str() const override {
      return pv_cmp::format("{}{}", this->g_.as_str(), this->s_.as_str());
    }
  };

} // namespace povu::types::pvst

#endif
