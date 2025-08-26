#ifndef POVU_GENOMICS_GRAPH_HPP
#define POVU_GENOMICS_GRAPH_HPP

#include <algorithm>
#include <deque>
#include <set>
#include <stack>
#include <string>
#include <sys/types.h>
#include <utility>
#include <vector>

#include "../common/compat.hpp"
#include "../common/types/core.hpp"
#include "../graph/types.hpp"
#include "../graph/pvst.hpp"
#include "../common/log.hpp"
#include "../graph/bidirected.hpp"

namespace povu::genomics::graph {
inline constexpr std::string_view MODULE = "povu::genomics::graph";

namespace pvst = povu::pvst;
namespace pgt = povu::types::graph;

// Maximum number of steps to take from flubble start to end
const pt::idx_t MAX_FLUBBLE_STEPS{20};

// direction for traversing a vertex in a bidirected graph
enum class dir_e { in, out };

const dir_e IN = dir_e::in;
const dir_e OUT = dir_e::out;

typedef pgt::id_or_t idx_or_t; // specifically for idx instead of id

/**
 * a collection of walks within a region of variation from start to end
 */
class RoV {
  std::vector<pgt::walk_t> walks_;
  const pvst::VertexBase *pvst_vtx;

public:
  // --------------
  // constructor(s)
  // --------------

  RoV(const pvst::VertexBase *v) : walks_(), pvst_vtx(v) {}

  // ---------
  // getter(s)
  // ---------

  pt::idx_t walk_count() const { return this->walks_.size(); }
  const pvst::VertexBase *get_pvst_vtx() const { return this->pvst_vtx; }
  const std::vector<pgt::walk_t> &get_walks() const { return this->walks_; }
  std::vector<pgt::walk_t> &get_walks_mut() { return this->walks_; }

  // ---------
  // setter(s)
  // ---------

  void set_walks(std::vector<pgt::walk_t> &&walks) { this->walks_ = walks; }

  // --------
  // other(s)
  // --------

  std::string as_str() const { return this->pvst_vtx->as_str(); }
};

void find_walks(const bd::VG &g, RoV &rov);
} // namespace povu::genomics::graph

#endif // POVU_GENOMICS_GRAPH_HPP
