#include "../graph/spanning_tree.hpp"
#include "../graph/tree.hpp"

namespace pst {
/**
  * @brief computes the PST from the spanning tree
  *
  * @param st the spanning tree
  * @return tree::Tree the PST
  */
tree::Tree compute_pst(spanning_tree::Tree &st);
} // namespace pst
