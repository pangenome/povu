#include "../graph/digraph.hpp"
#include <string>
#include "../core/core.hpp"
#include "../graph/bidirected.hpp"

namespace io {
  digraph::DiGraph gfa_to_digraph(const char* filename);

  bidirected::VariationGraph gfa_to_vg(const char* filename);
} // namespace io
