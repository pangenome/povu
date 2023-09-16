#include "../graph/digraph.hpp"
#include <string>
#include "../core/core.hpp"

namespace io {
  digraph::DiGraph gfa_to_digraph(const char* filename);
} // namespace io
