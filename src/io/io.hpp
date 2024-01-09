#include "../graph/digraph.hpp"
#include <string>
#include "../core/core.hpp"
#include "../graph/bidirected.hpp"

namespace io::from_gfa {
  digraph::DiGraph to_digraph(const char* filename);
  bidirected::VariationGraph to_vg(const char* filename, const core::config& app_config);
}; // namespace io::from_gfa
