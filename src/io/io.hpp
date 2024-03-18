#include "../core/core.hpp"
#include "../graph/bidirected.hpp"

namespace io::from_gfa {
  bidirected::VariationGraph to_vg(const char* filename, const core::config& app_config);
}; // namespace io::from_gfa
