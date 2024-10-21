#ifndef POVU_HPP
#define POVU_HPP

#include "./cli/app.hpp"
#include "./graph/graph.hpp"
#include "./common/types.hpp"
#include "./graph/flubble_tree.hpp"

namespace povu::bin {
void deconstruct(const povu::graph::Graph& g, std::size_t component_id, const core::config& app_config);
}

namespace povu::lib {
namespace pvtr = povu::tree;
namespace pgt = povu::graph_types;

std::vector<pgt::flubble> deconstruct_to_enum(const povu::graph::Graph& g, std::size_t component_id, const core::config& app_config);
pvtr::Tree<pgt::flubble> deconstruct_to_ft(const povu::graph::Graph& g, std::size_t component_id, const core::config& app_config);
}

#endif
