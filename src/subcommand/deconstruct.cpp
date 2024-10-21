#include "../algorithms/algorithms.hpp"
#include "../common/types.hpp"
#include "../graph/biedged.hpp"
#include "../graph/flubble_tree.hpp"
#include "../graph/spanning_tree.hpp"
#include "../io/io.hpp"
#include "../povu.hpp"

namespace povu::graph_ops {
namespace pst = povu::spanning_tree;

/**
 * @brief
 *
*/
pst::Tree biedge_and_cycle_equiv(const povu::graph::Graph& g, std::size_t component_id, const core::config& app_config)  {
  std::string fn_name = std::format("[povu::subcommand::{}]", __func__);

  // convert the bidirected variation graph into a biedged variation graph
  if (app_config.verbosity() > 2) { std::cerr << std::format("{} Bi-edging {}\n", fn_name, component_id); }
  biedged::BVariationGraph bg(g); // will add dummy vertices

  if (app_config.print_dot() && app_config.verbosity() > 4 ) { std::cout << "\n\n" << "Biedged" << "\n\n";
    bg.print_dot();
  }

  if (app_config.verbosity() > 2) { std::cerr << std::format("{} Generating spanning tree {}\n", fn_name, component_id); }
  pst::Tree st = bg.compute_spanning_tree();

  if (app_config.print_dot() && app_config.verbosity() > 4) { std::cout << "\n\n" << "Spanning Tree " << component_id << "\n\n";
    st.print_dot();
  }

  if (app_config.verbosity() > 2) { std::cerr << std::format("{} Computing Cycle Equivalence {}\n", fn_name, component_id); }
  povu::algorithms::eulerian_cycle_equiv(st);

  if (app_config.print_dot() && app_config.verbosity() > 4) { std::cout << "\n\n" << "Updated Spanning Tree " << component_id << "\n\n";
    st.print_dot();
  }

  return st;
}
}

namespace povu::bin {

namespace pst = povu::spanning_tree;
namespace pgt = povu::graph_types;
namespace pvtr = povu::tree;


void deconstruct(const povu::graph::Graph& g, std::size_t component_id, const core::config& app_config) {
  std::string fn_name = std::format("[povu::subcommand::{}]", __func__);

  pst::Tree st = povu::graph_ops::biedge_and_cycle_equiv(g, component_id, app_config);
  pvtr::Tree<pgt::flubble> flubble_tree = povu::graph::flubble_tree::st_to_ft(st);
  povu::io::bub::write_bub(flubble_tree, std::to_string(component_id), app_config);
}

} // namespace povu::bin

namespace povu::lib {

namespace pst = povu::spanning_tree;
namespace pgt = povu::graph_types;
namespace pvtr = povu::tree;


pvtr::Tree<pgt::flubble> deconstruct_to_ft(const povu::graph::Graph& g, std::size_t component_id, const core::config& app_config) {
  std::string fn_name = std::format("[povu::subcommand::{}]", __func__);

  pst::Tree st = povu::graph_ops::biedge_and_cycle_equiv(g, component_id, app_config);
  pvtr::Tree<pgt::flubble> flubble_tree = povu::graph::flubble_tree::st_to_ft(st);

  return flubble_tree;
}


std::vector<pgt::flubble> deconstruct_to_enum(const povu::graph::Graph& g, std::size_t component_id, const core::config& app_config) {
  std::string fn_name = std::format("[povu::subcommand::{}]", __func__);

  pst::Tree st = povu::graph_ops::biedge_and_cycle_equiv(g, component_id, app_config);
  std::vector<pgt::flubble> flubbles = povu::graph::flubble_tree::enumerate(st);
  povu::io::generic::write_txt(flubbles, std::to_string(component_id), app_config);

  return flubbles;
}
} // namespace povu::lib
