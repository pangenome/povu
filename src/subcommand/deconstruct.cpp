#include "./subcommands.hpp"

namespace povu::graph_ops {
namespace pst = povu::spanning_tree;
namespace pbd = povu::bidirected;

/**
 * @brief
 *
*/
pst::Tree biedge_and_cycle_equiv(const pbd::VG& g, std::size_t component_id, const core::config& app_config)  {
  std::string fn_name = std::format("[povu::subcommand::{}]", __func__);

  // g.print_dot();

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
  //povu::algorithms::eulerian_cycle_equiv(st);
  povu::algorithms::simple_cycle_equiv(st);

  if (app_config.print_dot() && app_config.verbosity() > 4) { std::cout << "\n\n" << "Updated Spanning Tree " << component_id << "\n\n";
    st.print_dot();
  }

  return st;
}

} // namespace povu::graph_ops

namespace povu::subcommands {

void deconstruct_component(const bd::VG& g, std::size_t component_id, const core::config& app_config) {
  std::string fn_name = std::format("[povu::deconstruct::{}]", __func__);

  pst::Tree st = povu::graph_ops::biedge_and_cycle_equiv(g, component_id, app_config);
  pvtr::Tree<pgt::flubble> flubble_tree = povu::graph::flubble_tree::st_to_ft(st);
  povu::io::bub::write_bub(flubble_tree, std::to_string(component_id), app_config);
}

/**
 * in this way the initial VG gors out of scope after the function read and componetize
 */
std::vector<bd::VG> get_components(const core::config &app_config) {
  std::string fn_name = std::format("[povu::deconstruct::{}]", __func__);

  std::chrono::duration<double> timeRefRead;
  auto t0 = pt::Time::now();

  /* Read the input gfa into a bidirected variation graph */
  if (app_config.verbosity() > 2)  { std::cerr << std::format ("{} Reading graph\n", fn_name); }
  bd::VG g = io::from_gfa::to_bd(app_config.get_input_gfa().c_str(), app_config);
  if (app_config.verbosity() > 1) {
    timeRefRead = pt::Time::now() - t0;
    povu::utils::report_time(std::cerr, fn_name, "read_gfa", timeRefRead);
    t0 = pt::Time::now();
  }

  if (app_config.verbosity() > 2)  { std::cerr << std::format("{} Finding components\n", fn_name); }
  std::vector<bd::VG> components =  bd::componetize(g, app_config);

  return components;
}


void do_deconstruct(const core::config &app_config) {
  std::string fn_name = std::format("[povu::deconstruct::{}]", __func__);
  std::vector<bd::VG> components =  get_components(app_config);

  if (app_config.verbosity() > 1) {
    std::cerr << std::format("{} Found {} components\n", fn_name, components.size());
  }

  /* Divide the vector into chunks for each thread */
  std::size_t chunk_size = components.size() / static_cast<std::size_t>(app_config.thread_count() );

  /* Create and launch threads */
  std::vector<std::thread> threads(app_config.thread_count());
  std::size_t start, end;
  for (unsigned int i {}; i < app_config.thread_count(); ++i) {
    start = i * chunk_size;
    end = (i == app_config.thread_count() - 1) ? components.size() : (i + 1) * chunk_size;

    threads[i] = std::thread([start, end,fn_name, app_config, &components] {
      for (std::size_t i{start}; i < end; i++) {

        std::size_t component_id {i + 1};

        if (app_config.verbosity()) {
          std::cerr << std::format("{} Handling component: {}\n", fn_name, component_id);
        }

        if (components[i].size() < 3) {
          if (app_config.verbosity() > 2) {
            std::cerr << std::format("{} Skipping component {} because it is too small. (size: {})\n", fn_name, component_id, components[i].size());
          }
          continue;
        }

        if (app_config.verbosity() > 3 && app_config.thread_count() == 1 && app_config.get_task() != core::task_t::info) {
          components[i].summary();
        }

        deconstruct_component(std::ref(components[i]), component_id, std::ref(app_config)); // Pass by reference
      }
    });
  }

  for (auto& thread : threads) { thread.join(); } // Wait for all threads to finish

  return;
}

} // namespace povu::subcommands
