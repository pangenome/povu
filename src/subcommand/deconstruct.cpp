#include "./subcommands.hpp"
#include <iostream>


namespace povu::subcommands {

void deconstruct_component(bd::VG *g, std::size_t component_id, const core::config& app_config) {
  std::string fn_name = std::format("[povu::deconstruct::{}]", __func__);

  g->untip();
  pst::Tree st { bd::compute_spanning_tree(*g) };
  delete g;
  povu::algorithms::simple_cycle_equiv(st); // find equivalence classes
  pvtr::Tree<pgt::flubble> flubble_tree = povu::graph::flubble_tree::st_to_ft(st);
  povu::io::bub::write_bub(flubble_tree, std::to_string(component_id), app_config);
}


/**
 * Read the input gfa into a bidirected variation graph
*/
bd::VG *get_vg(const core::config &app_config) {
  std::string fn_name = std::format("[povu::deconstruct::{}]", __func__);
  std::size_t ll = app_config.verbosity(); // log level

  std::chrono::duration<double> timeRefRead;
  auto t0 = pt::Time::now();

  if (ll > 2) std::cerr << std::format("{} Reading graph\n", fn_name);
  bd::VG *g = io::from_gfa::to_bd(app_config.get_input_gfa().c_str(), app_config);

  if (ll > 1) {
    timeRefRead = pt::Time::now() - t0;
    povu::utils::report_time(std::cerr, fn_name, "read_gfa", timeRefRead);
    t0 = pt::Time::now();
  }

  return g;
}

void do_deconstruct(const core::config &app_config) {
  std::string fn_name = std::format("[povu::deconstruct::{}]", __func__);
  std::size_t ll = app_config.verbosity(); // ll for log level, to avoid long names. good idea?

  bd::VG *g = get_vg(app_config);

  #ifdef DEBUG
  g->summary();
  #endif

  if (ll > 1) std::cerr << std::format("{} Finding components\n", fn_name);
  std::vector<bd::VG *> components = bd::componetize(*g);

  delete g;

  if (ll > 1) std::cerr << std::format("{} Found {} components\n", fn_name, components.size());

  /* Divide the number of components into chunks for each thread */
  unsigned int total_threads = std::thread::hardware_concurrency();
  std::size_t conf_num_threads = static_cast<std::size_t>(app_config.thread_count());
  unsigned int num_threads = (conf_num_threads > total_threads) ? total_threads : conf_num_threads;
  std::size_t chunk_size = components.size() / num_threads;

  /* Create and launch threads */
  std::vector<std::thread> threads(num_threads);
  std::size_t start, end;
  for (unsigned int i {}; i < num_threads; ++i) {
    start = i * chunk_size;
    end = (i == num_threads - 1) ? components.size() : (i + 1) * chunk_size;

    threads[i] = std::thread([start, end, fn_name, num_threads, app_config, &components] {
      for (std::size_t i{start}; i < end; i++) {

        std::size_t component_id {i + 1};

        if (app_config.verbosity()) {
          std::cerr << std::format("{} Handling component: {}\n", fn_name, component_id);
        }

        if (components[i]->vtx_count() < 3) {
          if (app_config.verbosity() > 2) {
            std::cerr << std::format("{} Skipping component {} because it is too small. (size: {})\n", fn_name, component_id, components[i]->vtx_count());
          }
          continue;
        }

        if (app_config.verbosity() > 3 && num_threads == 1) {
          components[i]->summary();
        }

        deconstruct_component(components[i], component_id, std::ref(app_config)); // Pass by reference
      }
    });
  }

  for (auto& thread : threads) { thread.join(); } // Wait for all threads to finish

  return;
}

} // namespace povu::subcommands
