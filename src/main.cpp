#include <cstddef>
#include <format>
#include <iostream>
#include <string>
#include <thread>
#include <vector>

#include "./cli/app.hpp"
#include "./cli/cli.hpp"
#include "./common/types.hpp"
#include "./common/utils.hpp"
#include "./graph/graph.hpp"
#include "./io/io.hpp"
#include "./povu.hpp"
#include "genomics/genomics.hpp"

namespace bd = povu::bidirected;
namespace pt = povu::types;
namespace pgt = povu::graph_types;

void do_info(const core::config &app_config) {
  std::string fn_name = std::format("[povu::main::{}]", __func__);


  // -----
  // read the input gfa into a bidirected variation graph
  // -----
  povu::graph::Graph g = io::from_gfa::to_pv_graph(app_config.get_input_gfa().c_str(), app_config);

  g.summary();
}


void do_call(const core::config& app_config) {
  std::string fn_name = std::format("[povu::main::{}]", __func__);

  std::chrono::duration<double> timeRefRead;
  auto t0 = pt::Time::now();


  // -----
  // read the input gfa into a bidirected variation graph
  // -----
  if (app_config.verbosity() > 2)  { std::cerr << std::format ("{} Reading graph\n", fn_name); }
  bd::VG bd_vg = io::from_gfa::to_bd(app_config.get_input_gfa().c_str(), app_config);

  if (app_config.verbosity() > 1) {
    timeRefRead = pt::Time::now() - t0;
    povu::utils::report_time(std::cerr, fn_name, "read_gfa", timeRefRead);
  }


  std::vector<std::filesystem::path> flubble_files = povu::io::generic::get_files(app_config.get_forest_dir(), ".flb");

  std::vector<pgt::flubble> canonical_flubbles;

  // TODO: do in parallel
  for (const auto& fp: flubble_files) {
    std::cerr << std::format("Reading flubble file: {}\n", fp.string());
    auto res = povu::io::bub::read_canonical_fl(fp.string());
    // append the results of the vector with these results
    canonical_flubbles.insert(canonical_flubbles.end(), res.begin(), res.end());
  }

  // ------
  // read from a flubble tree in flb in format
  // -----
  povu::genomics::call_variants(canonical_flubbles, bd_vg, app_config);

  return;
}


void do_deconstruct(const core::config &app_config) {
  std::string fn_name = std::format("[povu::main::{}]", __func__);

  std::chrono::duration<double> timeRefRead;
  auto t0 = pt::Time::now();

  // -----
  // read the input gfa into a bidirected variation graph
  // -----
  if (app_config.verbosity() > 2)  { std::cerr << std::format ("{} Reading graph\n", fn_name); }
  povu::graph::Graph g =  io::from_gfa::to_pv_graph(app_config.get_input_gfa().c_str(), app_config);

  if (app_config.verbosity() > 1) {
    timeRefRead = pt::Time::now() - t0;
    povu::utils::report_time(std::cerr, fn_name, "read_gfa", timeRefRead);
  }

  // -----
  //
  // -----
  if (app_config.verbosity() > 2)  { std::cerr << std::format("{} Finding components\n", fn_name); }
  std::vector<povu::graph::Graph> components =  povu::graph::componetize(g, app_config);

  if (app_config.verbosity() > 1) {
    std::cerr << std::format("{} Found {} components\n", fn_name, components.size());
  }

  // -----
  // Divide the vector into chunks for each thread//
  // -----
  std::size_t chunk_size = components.size() / static_cast<std::size_t>(app_config.thread_count() );

  // Create and launch threads
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

        povu::bin::deconstruct(std::ref(components[i]), component_id, std::ref(app_config)); // Pass by reference
      }
    });
  }

  // Wait for all threads to finish
  for (auto& thread : threads) { thread.join(); }

  return;
}


/**
 * @brief main function
 *
 * @param argc
 * @param argv
 * @return int
 */
int main(int argc, char *argv[]) {
  std::string fn_name = std::format("[povu::main::{}]", __func__);

  core::config app_config;
  cli::cli(argc, argv, app_config);

  if (app_config.verbosity()) { app_config.dbg_print(); }

  switch (app_config.get_task()) {
    case core::task_t::deconstruct:
      do_deconstruct(app_config);
      break;
    case core::task_t::call:
      do_call(app_config);
      break;
    case core::task_t::info:
      do_info(app_config);
      break;
    default:
      std::cerr << std::format("{} Task not recognized\n", fn_name);
      break;
  }

  return 0;
}
