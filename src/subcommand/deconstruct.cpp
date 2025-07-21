#include "./deconstruct.hpp"


namespace povu::subcommands::deconstruct {

void deconstruct_component(bd::VG *g, std::size_t component_id,
                           const core::config &app_config) {
  std::string fn_name = std::format("[povu::deconstruct::{}]", __func__);


#ifdef DEBUG
  if (app_config.verbosity() > 4) {
    std::cerr << "\n";
    g->print_dot(std::cerr);
    std::cerr << "\n";
  }
#endif

  pst::Tree st { bd::compute_spanning_tree(*g) };
  delete g;

  ptu::tree_meta tm = ptu::gen_tree_meta(st);

#ifdef DEBUG
  if (app_config.verbosity() > 4) {
    std::cerr << "\n";
    st.print_dot(std::cerr);
    std::cerr << "\n";
  }
#endif

  pvtr::Tree flubble_tree = pfl::find_flubbles(st, app_config);
  povu::tiny::find_tiny(st, flubble_tree, tm);
  povu::parallel::find_parallel(st, flubble_tree, tm);

  if (app_config.find_hubbles()) {
    povu::concealed::find_concealed(st, flubble_tree, tm);
    povu::midi::find_midi(st, flubble_tree, tm);
    povu::smothered::find_smothered(st, flubble_tree, tm);
  }

  povu::io::pvst::write_bub(flubble_tree, std::to_string(component_id), app_config);

  return;
}

std::pair<uint32_t, uint32_t> thread_count(const core::config &app_config, std::size_t item_count) {
  /* Divide the number of components into chunks for each thread */
  unsigned int total_threads = std::thread::hardware_concurrency();

  // std::cerr << std::format("{} Max threads: {}\n", fn_name, total_threads);
  std::size_t conf_num_threads = static_cast<std::size_t>(app_config.thread_count());
  uint32_t num_threads = (conf_num_threads > total_threads) ? total_threads : conf_num_threads;
  uint32_t chunk_size = item_count / num_threads;

  return std::make_pair(num_threads, chunk_size);
}

void do_deconstruct(const core::config &app_config) {
  std::string fn_name = std::format("[povu::deconstruct::{}]", __func__);
  std::size_t ll = app_config.verbosity(); // ll for log level, to avoid long names. good idea?

  bd::VG *g = get_vg(app_config);

  if (ll > 1) std::cerr << std::format("{} Finding components\n", fn_name);
  std::vector<bd::VG *> components = bd::componetize(*g);

  delete g;

  if (ll > 1) std::cerr << std::format("{} Found {} components\n", fn_name, components.size());

  auto [num_threads, chunk_size] = thread_count(app_config, components.size());

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
          components[i]->summary(false);
        }


        deconstruct_component(components[i], component_id, std::ref(app_config)); // Pass by reference
      }
    });
  }

  for (auto& thread : threads) { thread.join(); } // Wait for all threads to finish

  return;
}

} // namespace povu::subcommands
