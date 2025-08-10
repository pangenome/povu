#include "./gfa2vcf.hpp"
#include "./common.hpp"
#include "./call.hpp"

#include <filesystem>
#include <iostream>
#include <thread>
#include <cstdlib>

#include "../../include/graph/bidirected.hpp"
#include "../../include/graph/spanning_tree.hpp"
#include "../../include/common/types/pvst.hpp"
#include "../../include/common/tree_utils.hpp"
#include "../../include/algorithms/flubbles.hpp"
#include "../../include/algorithms/tiny.hpp"
#include "../../include/algorithms/parallel.hpp"
#include "../../include/algorithms/concealed.hpp"
#include "../../include/algorithms/midi.hpp"
#include "../../include/algorithms/smothered.hpp"
#include "../../include/genomics/variants.hpp"
#include "../io/pvst.hpp"
#include "../io/to_vcf.hpp"
#include "../io/common.hpp"

namespace fs = std::filesystem;
namespace bd = povu::bidirected;
namespace pvtr = povu::types::pvst;
namespace pcs = povu::subcommands::common;
namespace pic = povu::io::common;
namespace pvt = povu::types;
namespace pt = povu::types;
namespace pg = povu::genomics;
namespace piv = povu::io::to_vcf;
namespace pst = povu::spanning_tree;
namespace ptu = povu::tree_utils;
namespace pfl = povu::flubbles;

namespace povu::subcommands::gfa2vcf {

// Helper function to get references (similar to call::get_refs)
pt::status_t get_refs(core::config &app_config) {
  if (app_config.get_refs_input_fmt() == core::input_format_e::file_path) {
    std::vector<std::string> refs;
    pic::read_lines_to_vec_str(app_config.get_references_txt(), &refs);
    app_config.set_reference_paths(std::move(refs));
    return 0;
  }
  else if (app_config.get_refs_input_fmt() == core::input_format_e::params) {
    // If path prefixes are provided, we need to read all paths from GFA and filter by prefix
    if (!app_config.get_path_prefixes().empty()) {
      std::vector<std::string> filtered_refs = call::filter_paths_by_prefix(app_config);
      app_config.set_reference_paths(std::move(filtered_refs));
    }
    // If explicit reference paths are provided, they're already set
    return 0;
  }

  return -1;
}

void do_gfa2vcf(const core::config &app_config) {
  std::string fn_name = std::format("[povu::subcommands::{}]", __func__);
  std::size_t ll = app_config.verbosity();
  
  // Create a temporary directory for the forest files
  char temp_template[] = "/tmp/povu_gfa2vcf_XXXXXX";
  char* temp_dir = mkdtemp(temp_template);
  if (temp_dir == nullptr) {
    std::cerr << fn_name << " Error: Could not create temporary directory" << std::endl;
    exit(EXIT_FAILURE);
  }
  std::string temp_dir_str(temp_dir);
  
  if (ll > 0) {
    std::cerr << fn_name << " Using temporary directory: " << temp_dir_str << std::endl;
  }
  
  // Step 1: Decompose - generate forest of PVST files
  if (ll > 0) {
    std::cerr << fn_name << " Step 1: Decomposing graph..." << std::endl;
  }
  
  bd::VG *g = pcs::get_vg(app_config);
  
  if (ll > 1) std::cerr << std::format("{} Finding components\n", fn_name);
  std::vector<bd::VG *> components = bd::componetize(*g);
  
  delete g;
  
  if (ll > 1) std::cerr << std::format("{} Found {} components\n", fn_name, components.size());
  
  // Process each component and generate PVST files
  for (std::size_t i = 0; i < components.size(); i++) {
    std::size_t component_id = i + 1;
    
    if (ll > 2) {
      std::cerr << std::format("{} Handling component: {}\n", fn_name, component_id);
    }
    
    if (components[i]->vtx_count() < 3) {
      if (ll > 2) {
        std::cerr << std::format("{} Skipping component {} because it is too small. (size: {})\n", 
                                fn_name, component_id, components[i]->vtx_count());
      }
      delete components[i];
      continue;
    }
    
    // Process the component
    pst::Tree st { bd::compute_spanning_tree(*components[i]) };
    delete components[i];
    
    ptu::tree_meta tm = ptu::gen_tree_meta(st);
    
    pvtr::Tree flubble_tree = pfl::find_flubbles(st, app_config);
    povu::tiny::find_tiny(st, flubble_tree, tm);
    povu::parallel::find_parallel(st, flubble_tree, tm);
    
    if (app_config.find_hubbles()) {
      povu::concealed::find_concealed(st, flubble_tree, tm);
      povu::midi::find_midi(st, flubble_tree, tm);
      povu::smothered::find_smothered(st, flubble_tree, tm);
    }
    
    // Write PVST to temp directory
    core::config temp_config = app_config;
    temp_config.set_output_dir(temp_dir_str);
    povu::io::pvst::write_pvst(flubble_tree, std::to_string(component_id), temp_config);
  }
  
  // Step 2: Call - generate VCF from PVST files
  if (ll > 0) {
    std::cerr << fn_name << " Step 2: Calling variants..." << std::endl;
  }
  
  // Read PVST files from temp directory
  std::vector<pvtr::Tree> pvsts;
  std::vector<fs::path> fps = pic::get_files(temp_dir_str, ".pvst");
  
  if (fps.empty()) {
    std::cerr << fn_name << " Warning: No PVST files generated in " << temp_dir_str << std::endl;
    fs::remove_all(temp_dir_str);
    return;
  }
  
  for (std::size_t i = 0; i < fps.size(); i++) {
    pvtr::Tree pvst = povu::io::pvst::read_pvst(fps[i].string());
    pvst.comp_heights();
    pvsts.push_back(std::move(pvst));
  }
  
  // Load graph and references for variant calling
  get_refs(const_cast<core::config&>(app_config));
  bd::VG *call_graph = pcs::get_vg(app_config);
  
  // Generate VCF records
  pvt::VcfRecIdx vcf_recs = pg::gen_vcf_rec_map(pvsts, *call_graph);
  
  // Write VCF to stdout
  piv::write_vcfs_to_stdout(vcf_recs, *call_graph, app_config);
  
  delete call_graph;
  
  // Clean up the temporary directory
  fs::remove_all(temp_dir_str);
  
  return;
}

} // namespace povu::subcommands::gfa2vcf