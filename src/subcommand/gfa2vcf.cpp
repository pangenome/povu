#include "./gfa2vcf.hpp"
#include "./common.hpp"
#include "./call.hpp"
#include "./deconstruct.hpp"

#include <filesystem>
#include <iostream>
#include <cstdlib>

#include "../../include/graph/bidirected.hpp"
#include "../../include/common/types/pvst.hpp"
#include "../../include/genomics/variants.hpp"
#include "../io/pvst.hpp"
#include "../io/to_vcf.hpp"
#include "../io/common.hpp"

namespace fs = std::filesystem;
namespace bd = povu::bidirected;
namespace pvtr = povu::tree;
namespace pcs = povu::subcommands::common;
namespace pic = povu::io::common;
namespace pvt = povu::types::genomics;
namespace pt = povu::types;
namespace pg = povu::variants;
namespace piv = povu::io::to_vcf;

namespace povu::subcommands::gfa2vcf {

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
  
  // Step 1: Run decompose to generate forest of PVST files
  if (ll > 0) {
    std::cerr << fn_name << " Step 1: Decomposing graph..." << std::endl;
  }
  
  // Create a modified config for decompose with the temp directory
  core::config decompose_config = app_config;
  decompose_config.set_task(core::task_e::deconstruct);
  decompose_config.set_output_dir(temp_dir_str);
  
  // Reuse the existing deconstruct function
  deconstruct::do_deconstruct(decompose_config);
  
  // Step 2: Run call to generate VCF and output to stdout
  if (ll > 0) {
    std::cerr << fn_name << " Step 2: Calling variants..." << std::endl;
  }
  
  // Read PVST files from temp directory
  std::vector<pvtr::Tree> pvsts;
  call::read_pvsts(decompose_config, pvsts);
  
  // Get references
  core::config call_config = app_config;
  call::get_refs(call_config);
  
  // Load graph for variant calling
  bd::VG *g = pcs::get_vg(call_config);
  
  // Generate VCF records
  pvt::VcfRecIdx vcf_recs = pg::gen_vcf_rec_map(pvsts, *g);
  
  // Write VCF to stdout
  piv::write_vcfs_to_stdout(vcf_recs, *g, call_config);
  
  delete g;
  
  // Clean up the temporary directory
  fs::remove_all(temp_dir_str);
  
  return;
}

} // namespace povu::subcommands::gfa2vcf