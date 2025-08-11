#include "./gfa2vcf.hpp"
#include "./call.hpp"
#include "./deconstruct.hpp"

#include <filesystem>
#include <iostream>
#include <cstdlib>

namespace fs = std::filesystem;

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
  
  // Create a config for decompose with the temp directory
  core::config decompose_config = app_config;
  decompose_config.set_task(core::task_e::deconstruct);
  decompose_config.set_output_dir(temp_dir_str);
  
  // Run decompose
  deconstruct::do_deconstruct(decompose_config);
  
  // Step 2: Run call to generate VCF
  if (ll > 0) {
    std::cerr << fn_name << " Step 2: Calling variants..." << std::endl;
  }
  
  // Create a config for call with stdout output and the temp forest directory
  core::config call_config = app_config;
  call_config.set_task(core::task_e::call);
  call_config.set_forest_dir(temp_dir_str);
  call_config.set_stdout_vcf(true);  // Always output to stdout for gfa2vcf
  
  // Run call (it will handle everything including stdout output)
  call::do_call(call_config);
  
  // Clean up the temporary directory
  fs::remove_all(temp_dir_str);
  
  return;
}

} // namespace povu::subcommands::gfa2vcf