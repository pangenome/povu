#include <cstdlib>
#include <format>
#include <fstream>
#include <cstddef>
#include <iostream>
#include <ostream>
#include <filesystem>

#include "io.hpp"

namespace povu::io::generic {
namespace fs = std::filesystem;

void write_txt(const std::vector<pgt::flubble>& flubbles, const std::string& base_name, const core::config& app_config) {
  std::string fn_name = std::format("[povu::io::generic::{}]", __func__);

  std::string bub_file_name = std::format("{}/{}.txt", std::string{app_config.get_output_dir()}, base_name); // file path and name
  std::ofstream bub_file(bub_file_name);

  for (auto const& [s, e] : flubbles) {
    bub_file << s << e << std::endl;
  }

  bub_file.close();
}

std::vector<fs::path> get_files(const std::string& dir_path, const std::string& ext) {

  // Check if the directory exists
  if (!fs::exists(dir_path)) {
    std::cerr << "Directory does not exist: " << dir_path << std::endl;
    exit(EXIT_FAILURE);
  }

  std::vector<fs::path> files;
  for (const auto& entry : fs::directory_iterator(dir_path)) {
    if (entry.path().extension() == ext) {
      files.push_back(entry.path().string());
    }
  }

  return files;
}

} // namespace io::generic
