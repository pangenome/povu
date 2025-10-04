#ifndef IO_HPP
#define IO_HPP

#include <filesystem>   // for path
#include <string>       // for string, basic_string, operator+
#include <vector>       // for vector
#include <stdexcept>    // for invalid_argument
#include <string_view>  // for string_view

namespace povu::io::common
{
inline constexpr std::string_view MODULE = "povu::io::common";
namespace fs = std::filesystem;

#define FILE_ERROR(name)                                                       \
	throw std::invalid_argument(std::string("Failed to open file ") + name);

/**
 * @brief Get the list of files in a dir with a given name
 *
 */
std::vector<fs::path> get_files(const std::string &dir_path,
				const std::string &ext);

void read_lines_to_vec_str(const std::string &fp, std::vector<std::string> *v);

void create_dir_if_not_exists(const fs::path &out_dir);
}; // namespace povu::io::common

#endif
