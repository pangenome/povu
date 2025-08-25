#include "./common.hpp"

namespace povu::io::common {



std::vector<fs::path> get_files(const std::string& dir_path, const std::string& ext) {
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

/**
 * @brief Get the size of a file
 * @param fp file path
 * @return size of the file in bytes
 */
std::size_t get_file_size(const std::string &fp) {
  std::streampos begin, end;
  std::ifstream f(fp);

  if (!f) {
    FILE_ERROR(fp);
  }

  begin = f.tellg();
  f.seekg(0, std::ios::end);
  end = f.tellg();
  f.close();

  return end - begin;
}

/**
 * read the entire file into a string
 * should be faster for small inputs to read the entire file into a string and
 * process it at once
 * will perform whitespace normalization/formatting
 */
void read_lines_to_vec_str(const std::string &fp, std::vector<std::string> *v) {
  std::ifstream f{fp};
  std::string temp;

  if (!f) {
    FILE_ERROR(fp);
  }

  while (f >> temp) {
    v->push_back(temp);
  }
}

void fp_to_vector(const std::string &fp, std::vector<std::string> *v) {
  std::size_t file_size = get_file_size(fp);

  v->reserve(file_size);
  read_lines_to_vec_str(fp, v);
  v->shrink_to_fit();
}

} // namespace io::generic
