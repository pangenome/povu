
#include "./io.hpp"


namespace povu::io::bub {
namespace pc = povu::constants;
namespace pu = povu::utils;

/**
 * @brief Get the size of a file
 * @param fp file path
 * @return size of the file in bytes
 */
std::size_t get_file_size(const std::string& fp) {
  std::streampos begin, end;
  std::ifstream f (fp);

  if (!f) { FILE_ERROR(fp); }

  begin = f.tellg();
  f.seekg (0, std::ios::end);
  end = f.tellg();
  f.close();

  return end-begin;
}


/**
 * read the entire file into a string
 * should be faster for small inputs to read the entire file into a string and
 * process it at once
 * will perform whitespace normalization/formatting
 */
void read_lines_to_vector_str(const std::string& fp, std::vector<std::string>* v) {
  std::ifstream f{fp};
  std::string temp;

  if (!f) { FILE_ERROR(fp); }

  // read each line into a string
  while (std::getline(f, temp)) {
    v->push_back(temp);
  }
}


void fp_to_vector (const std::string& fp, std::vector<std::string>* v) {
  std::size_t file_size = get_file_size(fp);

  v->reserve(file_size);
  read_lines_to_vector_str(fp, v);
  v->shrink_to_fit();
}


std::vector<pgt::flubble> read_canonical_fl(const std::string& fp) {
  std::vector<pgt::flubble> canonical_fl;

  std::vector<std::string> lines;
  fp_to_vector(fp, &lines);

  std::vector<std::string> tokens;

  const std::size_t FL_COLS {3}; // number of columns in a .fl file

  for (const std::string& line : lines) {

    pu::split(line, pc::COL_SEP, &tokens);

    if (tokens.size() != FL_COLS) {
      std::cerr << std::format("ERROR: invalid number of columns. Expected {}, got {} in file {}\n", FL_COLS, tokens.size(), fp);
      std::exit(1);
    }

    // if it is a dummy or not a leaf
    if (tokens[1] == std::string(1, pc::NO_VALUE) || tokens[2] != std::string(1, pc::NO_VALUE)) {
      tokens.clear();
      continue;
    }

    pgt::flubble fl (tokens[1]) ;
    canonical_fl.push_back(fl);

    tokens.clear();
  }

  return canonical_fl;
}


void write_bub(const pvtr::Tree<pgt::flubble>& bt,
               const std::string& base_name,
               const core::config& app_config) {
  // TODO: combine and pass as single arg
  std::string bub_file_name = std::format("{}/{}.flb", std::string{app_config.get_output_dir()}, base_name); // file path and name
  std::ofstream bub_file(bub_file_name);

  if (!bub_file.is_open()) {
    std::cerr << "ERROR: could not open file " << bub_file_name << "\n";
    std::exit(1);
  }

  bub_file << constants::PVST_HEADER_SYMBOL << pc::COL_SEP << "0.1" << pc::COL_SEP << "." << pc::COL_SEP << "." << "\n";

  for (std::size_t i {}; i < bt.size(); ++i) {

    //std::cerr << "i: " << i << std::endl;

    bub_file << pc::PVST_FLUBBLE_SYMBOL << pc::COL_SEP;
    //bub_file << pc::COL_SEP;

    // vertex
    {
      bub_file << i;
    }
    bub_file << pc::COL_SEP;

    // range
    {
      std::optional<pgt::flubble> data = bt.get_vertex(i).get_data();
      if (data.has_value()) {
        auto [s, e] = data.value();

        bub_file << std::format("{},{}", s.as_str(), e.as_str());
      }
      else {
        bub_file << pc::NO_VALUE;
      }

    }
    bub_file << pc::COL_SEP;


    // children
    if (bt.is_leaf(i)) {
      bub_file << pc::NO_VALUE << "\n";
    }
    else {
      pu::print_with_comma(bub_file, bt.get_children(i), ',');
      bub_file << "\n";
    }
  }

  bub_file.close();
}
} // namespace bub
