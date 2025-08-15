#include "./from_pvst.hpp"
#include <sstream>

namespace povu::io::from_pvst {

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

// Split on commas, trim whitespace, parse each token as an integer.
// Ignores empty tokens (e.g. a trailing comma).
std::vector<pt::idx_t> split_numbers(std::string_view s) {
  std::vector<pt::idx_t> result;
  size_t pos = 0;

  while (pos < s.size()) {
    // find the next comma (or end)
    size_t comma = s.find(',', pos);

    // token = [pos, comma)
    auto token = s.substr(pos, comma - pos);

    // trim leading/trailing spaces
    size_t first = 0;
    while (first < token.size() && std::isspace(token[first])) ++first;
    size_t last = token.size();
    while (last > first && std::isspace(token[last-1])) --last;

    if (last > first) {
      // parse token[first..last)
      // use strtol on a null-terminated buffer
      std::string   buf(token.substr(first, last-first));
      char const*   ptr = buf.c_str();
      char*         end = nullptr;
      long          val = std::strtol(ptr, &end, 10);
      if (end != ptr) {
        result.push_back(static_cast<int>(val));
      }
      // else: token wasn’t a valid number—skip it
    }

    if (comma == std::string_view::npos) {
      break; // no more commas
    }

    // move to the next token
    pos = comma + 1;    // skip the comma
  }

  return result;
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

// TODO: [c] CLEANUP move to utils
/**
 * @brief split s based on > and < signs and using s.substr
 *        s is in the form of >1<2 or >1>2 or <1<2 or <1>2
 */
std::pair<pgt::id_or_t, pgt::id_or_t> str_to_id_or_t(const std::string &s) {
  std::string fn_name = pv_cmp::format("[povu::main::{}]", __func__);

  // find the first > or <
  auto first = s.find_first_of("><");
  auto last = s.find_last_of("><");

  // substring based on first and last occurences and store them as size_t
  pgt::id_or_t srt, end;

  srt.v_id = std::stoull(s.substr(first + 1, last - first - 1));
  srt.orientation = s[first] == '>' ? pgt::or_e::forward : pgt::or_e::reverse;

  end.v_id = std::stoull(s.substr(last + 1, s.size() - last - 1));
  end.orientation = s[last] == '>' ? pgt::or_e::forward : pgt::or_e::reverse;

  return {srt, end};
}

void fp_to_vector (const std::string& fp, std::vector<std::string>* v) {
  std::size_t file_size = get_file_size(fp);

  v->reserve(file_size);
  read_lines_to_vector_str(fp, v);
  v->shrink_to_fit();
}

pvtr::Tree read_pvst(const std::string &fp) {
  std::string fn_name = pv_cmp::format("[povu::main::{}]", __func__);

  pvtr::Tree pvst;

  // lines in the PVST
  std::vector<std::string> lines;
  fp_to_vector(fp, &lines);

  std::vector<std::string> tokens;

  // a map from line in the .pvst file to the vertex index in the PVST
  std::map<pt::idx_t, pt::idx_t> line_idx_to_pvst_idx;

  // map file vertex index to pvst vertex index
  std::map<pt::idx_t, pt::idx_t> file_v_idx_to_pvst_idx;

  // number of columns in a .pvst file
  const std::size_t PVST_COLS {4};

  for (pt::idx_t line_idx {} ; line_idx < lines.size(); line_idx++) {
    const std::string &line = lines[line_idx];

    pu::split(line, pc::COL_SEP, &tokens);

    if (tokens.size() != PVST_COLS) {
      std::stringstream err_msg;
      err_msg << "invalid number of columns. "
                << "File " << fp
                << ", line " << line_idx
                << ". Expected " << PVST_COLS
                << ", got " << tokens.size()
                << '\n';
      ERR("{}", err_msg.str());
      std::exit(1);
    }

    char typ = tokens[0][0];
    // id of the vertex in the file
    // assume that the value here will always be numeric
    pt::idx_t id = std::stoul(tokens[1]);

    pt::idx_t v_idx{pc::INVALID_IDX};
    const std::string &pvst_label = tokens[2];

    switch (typ) {
    case pc::PVST_DUMMY_SYMBOL: {
      pvst::Dummy root_v;
      v_idx = pvst.add_vertex(root_v);
      pvst.set_root_idx(v_idx);
      file_v_idx_to_pvst_idx[id] = v_idx;
      break;
    }
    case pc::PVST_FLUBBLE_SYMBOL: {
      auto [a, z] = str_to_id_or_t(pvst_label);
      pvst::Flubble v(pvst::vt_e::flubble, a, z);
      v_idx = pvst.add_vertex(v);
      file_v_idx_to_pvst_idx[id] = v_idx;
      break;
    }
    case pc::PVST_TINY_SYMBOL: {
      auto [a, z] = str_to_id_or_t(pvst_label);
      pvst::Flubble v(pvst::vt_e::tiny, a, z);
      v_idx = pvst.add_vertex(v);
      file_v_idx_to_pvst_idx[id] = v_idx;
      break;
    }
    case pc::PVST_OVERLAP_SYMBOL: {
      auto [a, z] = str_to_id_or_t(pvst_label);
      pvst::Flubble v(pvst::vt_e::parallel, a, z);
      v_idx = pvst.add_vertex(v);
      file_v_idx_to_pvst_idx[id] = v_idx;
      break;
    }
      // case pc::PVST_MIDI_SYMBOL: {
      //   auto [g, s] = str_to_id_or_t(pvst_label);
      //   pvst::MidiBubble v(g, s);
      //   v_idx = pvst.add_vertex(v);
      //   file_v_idx_to_pvst_idx[id] = v_idx;
      //   break;
      // }
      // case pc::PVST_CONCEALED_SYMBOL: {
      //   auto [f, s] = str_to_id_or_t(pvst_label);
      //   pvst::Concealed v(f, s);
      //   v_idx = pvst.add_vertex(v);
      //   file_v_idx_to_pvst_idx[id] = v_idx;
      //   break;
      // }
      // case pc::PVST_SMOTHERED_SYMBOL: {
      //   auto [f, s] = str_to_id_or_t(pvst_label);
      //   pvst::Smothered v(f, s);
      //   v_idx = pvst.add_vertex(v);
      //   file_v_idx_to_pvst_idx[id] = v_idx;
      //   break;
      // }
    default:
      break;
    }

    if (v_idx != pc::INVALID_IDX) {
      line_idx_to_pvst_idx[line_idx] = v_idx;
    }

    tokens.clear();
  }


  for (pt::idx_t line_idx{}; line_idx < lines.size(); line_idx++) {
    tokens.clear();
    // TODO: [B] PERFORMANCE use stringview instead of creating a new string
    const std::string &line = lines[line_idx];

    pu::split(line, pc::COL_SEP, &tokens);

    const std::string &ch = tokens[3];

    if (ch == ".") { // no children
      continue;
    }

    pt::idx_t p_pvst_idx = line_idx_to_pvst_idx[line_idx];
    for (const pt::idx_t &c_idx : split_numbers(ch)) {
      //pt::idx_t ch_pvst_idx = line_idx_to_pvst_idx[c_idx];
      pt::idx_t ch_pvst_idx = file_v_idx_to_pvst_idx[c_idx];
      pvst.add_edge(p_pvst_idx, ch_pvst_idx);
    }
  }

  return pvst;
}

} // namespace povu::io::from_pvst}
