#include <algorithm> // for find
#include <array>
#include <cctype>  // for isspace
#include <cstdlib> // for exit, size_t, strtol
#include <fstream> // for basic_ifstream, basic_ostream
#include <map>	   // for map
#include <memory>  // for make_unique
#include <sstream> // for basic_stringstream
#include <string>  // for basic_string, char_traits, string
#include <string_view>
#include <utility> // for get, pair
#include <vector>  // for vector

#include <log/log.h>		 // for log_fatal, log_fatal
#include <quilt/graph_types.hpp> // for v_end_e, side_n_id_t, side_n_idx_t
#include <quilt/types.hpp>	 // for qt
#include <quilt/utils.hpp>

#include "mto/common.hpp" // for FILE_ERROR
#include "mto/from_pvst.hpp"

namespace mto::from_pvst
{
namespace pgt = quilt::types::graph;

const std::vector<std::string> PVST_SUPPORTED_VERSIONS = {"0.0.3"};

/**
 * @brief Get the size of a file
 * @param fp file path
 * @return size of the file in bytes
 */
std::size_t get_file_size(const std::string &fp)
{
	std::streampos begin, end;
	std::ifstream f(fp);

	if (!f)
		FILE_ERROR(fp);

	begin = f.tellg();
	f.seekg(0, std::ios::end);
	end = f.tellg();
	f.close();

	return end - begin;
}

template <typename Iterable, typename U>
bool in_vector(const Iterable &c, const U &e)
{
	return !(std::find(c.begin(), c.end(), e) == c.end());
}

void check_pvst_version_support(const std::string &fp,
				const std::string &pvst_version)
{
	if (!in_vector(PVST_SUPPORTED_VERSIONS, pvst_version)) {
		std::string err = qs::format(
			"Unsupported PVST version in file {}, got {}. "
			"Supported "
			"versions are: {}.",
			fp, pvst_version,
			pu::concat_with(PVST_SUPPORTED_VERSIONS, ','));
		log_fatal("%s", err.c_str());

		std::exit(EXIT_FAILURE);
	}

	return;
}

// Split on commas, trim whitespace, parse each token as an integer.
// Ignores empty tokens (e.g. a trailing comma).
std::vector<qt::idx_t> split_numbers(std::string_view s)
{
	std::vector<qt::idx_t> result;
	size_t pos = 0;

	while (pos < s.size()) {
		// find the next comma (or end)
		size_t comma = s.find(',', pos);

		// token = [pos, comma)
		auto token = s.substr(pos, comma - pos);

		// trim leading/trailing spaces
		size_t first = 0;
		while (first < token.size() && std::isspace(token[first]))
			++first;
		size_t last = token.size();
		while (last > first && std::isspace(token[last - 1]))
			--last;

		if (last > first) {
			// parse token[first..last)
			// use strtol on a null-terminated buffer
			std::string buf(token.substr(first, last - first));
			char const *ptr = buf.c_str();
			char *end = nullptr;
			long val = std::strtol(ptr, &end, 10);
			if (end != ptr)
				result.push_back(static_cast<int>(val));

			// else: token wasn’t a valid number—skip it
		}

		if (comma == std::string_view::npos)
			break; // no more commas

		// move to the next token
		pos = comma + 1; // skip the comma
	}

	return result;
}

// TODO: [c] CLEANUP move to utils
/**
 * @brief split s based on > and < signs and using s.substr
 *        s is in the form of >1<2 or >1>2 or <1<2 or <1>2
 */
std::pair<pgt::id_or_t, pgt::id_or_t> str_to_id_or_t(const std::string &s)
{
	// find the first > or <
	auto first = s.find_first_of("><");
	auto last = s.find_last_of("><");

	// substring based on first and last occurences and store them as size_t
	pgt::id_or_t srt, end;

	srt.v_id = std::stoull(s.substr(first + 1, last - first - 1));
	srt.orientation =
		s[first] == '>' ? pgt::or_e::forward : pgt::or_e::reverse;

	end.v_id = std::stoull(s.substr(last + 1, s.size() - last - 1));
	end.orientation =
		s[last] == '>' ? pgt::or_e::forward : pgt::or_e::reverse;

	return {srt, end};
}

pvst::route_params_t
tokens_to_route_params(const std::vector<std::string> &tokens)
{
	const std::string &pvst_label = tokens[2];
	auto [l, r] = str_to_id_or_t(pvst_label);

	const char route_char = tokens[4][0];
	pvst::route_e route =
		route_char == 'L' ? pvst::route_e::s2e : pvst::route_e::e2s;
	return pvst::route_params_t{l, r, route};
}

pvst::Tree read_pvst(const std::string &fp)
{
	// bool dbg = "frst_dir/9.pvst" == fp ? true : false;

	// std::cerr << "Reading PVST file: " << fp << "\n";

	pvst::Tree pvst;

	// lines in the PVST
	std::vector<std::string> lines;
	mto::common::fp_to_vector(fp, &lines);

	std::vector<std::string> tokens;

	// a map from line in the .pvst file to the vertex index in the PVST
	std::map<qt::idx_t, qt::idx_t> line_idx_to_pvst_idx;

	// map file vertex index to pvst vertex index
	std::map<qt::idx_t, qt::idx_t> file_v_idx_to_pvst_idx;

	// qt::u32 x;

	for (qt::idx_t line_idx{}; line_idx < lines.size(); line_idx++) {
		const std::string &line = lines[line_idx];

		pu::split(line, pc::COL_SEP, &tokens);

		if (tokens.size() != pc::EXPECTED_PVST_COL_NUMS) {
			std::stringstream err_msg;
			err_msg << "invalid number of columns. "
				<< "File " << fp << ", line " << line_idx
				<< ". Expected " << pc::EXPECTED_PVST_COL_NUMS
				<< ", got " << tokens.size() << '\n';
			log_fatal("%s", err_msg.str().c_str());
			std::exit(EXIT_FAILURE);
		}

		char typ = tokens[0][0];
		// id of the vertex in the file
		// assume that the value here will always be numeric
		qt::idx_t id = std::stoul(tokens[1]);

		qt::idx_t v_idx{pc::INVALID_IDX};
		// const std::string &pvst_label = tokens[2];

		switch (typ) {
		case pc::PVST_DUMMY_SYMBOL: {
			pvst::Dummy root_v;
			v_idx = pvst.add_vertex(root_v);
			pvst.set_root_idx(v_idx);
			file_v_idx_to_pvst_idx[id] = v_idx;
			break;
		}
		case pc::PVST_FLUBBLE_SYMBOL: {
			auto rp = tokens_to_route_params(tokens);
			pvst::Flubble v =
				pvst::Flubble::parse(pvst::vf_e::flubble, rp);
			v_idx = pvst.add_vertex(v);
			file_v_idx_to_pvst_idx[id] = v_idx;
			break;
		}
		case pc::PVST_TINY_SYMBOL: {
			auto rp = tokens_to_route_params(tokens);
			pvst::Flubble v =
				pvst::Flubble::parse(pvst::vf_e::tiny, rp);

			v_idx = pvst.add_vertex(v);
			file_v_idx_to_pvst_idx[id] = v_idx;
			break;
		}
		case pc::PVST_OVERLAP_SYMBOL: {
			auto rp = tokens_to_route_params(tokens);
			pvst::Flubble v =
				pvst::Flubble::parse(pvst::vf_e::parallel, rp);
			v_idx = pvst.add_vertex(v);
			file_v_idx_to_pvst_idx[id] = v_idx;
			break;
		}
		case pc::PVST_MIDI_SYMBOL: {
			auto rp = tokens_to_route_params(tokens);
			pvst::MidiBubble v = pvst::MidiBubble::parse(rp);
			v_idx = pvst.add_vertex(v);
			file_v_idx_to_pvst_idx[id] = v_idx;
			break;
		}
		case pc::PVST_CONCEALED_SYMBOL: {
			auto rp = tokens_to_route_params(tokens);
			pvst::Concealed v = pvst::Concealed::parse(rp);
			v_idx = pvst.add_vertex(v);
			file_v_idx_to_pvst_idx[id] = v_idx;
			break;
		}
		case pc::PVST_SMOTHERED_SYMBOL: {
			auto rp = tokens_to_route_params(tokens);
			pvst::Smothered v = pvst::Smothered::parse(rp);
			v_idx = pvst.add_vertex(v);
			file_v_idx_to_pvst_idx[id] = v_idx;
			break;
		}
		case pc::PVST_HEADER_SYMBOL: {
			check_pvst_version_support(
				fp, tokens[1]); // may cause the prog to exit
			break;
		}
		default: {
			std::string err = qs::format(
				"Unknown vertex type in PVST file {}: L:{}", fp,
				line_idx + 1);
			log_error("%s", err.c_str());
			break;
		}
		}

		if (v_idx != pc::INVALID_IDX)
			line_idx_to_pvst_idx[line_idx] = v_idx;

		tokens.clear();
	}

	for (qt::idx_t line_idx{}; line_idx < lines.size(); line_idx++) {
		tokens.clear();
		// TODO: [B] PERFORMANCE use stringview instead of creating a
		// new string
		const std::string &line = lines[line_idx];

		pu::split(line, pc::COL_SEP, &tokens);

		const std::string &ch = tokens[3];

		if (ch == ".") // no children
			continue;

		qt::idx_t p_pvst_idx = line_idx_to_pvst_idx[line_idx];
		for (const qt::idx_t &c_idx : split_numbers(ch)) {
			// qt::idx_t ch_pvst_idx = line_idx_to_pvst_idx[c_idx];
			qt::idx_t ch_pvst_idx = file_v_idx_to_pvst_idx[c_idx];

			pvst.add_edge(p_pvst_idx, ch_pvst_idx);
		}
	}

	return pvst;
}

} // namespace mto::from_pvst
