#include <cstdlib>
#include <fstream> // for std::ifstream
#include <liteseq/gfa.h>
#include <ostream>
#include <string>

#include "../../include/common/log.hpp"
#include "../common/compat.hpp"
#include "./to_pvst.hpp"

namespace povu::io::to_pvst
{
namespace pc = povu::constants;
namespace pu = povu::utils;

inline void write_header_line(std::ofstream &bub_file) noexcept
{
	bub_file << pc::PVST_HEADER_SYMBOL << pc::COL_SEP << pc::PVST_VERSION
		 << pc::COL_SEP << pc::NO_VALUE << pc::COL_SEP << pc::NO_VALUE
		 << pc::COL_SEP << pc::NO_VALUE << "\n";
}

void write_pvst(const pvst::Tree &bt, const std::string &base_name,
		const core::config &app_config)
{
	// TODO: combine and pass as single arg
	std::string bub_file_name = pv_cmp::format(
		"{}/{}.pvst", std::string{app_config.get_output_dir()},
		base_name); // file path and name
	std::ofstream bub_file(bub_file_name);

	if (!bub_file.is_open()) {
		ERR("Could not open file {}", bub_file_name);
		std::exit(EXIT_FAILURE);
	}

	// writer header line
	// ------------------
	write_header_line(bub_file);

	// write the rest of the PVST
	// --------------------------
	for (std::size_t i{}; i < bt.vtx_count(); ++i) {
		const pvst::VertexBase &v = bt.get_vertex(i);
		switch (v.get_fam()) { // line identifier
		case pvst::vt_e::concealed:
			bub_file << pc::PVST_CONCEALED_SYMBOL << pc::COL_SEP;
			break;
		case pvst::vt_e::flubble:
			bub_file << pc::PVST_FLUBBLE_SYMBOL << pc::COL_SEP;
			break;
		case pvst::vt_e::dummy:
			bub_file << pc::PVST_DUMMY_SYMBOL << pc::COL_SEP;
			break;
		case pvst::vt_e::tiny:
			bub_file << pc::PVST_TINY_SYMBOL << pc::COL_SEP;
			break;
		case pvst::vt_e::parallel:
			bub_file << pc::PVST_OVERLAP_SYMBOL << pc::COL_SEP;
			break;
		case pvst::vt_e::smothered:
			bub_file << pc::PVST_SMOTHERED_SYMBOL << pc::COL_SEP;
			break;
		case pvst::vt_e::midi:
			bub_file << pc::PVST_MIDI_SYMBOL << pc::COL_SEP;
			break;
		default:
			ERR("Unknown vertex type in write_bub: {}", v.as_str());
			std::exit(EXIT_FAILURE);
		}

		// vertex idx
		bub_file << i << pc::COL_SEP;

		// vertex as str
		bub_file << v.as_str() << pc::COL_SEP;

		// children
		if (bt.is_leaf(i))
			bub_file << pc::NO_VALUE;
		else
			pu::print_with_comma(bub_file, bt.get_children(i), ',');

		bub_file << pc::COL_SEP;

		// route
		if (std::optional<pvst::route_params_t> rp =
			    v.get_route_params();
		    rp) {
			const auto &[_, __, rt] = *rp;
			bub_file << to_str(rt);
		}
		else {
			bub_file << pc::NO_VALUE;
		}

		bub_file << "\n";
	}

	bub_file.close();
}
} // namespace povu::io::to_pvst
