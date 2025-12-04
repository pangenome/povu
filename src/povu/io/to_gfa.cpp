#include "povu/io/to_gfa.hpp"
#include "povu/common/constants.hpp"
#include "povu/common/core.hpp"
#include "povu/graph/bidirected.hpp"
#include "povu/graph/types.hpp"
#include <filesystem>
#include <fstream> // for basic_ofstream, operator<<, bas...

namespace povu::io::to_gfa
{
void write_gfa(const bd::VG &g, const std::filesystem::path &fp)
{

	std::ofstream os(fp);

	/* helper fns */
	// map v end left and right to dot west and east for rectangular
	// vertices
	auto v_end_to_gfa = [](ptg::v_end_e e) -> std::string
	{
		return e == ptg::v_end_e::r ? "+" : "-";
	};

	auto v_end_to_gfa2 = [](ptg::v_end_e e) -> std::string
	{
		return e == ptg::v_end_e::l ? "+" : "-";
	};

	/* header */
	os << "H" << "\t" << "VN:Z:1.0\n";

	/* vertices */
	for (size_t v_idx{}; v_idx < g.vtx_count(); ++v_idx) {
		// const bd::Vertex &v = g.get_vertex_by_idx(v_idx);

		os << "S" << "\t" << g.v_idx_to_id(v_idx) << "\t" << "A"
		   << "\n";
	}

	/* edges */
	for (pt::u32 e_idx{}; e_idx < g.edge_count(); ++e_idx) {
		const bd::Edge &e = g.get_edge(e_idx);

		pt::idx_t v1_idx = e.get_v1_idx();
		std::string v1_e = v_end_to_gfa(e.get_v1_end());
		pt::idx_t v2_idx = e.get_v2_idx();
		std::string v2_e = v_end_to_gfa2(e.get_v2_end());

		os << "L" << "\t" << g.v_idx_to_id(v1_idx) << "\t" << v1_e
		   << "\t" << g.v_idx_to_id(v2_idx) << "\t" << v2_e << "\t"
		   << "0M"
		   << "\n";

		// os << pv_cmp::format("L{}\t{}\t{}\t{}\n", v1_idx, v1_e,
		// v2_idx,		     v2_e);
	}

	/* footer */
	// os << "}" << std::endl;
}
} // namespace povu::io::to_gfa
