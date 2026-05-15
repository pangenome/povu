#include "mto/to_gfa.hpp"

#include <filesystem> // for path
#include <fstream>    // for basic_ofstream, operator<<, basic_ostream

#include "povu/common/core.hpp"
#include "povu/common/log.hpp"
#include "povu/graph/bidirected.hpp"
#include "povu/graph/types.hpp"

namespace mto::to_gfa
{
void write_gfa(const bd::VG &g, const std::filesystem::path &fp)
{
	std::ofstream os(fp);
	if (!os.is_open())
		PL_ERR("Could not open file {} for writing", fp.string());

	auto vxt_pair_to_edge_pair = [&](const bd::Edge &e)
		-> std::tuple<pt::idx_t, std::string, pt::idx_t, std::string>
	{
		pt::idx_t v1_idx = e.get_v1_idx();
		pt::idx_t v2_idx = e.get_v2_idx();
		if (v1_idx == v2_idx) // self-loop
			return std::make_tuple(v1_idx, "+", v2_idx, "+");

		std::string v1_e =
			e.get_v1_end() == ptg::v_end_e::r ? "+" : "-";
		std::string v2_e =
			e.get_v2_end() == ptg::v_end_e::l ? "+" : "-";

		return std::make_tuple(v1_idx, v1_e, v2_idx, v2_e);
	};

	/* header */
	os << pv_cmp::format("H\tVN:Z:1.0\n");

	/* vertices */
	// A is a placeholder sequence
	for (size_t v_idx{}; v_idx < g.vtx_count(); ++v_idx)
		os << pv_cmp::format("S\t{}\tA\n", g.v_idx_to_id(v_idx));

	/* edges */
	for (pt::u32 e_idx{}; e_idx < g.edge_count(); ++e_idx) {
		const bd::Edge &e = g.get_edge(e_idx);

		auto [v1_idx, v1_e, v2_idx, v2_e] = vxt_pair_to_edge_pair(e);

		os << pv_cmp::format("L\t{}\t{}\t{}\t{}\t0M\n",
				     g.v_idx_to_id(v1_idx), v1_e,
				     g.v_idx_to_id(v2_idx), v2_e);
	}
}
} // namespace mto::to_gfa
