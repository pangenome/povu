#include <chrono>
#include <filesystem>
#include <fstream>
#include <memory>
#include <set>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "ita/genomics/genomics.hpp"
#include "mto/from_gfa.hpp"
#include "povu/common/app.hpp"
#include "povu/common/bounded_queue.hpp"
#include "povu/common/constants.hpp"
#include "povu/graph/pvst.hpp"
#include "povu/graph/types.hpp"

namespace
{
std::filesystem::path write_sne_subr_fixture()
{
	const std::filesystem::path fp =
		std::filesystem::temp_directory_path() /
		("povu_sne_subr_chunk_gate_" +
		 std::to_string(std::chrono::steady_clock::now()
					.time_since_epoch()
					.count()) +
		 ".gfa");

	std::ofstream out(fp);
	out << "H\tVN:Z:1.0\n";
	out << "S\t1\tA\n";
	out << "S\t2\tC\n";
	out << "S\t3\tG\n";
	out << "S\t4\tT\n";
	out << "S\t5\tA\n";
	out << "L\t1\t+\t2\t+\t0M\n";
	out << "L\t2\t+\t3\t+\t0M\n";
	out << "L\t3\t+\t4\t+\t0M\n";
	out << "L\t4\t+\t5\t+\t0M\n";
	out << "L\t5\t-\t4\t-\t0M\n";
	out << "L\t4\t-\t3\t-\t0M\n";
	out << "L\t3\t-\t2\t-\t0M\n";
	out << "L\t2\t-\t1\t-\t0M\n";
	out << "P\tref\t1+,2+,3+,4+,5+\t*\n";
	out << "P\talt\t5-,4-,3-,2-,1-\t*\n";
	return fp;
}

core::config sne_subr_config(const std::filesystem::path &gfa_fp)
{
	core::config app_config;
	app_config.set_input_gfa(gfa_fp.string());
	app_config.set_inc_vtx_labels(true);
	app_config.set_inc_refs(true);
	app_config.set_chunk_size(100);
	app_config.set_queue_len(4);
	app_config.set_thread_count(1);
	return app_config;
}

pvst::Tree sne_subr_pvst()
{
	namespace pgt = povu::types::graph;

	pvst::Tree pvst;
	const pt::idx_t root_idx = pvst.add_vertex(pvst::Dummy{});
	pvst.set_root_idx(root_idx);

	pvst::Flubble flubble = pvst::Flubble::create(
		pgt::id_or_t{1, pgt::or_e::forward},
		pgt::id_or_t{5, pgt::or_e::forward}, pc::INVALID_IDX,
		pc::INVALID_IDX);
	const pt::idx_t flubble_idx = pvst.add_vertex(flubble);
	pvst.add_edge(root_idx, flubble_idx);
	pvst.comp_heights();

	return pvst;
}
} // namespace

TEST(GenomicsTest, RunsSneForFinalPartialChunkAndEmitsSingleSubr)
{
	const std::filesystem::path gfa_fp = write_sne_subr_fixture();
	core::config app_config = sne_subr_config(gfa_fp);
	std::unique_ptr<bd::VG> graph(mto::from_gfa::to_bd(app_config));

	std::vector<pvst::Tree> pvsts;
	pvsts.emplace_back(sne_subr_pvst());

	pbq::bounded_queue<iv::VcfRecIdx> q(app_config.get_queue_len());
	ig::gen_vcf_rec_map(pvsts, *graph, std::set<pt::id_t>{0}, q,
			    app_config);

	std::size_t subr_count = 0;
	while (auto opt_rec_idx = q.pop()) {
		for (auto &[_, recs] : opt_rec_idx->get_recs_mut()) {
			for (const iv::VcfRec &rec : recs) {
				if (rec.get_var_type() == ir::var_type_e::subr)
					++subr_count;
			}
		}
	}

	std::filesystem::remove(gfa_fp);
	EXPECT_EQ(subr_count, 1U);
}
