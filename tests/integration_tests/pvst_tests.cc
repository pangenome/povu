#include <gtest/gtest.h>
#include <vector>

#include "../../app/cli/app.hpp"
#include "../../include/algorithms/flubbles.hpp"
#include "../../include/graph/bidirected.hpp"
#include "../../include/graph/pvst.hpp"
#include "../../include/graph/spanning_tree.hpp"

namespace pfl = povu::flubbles;

bd::VG *create_test_vg()
{
	bd::VG *vg = new bd::VG(7, 10, 0); // 7 vertices, 10 edges, 0 references

	// Add vertices
	vg->add_vertex(1, "AAT");
	vg->add_vertex(2, "GTC");
	vg->add_vertex(3, "GTG");
	vg->add_vertex(4, "TA");
	vg->add_vertex(5, "AA");
	vg->add_vertex(6, "TTG");
	vg->add_vertex(7, "C");

	// Add edges
	vg->add_edge(1, bd::v_end_e::r, 3, bd::v_end_e::l);

	vg->add_edge(1, bd::v_end_e::r, 4, bd::v_end_e::l);
	vg->add_edge(2, bd::v_end_e::r, 4, bd::v_end_e::l);
	vg->add_edge(3, bd::v_end_e::r, 4, bd::v_end_e::l);
	vg->add_edge(3, bd::v_end_e::l, 2, bd::v_end_e::l);
	vg->add_edge(4, bd::v_end_e::r, 5, bd::v_end_e::l);
	vg->add_edge(4, bd::v_end_e::r, 6, bd::v_end_e::l);
	vg->add_edge(4, bd::v_end_e::l, 7, bd::v_end_e::l);
	vg->add_edge(5, bd::v_end_e::r, 6, bd::v_end_e::l);
	vg->add_edge(6, bd::v_end_e::r, 7, bd::v_end_e::l);

	return vg;
}

core::config create_test_config()
{
	core::config app_config;
	app_config.set_verbosity(0);
	app_config.set_thread_count(1);
	return app_config;
}

pvst::Tree decompose(bd::VG *g, const core::config &conf)
{
	pst::Tree st = pst::Tree::from_bd(*g);
	pvst::Tree pvst = pfl::find_flubbles(st, conf);
	return pvst;
}

pvst::Tree decompose_and_cleanup()
{
	bd::VG *vg = create_test_vg();
	std::vector<bd::VG *> components = bd::VG::componetize(*vg);
	delete vg;
	core::config conf = create_test_config();
	pvst::Tree pvst = decompose(components.at(0), conf);
	for (auto vg_component : components) {
		delete vg_component;
	}
	return pvst;
}

TEST(PVSTTest, VertexCount)
{
	pvst::Tree pvst = decompose_and_cleanup();

	const pt::idx_t EXPECTED_PVST_VTX_COUNT = 3;
	EXPECT_EQ(pvst.vtx_count(), EXPECTED_PVST_VTX_COUNT);
}

TEST(PVSTTest, HasVertices)
{
	pvst::Tree pvst = decompose_and_cleanup();

	std::set<std::string> expected_labels = {".", ">1>7", ">4>6"};

	for (pt::idx_t i = 0; i < pvst.vtx_count(); ++i) {
		const pvst::VertexBase &v = pvst.get_vertex(i);
		EXPECT_TRUE(expected_labels.find(v.as_str()) !=
			    expected_labels.end());
	}
}

TEST(PVSTTest, VertexHierarchy)
{
	pvst::Tree pvst = decompose_and_cleanup();

	std::set<std::string> expected_labels = {".", ">1>7", ">4>6"};

	for (pt::idx_t i = 0; i < pvst.vtx_count(); ++i) {
		const pvst::VertexBase &v = pvst.get_vertex(i);

		const std::vector<pt::idx_t> &children = pvst.get_children(i);

		if (v.as_str() == ".") {
			EXPECT_EQ(pvst.root_idx(), i);
			EXPECT_EQ(v.get_fam(), pvst::vf_e::dummy);
			EXPECT_EQ(children.size(), 1);
			EXPECT_EQ(pvst.get_vertex(children.front()).as_str(),
				  ">1>7");
			continue; // root has no parent
		}

		const pvst::VertexBase &parent = pvst.get_parent(i);

		if (v.as_str() == ">1>7") {
			EXPECT_EQ(v.get_fam(), pvst::vf_e::flubble);
			EXPECT_EQ(parent.as_str(), ".");
			EXPECT_EQ(children.size(), 1);
			EXPECT_EQ(pvst.get_vertex(children.front()).as_str(),
				  ">4>6");
		}
		else if (v.as_str() == ">4>6") {
			EXPECT_EQ(v.get_fam(), pvst::vf_e::flubble);
			EXPECT_EQ(parent.as_str(), ">1>7");
			EXPECT_TRUE(children.empty());
		}
		else {
			FAIL() << "Unexpected vertex label: " << v.as_str();
		}
	}
}
