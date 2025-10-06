#include <gtest/gtest.h>
// #include <vector>
// #include "povu/algorithms/flubbles.hpp"
// #include "povu/common/app.hpp"

#include "povu/graph/bidirected.hpp"
#include "povu/graph/pvst.hpp"
#include "povu/graph/spanning_tree.hpp"

namespace povu::unit_tests_spanning_tree
{
// namespace pfl = povu::flubbles;

bd::VG *create_test_vg()
{
	bd::VG *vg = new bd::VG(7, 9, 0); // 7 vertices, 10 edges, 0 references

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
	vg->add_edge(4, bd::v_end_e::r, 5, bd::v_end_e::l);
	vg->add_edge(4, bd::v_end_e::r, 6, bd::v_end_e::l);
	vg->add_edge(4, bd::v_end_e::l, 7, bd::v_end_e::r);
	vg->add_edge(5, bd::v_end_e::r, 7, bd::v_end_e::l);
	vg->add_edge(6, bd::v_end_e::r, 7, bd::v_end_e::l);

	return vg;
}

TEST(SpanningTreeTest, ConstructSpanningTree)
{
	bd::VG *vg = create_test_vg();
	pst::Tree st = pst::Tree::from_bd(*vg);
	delete vg;
}

} // namespace povu::unit_tests_spanning_tree
