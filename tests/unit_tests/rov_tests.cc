#include <gtest/gtest.h>
#include <vector>

#include "povu/genomics/graph.hpp"
#include "povu/graph/types.hpp"
#include "povu/variation/rov.hpp"

namespace povu::unit_tests_rov
{
namespace pgg = povu::genomics::graph;

constexpr pvr::var_type_e ins = pvr::var_type_e::ins;
constexpr pvr::var_type_e del = pvr::var_type_e::del;
constexpr pvr::var_type_e sub = pvr::var_type_e::sub;

constexpr ptg::or_e fwd = ptg::or_e::forward;
constexpr ptg::or_e f = ptg::or_e::forward;
constexpr ptg::or_e rev = ptg::or_e::reverse;
constexpr ptg::or_e r = ptg::or_e::reverse;

TEST(FindHidden, Substitution)
{
	pvr::RoV r{nullptr};
	ptg::walk_t w1{{1, fwd}, {2, fwd}, {4, fwd}};
	ptg::walk_t w2{{1, fwd}, {3, fwd}, {4, fwd}};
	r.set_walks({w1, w2});

	// pvr::find_hidden(r);
}

TEST(FindHidden, SubNarrowDown)
{
	pvr::RoV r{nullptr};
	ptg::walk_t w1{{3, fwd}, {4, fwd}, {5, fwd}, {8, fwd}};
	ptg::walk_t w2{{3, fwd}, {4, fwd}, {7, fwd}, {8, fwd}};
	r.set_walks({w1, w2});

	std::vector<pvr::raw_variant> exp_rvs{
		{{2, 1}, {2, 1}, sub},
	};
	const pt::u32 EXP_VAR_COUNT = exp_rvs.size();

	// pvr::find_hidden(r);

	// const std::vector<pvr::raw_variant> &pv =
	//	r.get_irreducibles().at(0).variants;

	// ASSERT_EQ(pv.size(), EXP_VAR_COUNT);

	// for (pt::u32 i{}; i < EXP_VAR_COUNT; i++) {
	//	const pvr::raw_variant &found_rv = pv[i];
	//	const pvr::raw_variant &exp_rv = exp_rvs[i];

	//	EXPECT_EQ(found_rv, exp_rv);
	// }
}

TEST(FindHidden, Insertion)
{
	pvr::RoV r{nullptr};
	ptg::walk_t w1{{1, fwd}, {2, fwd}, {3, fwd}};
	ptg::walk_t w2{{1, fwd}, {3, fwd}};
	r.set_walks({w1, w2});

	// pvr::find_hidden(r);
}

TEST(FindHidden, Deletion)
{
	pvr::RoV r{nullptr};
	ptg::walk_t w1{{1, fwd}, {3, fwd}};
	ptg::walk_t w2{{1, fwd}, {2, fwd}, {3, fwd}};
	r.set_walks({w1, w2});

	// pvr::find_hidden(r);

	pvr::RoV r2{nullptr};
	w1 = {{8, f}, {14, f}};
	w2 = {{8, f}, {11, f}, {12, f}, {13, f}, {14, f}};
	r2.set_walks({w1, w2});

	// pvr::find_hidden(r2);
}

TEST(FindHidden, NarrowDown)
{
	pvr::RoV r{nullptr};
	ptg::walk_t w1{{1, fwd}, {2, fwd}, {3, fwd}, {5, fwd}};
	ptg::walk_t w2{{1, fwd}, {2, fwd}, {4, fwd}, {5, fwd}};
	r.set_walks({w1, w2});

	// pvr::find_hidden(r);
}

TEST(FindHidden, IndelsAndNarrowDown)
{
	pvr::RoV r{nullptr};
	ptg::walk_t w1{{1, fwd}, {2, fwd}, {3, fwd}, {5, fwd}};
	ptg::walk_t w2{{1, fwd}, {3, fwd}, {5, fwd}};
	r.set_walks({w1, w2});

	// pvr::find_hidden(r);
}

TEST(FindHidden, IndelSub)
{
	pvr::RoV r{nullptr};
	ptg::walk_t w1{{1, f}, {2, f}, {4, f}, {5, f}, {6, f}, {7, f}};
	ptg::walk_t w2{{1, f}, {3, f}, {7, f}};
	r.set_walks({w1, w2});

	// pvr::find_hidden(r);

	pvr::RoV r2{nullptr};
	w1 = {{1, f}, {2, f}, {3, f}, {7, f}, {8, f}};
	w2 = {{1, f}, {3, f}, {8, f}};
	r2.set_walks({w1, w2});

	// pvr::find_hidden(r2);
}

TEST(FindHidden, IndelSubNarrowDown)
{
	pvr::RoV r{nullptr};
	ptg::walk_t w1{{97, f}, {98, f}, {2, f}, {4, f},  {5, f},  {6, f},
		       {7, f},	{8, f},	 {9, f}, {12, f}, {13, f}, {15, f}};

	ptg::walk_t w2{{97, f}, {98, f}, {8, f},  {9, f},
		       {11, f}, {12, f}, {13, f}, {15, f}};
	r.set_walks({w1, w2});

	// pvr::find_hidden(r);
}

TEST(FindHidden, IndelsAndMismatchNarrowDown)
{
	pvr::RoV r{nullptr};
	ptg::walk_t w1{{1, fwd}, {2, fwd}, {4, fwd}, {5, fwd}};
	ptg::walk_t w2{{1, fwd}, {3, fwd}, {5, fwd}};
	r.set_walks({w1, w2});

	// pvr::find_hidden(r);
}
} // namespace povu::unit_tests_rov
