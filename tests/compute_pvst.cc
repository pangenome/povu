#include <gtest/gtest.h>
#include <iostream>

#include "../src/pvst/pvst.hpp"
#include "../src/core/core.hpp"

// global config (test config)
core::config test_config;

TEST(PVSTTest, AdjacentFlubbles) {
  std::vector<std::pair<std::size_t, std::size_t>> v = {
		{2, 0},
		{4, 1},
		{6, 2},
		{8, 0},
		{10, 3},
		{12, 4},
		{14, 0},
		{15, 5},
		{17, 6},
		{19, 0},
  };

  tree::Tree t = pvst::compute_pvst(v, test_config);
  t.print_dot(true);
  EXPECT_EQ(1, 1);
}

TEST(PVSTTest, NestedFlubbles) {

}

//
TEST(PVSTTest, Redunduncies) {
  std::vector<std::pair<std::size_t, std::size_t>> v = {
	{2, 7},
	{4, 8},
	{6, 9},
	{8, 9},
	{10, 10},
	{11, 0},
	{14, 9},
	{16, 8},
	{17, 3},
	{19, 3},
	{22, 2},
	{23, 4},
	{26, 1},
	{27, 3},
	{29, 3},
	{32, 7},
	{34, 7},
  };

  tree::Tree t = pvst::compute_pvst(v, test_config);
    t.print_dot(true);
  EXPECT_EQ(1, 1);
}
