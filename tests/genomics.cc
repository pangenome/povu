#include "../src/core/core.hpp"
#include "../src/genomics/genomics.hpp"


TEST(GenomicsTest, ExtractCanonicalFlubbles) {
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
  std::vector<std::pair<std::size_t, std::size_t>> cfl =
	genomics::extract_canonical_flubbles(t);

  // print cfl
  for (auto flubble : cfl) {
	std::cout << flubble.first << " " << flubble.second << std::endl;
	}
  EXPECT_EQ(1, 2);
}
