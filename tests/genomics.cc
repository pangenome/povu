#include "../src/common/common.hpp"
#include "../src/genomics/genomics.hpp"


TEST(GenomicsTest, ExtractCanonicalFlubbles) {
  std::vector<graph_types::eq_n_id_t> v = {
    {7, 2},
    {8, 4},
    {9, 6},
    {9, 8},
    {10, 10},
    {0, 11},
    {9, 14},
    {8, 16},
    {3, 17},
    {3, 19},
    {2, 22},
    {4, 23},
    {1, 26},
    {3, 27},
    {3, 29},
    {7, 32},
    {7, 34},
    };


  //tree::Tree t = pvst::compute_pvst(v, test_config);
  //std::vector<std::pair<std::size_t, std::size_t>> cfl = genomics::extract_canonical_flubbles(t);

  // print cfl
  //for (auto flubble : cfl) {
  //std::cout << flubble.first << " " << flubble.second << std::endl;
  //}
  //EXPECT_EQ(1, 2);
}
