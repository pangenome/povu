#ifndef PV_ALGOS_HPP
#define PV_ALGOS_HPP

//#include <vector>

//#include "../common/types.hpp"
#include "../graph/spanning_tree.hpp"
//#include "../graph/tree.hpp"
#include "WFAligner.hpp"

namespace povu::algorithms {
namespace pst = povu::spanning_tree;


/**
 * @brief
 */
void eulerian_cycle_equiv(pst::Tree &t);

} // namespace povu::algorithms

namespace povu::align {

std::string wfa2(wfa::WFAlignerGapAffine& aligner, const std::string& query, const std::string& text);
} // namespace align



#endif
