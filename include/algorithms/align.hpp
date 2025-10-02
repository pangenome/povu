#ifndef PV_ALN_HPP
#define PV_ALN_HPP

#include "WFAligner.hpp"

namespace povu::align
{
std::string wfa2(wfa::WFAlignerGapAffine &aligner, const std::string &query,
		 const std::string &text);
} // namespace povu::align

#endif // PV_ALN_HPP
