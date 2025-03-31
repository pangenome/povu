#include <string>
#include "WFAligner.hpp"

#include "./align.hpp"

namespace povu::align {

std::string wfa2(wfa::WFAlignerGapAffine& aligner, const std::string& query, const std::string& text) {
  aligner.alignEnd2End(query, text); // Align
  return aligner.getAlignment(); // edit transcript
}

} // namespace align
