#include "./pvst.hpp"

namespace povu::types::pvst {
traversal_params_t null_tp() {
  return traversal_params_t{pgt::id_or_t{pc::INVALID_ID, pgt::or_e::forward},
                            pgt::id_or_t{pc::INVALID_ID, pgt::or_e::forward},
                            false, true};
}

std::ostream &operator<<(std::ostream &os, vt_e t) {
  switch (t) {
  case vt_e::dummy:
    os << "dummy";
    break;
  case vt_e::flubble:
    os << "flubble";
    break;
  case vt_e::tiny:
    os << "tiny";
    break;
  case vt_e::parallel:
    os << "parallel";
    break;
  case vt_e::slubble:
    os << "slubble";
    break;
  case vt_e::smothered:
    os << "smothered";
    break;
  case vt_e::midi:
    os << "midi";
    break;
  }
  return os;
}
} // namespace povu::types::pvst
