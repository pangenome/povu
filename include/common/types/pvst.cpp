#include "./pvst.hpp"

namespace povu::types::pvst {
traversal_params_t null_tp() {
  return traversal_params_t{pgt::id_or_t{pc::INVALID_ID, pgt::or_e::forward},
                            pgt::id_or_t{pc::INVALID_ID, pgt::or_e::forward},
                            false, true};
}
}
