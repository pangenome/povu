#include <vector>

#include "../graph/bidirected.hpp"



namespace genomics {
using namespace common::typedefs;
using namespace graph_types;

// TODO: make use of this or delete
enum variant_type {
    SNP,
    DEL,
    INS,
    INV,
    DUP,
    CNV,
    BND
};

// TODO which version of VCF is best?
enum output_format {
    VCF, //  currently outputs v4.2
    PAF, // not yet supported
};

void call_variants(const std::vector<canonical_sese>& canonical_flubbles, const bidirected::VariationGraph& bd_vg, const core::config& app_config);

} // namespace genomics
