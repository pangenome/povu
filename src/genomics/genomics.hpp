#include <vector>

#include "../core/core.hpp"
#include "../pvst/pvst.hpp"
#include "../graph/tree.hpp"
#include "../graph/digraph.hpp"
#include "../graph/bidirected.hpp"
#include "../graph/spanning_tree.hpp"

namespace graph_operations {

using namespace graph_types;

struct node {
    //std::vector<std::size_t> ids;
    std::set<std::pair<std::size_t, std::size_t>> ids;
    std::size_t id;
    std::set<std::size_t> children;
    std::size_t parent {core::constants::UNDEFINED_SIZE_T};
};


}

namespace genomics {
typedef std::pair<bidirected::VertexEnd, id_t> side_n_id_t; // TODO: replace with struct
typedef std::vector<side_n_id_t> subpath_t;
typedef std::vector<subpath_t> subpaths_t;

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

void call_variants(const tree::Tree& pvst_, const bidirected::VariationGraph& bd_vg, const core::config& app_config);

/**
  * @brief Extracts the canonical flubbles from the PVST.
  *
  * @param pvst_ The PVST.
  * @return A vector of pairs of vertex ids that represent the canonical flubbles.
 */
std::vector<std::pair<std::size_t, std::size_t>> extract_canonical_flubbles(const tree::Tree& pvst_);
} // namespace genomics

namespace vcf {

const std::string UNDEFINED_PATH_LABEL = "undefined";
const std::size_t UNDEFINED_PATH_ID = core::constants::INVALID_ID;
const std::size_t UNDEFINED_PATH_POS = core::constants::INVALID_ID;

struct vcf_record {
  std::string chrom;
  std::size_t pos;
  //std::string id;
  //std::string ref;
  //std::string alt;
  //std::string qual;
  //std::string filter;
  //std::string info;
  //std::string format;
  std::string ref;
  std::vector<std::string> alt;
};

std::map<std::size_t, std::vector<vcf::vcf_record>>
gen_vcf_records(
  const bidirected::VariationGraph& bd_vg,
  const std::vector<std::vector<std::set<std::size_t>>>& haplotypes,
  const std::vector<std::vector<std::vector<bidirected::side_n_id_t>>>& all_paths,
  const core::config& app_config);

void write_vcfs(const std::map<std::size_t, std::vector<vcf_record>>& vcf_records,
                const bidirected::VariationGraph& bd_vg, const core::config& app_config);
bool print_vcf(std::vector<std::vector<std::size_t>> paths,
               digraph::DiGraph dg,
               std::string ref_path_name,
               std::size_t ref_path_id);

void print_vcf_header(std::string ref_name);
} // namespace vcf
