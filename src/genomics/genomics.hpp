#include <vector>

#include "../pvst/pvst.hpp"
#include "../graph/tree.hpp"
#include "../graph/digraph.hpp"
#include "../core/core.hpp"
#include "../graph/bidirected.hpp"


namespace genomics {

typedef std::size_t id_t; //
typedef id_t vertex_id_t; //
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
  
enum output_format {
	VCF, // which version of VCF?
	PAF, // not yet supported
};


//void call_variants(tree::Tree pvst_, digraph::DiGraph dg, core::config app_config);

void call_variants(const tree::Tree& pvst_, const bidirected::VariationGraph& bd_vg, const core::config& app_config);
std::vector<std::pair<std::size_t, std::size_t>> extract_canonical_flubbles(const tree::Tree& pvst_);
  
} // namespace genomics

namespace vcf {

const std::string UNDEFINED_PATH_LABEL = "undefined";
const std::size_t UNDEFINED_PATH_ID = core::constants::SIZE_T_MAX;
const std::size_t UNDEFINED_PATH_POS = core::constants::SIZE_T_MAX;

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

void write_vcfs(const std::map<std::size_t, std::vector<vcf_record>>& vcf_records,
				const bidirected::VariationGraph& bd_vg);
bool print_vcf(std::vector<std::vector<std::size_t>> paths,
			   digraph::DiGraph dg,
			   std::string ref_path_name,
			   std::size_t ref_path_id);

void print_vcf_header(std::string ref_name);

} // namespace vcf
