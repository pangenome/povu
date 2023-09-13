#include <vector>

#include "../pvst/pvst.hpp"
#include "../graph/tree.hpp"
#include "../graph/digraph.hpp"
namespace genomics {
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
  
  
void call_variants(tree::Tree pvst_, digraph::DiGraph dg);


} // namespace genomics

namespace vcf {
void print_vcf_variants(std::vector<std::vector<std::size_t>> paths, digraph::DiGraph dg);
void print_vcf_header();

} // namespace vcf
