#include <fstream>
#include <cstddef>
#include <iostream>
#include <stack>
#include <queue>
#include <string>
#include <map>
#include <set>
#include <unordered_set>
#include <vector>
#include <utility>
#include <format>
#include <numeric>

#include "../pvst/pvst.hpp"
#include "../graph/tree.hpp"
#include "../graph/digraph.hpp"
#include "../core/constants.hpp"
#include "../core/core.hpp"
#include "./genomics.hpp"
#include "../core/utils.hpp"

namespace vcf {

void write_vcf(const std::string& ref_name, const std::vector<vcf_record>& vcf_records) {
  	std::string vcf_file_name = std::format("{}/{}.vcf", ".", ref_name);
	std::ofstream vcf_file(vcf_file_name);

	if (!vcf_file.is_open()) {
	  std::cerr << "ERROR: could not open file " << vcf_file_name << "\n";
	  std::exit(1);
	}
	
	// write the header
	vcf_file << "##fileformat=VCFv4.2\n";
	vcf_file << "##fileDate=" << utils::today() << std::endl;
	vcf_file << "##source=povu\n";
	vcf_file << "##reference=" << ref_name << "\n";
	vcf_file << "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">\n";
	vcf_file << "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">\n";
	vcf_file << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
	vcf_file << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n";

	
	// write the records
	for (const vcf_record& vcf_rec : vcf_records) {
	  vcf_file << vcf_rec.chrom << "\t" << vcf_rec.pos << "\t.\t" << vcf_rec.ref
			   << "\t" << utils::concat_with(vcf_rec.alt, ',') << "\t.\t.\t.\tGT\t0/1\n";
	}

	vcf_file.close();
}
  
void write_vcfs(const std::map<std::size_t,
				std::vector<vcf_record>>& vcf_records,
				const bidirected::VariationGraph& bd_vg) {

  // this map is redundant
  std::map<std::size_t, std::string> path_id_name_map; //  id to path name
  for (auto p : bd_vg.get_paths()) {
	path_id_name_map[p.id] = p.name;
  }

 
  path_id_name_map[UNDEFINED_PATH_ID] = UNDEFINED_PATH_LABEL;

  for (auto& [ref_id, vcf_recs]: vcf_records) {

	std::string ref_name = path_id_name_map[ref_id];
   std::cout << "writing vcf for " << ref_name << "\n";
	write_vcf(ref_name, vcf_recs);
  }
}

  
} // namespace vcf
