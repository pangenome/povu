#include "../cli/app.hpp"
#include "../graph/bidirected.hpp"

namespace io::from_gfa {
  bidirected::VariationGraph to_vg(const char* filename, const core::config& app_config);
}; // namespace io::from_gfa

namespace vcf {
namespace constants = common::constants;

const std::string UNDEFINED_PATH_LABEL = "undefined";
const std::size_t UNDEFINED_PATH_ID {constants::INVALID_ID};
const std::size_t UNDEFINED_PATH_POS {constants::INVALID_ID};

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
  const std::vector<std::vector<std::vector<bidirected::id_n_orientation_t>>>& all_paths,
  const core::config& app_config);

void write_vcfs(const std::map<std::size_t, std::vector<vcf_record>>& vcf_records,
                const bidirected::VariationGraph& bd_vg, const core::config& app_config);
} // namespace vcf
