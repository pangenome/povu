#include <format>
#include <fstream>
#include <cstddef>
#include <iostream>
#include <map>
#include <ostream>
#include <set>
#include <string>
#include <vector>

#include "../cli/app.hpp"
#include "../common/utils.hpp"
#include "../graph/bidirected.hpp"
#include "./io.hpp"

namespace povu::io::vcf {
namespace pu = povu::utils;
namespace pc = povu::constants;
namespace pgt = povu::graph_types;
namespace bd = povu::bidirected;

  // TODO remove
/**
 * @brief Given a bidirected::VariationGraph and a subpath, return the string sequence
 */
std::string path_to_seq(const bd::VariationGraph& bd_vg,
                        const std::vector<pgt::id_n_orientation_t>& p) {
  std::string seq = "";
  for (auto [id, orientation]: p) {
    bd::Vertex const& v = bd_vg.get_vertex(id);
    seq += orientation == pgt::orientation_t::forward ? v.get_label() : v.get_rc_label();
  }

  return seq;
}


void write_vcf(const std::string& ref_name,
               const std::vector<vcf_record>& vcf_records,
               const core::config& app_config) {
  std::string fn_name = std::format("[povu::io::vcf::{}]", __func__);
  std::string vcf_file_name = std::format("{}/{}.vcf", std::string(app_config.get_output_dir()), ref_name);
  std::ofstream vcf_file(vcf_file_name);

    if (!vcf_file.is_open()) {
      std::cerr << std::format("{} ERROR: could not open file {}\n", fn_name, vcf_file_name);
      std::exit(1);
    }

    // write the header
    vcf_file << "##fileformat=VCFv4.2\n";
    vcf_file << "##fileDate=" << pu::today() << std::endl;
    vcf_file << "##source=povu\n";
    vcf_file << "##reference=" << ref_name << "\n";
    vcf_file << "##INFO=<ID=AT,Number=R,Type=String,Description=\"Allele Traversal as path in graph\">\n";
    vcf_file << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\n";

    // write the records
    for (const vcf_record& vcf_rec : vcf_records) {
      vcf_file << app_config.get_chrom()
               << "\t" << (vcf_rec.pos == pc::UNDEFINED_PATH_POS ? std::to_string(-1) : std::to_string(vcf_rec.pos))
               << "\t" << vcf_rec.id
               << "\t" << (vcf_rec.ref.empty() ? "." : vcf_rec.ref)
               << "\t" << (pu::concat_with(vcf_rec.alt, ',') == "" ? "." : pu::concat_with(vcf_rec.alt, ','))
               << "\t" << 60 // qual
               << "\t." // filter
               << "\t" << vcf_rec.format // info
               << "\tGT" // format
               << "\n";
    }

    vcf_file.close();
}

void write_vcfs(const std::map<std::size_t,
                std::vector<vcf_record>>& vcf_records,
                const bd::VariationGraph& bd_vg,
                const core::config& app_config) {
  std::string fn_name = std::format("[povu::vcf::{}]", __func__);

  // this map is redundant
  std::map<std::size_t, std::string> path_id_name_map; //  id to path name
  for (auto p : bd_vg.get_paths()) {
    path_id_name_map[p.id] = p.name;
  }

  path_id_name_map[pc::UNDEFINED_PATH_ID] = pc::UNDEFINED_PATH_LABEL;

  for (auto& [ref_id, vcf_recs]: vcf_records) {
    std::string ref_name = path_id_name_map[ref_id];

    if (app_config.verbosity() > 0) {
      std::cerr << std::format("{} writing vcf for {}\n", fn_name, ref_name);
    }

    write_vcf(ref_name, vcf_recs, app_config);
  }
}


} // namespace io::vcf
