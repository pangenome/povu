#include <fstream>
#include <cstddef>
#include <iostream>
#include <ostream>
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

// #include "../pvst/pvst.hpp"
// #include "../graph/tree.hpp"
//#include "../graph/digraph.hpp"
//#include "../core/constants.hpp"
#include "../cli/app.hpp"
#include "./io.hpp"
#include "../common/utils.hpp"
#include "../graph/bidirected.hpp"

namespace vcf {

// TODO: move to utils?
/**
 * @brief Given a bidirected::VariationGraph and a subpath, return the string sequence
 */
std::string path_to_seq(const bidirected::VariationGraph& bd_vg,
                        const std::vector<bidirected::side_n_id_t>& p) {
  std::string seq = "";
  for (auto [side, id]: p) {
    bidirected::Vertex const& v = bd_vg.get_vertex(id);
    seq += side == bidirected::VertexEnd::l ? v.get_label() : v.get_rc_label();
  }

  return seq;
}


void write_vcf(const std::string& ref_name,
               const std::vector<vcf_record>& vcf_records,
               const core::config& app_config) {
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

      vcf_file << app_config.get_chrom()
               << "\t" << (vcf_rec.pos == UNDEFINED_PATH_POS ? std::to_string(-1) : std::to_string(vcf_rec.pos))
               << "\t."
               << "\t" << vcf_rec.ref
               << "\t" << utils::concat_with(vcf_rec.alt, ',') << "\t.\t.\t.\tGT\t0/1\n";
    }

    vcf_file.close();
}

void write_vcfs(const std::map<std::size_t,
                std::vector<vcf_record>>& vcf_records,
                const bidirected::VariationGraph& bd_vg,
                const core::config& app_config) {
  std::string fn_name = std::format("[povu::vcf::{}]", __func__);

  // this map is redundant
  std::map<std::size_t, std::string> path_id_name_map; //  id to path name
  for (auto p : bd_vg.get_paths()) {
    path_id_name_map[p.id] = p.name;
  }

  path_id_name_map[UNDEFINED_PATH_ID] = UNDEFINED_PATH_LABEL;

  for (auto& [ref_id, vcf_recs]: vcf_records) {
    std::string ref_name = path_id_name_map[ref_id];

    if (app_config.verbosity() > 0) {
      std::cerr << fn_name << " writing vcf for " << ref_name << "\n";
    }

    write_vcf(ref_name, vcf_recs, app_config);
  }
}

// TODO: rewrite to reduce complexity
/**
 * @brief
 */
std::map<std::size_t, std::vector<vcf::vcf_record>> gen_vcf_records(
  const bidirected::VariationGraph& bd_vg,
  const std::vector<std::vector<std::set<std::size_t>>>& haplotypes,
  const std::vector<std::vector<std::vector<bidirected::side_n_id_t>>>& all_paths,
  const core::config& app_config
  ) {

  // all our VCF files will be in this map
  std::map<std::size_t, std::vector<vcf::vcf_record>> vcf_records;

  std::vector<std::string> reference_paths = app_config.get_reference_paths();

  utils::TwoWayMap<std::size_t, std::string> paths_map; // path id to its name

  for (auto p : bd_vg.get_paths()) {
    paths_map.insert(p.id, p.name);
  }

  if (app_config.gen_undefined_vcf()) {
    // add the undefined reference for flubbles that don't have the a reference
    paths_map.insert(vcf::UNDEFINED_PATH_ID, vcf::UNDEFINED_PATH_LABEL);
    reference_paths.push_back(vcf::UNDEFINED_PATH_LABEL);
  }

  //std::size_t num_flubbles = haplotypes.size();
  std::size_t num_paths{};
  // for each flubble
  for (std::size_t fl_idx{}; fl_idx < all_paths.size(); ++fl_idx) {
    num_paths = all_paths[fl_idx].size();

    // does this flubble contain a reference path?
    bool has_ref{false};

    // while being strand aware
    // spell all the sequences in the flubble as a vector of strings
    // where each string is a path
    std::vector<std::string> path_seqs;
    for (std::size_t path_idx{}; path_idx < num_paths; ++path_idx) {
      const std::vector<bidirected::side_n_id_t>& p = all_paths[fl_idx][path_idx];
      path_seqs.push_back(path_to_seq(bd_vg, p));
    }

    // for each path in the flubble
    // if the path is not a reference path
    // then it is a variant
    // and we need to add it to the vcf
    // we look at the haplotypes that pass through the path
    for (std::size_t path_idx{}; path_idx < num_paths; ++path_idx) {

      // get the haplotypes that pass through the path at path_idx
      std::set<std::size_t> path_haplotypes = haplotypes[fl_idx][path_idx];

      // not every path will be supported by at least one haplotype
      if (path_haplotypes.empty()) {
        //std::cerr << "ERROR: path " << path_idx << " in flubble " << fl_idx << " has no haplotypes\n";
        continue;
      }

      // because we want a VCF for each reference path then
      // for each reference path we check if it is in the set of haplotypes
      // for that path
      // if it is then remove it from the path seqs and the rest are variants
      // we then populate the vcf records map which contains
      // the vcf records for each reference path
      for (const std::string& ref_name : reference_paths) {
        std::size_t ref_id = paths_map.get_key(ref_name);

        if (path_haplotypes.count(ref_id)) {

          has_ref = true;

          // the path contains the reference haplotype ref_id
          vcf::vcf_record vcf_rec;

          // TODO: this is a hack, a path may differ internally by a little
          // use the intersection to determine this

          // get the position of that path in the reference whose id is ref_id
          auto [_, v_id] = all_paths[fl_idx][path_idx].front();
          bidirected::Vertex const& v = bd_vg.get_vertex(v_id);

          // TODO: can this loop be avoided?
          for (const bidirected::PathInfo& vertex_path : v.get_paths()) {

            if (vertex_path.path_id == ref_id) {
              vcf_rec.pos = vertex_path.step_index;
            }
          }
          vcf_rec.ref = path_seqs[path_idx];
          vcf_rec.alt = utils::immutable_erase(path_seqs, path_idx);

          vcf_records[ref_id].push_back(vcf_rec);
        }
      }
    }

    if (app_config.gen_undefined_vcf() && !has_ref) {
      vcf::vcf_record vcf_rec;
      vcf_rec.pos = vcf::UNDEFINED_PATH_POS;
      vcf_rec.ref = "";
      vcf_rec.alt = path_seqs;

      vcf_records[vcf::UNDEFINED_PATH_ID].push_back(vcf_rec);
    }
  }

  return vcf_records;
}


} // namespace vcf
