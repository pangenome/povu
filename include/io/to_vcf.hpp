#ifndef PV_IO_VCF_HPP
#define PV_IO_VCF_HPP

#include <cstddef> // std::size_t
#include <fstream> // std::ofstream
#include <vector>

#include "../../include/common/types/compat.hpp"
#include "../../include/common/types/genomics.hpp"
#include "../../include/common/types/graph.hpp"
#include "../../include/common/utils.hpp"
#include "../../include/graph/bidirected.hpp"
#include "../../app/cli/app.hpp"

namespace povu::io::to_vcf {
inline constexpr std::string_view MODULE = "povu::io::to_vcf";
namespace pvt = povu::types::genomics;
namespace bd = povu::bidirected;
namespace pu = povu::utils;
namespace pgt = povu::types::graph;
namespace pt = povu::types;

void write_vcfs(const pvt::VcfRecIdx &vcf_recs, const bd::VG &g, const core::config &app_config);
void write_combined_vcf_to_stdout(const pvt::VcfRecIdx &vcf_recs, const bd::VG &g, const core::config &app_config);
} // namespace povu::io::vcf


#endif // PV_IO_VCF_HPP
