#ifndef PV_IO_VCF_HPP
#define PV_IO_VCF_HPP

#include "../../include/common/genomics.hpp"
#include "../../include/graph/bidirected.hpp"

namespace povu::io::to_vcf {
namespace pvt = povu::types::genomics;
namespace bd = povu::bidirected;

void write_vcfs(pvt::VcfRecIdx &vcf_recs,const bd::VG &g, const core::config &app_config);

/*
using povu::genomics::vcf::vcf_record;



void write_vcfs(const std::map<std::size_t,
              std::vector<vcf_record>>& vcf_records,
              const bidirected::VariationGraph& bd_vg,
              const core::config& app_config);
*/
} // namespace povu::io::vcf


#endif // PV_IO_VCF_HPP
