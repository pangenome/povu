#ifndef POVU_VALIDATE_VCF_HPP
#define POVU_VALIDATE_VCF_HPP

#include "povu/common/app.hpp"	     // for config
#include "povu/graph/bidirected.hpp" // for bidirected
#include "povu/io/from_vcf.hpp"	     // for read_vcf

namespace povu::validate::vcf
{
constexpr std::string_view MODULE = "povu::validate::vcf";

void validate_vcf_records(const bd::VG &g,
			  const povu::io::from_vcf::VCFile &vcf_file,
			  const core::config &app_config);
} // namespace povu::validate::vcf
#endif // POVU_VALIDATE_VCF_HPP
