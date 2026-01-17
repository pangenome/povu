#ifndef ZIEN_VALIDATE_HPP
#define ZIEN_VALIDATE_HPP

#include "mto/from_vcf.hpp" // for read_vcf

#include "povu/common/app.hpp"	     // for config
#include "povu/graph/bidirected.hpp" // for bd::VG

namespace zien::validate
{
constexpr std::string_view MODULE = "povu::validate::vcf";

std::vector<pt::u32> validate_vcf_records(const bd::VG &g,
					  const mto::from_vcf::VCFile &vcf_file,
					  const core::config &app_config,
					  bool output_to_file = true);
} // namespace zien::validate
#endif // ZIEN_VALIDATE_HPP

namespace zv = zien::validate;
