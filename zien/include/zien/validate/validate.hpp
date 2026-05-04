#ifndef ZIEN_VALIDATE_HPP
#define ZIEN_VALIDATE_HPP

#include <mto/from_vcf.hpp>	    // for read_vcf
#include <quilt/app.hpp>	    // for config
#include <oza/graph/bidirected.hpp> // for bd::VG
#include <quilt/types.hpp>	    // for qt

namespace zien::validate
{
constexpr std::string_view MODULE = "povu::validate::vcf";

std::vector<qt::u32> validate_vcf_records(const bd::VG &g,
					  const mto::from_vcf::VCFile &vcf_file,
					  const core::config &app_config,
					  bool output_to_file = true);
} // namespace zien::validate
#endif // ZIEN_VALIDATE_HPP

namespace zv = zien::validate;
