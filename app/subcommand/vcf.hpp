#ifndef PV_SUBCOMMANDS_VCF_HPP
#define PV_SUBCOMMANDS_VCF_HPP

#include <string_view> // for string_view

#include "povu/common/app.hpp"

namespace povu::subcommands::vcf
{
constexpr std::string_view MODULE = "povu::subcommands::vcf";

void do_vcf(const core::config &app_config);
} // namespace povu::subcommands::vcf

#endif
