#ifndef POVU_GFA2VCF_HPP
#define POVU_GFA2VCF_HPP

#include <string>
#include <vector>

#include "../cli/app.hpp"

namespace povu::subcommands::gfa2vcf {

void do_gfa2vcf(const core::config &app_config);

} // namespace povu::subcommands::gfa2vcf

#endif // POVU_GFA2VCF_HPP