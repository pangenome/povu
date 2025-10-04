#ifndef POVU_GFA2VCF_HPP
#define POVU_GFA2VCF_HPP

#include <chrono>               // for filesystem

#include "povu/common/app.hpp"  // for config

namespace povu::subcommands::gfa2vcf
{

namespace fs = std::filesystem;

void do_gfa2vcf(const core::config &app_config);

} // namespace povu::subcommands::gfa2vcf

#endif // POVU_GFA2VCF_HPP
