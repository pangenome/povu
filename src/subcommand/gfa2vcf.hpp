#ifndef POVU_GFA2VCF_HPP
#define POVU_GFA2VCF_HPP

#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <string>
#include <vector>

#include "./call.hpp"
#include "./deconstruct.hpp"
#include "../../include/common/types/compat.hpp"
#include "../cli/app.hpp"



namespace povu::subcommands::gfa2vcf {

namespace fs = std::filesystem;

void do_gfa2vcf(const core::config &app_config);

} // namespace povu::subcommands::gfa2vcf

#endif // POVU_GFA2VCF_HPP
