#include <iostream>

#include "./subcommand/call.hpp"
#include "./subcommand/decompose.hpp"
#include "./subcommand/gfa2vcf.hpp"
#include "./subcommand/info.hpp"

namespace pv = povu::subcommands;

constexpr std::string_view MODULE = "povu::main";

int main(int argc, char *argv[]) {

  core::config app_config;
  cli::cli(argc, argv, app_config);

  if (app_config.verbosity()) { app_config.dbg_print(); }

  switch (app_config.get_task()) {
  case core::task_e::deconstruct:
    pv::decompose::do_decompose(app_config);
    break;
  case core::task_e::call:
    pv::call::do_call(app_config);
    break;
  case core::task_e::gfa2vcf:
    pv::gfa2vcf::do_gfa2vcf(app_config);
    break;
  case core::task_e::info:
    pv::info::do_info(app_config);
    break;
  default:
    // the help text handles this case
    break;
  }

  return 0;
}
