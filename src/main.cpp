#include <iostream>

#include "../include/common/types/compat.hpp"

#include "./subcommand/call.hpp"
#include "./subcommand/deconstruct.hpp"
#include "./subcommand/gfa2vcf.hpp"
#include "./subcommand/info.hpp"



namespace pv = povu::subcommands;

int main(int argc, char *argv[]) {
  std::string fn_name = pv_cmp::format("[povu::main::{}]", __func__);

  core::config app_config;
  cli::cli(argc, argv, app_config);

  if (app_config.verbosity()) { app_config.dbg_print(); }

  switch (app_config.get_task()) {
  case core::task_e::deconstruct:
    pv::deconstruct::do_deconstruct(app_config);
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
    std::cerr << pv_cmp::format("{} Task not recognized\n", fn_name);
    break;
  }

  return 0;
}
