// #include <string_view> // for string_view

#include "./cli/cli.hpp"	      // for cli
#include "./subcommand/call.hpp"      // for do_call
#include "./subcommand/decompose.hpp" // for do_decompose
#include "./subcommand/gfa2vcf.hpp"   // for do_gfa2vcf
#include "./subcommand/info.hpp"      // for do_info
#include "povu/common/app.hpp"	      // for task_e, config

// constexpr std::string_view MODULE = "povu::main";
namespace pv = povu::subcommands;

int main(int argc, char *argv[])
{
	core::config app_config;
	cli::cli(argc, argv, app_config);

	if (app_config.verbosity())
		app_config.dbg_print();

	switch (app_config.get_task()) {
	case core::task_e::decompose:
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
	default: // the help text handles this case
		break;
	}

	return 0;
}
