#include "povu/call.hpp"       // for do_call
#include "povu/cli.hpp"	       // for cli
#include "povu/common/app.hpp" // for task_e, config
#include "povu/decompose.hpp"  // for do_decompose
#include "povu/gfa2vcf.hpp"    // for do_gfa2vcf
#include "povu/info.hpp"       // for do_info
#include "povu/prune.hpp"      // for do_prune
#include "povu/vcf.hpp"	       // for do_vcf
#include "povu/view.hpp"       // for do_view

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
	case core::task_e::prune:
		pv::prune::do_prune(app_config);
		break;
	case core::task_e::vcf:
		pv::vcf::do_vcf(app_config);
		break;
	case core::task_e::view:
		pv::view::do_view(app_config);
		break;
	default:
		break;
	}

	return 0;
}
