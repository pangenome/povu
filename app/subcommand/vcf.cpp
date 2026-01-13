#include "./vcf.hpp"

#include <atomic> // for atomic
#include <thread> // for thread
#include <vector> // for vector

#include "povu/common/app.hpp"	      // for config
#include "povu/graph/bidirected.hpp"  // for bidirected
#include "povu/graph/pvst.hpp"	      // for pvst
#include "povu/io/from_gfa.hpp"	      // for to_bd
#include "povu/io/from_vcf.hpp"	      // for read_vcf
#include "povu/io/to_gfa.hpp"	      // for write_gfa
				      //
#include "zien/tui/tui.hpp"	      // for view
#include "zien/validate/validate.hpp" // for validate_vcf_records

namespace povu::subcommands::vcf
{

std::pair<bd::VG *, povu::io::from_vcf::VCFile>
data_loader(const core::config &app_config)
{

	povu::io::from_vcf::VCFile vcf_file;
	bd::VG *g = povu::io::from_gfa::to_bd(app_config);

	pt::u32 ll = app_config.verbosity(); // ll for log level
	const core::vcf_subcommand &vcf_opts = app_config.get_vcf_subcommand();

	// Load VCF
	povu::io::from_vcf::read_vcf(vcf_opts.get_input_vcf(), ll, vcf_file);

	return {g, vcf_file};
}

void handle_tui(const core::config &app_config)
{

	// 1. Initialize ncurses ONCE for this entire TUI session
	zt::NcursesGuard guard;

	std::atomic<bool> is_loading(true);

	// We need these to persist after the thread finishes
	bd::VG *g = nullptr;
	povu::io::from_vcf::VCFile vcf_file;
	std::vector<pt::u32> invalid_recs;

	// 1. Launch the loading thread
	std::thread loader(
		[&]()
		{
			auto res = data_loader(app_config);

			g = res.first;
			vcf_file = res.second;

			invalid_recs = zv::validate_vcf_records(
				*g, vcf_file, app_config, false);

			// 2. Signal that we are done
			is_loading.store(false);
		});

	// 3. Run the spinner on the MAIN thread
	// Note: Ncurses must be initialized inside show_loading_spinner or
	// before it
	zt::show_loading_spinner(is_loading);

	// 4. Wait for the loader thread to safely finish
	if (loader.joinable())
		loader.join();

	is_loading.store(false); // stop loading spinner

	zt::view(*g, vcf_file, invalid_recs);

	delete g;
}

void handle_report(const core::config &app_config)
{
	auto [g, vcf_file] = data_loader(app_config);
	zv::validate_vcf_records(*g, vcf_file, app_config);
	delete g;
}

void do_vcf(const core::config &app_config)
{

	const core::vcf_subcommand &vcf_opts = app_config.get_vcf_subcommand();

	switch (vcf_opts.get_vcf_options()) {
	case core::vcf_options::report:
		handle_report(app_config);
		break;
	case core::vcf_options::tui:
		handle_tui(app_config);
		break;
	}
}

} // namespace povu::subcommands::vcf
