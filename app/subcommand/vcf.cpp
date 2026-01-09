#include "./vcf.hpp"

#include "povu/common/app.hpp"	     // for config
#include "povu/graph/bidirected.hpp" // for bidirected
#include "povu/graph/pvst.hpp"	     // for pvst
#include "povu/io/from_gfa.hpp"	     // for to_bd
#include "povu/io/from_vcf.hpp"	     // for read_vcf
#include "povu/io/to_gfa.hpp"	     // for write_gfa
#include "povu/validate/vcf.hpp"     // for validate_vcf_records

namespace povu::subcommands::vcf
{
void do_vcf(const core::config &app_config)
{
	pt::u32 ll = app_config.verbosity();		   // ll for log level
	bd::VG *g = povu::io::from_gfa::to_bd(app_config); // read graph

	const core::vcf_subcommand &vcf_opts = app_config.get_vcf_subcommand();
	povu::io::from_vcf::VCFile vcf_file;
	switch (vcf_opts.get_vcf_options()) {
	case core::vcf_options::verify:
		povu::io::from_vcf::read_vcf(vcf_opts.get_input_vcf(), ll,
					     vcf_file);
		povu::validate::vcf::validate_vcf_records(*g, vcf_file,
							  app_config);
		break;
	case core::vcf_options::compare:
		break;
	}

	delete g;

	return;
}

} // namespace povu::subcommands::vcf
