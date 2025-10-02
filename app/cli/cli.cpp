#include "./cli.hpp"

namespace cli
{

struct decomopose_opts {
	args::Group decompose;
	args::Flag hairpins;
	args::Flag subflubbles;

	explicit decomopose_opts(args::Subparser &p)
	    : decompose(p, "Decompose options",
			args::Group::Validators::DontCare),
	      hairpins(decompose, "hairpins",
		       "Find hairpins in the variation graph [default: false]",
		       {'h', "hairpins"}),
	      subflubbles(decompose, "subfubbles",
			  "Find subflubbles in the variation graph [default: "
			  "false]",
			  {'s', "subflubbles"})
	{}
};

struct streaming_opts {
	args::Group streaming;
	args::ValueFlag<std::size_t> chunk_size;
	args::ValueFlag<std::size_t> queue_length;

	explicit streaming_opts(args::Subparser &p)
	    : streaming(p, "Streaming options",
			args::Group::Validators::DontCare),
	      chunk_size(streaming, "chunk_size",
			 "Number of variants to process in each chunk "
			 "[default: 100]",
			 {'c', "chunk-size"}),
	      queue_length(streaming, "queue_length",
			   "Number of chunks to buffer [default: 4]",
			   {'q', "queue-length"})
	{}
};

struct output_opts {
	args::Group outsel;
	args::ValueFlag<std::string> output_dir;
	args::Flag stdout_vcf;

	explicit output_opts(args::Subparser &p)
	    : outsel(p, "Output destination (choose exactly one)",
		     args::Group::Validators::Xor),
	      output_dir(outsel, "output_dir", "Output directory [default: .]",
			 {'o', "output-dir"}),
	      stdout_vcf(outsel, "stdout_vcf",
			 "Output single VCF to stdout instead of separate "
			 "files [default: false]",
			 {"stdout"})
	{}
};

// Holds the three "reference source" options and enforces XOR.
struct reference_opts {
	args::Group refsel;
	args::ValueFlag<std::string> prefix_list;
	args::ValueFlagList<std::string> path_prefixes;
	args::PositionalList<std::string> refs_positional;

	/* One of ref_list, path_prefixes, or list of references must be
	 * set—never multiple, and never none */
	explicit reference_opts(args::Subparser &p)
	    : refsel(p, "Reference source (choose exactly one)",
		     args::Group::Validators::Xor),
	      prefix_list(refsel, "prefix_list",
			  "path to file containing reference name prefixes "
			  "[optional]",
			  {'r', "prefix-list"}),
	      path_prefixes(refsel, "path_prefix",
			    "All paths beginning with NAME used as reference "
			    "(multiple allowed) [optional]",
			    {'P', "path-prefix"}),
	      refs_positional(
		      refsel, "refs",
		      "list of refs to use as reference haplotypes [optional]")
	{}
};

void populate_ref_ops(reference_opts &ref_opts, core::config &app_config)
{
	if (ref_opts.prefix_list) {
		app_config.set_ref_input_format(
			core::input_format_e::file_path);
		app_config.set_reference_txt_path(
			std::move(args::get(ref_opts.prefix_list)));
		return;
	}

	app_config.set_ref_input_format(core::input_format_e::params);
	auto prefixes = ref_opts.path_prefixes
				? args::get(ref_opts.path_prefixes)
				: args::get(ref_opts.refs_positional);

	for (auto &&p : prefixes) {
		app_config.add_ref_name_prefix(p);
	}
}

void call_handler(args::Subparser &parser, core::config &app_config)
{
	args::Group arguments("arguments");
	args::ValueFlag<std::string> input_gfa(
		parser, "gfa", "path to input gfa [required]",
		{'i', "input-gfa"}, args::Options::Required);
	args::ValueFlag<std::string> forest_dir(
		parser, "forest_dir",
		"dir containing flubble forest [default: .]",
		{'f', "forest-dir"});
	streaming_opts stream_opts(parser);
	output_opts out_opts(parser);
	reference_opts ref_opts(parser);

	parser.Parse();

	app_config.set_task(core::task_e::call);

	// set mandatory call opts for graph
	app_config.set_inc_vtx_labels(true);
	app_config.set_inc_refs(true);

	// input gfa is already a c_str
	app_config.set_input_gfa(args::get(input_gfa));

	if (forest_dir) {
		app_config.set_forest_dir(args::get(forest_dir));
	}

	{ // output options
		if (out_opts.output_dir) {
			app_config.set_output_dir(
				args::get(out_opts.output_dir));
		}
		else if (out_opts.stdout_vcf) {
			app_config.set_stdout_vcf(true);
		}
	}

	{ // streaming options
		if (stream_opts.chunk_size) {
			app_config.set_chunk_size(
				args::get(stream_opts.chunk_size));
		}

		if (stream_opts.queue_length) {
			app_config.set_queue_len(
				args::get(stream_opts.queue_length));
		}
	}

	// ref handling
	populate_ref_ops(ref_opts, app_config);
}

void gfa2vcf_handler(args::Subparser &parser, core::config &app_config)
{
	args::Group arguments("arguments");
	args::ValueFlag<std::string> input_gfa(
		parser, "gfa", "path to input gfa [required]",
		{'i', "input-gfa"}, args::Options::Required);
	decomopose_opts decomp_opts(parser);
	streaming_opts stream_opts(parser);
	output_opts out_opts(parser);
	reference_opts ref_opts(parser);

	parser.Parse();

	// set mandatory call opts for graph
	app_config.set_inc_vtx_labels(true);
	app_config.set_inc_refs(true);

	app_config.set_task(core::task_e::gfa2vcf);
	app_config.set_input_gfa(args::get(input_gfa));

	{ // set decompose options
		if (decomp_opts.hairpins) {
			app_config.set_hairpins(true);
		}

		if (decomp_opts.subflubbles) {
			app_config.set_subflubbles(true);
		}
	}

	{ // output options
		if (out_opts.output_dir) {
			app_config.set_output_dir(
				args::get(out_opts.output_dir));
		}
		else if (out_opts.stdout_vcf) {
			app_config.set_stdout_vcf(true);
		}
	}

	{ // streaming options
		if (stream_opts.chunk_size) {
			app_config.set_chunk_size(
				args::get(stream_opts.chunk_size));
		}

		if (stream_opts.queue_length) {
			app_config.set_queue_len(
				args::get(stream_opts.queue_length));
		}
	}

	// ref handling
	populate_ref_ops(ref_opts, app_config);
}

void info_handler(args::Subparser &parser, core::config &app_config)
{
	args::Group arguments("arguments");
	args::ValueFlag<std::string> input_gfa(
		parser, "gfa", "path to input gfa [required]",
		{'i', "input-gfa"}, args::Options::Required);
	args::Flag tips(parser, "tips", "print the tips", {'t', "print_tips"});

	parser.Parse();
	app_config.set_task(core::task_e::info);

	if (tips) {
		app_config.set_print_tips(true);
	}

	// input gfa is already a c_str
	app_config.set_input_gfa(args::get(input_gfa));
}

void decompose_handler(args::Subparser &parser, core::config &app_config)
{
	args::Group arguments("arguments");
	args::ValueFlag<std::string> input_gfa(
		parser, "gfa", "path to input gfa [required]",
		{'i', "input-gfa"}, args::Options::Required);
	args::ValueFlag<std::string> output_dir(parser, "output_dir",
						"Output directory [default: .]",
						{'o', "output-dir"});
	decomopose_opts decomp_opts(parser);

	parser.Parse();
	app_config.set_task(core::task_e::decompose);

	{ // set decompose options
		if (decomp_opts.hairpins) {
			app_config.set_hairpins(true);
		}

		if (decomp_opts.subflubbles) {
			app_config.set_subflubbles(true);
		}
	}

	// input gfa is already a c_str
	app_config.set_input_gfa(args::get(input_gfa));

	if (output_dir) {
		app_config.set_output_dir(args::get(output_dir));
	}
}

int cli(int argc, char **argv, core::config &app_config)
{

	args::ArgumentParser p(
		"Explore genomic variation in a variation graph");
	args::Group commands(p, "commands");

	args::Command gfa2vcf(commands, "gfa2vcf",
			      "Convert GFA to VCF (decompose + call)",
			      [&](args::Subparser &parser)
			      { gfa2vcf_handler(parser, app_config); });
	args::Command decompose(commands, "decompose",
				"Find regions of variation",
				[&](args::Subparser &parser)
				{ decompose_handler(parser, app_config); });
	args::Command call(commands, "call",
			   "Generate a VCF from regions of variation",
			   [&](args::Subparser &parser)
			   { call_handler(parser, app_config); });
	args::Command info(commands, "info",
			   "Print graph information [uses 1 thread]",
			   [&](args::Subparser &parser)
			   { info_handler(parser, app_config); });

	args::Group arguments(p, "arguments", args::Group::Validators::DontCare,
			      args::Options::Global);
	args::Flag version(arguments, "version", "The current version of povu",
			   {"version"});
	args::ValueFlag<int> verbosity(arguments, "verbosity",
				       "Level of output [default: 0]",
				       {'v', "verbosity"});
	args::ValueFlag<int> thread_count(
		arguments, "threads", "Number of threads to use [default: 1]",
		{'t', "threads"});
	args::Flag progress(arguments, "progress", "Show progress bars",
			    {"progress"});
	args::HelpFlag h(arguments, "help", "help", {'h', "help"});

	try {
		p.ParseCLI(argc, argv);
	}
	catch (args::Help &_) {
		std::cout << p;
	}
	catch (args::Error &e) {
		// only run this if the user is not requesting to print the
		// version
		if (!version) {
			std::cerr << e.what() << std::endl << p;
			return 1;
		}
	}

	if (version) {
		std::cout << VERSION << std::endl;
		std::exit(EXIT_SUCCESS);
	}

	if (args::get(verbosity)) {
		app_config.set_verbosity(args::get(verbosity));
	}

	if (thread_count) {
		app_config.set_thread_count(args::get(thread_count));
	}

	if (progress) {
		app_config.set_progress(true);
	}

	return 0;
}
} // namespace cli
