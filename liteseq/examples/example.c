#include "liteseq/refs.h"
#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>

#include <liteseq/gfa.h>
#include <liteseq/types.h>
#include <liteseq/version.h>
#include <log.h>
#include <stdlib.h>

gfa_config gen_config(const char *fp)
{
	gfa_config conf_a = {
		.fp = fp,
		.inc_vtx_labels = false,
		.inc_refs = false,
	};

	gfa_config conf_b = {
		.fp = fp,
		.inc_vtx_labels = true,
		.inc_refs = true,
	};

	(void)conf_a; // to avoid unused variable warning
	(void)conf_b; // to avoid unused variable warning

	gfa_config conf_c = {
		.fp = fp,
		.inc_vtx_labels = true,
		.inc_refs = true,
	};

	return conf_c;
}

int main(int argc, char *argv[])
{
	// Check if no arguments or '-h' is passed
	if (argc != 2 || strcmp(argv[1], "-h") == 0) {
		fprintf(stderr, "LiteSeq version %s\n", LITESEQ_VERSION_STRING);
		fprintf(stderr, "Usage: %s <file_path>\n", argv[0]);
		exit(ERROR_CODE_INVALID_ARGUMENT);
	}

	const char *fp = argv[1];
	gfa_config conf = gen_config(fp);
	gfa_props *gfa = gfa_new(&conf);

	if (gfa->status != 0) {
		log_fatal("GFA parsing failed. Status: %d", gfa->status);
		gfa_free(gfa);
		exit(EXIT_FAILURE);
	}

	log_info("GFA parsed successfully");
	log_info("GFA version %s", to_string_gfa_version(gfa->version));
	log_info("vtx count %u", gfa->s_line_count);

	/* confirm the refs parsed */
	if (gfa->inc_refs) {
		printf("refs:\n");
		for (size_t i = 0; i < gfa->ref_count; i++) {
			printf("%s\n", get_tag(gfa->refs[i]));
		}
	} else {
		printf("No refs parsed\n");
	}

	log_info("Freeing GFA resources");

	gfa_free(gfa);

	return EXIT_SUCCESS;
}
