#if defined(__linux__)
#define _GNU_SOURCE // Enable GNU extensions for Linux-specific behavior
#elif defined(_WIN32)
// Provide portable macros for missing functionality on Windows
#define strtok_r strtok_s
#define strndup _strndup
#define strdup _strdup
#elif defined(__APPLE__)
// No additional definitions required for macOS
#else
#error "Platform not supported"
#endif

#include "../include/liteseq/gfa.h"
#include "../include/liteseq/types.h"
#include "../src/internal/lq_io.h"
#include "../src/internal/lq_utils.h"

#include "./gfa_l.h"
#include "./gfa_s.h"
#include "./refs/ref_impl.h"

#include <log.h>

#include <limits.h>
#include <pthread.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
/*
 * Utility Functions
 * -----------------
 */

#define H_LINE_VERSION_IDX 1 // the index of the version token in the H line

DEFINE_ENUM_AND_STRING(gfa_version, GFA_VERSION_ITEMS)

/*
 * GFA Parse Functions
 * -------------------
 */

status_t set_version(const char *h_line, const char *newline, gfa_props *g)
{
	char *tokens[EXPECTED_H_LINE_TOKENS] = {NULL};
	struct split_str_params p = {
		.str = h_line,
		.up_to = newline,
		.delimiter = TAB_CHAR,
		.fallbacks = "",
		.fallback_chars_count = 0,
		.max_splits = EXPECTED_H_LINE_TOKENS,

		.tokens_found = 0,
		.tokens = tokens,
		.end = NULL,
	};

	status_t res = split_str(&p);
	if (res != SUCCESS) {
		log_fatal("Could not parse H line");
		return res;
	}

	char *version_str = tokens[H_LINE_VERSION_IDX];
	char *nl = strchr(version_str, NEWLINE);
	g->version = from_string_gfa_version(version_str);
	if (g->version == gfa_version_INVALID) {
		log_fatal("Unsupported GFA version: %s", version_str);
		return ERROR_CODE_INVALID_ARGUMENT;
	}

	tokens_free(tokens, EXPECTED_H_LINE_TOKENS);

	return SUCCESS;
}

/**
 * @brief used to extract the v id from an S line
 *
 * @param [in] str the start of the S line
 * @return the vertex id
 */
u32 get_num_vid(const char *str, u32 linum)
{
	char *start = NULL;
	char *end = NULL;
	int count = 0;
	do {
		start = end;
		end = strchr(str, TAB_CHAR);
		if (end == NULL)
			log_fatal("Badly formatted S Line on line %u", linum);
		count++;
	} while (count < 2);

	u32 num = strtoull(start, &end, 10);

	return num;
}

void set_v_id_bounds(const char *curr_char, gfa_props *g, idx_t linum)
{
	u32 curr_v_id = get_num_vid(curr_char, linum);
	if (curr_v_id > g->max_v_id)
		g->max_v_id = curr_v_id;
	if (curr_v_id < g->min_v_id)
		g->min_v_id = curr_v_id;
}

/**
 * Process memory mapped file line by line and count the number of S, L, and P
 * lines
 *
 * @param [in] gfa_props gfa metadata
 * @return 0 on success, -1 on failure
 */
status_t analyse_gfa_structure(gfa_props *g)
{
	idx_t curr_line = 0;
	g->s_line_count = 0;
	g->l_line_count = 0;
	g->p_line_count = 0;
	g->w_line_count = 0;

	char *curr_char = g->start;
	while (curr_char < g->end) {
		// Find the next newline
		char *newline = memchr(curr_char, NEWLINE, g->end - curr_char);
		// If no newline is found, process the remainder of the file
		if (!newline) {
			newline = g->end;
		}

		switch (curr_char[0]) {
		case GFA_S_LINE:
			g->s_line_count++;
			set_v_id_bounds(curr_char, g, curr_line);
			break;
		case GFA_L_LINE:
			g->l_line_count++;
			break;
		case GFA_P_LINE:
			g->p_line_count++;
			break;
		case GFA_W_LINE:
			g->w_line_count++;
			break;
		case GFA_H_LINE:
			if (set_version(curr_char, newline, g) != 0) {
				fprintf(stderr, "[liteseq::gfa] Failed to set "
						"GFA version\n");
				return -1;
			}
			break;
		default: // unsupported line type
			fprintf(stderr,
				"[liteseq::gfa] Unsupported line type: [%c] on "
				"line: [%d]\n",
				curr_char[0], curr_line);
			return -2;
		}

		// Move to the next line
		curr_char = newline + 1; // Skip the newline character
		curr_line++;
	}

	// assumes at least one vertex found
	// TODO: [c] is that a safe assumption?
	g->vtx_arr_size = g->max_v_id + 1;

	return 0;
}

status_t index_lines(gfa_props *gfa)
{
	idx_t s_idx = 0;
	idx_t l_idx = 0;
	idx_t p_idx = 0;
	idx_t w_idx = 0;
	idx_t line_count = 0;	      // reset line count
	char *curr_char = gfa->start; // reset the current character

	char *newline;
	line curr_line;
	idx_t line_length;

	/*  allocate ... */
	gfa->s_lines = (line *)malloc(gfa->s_line_count * sizeof(line));
	gfa->l_lines = (line *)malloc(gfa->l_line_count * sizeof(line));
	gfa->p_lines = (line *)malloc(gfa->p_line_count * sizeof(line));
	gfa->w_lines = (line *)malloc(gfa->w_line_count * sizeof(line));

	if (!gfa->s_lines || !gfa->l_lines || !gfa->p_lines) {
		perror("Failed to allocate memory for line indices");
		return -1;
	}

	while (curr_char < gfa->end) {
		// Find the next newline
		newline = memchr(curr_char, NEWLINE, gfa->end - curr_char);
		// If no newline is found, process the remainder of the file
		if (!newline) {
			newline = gfa->end;
		}

		line_length = (idx_t)(newline - curr_char);
		curr_line = (line){.start = curr_char,
				   .line_idx = line_count++,
				   .len = line_length};

		switch (curr_char[0]) {
		case GFA_S_LINE:
			gfa->s_lines[s_idx++] = curr_line;
			break;
		case GFA_L_LINE:
			gfa->l_lines[l_idx++] = curr_line;
			break;
		case GFA_W_LINE:
			gfa->w_lines[w_idx++] = curr_line;
			break;
		case GFA_P_LINE:
			gfa->p_lines[p_idx++] = curr_line;
			break;
		case GFA_H_LINE:
			line_count++;
			break;
		default: // unsupported line type
			fprintf(stderr,
				"Unsupported line type: [%c] on line: [%d]\n",
				curr_char[0], line_count);
			return -1;
		}

		// Move to the next line
		curr_char = newline + 1; // Skip the newline character
	}

	return 0;
}

status_t set_ref_loci(gfa_props *gfa)
{
	vtx **vs = gfa->v;
	if (vs == NULL)
		return ERROR_CODE_INVALID_ARGUMENT;

	for (int i = 0; i < gfa->ref_count; i++) {
		struct ref *r = gfa->refs[i];
		struct ref_walk *rw = r->walk;
		idx_t pos = 1; // DNA is 1 indexed
		for (int j = 0; j < rw->step_count; j++) {
			rw->loci[j] = pos;
			idx_t v_id = rw->v_ids[j];

			if (vs[v_id] == NULL) {
				continue;
			}
			idx_t l = pos += strlen(vs[v_id]->seq);
		}
		set_hap_len(r, pos - 1);
	}

	return SUCCESS;
}

status_t populate_gfa(gfa_props *gfa)
{
	pthread_t thread_s, thread_l, thread_p;

	struct s_thread_meta s_meta = {
		.vertices = gfa->v,
		.s_lines = gfa->s_lines,
		.s_line_count = gfa->s_line_count,
		.inc_vtx_labels = gfa->inc_vtx_labels,
	};

	struct l_thread_meta l_meta = {
		.edges = gfa->e,
		.l_lines = gfa->l_lines,
		.l_line_count = gfa->l_line_count,
	};

	struct ref_thread_data ref_meta = {
		.refs = gfa->refs,
		.p_lines = gfa->p_lines,
		.w_lines = gfa->w_lines,
		.p_line_count = gfa->p_line_count,
		.w_line_count = gfa->w_line_count,
	};

	if (pthread_create(&thread_s, NULL, t_handle_s, (void *)&s_meta) != 0)
		return FAILURE; // Failed to create thread for S lines

	// res = pop_l(&thread_l, &l_meta);
	if (pthread_create(&thread_l, NULL, t_handle_l, (void *)&l_meta) != 0) {
		return FAILURE; // Failed to create thread for L lines
	}

	if (gfa->inc_refs) {
		if (pthread_create(&thread_p, NULL, t_handle_p,
				   (void *)&ref_meta) != 0)
			return FAILURE; // Failed to create thread for P lines
	}

	// TODO: what if one of the threads fails and returns early
	/* Wait for threads to finish */
	pthread_join(thread_s, NULL);
	pthread_join(thread_l, NULL);
	if (gfa->inc_refs)
		pthread_join(thread_p, NULL);

	if (gfa->inc_refs && gfa->inc_vtx_labels) {
		status_t res = set_ref_loci(gfa);
		if (res != SUCCESS) {
			log_fatal("Failed to set reference loci");
			return res;
		}
	}

	return SUCCESS;
}

/**
 * pre-allocate memory for vertices, edges, and maybe references
 */
status_t preallocate_gfa(gfa_props *p)
{
	p->v = calloc(p->vtx_arr_size, sizeof(vtx *));
	// p->v = malloc(sizeof(vtx) * p->vtx_arr_size);
	if (!p->v)
		return ERROR_CODE_OUT_OF_MEMORY;
	// init with NULLs makes freeing more straightforward
	// memset assumes seq pointers are at offset 0 in the vtx structure
	// memset(p->v, 0, sizeof(vtx) * p->vtx_arr_size);
	// if (p->inc_vtx_labels) {}

	p->e = malloc(p->l_line_count * sizeof(edge));
	if (!p->e)
		return ERROR_CODE_OUT_OF_MEMORY;

	if (p->inc_refs) { // Initialize references (paths)
		p->ref_count = p->p_line_count + p->w_line_count;
		p->refs = malloc(sizeof(struct ref *) * p->ref_count);
		if (!p->refs)
			return ERROR_CODE_OUT_OF_MEMORY;
	}

	return SUCCESS;
}

gfa_props *init_gfa(const gfa_config *conf)
{
	gfa_props *p = (gfa_props *)malloc(sizeof(gfa_props));

	p->fp = conf->fp;
	p->inc_vtx_labels = conf->inc_vtx_labels;
	p->inc_refs = conf->inc_refs;

	p->start = NULL;
	p->end = NULL;
	p->l_lines = NULL;
	p->s_lines = NULL;
	p->p_lines = NULL;
	p->w_lines = NULL;

	p->s_line_count = 0;
	p->l_line_count = 0;
	p->p_line_count = 0;
	p->w_line_count = 0;

	p->ref_count = 0;

	p->min_v_id = UINT32_MAX;
	p->max_v_id = 0;

	p->v = NULL;
	p->e = NULL;
	p->refs = NULL;

	p->file_size = 0;
	p->status = -1;

	return p;
}

gfa_props *gfa_new(const gfa_config *conf)
{
	gfa_props *p = init_gfa(conf); // set up the config
	if (p == NULL) {
		log_fatal("init gfa failed");
		return NULL;
	}

	char *mapped;	      // pointer to the start of the memory mapped file
	char *end;	      // pointer to the end of the memory mapped file
	size_t file_size = 0; // size of the memory mapped file

	p->status = -1; // status of a given operation

	open_mmap(p->fp, &mapped, &file_size);
	if (mapped == NULL) { // Failed to mmap file
		return p;
	}
	end = mapped + file_size;

	p->start = mapped;
	p->end = end;
	p->file_size = file_size;

	p->status = analyse_gfa_structure(p);
	if (p->status != 0) {
		fprintf(stderr, "Error: GFA file structure analysis failed\n");
		return p;
	}

	if (p->s_line_count == 0 && p->l_line_count == 0 &&
	    p->p_line_count == 0) {
		fprintf(stderr, "Error: GFA has no vertices edges or paths\n");
		return p;
	}

	index_lines(p);
	preallocate_gfa(p);
	p->status = populate_gfa(p);
	p->status = 0;

	return p;
}

void gfa_free(gfa_props *gfa)
{
	close_mmap(gfa->start, gfa->file_size);

	if (gfa->s_lines)
		free(gfa->s_lines);

	if (gfa->l_lines)
		free(gfa->l_lines);

	if (gfa->p_lines)
		free(gfa->p_lines);

	if (gfa->w_lines)
		free(gfa->w_lines);

	if (gfa->e)
		free(gfa->e);

	if (gfa->inc_vtx_labels)
		for (idx_t i = 0; i < (gfa->max_v_id + 1); i++)
			if (gfa->v[i] != NULL && gfa->v[i]->seq != NULL)
				free(gfa->v[i]->seq);

	for (idx_t i = 0; i < gfa->vtx_arr_size; i++)
		if (gfa->v[i] != NULL)
			free(gfa->v[i]);
	if (gfa->v)
		free(gfa->v);

	if (gfa->refs) {
		for (idx_t i = 0; i < gfa->ref_count; i++)
			destroy_ref(&(gfa->refs[i]));
		if (gfa->refs)
			free(gfa->refs);
	}

	free(gfa);
}
