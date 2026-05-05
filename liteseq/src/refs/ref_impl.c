#ifdef __linux__
#define _GNU_SOURCE // Enable GNU extensions on Linux
#elif defined(_WIN32)
#define strtok_r strtok_s // Map strtok_r to strtok_s on Windows
#define strndup _strndup  // Map strndup to _strndup
#define strdup _strdup	  // Map strdup to _strdup
#elif defined(__APPLE__)
/* On macOS, strndup and strdup exist, but no mappings are needed. */
#else
#error "Platform not supported"
#endif

#include <assert.h>
#include <stdio.h>
#include <string.h>

#include <log/log.h>

#include "../../include/liteseq/refs.h"
#include "../../include/liteseq/types.h"
#include "../../src/internal/lq_utils.h"

#include "./ref_impl.h"
#include "./ref_name.h"
#include "./ref_walk.h"

#define PANSN_ID_PARTS_COUNT 3

#define READ_P_LINE_TOKENS 3 // the number of tokens expected in a P line
#define READ_W_LINE_TOKENS 7 // the number of tokens expected in a P line

// the column of the W line that contains the walk data
#define W_LINE_WALK_COL 6

#define P_LINE_NAME_COL 1
#define P_LINE_WALK_COL 2

#define DEFAULT_HAP_LEN 0

// how PanSN fields are organised in a W line
#define PANSN_SAMPLE_COL 1
#define PANSN_HAP_ID_COL 2
#define PANSN_CONTIG_NAME_COL 3

struct ref *alloc_ref(enum gfa_line_prefix line_prefix,
		      struct ref_walk **r_walk, struct ref_id **id)
{
	struct ref *r = malloc(sizeof(struct ref));
	if (!r)
		return NULL;

	r->line_prefix = line_prefix;
	r->walk = *r_walk;
	r->id = *id;

	return r;
}

void destroy_ref(struct ref **r)
{
	if (r == NULL || *r == NULL) {
		return;
	}

	destroy_ref_id(&(*r)->id);
	destroy_ref_walk(&(*r)->walk);

	if (*r != NULL) {
		free(*r);
		*r = NULL;
	}
}

struct ref *get_ref(gfa_props *gfa, idx_t ref_idx)
{
	if (!gfa || ref_idx >= gfa->ref_count) {
		return NULL;
	}

	return gfa->refs[ref_idx];
}

const char *get_tag(const struct ref *r)
{
	if (!r)
		return NULL;

	return r->id->tag;
}

id_t get_hap_id(const struct ref *r)
{
	if (get_ref_id_type(r) == REF_ID_RAW)
		return NULL_ID;

	return r->id->value.id_value->hap_id;
}

const char *get_contig_name(const struct ref *r)
{
	if (get_ref_id_type(r) == REF_ID_RAW)
		return NULL;

	return r->id->value.id_value->contig_name;
}

const char *get_sample_name(const struct ref *r)
{
	if (!r)
		return NULL;

	switch (r->id->type) {
	case REF_ID_PANSN:
		return r->id->value.id_value->sample_name;
	case REF_ID_RAW:
		return get_tag(r);
	}

	return NULL;
}

idx_t get_step_count(const struct ref *r)
{
	return r->walk->step_count;
}

const id_t *get_walk_v_ids(const struct ref *r)
{
	return r->walk->v_ids;
}

const enum strand *get_walk_strands(const struct ref *r)
{
	return r->walk->strands;
}

idx_t get_hap_len(const struct ref *r)
{
	return r->walk->hap_len;
}

enum ref_id_type get_ref_id_type(const struct ref *r)
{
	return r->id->type;
}

enum gfa_line_prefix get_line_prefix(const struct ref *r)
{
	if (!r)
		return (enum gfa_line_prefix) - 1;

	return r->line_prefix;
}

status_t set_hap_len(struct ref *r, idx_t hap_len)
{
	if (!r)
		return ERROR_CODE_INVALID_ARGUMENT;

	r->walk->hap_len = hap_len;

	return SUCCESS;
}

static inline bool is_step_sep(enum gfa_line_prefix line_prefix, const char c)
{
	switch (line_prefix) {
	case P_LINE:
		return c == P_LINE_FORWARD_SYMBOL || c == P_LINE_REVERSE_SYMBOL;
	case W_LINE:
		return c == W_LINE_FORWARD_SYMBOL || c == W_LINE_REVERSE_SYMBOL;
	default:
		log_fatal("Invalid line prefix in is_step_sep");
		exit(1);
	}
}

/**
 * @brief Count the number of steps in a P line path; it is the number
 * of commas
 * + 1
 *
 * @param [in] str the path string
 * @return the number of steps in the path
 */
idx_t count_steps(enum gfa_line_prefix line_prefix, const char *str)
{
	idx_t steps = 0;
	for (; *str; str++)
		if (is_step_sep(line_prefix, *str))
			steps++;

	return steps;
}

struct line_metadata {
	idx_t required_tokens;
	idx_t id_token_count;
	const idx_t *id_token_indices;
	idx_t data_col_index;
	enum gfa_line_prefix line_prefix;
	status_t (*parse_data_line)(const char *data_str, struct ref_walk **w);
};

// Define metadata for W_LINE and P_LINE
static const struct line_metadata metadata[] = {
	[W_LINE] = {.required_tokens = READ_W_LINE_TOKENS,
		    .id_token_count = W_LINE_ID_TOKEN_COUNT,
		    .id_token_indices =
			    (const idx_t[]){PANSN_SAMPLE_COL, PANSN_HAP_ID_COL,
					    PANSN_CONTIG_NAME_COL},
		    .data_col_index = W_LINE_WALK_COL,
		    .line_prefix = W_LINE,
		    .parse_data_line = parse_data_line_w},
	[P_LINE] = {.required_tokens = READ_P_LINE_TOKENS,
		    .id_token_count = P_LINE_ID_TOKEN_COUNT,
		    .id_token_indices = (const idx_t[]){P_LINE_NAME_COL},
		    .data_col_index = P_LINE_WALK_COL,
		    .line_prefix = P_LINE,
		    .parse_data_line = parse_data_line_p}};

// for testing
const struct line_metadata *get_line_metadata(enum gfa_line_prefix prefix)
{
	if (prefix == W_LINE || prefix == P_LINE)
		return &metadata[prefix];

	return NULL;
}

// Consolidated line parsing logic using metadata
struct ref *parse_line_generic(const char *line, u32 len,
			       const struct line_metadata *meta)
{
	char *tokens[MAX_TOKENS] = {NULL};
	struct split_str_params p = {
		// input
		.str = line,
		.up_to = line + len,
		.delimiter = TAB_CHAR,
		.max_splits = meta->required_tokens,
		.fallbacks = (const char[]){NEWLINE, NULL_CHAR},
		.fallback_chars_count = 2,

		// output
		.tokens_found = 0,
		.tokens = tokens,
	};

	// log_info("%d", len);

	status_t s = split_str(&p);
	if (p.tokens_found < meta->required_tokens) {
		log_fatal("Failed to split %d-line. Found %u tokens.",
			  p.delimiter, p.tokens_found);
		return NULL;
	}

	// Extract ID tokens
	char *id_tokens[meta->id_token_count];
	for (idx_t i = 0; i < meta->id_token_count; i++) {
		id_tokens[i] = tokens[meta->id_token_indices[i]];
	}

	const char **tok = (const char **)id_tokens;
	struct ref_id *r_id = alloc_ref_id(tok, meta->id_token_count);
	if (!r_id) {
		log_fatal("Failed to allocate ref_id for %d-line.",
			  meta->line_prefix);
		return NULL;
	}

	// Parse the data string
	const char *data_str = tokens[meta->data_col_index];
	idx_t step_count = count_steps(meta->line_prefix, data_str);
	struct ref_walk *w = alloc_ref_walk(step_count);
	if (!w) {
		log_fatal("Failed to allocate walk for %d-line.",
			  meta->line_prefix);
		destroy_ref_id(&r_id);
		exit(1);
	}

	status_t res = meta->parse_data_line(data_str, &w);
	if (res != SUCCESS) {
		log_error("Failed to parse data for %d-line.",
			  meta->line_prefix);
		return NULL;
	}

	for (size_t i = 0; i < MAX_TOKENS; i++) {
		if (tokens[i] != NULL)
			free(tokens[i]);
		tokens[i] = NULL;
	}
	/* struct ref_walk *w = alloc_ref_walk(0); */

	return alloc_ref(meta->line_prefix, &w, &r_id);
}

struct ref *parse_ref_line(enum gfa_line_prefix prefix, const char *line,
			   u32 len)
{
	switch (prefix) {
	case (P_LINE):
		return parse_line_generic(line, len, &metadata[P_LINE]);
	case (W_LINE):
		return parse_line_generic(line, len, &metadata[W_LINE]);
	default:
		log_error("%s Unsupported line prefix.");
		return NULL;
	}
}

/**
 * @brief a wrapper function for handle_p_lines
 */
void *t_handle_p(void *ref_metadata)
{
	struct ref_thread_data *data = (struct ref_thread_data *)ref_metadata;
	line *pl = data->p_lines;
	line *wl = data->w_lines;
	idx_t p_line_count = data->p_line_count;
	idx_t w_line_count = data->w_line_count;
	struct ref **refs = data->refs;
	idx_t ref_idx = 0;

	for (idx_t i = 0; i < p_line_count; i++)
		refs[ref_idx++] =
			parse_ref_line(P_LINE, pl[i].start, pl[i].len);

	for (idx_t i = 0; i < w_line_count; i++)
		refs[ref_idx++] =
			parse_ref_line(W_LINE, wl[i].start, wl[i].len);

	return NULL;
}
