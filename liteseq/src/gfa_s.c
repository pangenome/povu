#include <pthread.h>

#include <log.h>
#include <stdlib.h>

#include "./gfa_s.h"

#include "../include/liteseq/gfa.h"
#include "../include/liteseq/types.h"
#include "../src/internal/lq_utils.h"

#define EXPECTED_S_LINE_TOKENS 3 // the number of tokens expected in a S line
#define S_LINE_TYPE_IDX 0 // the index of the line type token in the S line
#define S_LINE_V_ID_IDX 1 // the index of the vertex ID token in the S line
#define S_LINE_SEQ_IDX 2  // the index of the sequence token in the S line

vtx *get_vtx(gfa_props *gfa, id_t v_id)
{
	return gfa->v[v_id];
}

status_t handle_s(const char *s_line, u32 line_len, char **tokens,
		  bool inc_vtx_labels, vtx **vertices)
{
	struct split_str_params p = {
		.str = s_line,
		.up_to = s_line + line_len,
		.delimiter = TAB_CHAR,
		.fallbacks = "",
		.fallback_chars_count = 0,
		.max_splits = 3,

		.tokens_found = 0,
		.tokens = tokens,
		.end = NULL,
	};

	status_t res = split_str(&p);
	if (res != SUCCESS) {
		log_fatal("Could not parse S line");
		return res;
	}

	vtx *v = malloc(sizeof(vtx));
	if (v == NULL) {
		log_fatal("Could not allocate memory for vertex");
		return FAILURE;
	}
	v->id = strtoul(tokens[S_LINE_V_ID_IDX], NULL, 10);
	v->seq = inc_vtx_labels ? tokens[S_LINE_SEQ_IDX] : NULL;
	/* vtx v = {.id = strtoul(tokens[S_LINE_V_ID_IDX], NULL, 10), */
	/*	 .seq = inc_vtx_labels ? tokens[S_LINE_SEQ_IDX] : NULL}; */
	vertices[v->id] = v;

	free(tokens[S_LINE_TYPE_IDX]); // free the line type token
	free(tokens[S_LINE_V_ID_IDX]); // free the vertex ID token

	// free the sequence token only if vertex labels are not included
	if (!inc_vtx_labels)
		free(tokens[S_LINE_SEQ_IDX]);

	return SUCCESS;
}

/**
 * @brief a wrapper function for handle_s_lines
 */
void *t_handle_s(void *s_meta)
{
	struct s_thread_meta *meta = (struct s_thread_meta *)s_meta;
	vtx **vtxs = meta->vertices;
	line *sl = meta->s_lines;
	idx_t line_count = meta->s_line_count;
	bool inc_vtx_labels = meta->inc_vtx_labels;

	// temporary storage for the tokens extracted from a given line
	char *tokens[EXPECTED_S_LINE_TOKENS] = {NULL};

	for (idx_t i = 0; i < line_count; i++)
		handle_s(sl[i].start, sl[i].len, tokens, inc_vtx_labels, vtxs);

	return NULL;
}
