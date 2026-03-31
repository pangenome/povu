#include <pthread.h>
#include <stdlib.h>

#include <log.h>

#include "../include/liteseq/gfa.h"
#include "../include/liteseq/types.h"
#include "../src/internal/lq_utils.h"
#include "./gfa_l.h"

#define L_LINE_TYPE_IDX 0      // the index of the line type token in the L line
#define L_LINE_V1_ID_IDX 1     //  first vertex ID token in the L line
#define L_LINE_V1_STRAND_IDX 2 // first vertex strand token in the L line
#define L_LINE_V2_ID_IDX 3     // second vertex ID token in the L line
#define L_LINE_V2_STRAND_IDX 4 // second vertex strand token in the L line

/**
 * In a self loop the source (src) and sink (snk) are the same value
 * A self loop can be in the forward, reverse, or mixed strand
 *
 * Examples:
 * L 1 + 1 + is a forward self loop
 * L 1 - 1 - is a reverse self loop
 *
 * L 1 + 1 -  and L 1 - 1 + are mixed self loops that aren't
 * representable in a bidirected graph without node duplication
 *
 * @param [in] l_line the line to parse
 * @param [in] line_length the length of the line
 * @param [in] idx the index of the line
 * @param [in] tokens the tokens to parse
 * @param [in] c the config
 * @param [out] e the edges to populate
 * @return 0 on success, -1 on failure
 */
status_t handle_l(const char *l_line, u32 line_len, size_t idx, char **tokens,
		  edge *edges)
{
	struct split_str_params p = {
		.str = l_line,
		.up_to = NULL,
		.delimiter = TAB_CHAR,
		.fallbacks = "",
		.fallback_chars_count = 0,
		.max_splits = EXPECTED_L_LINE_TOKENS,

		.tokens_found = 0,
		.tokens = tokens,
		.end = NULL,
	};

	status_t res = split_str(&p);

	// Parse vertex IDs
	size_t v1_id = strtoul(tokens[L_LINE_V1_ID_IDX], NULL, 10);
	size_t v2_id = strtoul(tokens[L_LINE_V2_ID_IDX], NULL, 10);

	// parse strand symbols. A strand symbol is either + or -
	char v1_strand_symbol = tokens[L_LINE_V1_STRAND_IDX][0];
	char v2_strand_symbol = tokens[L_LINE_V2_STRAND_IDX][0];

	// Determine vertex sides based on strand symbols
	vtx_side_e v1_side, v2_side;
	if (unlikely(v1_id == v2_id)) { // check for self loop
		if (v1_strand_symbol != v2_strand_symbol) {
			fprintf(stderr,
				"Error: Invalid self loop: %ld %c and %ld %c\n",
				v1_id, v1_strand_symbol, v2_id,
				v2_strand_symbol);
			return -1;
		} else if (v1_strand_symbol == v2_strand_symbol) {
			v1_side = LEFT;
			v2_side = RIGHT;
		}
	} else {
		v1_side = (v1_strand_symbol == '+') ? RIGHT : LEFT;
		v2_side = (v2_strand_symbol == '+') ? LEFT : RIGHT;
	}

	// populate the edge
	edges[idx] = (edge){.v1_id = v1_id,
			    .v1_side = v1_side,
			    .v2_id = v2_id,
			    .v2_side = v2_side};

	tokens_free(tokens, EXPECTED_L_LINE_TOKENS);

	return 0;
}

/**
 * @brief a wrapper function for handle_l_lines
 */
void *t_handle_l(void *l_meta)
{
	struct l_thread_meta *meta = (struct l_thread_meta *)l_meta;
	edge *edges = meta->edges;
	line *ll = meta->l_lines;
	idx_t line_count = meta->l_line_count;

	// temporary storage for the tokens extracted from a given line
	char *tokens[EXPECTED_L_LINE_TOKENS] = {NULL};

	for (idx_t i = 0; i < line_count; i++)
		handle_l(ll[i].start, ll[i].len, i, tokens, edges);

	return NULL;
}
