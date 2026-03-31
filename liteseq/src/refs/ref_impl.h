#ifndef LQ_REF_IMPL_H
#define LQ_REF_IMPL_H

#include <pthread.h>

#include "../../include/liteseq/refs.h"
#include "../../include/liteseq/types.h"
#include "../include/liteseq/gfa.h"

#ifdef __cplusplus
extern "C" { // Ensure the function has C linkage
namespace liteseq
{
#endif

#define P_LINE_ID_TOKEN_COUNT 1
#define W_LINE_ID_TOKEN_COUNT 3

struct ref_thread_data {
	struct ref **refs;
	line *p_lines; // metadata for a P line
	line *w_lines; // metadata for a W line
	idx_t p_line_count;
	idx_t w_line_count;
};

void *t_handle_p(void *ref_metadata);

// fns I want to expose only for testing
#ifdef TESTING

struct ref *alloc_ref(enum gfa_line_prefix line_prefix,
		      struct ref_walk **r_walk, struct ref_id **id);

/**
 * Function to retrieve metadata for a given line prefix.
 * This function provides access to line-specific metadata, which includes
 * required token counts, token indices, and parsing logic. It is primarily
 * used for testing purposes or when controlled access to metadata is needed.
 *
 * @param prefix The line prefix for which metadata is requested (e.g., W_LINE,
 * P_LINE).
 * @return A pointer to the corresponding metadata, or NULL if the prefix is
 * unsupported.
 */
const struct line_metadata *get_line_metadata(enum gfa_line_prefix prefix);
#endif

#ifdef __cplusplus
} // namespace liteseq
} // extern "C"
#endif

#endif // LQ_REF_IMPL_H
