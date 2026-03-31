#ifndef LQ_REF_WALK_H
#define LQ_REF_WALK_H

#include "../../include/liteseq/refs.h"
#include "../../include/liteseq/types.h"

#ifdef __cplusplus
extern "C" { // Ensure the function has C linkage
namespace liteseq
{
#endif

void destroy_ref_walk(struct ref_walk **w);
struct ref_walk *alloc_ref_walk(idx_t step_count);
status_t parse_data_line_w(const char *str, struct ref_walk **empty_r_walk);
status_t parse_data_line_p(const char *str, struct ref_walk **empty_r_walk);

#ifdef TESTING
struct ref_walk *alloc_ref_walk(idx_t step_count);
void destroy_ref_walk(struct ref_walk **r_walk);

idx_t count_steps(enum gfa_line_prefix line_type, const char *str);
#endif // TESTING

#ifdef __cplusplus
} // namespace liteseq
} // extern "C"
#endif

#endif // LQ_REF_WALK_H
