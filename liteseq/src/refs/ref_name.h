#ifndef LQ_REF_NAME_H
#define LQ_REF_NAME_H

#include "../../include/liteseq/refs.h"
#include "../../include/liteseq/types.h"

#ifdef __cplusplus
extern "C" { // Ensure the function has C linkage
namespace liteseq
{
#endif

void destroy_ref_id(struct ref_id **r_id);
struct ref_id *alloc_ref_id(const char **tokens, idx_t token_count);

#ifdef TESTING
struct pansn *try_extract_pansn_from_str(const char *name, const char delim);
void destroy_pansn(struct pansn **pn);
struct pansn *alloc_pansn(const char **tokens);

#endif // TESTING

#ifdef __cplusplus
} // namespace liteseq
} // extern "C"
#endif

#endif // LQ_REF_NAME_H
