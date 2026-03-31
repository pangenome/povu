#ifndef LQ_REFS_H
#define LQ_REFS_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h> // uint8_t, uint32_t

#include "./types.h"

#ifdef __cplusplus
extern "C" { // Ensure the function has C linkage
namespace liteseq
{
#endif

#define P_LINE_FORWARD_SYMBOL '+'
#define P_LINE_REVERSE_SYMBOL '-'
#define W_LINE_FORWARD_SYMBOL '>'
#define W_LINE_REVERSE_SYMBOL '<'

/* ref name related types */
struct pansn {
	char *sample_name;
	id_t hap_id;
	char *contig_name;
};

enum ref_id_type {
	REF_ID_PANSN,
	REF_ID_RAW,
};

/**
 * @brief: A union to hold either a parsed PanSN structure or a raw string.
 *
 * Note to future implementers:
 * If adding new types to this union, ensure that the new type implements
 * the necessary operations or mechanisms for tagging. Proper tagging support
 * must be maintained for all reference ID representations.
 */
union ref_id_value {
	struct pansn *id_value;
	char *raw;
};

struct ref_id {
	enum ref_id_type type;
	union ref_id_value value;
	char *tag;
};

/* ref walk related types */
struct ref_walk {
	// the actual walk
	enum strand *strands;
	id_t *v_ids;
	idx_t *loci; // the locus of the first base of each step
	// walk metadata
	idx_t step_count; // the number of steps
	idx_t hap_len;	  // the length of the haplotype in bases
};

/* ref itself */
struct ref {
	enum gfa_line_prefix line_prefix;
	struct ref_walk *walk;
	struct ref_id *id;
};
struct ref *parse_ref_line(enum gfa_line_prefix line_type, const char *line,
			   u32 len);
void destroy_ref(struct ref **r);

/*
 * ------------------
 * Accessor functions
 * ------------------
 */

/* ref name related */
const char *get_tag(const struct ref *r);
const char *get_sample_name(const struct ref *r);
idx_t get_hap_id(const struct ref *r);
const char *get_contig_name(const struct ref *r);
enum gfa_line_prefix get_line_prefix(const struct ref *r);
enum ref_id_type get_ref_id_type(const struct ref *r);

/* ref walk related */
idx_t get_hap_len(const struct ref *r);
status_t set_hap_len(struct ref *r, idx_t hap_len);
idx_t get_step_count(const struct ref *r);
const id_t *get_walk_v_ids(const struct ref *r);
const enum strand *get_walk_strands(const struct ref *r);

#ifdef __cplusplus
} // namespace liteseq
} // extern "C"
#endif

#endif /* LQ_REFS_H */
