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

#include <log.h>
#include <stdio.h> // for snprintf
#include <stdlib.h>
#include <string.h> // for strlen

#include "../../include/liteseq/refs.h"
#include "../internal/lq_utils.h"

const char DELIM = HASH_CHAR;

/* Define constants globally or within the file */
#define PANSN_MAX_TOKENS 3
#define PANSN_PARSER_BUF_SIZE PANSN_MAX_TOKENS + 1

// how PanSN fields are organised in a the tokens passed to alloc_pansn
// different from the ones in the W line which are offset +1
#define PANSN_SAMPLE_COL 0
#define PANSN_HAP_ID_COL 1
#define PANSN_CONTIG_NAME_COL 2

void destroy_pansn(struct pansn **pn)
{
	if (pn == NULL || *pn == NULL)
		return;

	if ((*pn)->sample_name != NULL) {
		free((*pn)->sample_name);
		(*pn)->sample_name = NULL;
	}

	if ((*pn)->contig_name != NULL) {
		free((*pn)->contig_name);
		(*pn)->contig_name = NULL;
	}

	if (*pn != NULL) {
		free(*pn);
		*pn = NULL;
	}
}

struct pansn *alloc_pansn(const char *tokens[PANSN_MAX_TOKENS])
{
	struct pansn *pn = malloc(sizeof(struct pansn));
	if (!pn)
		return NULL;

	const char *sn = tokens[PANSN_SAMPLE_COL];
	const char *h = tokens[PANSN_HAP_ID_COL];
	const char *cn = tokens[PANSN_CONTIG_NAME_COL];

	pn->sample_name = sn ? strdup(sn) : NULL;
	pn->hap_id = h ? atol(h) : NULL_ID;
	pn->contig_name = cn ? strdup(cn) : NULL;

	if (!pn->sample_name || !pn->contig_name || pn->hap_id == NULL_ID) {
		destroy_pansn(&pn);
		return NULL;
	}

	return pn;
}

char *alloc_pansn_tag(const struct pansn *pn)
{
	// idx_t hap_id_len = (idx_t)log10(pn->hap_id);
	idx_t hap_id_len = count_digits(pn->hap_id);

	char hap_id_str[MAX_DIGITS];
	snprintf(hap_id_str, sizeof hap_id_str, "%u", pn->hap_id);

	size_t need =
		strlen(pn->sample_name) + hap_id_len + strlen(pn->contig_name);
	need += 2; // +2 for the two '#' characters
	need += 1; // +1 for the null terminator

	char *tag = (char *)malloc(need);
	if (!tag)
		return NULL;

	sprintf(tag, "%s#%u#%s", pn->sample_name, pn->hap_id, pn->contig_name);

	return tag;
}

struct pansn *free_pansn_buf(char **out_tokens)
{
	for (idx_t i = 0; i < PANSN_MAX_TOKENS; i++) {
		if (out_tokens[i] != NULL) {
			free(out_tokens[i]);
			out_tokens[i] = NULL;
		}
	}
	return NULL;
}

/* used for PanSN in P lines */
struct pansn *try_extract_pansn_from_str(const char *name, const char delim)
{
	char *out_tokens[PANSN_MAX_TOKENS] = {NULL};
	struct split_str_params p = {
		// input
		.str = name,
		.up_to = NULL,
		.delimiter = delim,
		.max_splits = PANSN_MAX_TOKENS,
		.fallbacks = (const char[]){NEWLINE, NULL_CHAR},
		.fallback_chars_count = 2,
		// output
		.tokens_found = 0,
		.tokens = out_tokens,
	};
	status_t res = split_str(&p);
	if (res != SUCCESS)
		return free_pansn_buf(out_tokens);

	if (p.tokens_found != PANSN_MAX_TOKENS) {
		return free_pansn_buf(out_tokens);
	}

	idx_t offset = 0;
	int sample_idx = PANSN_SAMPLE_COL - offset;
	int contig_idx = PANSN_CONTIG_NAME_COL - offset;
	int hap_id_idx = PANSN_HAP_ID_COL - offset;

	for (int i = 0; i < PANSN_MAX_TOKENS; i++) {
		if (strlen(out_tokens[i]) == 0) {
			return free_pansn_buf(out_tokens);
		}
	}

	char *endptr;

	strtoul(out_tokens[hap_id_idx], &endptr, 10);
	if (*endptr != NULL_CHAR) {
		return free_pansn_buf(out_tokens);
	}

	if (strchr(out_tokens[contig_idx], delim) != NULL) {
		return free_pansn_buf(out_tokens);
	}

	const char **o = (const char **)out_tokens;
	struct pansn *pn = alloc_pansn(o);
	free_pansn_buf(out_tokens);

	return pn;
}

struct pansn *try_create_pansn(const char **id_tokens, idx_t token_count,
			       char delim)
{
	if (token_count == 1) {
		const char *ref_name = id_tokens[PANSN_SAMPLE_COL];
		return try_extract_pansn_from_str(ref_name, DELIM);
	} else if (token_count == 3) {
		return alloc_pansn(id_tokens);
	}

	return NULL;
}

struct ref_id *alloc_ref_id(const char **id_tokens, idx_t token_count)
{
	struct ref_id *r_id = malloc(sizeof(struct ref_id));
	if (!r_id)
		return NULL;

	struct pansn *pn = try_create_pansn(id_tokens, token_count, DELIM);

	if (pn) {
		r_id->type = REF_ID_PANSN;
		r_id->value.id_value = pn;
		// create the tag
		r_id->tag = alloc_pansn_tag(pn);
		if (!r_id->tag) {
			free(r_id);
			destroy_pansn(&pn);
			return NULL;
		}
	} else {
		r_id->type = REF_ID_RAW;
		r_id->value.raw = strdup(id_tokens[PANSN_SAMPLE_COL]);
		if (!r_id->value.raw) {
			free(r_id);
			return NULL;
		}
		// for raw ids, tag is the raw string
		r_id->tag = r_id->value.raw;
	}

	return r_id;
}

void destroy_ref_id(struct ref_id **r_id)
{
	if (r_id == NULL || *r_id == NULL)
		return;

	if ((*r_id)->type == REF_ID_PANSN) {
		if ((*r_id)->value.id_value->sample_name != NULL) {
			free((*r_id)->value.id_value->sample_name);
			(*r_id)->value.id_value->sample_name = NULL;
		}
		if ((*r_id)->value.id_value->contig_name != NULL) {
			free((*r_id)->value.id_value->contig_name);
			(*r_id)->value.id_value->contig_name = NULL;
		}
		if ((*r_id)->tag != NULL) {
			free((*r_id)->tag);
			(*r_id)->tag = NULL;
		}

		if ((*r_id)->value.id_value != NULL) {
			free((*r_id)->value.id_value);
			(*r_id)->value.id_value = NULL;
		}

	} else if ((*r_id)->type == REF_ID_RAW) {
		if ((*r_id)->value.raw != NULL) {
			free((*r_id)->value.raw);
			(*r_id)->value.raw = NULL;
		}
	}

	if (*r_id != NULL) {
		free(*r_id);
		*r_id = NULL;
	}
}
