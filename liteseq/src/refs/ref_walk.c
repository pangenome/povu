#include <log.h>
#include <stdlib.h>
#include <string.h> // for memset

#include "../../include/liteseq/refs.h"

void destroy_ref_walk(struct ref_walk **w)
{

	if (w == NULL || *w == NULL)
		return;

	if ((*w)->strands != NULL) {
		free((*w)->strands);
		(*w)->strands = NULL;
	}

	if ((*w)->v_ids != NULL) {
		free((*w)->v_ids);
		(*w)->v_ids = NULL;
	}

	if ((*w)->loci != NULL) {
		free((*w)->loci);
		(*w)->loci = NULL;
	}

	if (*w != NULL) {
		free(*w);
		*w = NULL;
	}
}

struct ref_walk *alloc_ref_walk(idx_t step_count)
{
	struct ref_walk *w = malloc(sizeof(struct ref_walk));
	if (!w)
		return NULL;

	w->strands = NULL;
	w->v_ids = NULL;
	w->loci = NULL;

	w->strands = malloc(sizeof(enum strand) * step_count);
	if (!w->strands) {
		destroy_ref_walk(&w);
		return NULL;
	}

	w->v_ids = malloc(sizeof(id_t) * step_count);
	if (!w->v_ids) {
		destroy_ref_walk(&w);
		return NULL;
	}

	w->loci = malloc(sizeof(idx_t) * step_count);
	if (!w->loci) {
		destroy_ref_walk(&w);
		return NULL;
	}

	w->step_count = step_count;
	w->hap_len = 0; // default to 0

	return w;
}

status_t parse_data_line_w(const char *str, struct ref_walk **empty_r_walk)
{
	if (str == NULL) {
		log_fatal("Input string is NULL");
		return ERROR_CODE_INVALID_ARGUMENT;
	}

	// we assume the payload is already allocated and empty
	// we just fill it
	struct ref_walk *w = *empty_r_walk;

	char str_num[MAX_DIGITS] = {0};
	size_t digit_pos = 0; // current position in str_num
	size_t v_id = 0;
	idx_t step_count = 0;
	enum strand curr_or;

	for (; *str; str++) {
		switch (*str) {
		case W_LINE_FORWARD_SYMBOL:
		case W_LINE_REVERSE_SYMBOL:
			// clang-format off
			if (step_count > 0) {
				id_t v_id = strtoul(str_num, NULL, 10);
				w->v_ids[step_count - 1] = strtoul(str_num, NULL, 10);
				w->strands[step_count - 1] = curr_or;
			}
			curr_or = (*str == W_LINE_FORWARD_SYMBOL) ? STRAND_FWD : STRAND_REV;
			// clang-format on

			step_count++;
			digit_pos = 0;
			memset(str_num, 0, sizeof(char) * MAX_DIGITS);
			break;
		default: // accumulate the str
			if (digit_pos >= MAX_DIGITS) {
				log_fatal("Vertex ID exceeds maximum length of "
					  "%d digits",
					  MAX_DIGITS);
				return ERROR_CODE_OUT_OF_BOUNDS;
			}
			str_num[digit_pos++] = *str;
		}
	}

	// Handle trailing token if valid
	if (digit_pos > 0) {
		w->v_ids[step_count - 1] = strtoul(str_num, NULL, 10);
		w->strands[step_count - 1] = curr_or;
	}

	return SUCCESS;
}

status_t parse_data_line_p(const char *str, struct ref_walk **empty_r_walk)
{
	if (str == NULL) {
		log_fatal("Input string is NULL");
		return ERROR_CODE_INVALID_ARGUMENT;
	}

	// we assume the payload is already allocated and empty
	// we just fill it
	struct ref_walk *w = *empty_r_walk;

	char str_num[MAX_DIGITS] = {0};
	size_t digit_pos = 0; // current position in str_num
	size_t v_id = 0;
	idx_t step_count = 0;

	for (; *str; str++) {
		switch (*str) {
		case P_LINE_FORWARD_SYMBOL:
			v_id = strtoul(str_num, NULL, 10);
			w->v_ids[step_count] = v_id;
			w->strands[step_count] = STRAND_FWD;
			break;
		case P_LINE_REVERSE_SYMBOL:
			v_id = strtoul(str_num, NULL, 10);
			w->v_ids[step_count] = v_id;
			w->strands[step_count] = STRAND_REV;
			break;
		case COMMA_CHAR:
			step_count++;
			digit_pos = 0;
			memset(str_num, 0, sizeof(char) * MAX_DIGITS);
			break;
		default: // accumulate the str
			str_num[digit_pos++] = *str;
			if (digit_pos >= MAX_DIGITS) {
				log_fatal("Vertex ID %s exceeds maximum length "
					  "of "
					  "%d digits",
					  str_num, MAX_DIGITS);
				return EXIT_FAILURE;
			}
		}
	}

	return SUCCESS;
}
