#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <log/log.h>
#include <math.h>

#include "../../include/liteseq/types.h"
#include "./lq_utils.h"

uint8_t encodeBase(char base)
{
	switch (base) {
	case 'A':
		return 0b000;
	case 'T':
		return 0b001;
	case 'C':
		return 0b010;
	case 'G':
		return 0b011;
	case 'N':
		return 0b100;
	case 'I':
		return 0b101; // Invalid base. Internal to far. May indicate end
			      // of sequence
	default:
		fprintf(stderr, "Invalid base character: %c\n", base);
		exit(EXIT_FAILURE); // Exit the program on invalid input
	}
}

#ifdef NDEBUG
inline bool in_alphabet(char c)
{
	if (c < 'A' || c > 'Z')
		return false; // Out of range for uppercase letters
	return (ALPHABET_MASK & (1 << (c - 'A'))) != 0;
}
#else
bool in_alphabet(char c)
{
	if (c < 'A' || c > 'Z')
		return false; // Out of range for lowercase letters
	return (ALPHABET_MASK & (1 << (c - 'A'))) != 0;
}
#endif

void validate_character(char c)
{
	if (!in_alphabet(c)) {
		fprintf(stderr, "character (%c), ASCII (%d)", c, c);
		perror("Invalid character in sequence\n");
		exit(EXIT_FAILURE); // Exit on invalid input
	}
}

struct match_result {
	const char *delim_ptr;
	int len;
	bool at_fallback;
};

#define NULL_MACTH_RES                                                         \
	(struct match_result)                                                  \
	{                                                                      \
		NULL, -1, false                                                \
	}

const char *strchr_bounded(const char *start, char delim)
{
	const char *curr = start;
	while (true) {
		if (*curr == delim || *curr == '\0' || *curr == '\n')
			return curr;
		curr++;
	}

	return NULL;
}

struct match_result find_delim(const char *start, char c,
			       const char *fallback_chars,
			       int fallback_chars_count)
{
	const char *fn = "[liteseq::utils::find_delim]";

	char *res;
	res = strchr_bounded(start, c); // Locate the position of the token
	/* res = strchr(start, c); // Locate the position of the token */
	if (res != NULL)
		return (struct match_result){res, (int)(res - start), false};

	if (!fallback_chars)
		return NULL_MACTH_RES;

	// If not found, look for any of the fallback characters
	for (int i = 0; i < fallback_chars_count; i++) {
		res = strchr(start, fallback_chars[i]);
		if (res != NULL)
			return (struct match_result){res, (int)(res - start),
						     true};
	}

	return NULL_MACTH_RES;
}

status_t split_str2(struct split_str_params *p)
{
	const char *fn = "[liteseq::utils::split_str]";

	const char *str = p->str;
	const char *up_to = p->up_to;
	const char c = p->delimiter;
	const char *fallbacks = p->fallbacks;
	idx_t fallback_chars_count = p->fallback_chars_count;
	idx_t max_tokens = p->max_splits;
	char **all_tokens = p->tokens;

	idx_t tokens_found = 0;
	int len = 0;

	bool limit_reached = false;

	int len_read = 0;
	int max_len = p->up_to - p->str;

	while (len != -1 && tokens_found < max_tokens) {
		struct match_result match_res =
			find_delim(str, c, fallbacks, fallback_chars_count);
		len = match_res.len;

		if (len == -1)
			break;

		len_read += len + 1;

		if (up_to != NULL && match_res.delim_ptr > up_to) {
			limit_reached = true;
			len = (int)(up_to - str);
		}

		if (len <= 0)
			break;

		char *tok = malloc(sizeof(char) * len + 1);
		if (tok == NULL) {
			log_error("%s Memory allocation failed\n", fn);
			return ERROR_CODE_FAILURE;
		}

		memcpy(tok, str, len);
		tok[len] = '\0';

		if (len_read > max_len) {
			log_error("%s Attempted to read past the up_to limit\n",
				  fn);
			return ERROR_CODE_FAILURE;
		} else if (max_len == 13371958) {
			printf("%s\n", tok);
			log_info("%d", len_read);
		}

		all_tokens[tokens_found++] = tok;

		str += len + 1;
		if (match_res.at_fallback || limit_reached)
			break;
	}

	p->tokens_found = tokens_found;
	p->end = str;

	return SUCCESS;
}

status_t split_str(struct split_str_params *p)
{
	const char *fn = "[liteseq::utils::split_str]";

	const char *str = p->str;
	const char *up_to = p->up_to;
	const char c = p->delimiter;
	const char *fallbacks = p->fallbacks;
	idx_t fallback_chars_count = p->fallback_chars_count;
	idx_t max_tokens = p->max_splits;
	char **all_tokens = p->tokens;

	idx_t tokens_found = 0;
	int len = 0;

	bool limit_reached = false;

	while (len != -1 && tokens_found < max_tokens) {
		struct match_result match_res =
			find_delim(str, c, fallbacks, fallback_chars_count);
		len = match_res.len;

		if (len == -1)
			break;

		if (up_to != NULL && match_res.delim_ptr > up_to) {
			limit_reached = true;
			len = (int)(up_to - str);
		}

		if (len <= 0)
			break;

		char *tok = malloc(sizeof(char) * len + 1);
		if (tok == NULL) {
			log_error("%s Memory allocation failed\n", fn);
			return ERROR_CODE_FAILURE;
		}

		memcpy(tok, str, len);
		tok[len] = '\0';

		all_tokens[tokens_found++] = tok;

		str += len + 1;
		if (match_res.at_fallback || limit_reached)
			break;
	}

	p->tokens_found = tokens_found;
	p->end = str;

	return SUCCESS;
}

idx_t count_digits(idx_t num)
{
	if (num == 0)
		return 1;
	return (idx_t)log10(num) + 1;
}

void tokens_free(char **tokens, u32 N)
{
	for (size_t i = 0; i < N && tokens[i] != NULL; i++) {
		free(tokens[i]);
		tokens[i] = NULL;
	}
	// TODO is a memset here faster or better suited in case of errors?
}
