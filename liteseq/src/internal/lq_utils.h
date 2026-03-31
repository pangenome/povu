#ifndef LQ_UTILS_H
#define LQ_UTILS_H

#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "../include/liteseq/types.h"

#ifdef __cplusplus
extern "C" {

namespace liteseq
{

#endif

/**
 * A 32 bit bitstring where each bit corresponds to a character in the alphabet
 * Bitstring for 'A' to 'Z' (32 bits) which are 26 but padded to 32
 * Map characters relative to 'A', so bit_position = char - 'A'.
 * 'A' → Bit 0, 'C' → Bit 2, 'G' → Bit 6, 'T' → Bit 19, 'N' → Bit 13.
 */
#define ALPHABET_MASK                                                          \
	((uint32_t)((1 << ('A' - 'A')) | /* 'A' */                             \
		    (1 << ('C' - 'A')) | /* 'C' */                             \
		    (1 << ('G' - 'A')) | /* 'G' */                             \
		    (1 << ('T' - 'A')) | /* 'T' */                             \
		    (1 << ('N' - 'A'))	 /* 'N' */                             \
		    ))

/**
 * mark a variable as unused
 */
#define UNUSED(x) (void)(x)

/**
 * Likely and unlikely macros to wrap around __builtin_expect
 */
#ifdef __GNUC__
#define likely(x) __builtin_expect(!!(x), 1)
#define unlikely(x) __builtin_expect(!!(x), 0)
#else
#define likely(x) (x)
#define unlikely(x) (x)
#endif

/**
 * Validates that the character is in the alphabet
 * if not it will print an error message and exit the program
 */
void validate_character(char c);

/**
   Function to check if a character is in the alphabet
*/
bool in_alphabet(char c);

/**
 * Encodes a base character to a 3-bit representation.
 * Map the character to a 3-bit encoding
 */
uint8_t encodeBase(char base);

struct split_str_params {
	// input, not mutated.
	const char *str;   // input string
	const char *up_to; // a pointer to a char in the string which we shall
			   // not read past. if up_to is NULL, the function will
			   // depend on max tokens to stop reading the string
	char delimiter;	   // primary delimiter for split
	const char *fallbacks;	    // fallback characters for split
	idx_t fallback_chars_count; // number of fallback characters in fallback
	idx_t max_splits;	    // the max number of tokens to extract

	// output
	idx_t tokens_found; // output, number of tokens found, mutated
	char **tokens;	    // output, only part changed by the function
	const char *end;    // output, pointer to the end of the last token
};

idx_t count_digits(idx_t num);
void tokens_free(char **tokens, u32 N);
/**
 * Tokenises a line into tokens based on a delimiter.
 * The tokens are allocated using malloc and should be freed by the caller.
 * Returns the number of tokens found, or -1 on error.
 */
status_t split_str(struct split_str_params *params);

status_t split_str2(struct split_str_params *params);

#ifdef __cplusplus

} // liteseq
} // extern "C"
#endif

#endif /* LQ_UTILS_H */
