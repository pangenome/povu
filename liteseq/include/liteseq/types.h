#ifndef LQ_TYPES_H
#define LQ_TYPES_H

#include <stddef.h>
#include <stdint.h> // uint8_t, uint32_t

#ifdef __cplusplus
extern "C" {
namespace liteseq
{
#endif

/*
 * =======
 * Aliases
 * =======
 */
typedef uint8_t byte;
typedef uint32_t idx_t;
typedef uint32_t id_t;
typedef uint32_t u32;
typedef int8_t status_t;

/*
 * ============================================
 * Constants for null or invalid sentinel values
 * ============================================
 */
#define NULL_ID UINT32_MAX
#define NULL_IDX UINT32_MAX
#define INVALID_LEN UINT32_MAX

/*
 * ========================
 * Common character values
 * ========================
 */
#define COMMA_CHAR ','
#define NEWLINE '\n'
#define TAB_CHAR '\t'
#define HASH_CHAR '#'
#define NULL_CHAR '\0'

/*
 * ========================================
 * Error codes used throughout the library.
 * ========================================
 */
typedef enum {
	ERROR_CODE_SUCCESS = 0,
	ERROR_CODE_FAILURE = -1,
	ERROR_CODE_NULL_POINTER = -2,
	ERROR_CODE_OUT_OF_MEMORY = -3,
	ERROR_CODE_OUT_OF_BOUNDS = -6,
	ERROR_CODE_INVALID_ARGUMENT = -7,
	ERROR_CODE_NOT_FOUND = -10,
	ERROR_CODE_NOT_IMPLEMENTED = -15,
	ERROR_CODE_UNKNOWN = -17
} ERROR_CODE;

#define SUCCESS ERROR_CODE_SUCCESS
#define FAILURE ERROR_CODE_FAILURE

/*
 * =========================
 * Parsing related constants
 * =========================
 */
#define BUFFER_SIZE (1024 * 1024) // Default buffer size (1 MB)
#define EXPECTED_HEADER_LENGTH 512
#define EXPECTED_SEQ_LENGTH (1 << 30) // 2^30 ≈ 10^9
#define MAX_COMPRESSED_RUN 31
#define MAX_DIGITS 12 // Maximum number of digits in a number
#define MAX_TOKENS 10 // Maximum tokens extracted from a line (GFA)

/*
 * ======================
 * GFA specific constants
 * ======================
 */
// GFA line types
#define GFA_H_LINE 'H' // Header
#define GFA_S_LINE 'S' // Segment
#define GFA_L_LINE 'L' // Link
#define GFA_P_LINE 'P' // Path
#define GFA_W_LINE 'W' // Walk (GFA >= 1.1)

// Expected token counts for different GFA line types
#define EXPECTED_P_LINE_TOKENS 3
#define EXPECTED_L_LINE_TOKENS 5
#define EXPECTED_H_LINE_TOKENS 2

/*
 * ============
 * Common enums
 * ============
 */
enum gfa_line_prefix {
	P_LINE,
	W_LINE,
};

enum strand {
	STRAND_FWD, // aka '+'
	STRAND_REV  // aka '-'
};

#ifdef __cplusplus
} // liteseq
} // extern "C"
#endif

#endif /* LQ_TYPES_H */
