#include <gtest/gtest.h>

#include <string>

#include "../src/internal/lq_utils.h"
#include <liteseq/types.h>

using namespace liteseq;

TEST(EncodeBase, ValidBases)
{
	EXPECT_EQ(encodeBase('A'), 0b000);
	EXPECT_EQ(encodeBase('T'), 0b001);
	EXPECT_EQ(encodeBase('C'), 0b010);
	EXPECT_EQ(encodeBase('G'), 0b011);
	EXPECT_EQ(encodeBase('N'), 0b100);
	EXPECT_EQ(encodeBase('I'), 0b101); // Invalid base. Internal to far. May
					   // indicate end of sequence
}

TEST(Tokenise, Basic)
{
	const char *input = "token1\ttoken2\ttoken3";
	const char delim = TAB_CHAR;
	const idx_t fallback_chars_count = 2;
	const char fallbacks[fallback_chars_count] = {NEWLINE, NULL_CHAR};
	// idx_t line_length = strlen(input);
	char *tokens[MAX_TOKENS] = {NULL};
	idx_t tokens_found = 0;
	const idx_t expected_tokens_found = 3;
	const char *out_tokens[] = {"token1", "token2", "token3"};

	struct split_str_params p = {
		input, // str: Must be the first member in the struct
		NULL,
		delim,		      // delimiter
		fallbacks,	      // fallbacks
		fallback_chars_count, // fallback_chars_count
		10,		      // max_tokens
		tokens_found,	      // tokens_found
		tokens,		      // tokens
		nullptr		      // end
	};

	status_t res = split_str(&p);

	ASSERT_EQ(res, SUCCESS);
	ASSERT_EQ(p.tokens_found, expected_tokens_found);

	for (idx_t i = 0; i < p.tokens_found; i++) {
		ASSERT_STREQ(p.tokens[i], out_tokens[i]);
	}

	for (size_t i = 0; i < MAX_TOKENS && tokens[i] != NULL; i++) {
		free(tokens[i]);
		tokens[i] = NULL;
	}
}

TEST(Tokenise, NoDelim)
{
	const char *input = "single_token_no_delim";
	const char delim = TAB_CHAR;
	const idx_t fallback_chars_count = 2;
	const char fallbacks[fallback_chars_count] = {NEWLINE, NULL_CHAR};
	// idx_t line_length = strlen(input);
	char *tokens[MAX_TOKENS] = {NULL};
	idx_t tokens_found = 0;
	const idx_t expected_tokens_found = 1;
	const char *out_tokens[] = {"single_token_no_delim"};

	struct split_str_params p = {
		input, // str: Must be the first member in the struct
		NULL,
		delim,		      // delimiter
		fallbacks,	      // fallbacks
		fallback_chars_count, // fallback_chars_count
		10,		      // max_tokens
		tokens_found,	      // tokens_found
		tokens,		      // tokens
		nullptr		      // end
	};

	status_t res = split_str(&p);

	ASSERT_EQ(res, SUCCESS);
	ASSERT_EQ(p.tokens_found, expected_tokens_found);

	for (idx_t i = 0; i < p.tokens_found; i++) {
		ASSERT_STREQ(p.tokens[i], out_tokens[i]);
	}

	for (size_t i = 0; i < MAX_TOKENS && tokens[i] != NULL; i++) {
		free(tokens[i]);
		tokens[i] = NULL;
	}
}

TEST(Tokenise, SLine)
{
	const char *input = "S\t1\tAT\nL\t1\t2\t+\t0M";
	auto input_str = std::string(input);
	size_t newline_pos = input_str.find_first_of(NEWLINE);
	const char delim = TAB_CHAR;
	const idx_t fallback_chars_count = 0;
	const char *fallbacks = "";
	const idx_t EXPECTED_S_LINE_TOKENS = 3;
	// idx_t line_length = strlen(input);
	char *tokens[EXPECTED_S_LINE_TOKENS] = {NULL};
	idx_t tokens_found = 0;
	const idx_t expected_tokens_found = EXPECTED_S_LINE_TOKENS;
	const char *out_tokens[EXPECTED_S_LINE_TOKENS] = {"S", "1", "AT"};

	struct split_str_params p = {
		input, // str: Must be the first member in the struct
		input + newline_pos,	// up_to
		delim,			// delimiter
		fallbacks,		// fallbacks
		fallback_chars_count,	// fallback_chars_count
		EXPECTED_S_LINE_TOKENS, // max_tokens
		tokens_found,		// tokens_found
		tokens,			// tokens
		nullptr			// end
	};

	status_t res = split_str(&p);
	ASSERT_EQ(res, SUCCESS);
	ASSERT_EQ(p.tokens_found, expected_tokens_found);
	for (idx_t i = 0; i < p.tokens_found; i++) {
		ASSERT_STREQ(p.tokens[i], out_tokens[i]);
	}
}
