#include <gtest/gtest.h>

#include <array>
#include <liteseq/refs.h>
#include <liteseq/types.h>

#include "../src/refs/ref_impl.h"
#include "../src/refs/ref_name.h"
#include "../src/refs/ref_walk.h"

using namespace liteseq;

const char DELIM = HASH_CHAR;
#define PANSN_MAX_TOKENS 3

// Pads an initializer_list to a fixed size with NULL_ID
template <typename T, std::size_t MAX_COLUMNS>
constexpr std::array<T, MAX_COLUMNS>
pad_to_max(const std::initializer_list<T> &input)
{
	std::array<T, MAX_COLUMNS> padded = {};
	auto it = input.begin();

	// Fill the array with input values
	for (std::size_t i = 0; i < input.size() && i < MAX_COLUMNS; ++i) {
		padded[i] = *it;
		++it;
	}

	// Pad the remainder with NULL_ID
	for (std::size_t i = input.size(); i < MAX_COLUMNS; ++i) {
		padded[i] = static_cast<T>(NULL_ID); // Ensure type consistency
	}
	return padded;
}

TEST(HapLen, SetsHapLen)
{
	const idx_t STEP_COUNT = 5;
	struct ref_walk *rw = alloc_ref_walk(STEP_COUNT);

	// set directly
	rw->hap_len = 42;

	const char *tokens[1] = {"sample#1#contig"};
	struct ref_id *r_id = alloc_ref_id(tokens, P_LINE_ID_TOKEN_COUNT);

	struct ref *r = alloc_ref(P_LINE, &rw, &r_id);
	ASSERT_EQ(get_hap_len(r), 42);

	idx_t res = set_hap_len(r, 100);
	ASSERT_EQ(res, SUCCESS);
	ASSERT_EQ(get_hap_len(r), 100);

	destroy_ref(&r);
	ASSERT_EQ(r, nullptr);
}

TEST(PanSNStruct, AllocAndFree)
{
	const char *tokens[PANSN_MAX_TOKENS] = {"sampleA", "0", "contig_5"};
	struct pansn *pn = alloc_pansn(tokens);
	ASSERT_NE(pn, nullptr);
	ASSERT_STREQ(pn->sample_name, "sampleA");
	ASSERT_EQ(pn->hap_id, 0);
	ASSERT_STREQ(pn->contig_name, "contig_5");

	destroy_pansn(&pn);
	ASSERT_EQ(pn, nullptr);
}

// Demonstrate some basic assertions.
TEST(ParsePanSNRefId, Valid)
{

	const idx_t N = 3;

	// input
	const char *names[] = {"chm13#0#Chr1", "chm13#1#ChrX",
			       "sampleA#2#contig_5"};

	// expected output
	const char *expected_samples[] = {"chm13", "chm13", "sampleA"};
	const char *expected_contigs[] = {"Chr1", "ChrX", "contig_5"};
	const id_t expected_hap_ids[] = {0, 1, 2};

	for (idx_t i = 0; i < N; i++) {
		// char *name = strdup(names[i]);
		pansn *pn = try_extract_pansn_from_str(names[i], DELIM);
		ASSERT_NE(pn, nullptr);
		ASSERT_STREQ(pn->sample_name, expected_samples[i]);
		ASSERT_STREQ(pn->contig_name, expected_contigs[i]);
		ASSERT_EQ(pn->hap_id, expected_hap_ids[i]);
		destroy_pansn(&pn);
	}
}

TEST(ParsePanSNRefId, Invalid)
{
	const idx_t N = 6;
	// input
	const char *names[] = {
		"chm13__LPA__tig00000001", "HG002__LPA__tig00000001",
		"HG002__LPA__tig00000005", "#1#contig",
		"invalid_name_no_delim",   "too#many#delims#here",
		"hap#not#number"};

	// expected output: all should fail to parse
	for (idx_t i = 0; i < N; i++) {
		char *name = strdup(names[i]);
		pansn *pn = try_extract_pansn_from_str(name, HASH_CHAR);
		ASSERT_EQ(pn, nullptr);
	}
}

TEST(AllocRefIdPLine, ValidPanSN)
{
	const char *name = "chm13#0#Chr1";
	const char *tokens[1] = {name};
	struct ref_id *r_id = alloc_ref_id(tokens, P_LINE_ID_TOKEN_COUNT);

	ASSERT_NE(r_id, nullptr);
	ASSERT_EQ(r_id->type, REF_ID_PANSN);
	ASSERT_STREQ(r_id->value.id_value->sample_name, "chm13");
	ASSERT_STREQ(r_id->value.id_value->contig_name, "Chr1");
	ASSERT_EQ(r_id->value.id_value->hap_id, 0);
	ASSERT_STREQ(r_id->tag, tokens[0]);

	destroy_ref_id(&r_id);
	ASSERT_EQ(r_id, nullptr);
}

TEST(AllocRefIdPLine, ValidRaw)
{
	const char *name = "chm13__LPA__tig00000001";
	const char *tokens[1] = {name};
	struct ref_id *r_id = alloc_ref_id(tokens, P_LINE_ID_TOKEN_COUNT);

	ASSERT_NE(r_id, nullptr);
	ASSERT_EQ(r_id->type, REF_ID_RAW);
	ASSERT_STREQ(r_id->value.raw, name);
	ASSERT_STREQ(r_id->tag, name);

	destroy_ref_id(&r_id);
	ASSERT_EQ(r_id, nullptr);
}

TEST(StepCount, W_LINE)
{
	const char *w_line_data_strs[] = {
		">1>4>5>6>7>9", ">1>2>4>5>6>8>9", ">1>14>235>634>7>9", "",
		">343",
	};

	const idx_t expected_step_counts[] = {6, 7, 6, 0, 1};

	const idx_t N = 5;
	idx_t res;
	for (idx_t i = 0; i < N; i++) {
		res = count_steps(W_LINE, w_line_data_strs[i]);
		ASSERT_EQ(res, expected_step_counts[i]);
	}
}

TEST(StepCount, P_LINE)
{
	const char *p_line_data_strs[] = {
		"3+,4+,6+,8+,9+",
		"612+,614+,615+,617+,619+,621+,623+,624+,626+,628+,630+,631+,"
		"633+,634+",
		"",
		"123+",
		"11+,13+,15+,16+,18+,20+,21+,22+,23+,25+,15+,16+,26+,27+,15+,"
		"16+,28+",
	};

	const idx_t expected_step_counts[] = {5, 14, 0, 1, 17};

	const idx_t N = 5;
	idx_t res;
	for (idx_t i = 0; i < N; i++) {
		res = count_steps(P_LINE, p_line_data_strs[i]);
		ASSERT_EQ(res, expected_step_counts[i]);
	}
}

TEST(ParseDataLine, P_LINE)
{
	const char *p_line_data_strs[] = {
		"3+,4+,6+,8+,9+",
		"612+,614+,615+,617+,619+,621+,623+,624+,626+,628+,630+,631+,"
		"633+,634+",
		"",
		"123+",
		"11+,13+,15+,16+,18+,20+,21+,22+,23+,25+,15+,16+,26+,27+,15+,"
		"16+,28+",
	};

	const idx_t expected_step_counts[] = {5, 14, 0, 1, 17};

	// #define MAX_COLUMNS 17 // Maximum row length
	constexpr idx_t MAX_COLUMNS = 17; // Maximum row length
	const std::array<std::array<id_t, MAX_COLUMNS>, 5> p_v_ids = {
		pad_to_max<id_t, MAX_COLUMNS>({3, 4, 6, 8, 9, 9}),
		pad_to_max<id_t, MAX_COLUMNS>({612, 614, 615, 617, 619, 621,
					       623, 624, 626, 628, 630, 631,
					       633, 634}),
		pad_to_max<id_t, MAX_COLUMNS>({}),
		pad_to_max<id_t, MAX_COLUMNS>({123}),
		pad_to_max<id_t, MAX_COLUMNS>({11, 13, 15, 16, 18, 20, 21, 22,
					       23, 25, 15, 16, 26, 27, 15, 16,
					       28}),
	};

	const idx_t N = 5;
	idx_t step_count;
	for (idx_t i = 0; i < N; i++) {
		const char *data_str = p_line_data_strs[i];
		step_count = count_steps(P_LINE, data_str);
		ASSERT_EQ(step_count, expected_step_counts[i]);
		struct ref_walk *rw = alloc_ref_walk(step_count);
		ASSERT_NE(rw, nullptr);

		status_t res = parse_data_line_p(data_str, &rw);
		ASSERT_EQ(SUCCESS, res);

		for (idx_t j = 0; j < rw->step_count; j++) {
			ASSERT_EQ(p_v_ids[i][j], rw->v_ids[j]);
		}
		destroy_ref_walk(&rw);
		ASSERT_EQ(rw, nullptr);
	}
}

TEST(ParseDataLine, W_LINE)
{
	const char *w_line_data_strs[] = {
		">1>4>5>6>7>9", ">1>2>4>5>6>8>9", ">1>14>235>634>7>9", "",
		">343",
	};

	const idx_t expected_step_counts[] = {6, 7, 6, 0, 1};

	constexpr idx_t MAX_COLUMNS = 7; // Maximum row length
	const std::array<std::array<id_t, MAX_COLUMNS>, 5> w_v_ids = {
		pad_to_max<id_t, MAX_COLUMNS>({1, 4, 5, 6, 7, 9}),
		pad_to_max<id_t, MAX_COLUMNS>({1, 2, 4, 5, 6, 8, 9}),
		pad_to_max<id_t, MAX_COLUMNS>({1, 14, 235, 634, 7, 9}),
		pad_to_max<id_t, MAX_COLUMNS>({}),
		pad_to_max<id_t, MAX_COLUMNS>({343}),
	};

	const idx_t N = 5;
	idx_t computed_step_count, expected_step_count;
	for (idx_t i = 0; i < N; i++) {
		const char *data_str = w_line_data_strs[i];
		computed_step_count = count_steps(W_LINE, data_str);
		expected_step_count = expected_step_counts[i];

		ASSERT_EQ(computed_step_count, expected_step_count);
		struct ref_walk *rw = alloc_ref_walk(expected_step_count);
		ASSERT_NE(rw, nullptr);

		status_t res = parse_data_line_w(data_str, &rw);
		ASSERT_EQ(SUCCESS, res);

		// const id_t *v_id_ptr = &w_v_ids[i][0];
		for (idx_t j = 0; j < rw->step_count; j++) {
			ASSERT_EQ(w_v_ids[i][j], rw->v_ids[j]);
		}
	}
}

TEST(AllocRef, Valid)
{
	const idx_t STEP_COUNT = 5;
	struct ref_walk *rw = alloc_ref_walk(STEP_COUNT);

	const char *tokens[1] = {"sample#1#contig"};
	struct ref_id *r_id = alloc_ref_id(tokens, P_LINE_ID_TOKEN_COUNT);

	struct ref *r = alloc_ref(P_LINE, &rw, &r_id);
	ASSERT_NE(r, nullptr);
	ASSERT_EQ(r->line_prefix, P_LINE);
	ASSERT_EQ(rw->step_count, STEP_COUNT);
	ASSERT_EQ(r->id->type, REF_ID_PANSN);
	ASSERT_STREQ(r->id->value.id_value->sample_name, "sample");
	ASSERT_STREQ(r->id->value.id_value->contig_name, "contig");
	ASSERT_EQ(r->id->value.id_value->hap_id, 1);
	ASSERT_STREQ(r->id->tag, "sample#1#contig");

	destroy_ref(&r);
	ASSERT_EQ(r, nullptr);
}

TEST(ParseWLines, Valid)
{
	const char *w_lines[] = {
		"W\tshort\t1\tchr\t0\t8\t>1>4>5<6>7>9",
		"W\talt\t0\tchr\t0\t9\t>1>2>4>5>6>8>9",
		"W\tshort\t2\tchr\t0\t8\t>1>2<5>6>7<9",
	};

	const id_t hap_ids[] = {1, 0, 2};
	const char *sample_names[] = {"short", "alt", "short"};
	const char *expected_tags[] = {"short#1#chr", "alt#0#chr",
				       "short#2#chr"};

	const id_t expected_v_ids[][7] = {
		{1, 4, 5, 6, 7, 9, NULL_ID},
		{1, 2, 4, 5, 6, 8, 9},
		{1, 2, 5, 6, 7, 9, NULL_ID},
	};

	// clang-format off
	const enum strand exp_strands[][7] = {
		{STRAND_FWD, STRAND_FWD, STRAND_FWD, STRAND_REV, STRAND_FWD, STRAND_FWD, STRAND_FWD},
		{STRAND_FWD, STRAND_FWD, STRAND_FWD, STRAND_FWD, STRAND_FWD, STRAND_FWD, STRAND_FWD},
		{STRAND_FWD, STRAND_FWD, STRAND_REV, STRAND_FWD, STRAND_FWD, STRAND_REV, STRAND_FWD},
	};
	// clang-format on

	const char *expected_contig_name = "chr";

	const idx_t N = 3;
	for (idx_t i = 0; i < N; i++) {
		char *w_line = strdup(w_lines[i]);
		u32 line_len = (u32)strlen(w_line);

		if (w_line == NULL) {
			fprintf(stderr, "Memory allocation failed\n");
			exit(EXIT_FAILURE);
		}

		ref *r = parse_ref_line(W_LINE, w_line, line_len);
		ASSERT_NE(r, nullptr);

		// ASSERT_EQ(get_line_prefix(r), W_LINE);
		ASSERT_EQ(get_ref_id_type(r), REF_ID_PANSN);
		ASSERT_STREQ(get_tag(r), expected_tags[i]);
		ASSERT_EQ(get_hap_id(r), hap_ids[i]);
		ASSERT_STREQ(get_contig_name(r), expected_contig_name);
		ASSERT_STREQ(get_sample_name(r), sample_names[i]);

		for (idx_t j = 0; j < get_step_count(r); j++) {
			ASSERT_EQ(get_walk_v_ids(r)[j], expected_v_ids[i][j]);
			ASSERT_EQ(get_walk_strands(r)[j], exp_strands[i][j]);
		}
		destroy_ref(&r);
		ASSERT_EQ(r, nullptr);
	}
}

TEST(ParsePLines, Valid)
{

	const char *p_lines[] = {
		"P\tchm13__LPA__tig00000001\t3+,4+,6-,8+,9+",
		"P\tchm13#0#tig00000001\t3+,4+,6-,8+,9+",
	};

	const char *expected_tags[] = {"chm13__LPA__tig00000001",
				       "chm13#0#tig00000001"};

	const id_t expected_v_ids[][7] = {
		{3, 4, 6, 8, 9, 9},
		{3, 4, 6, 8, 9, 9},
	};

	const enum ref_id_type expected_ref_name_types[] = {REF_ID_RAW,
							    REF_ID_PANSN};

	const char *expected_sample_names[] = {"chm13__LPA__tig00000001",
					       "chm13"};

	const enum strand exp_strands[][7] = {
		{STRAND_FWD, STRAND_FWD, STRAND_REV, STRAND_FWD, STRAND_FWD},
		{STRAND_FWD, STRAND_FWD, STRAND_REV, STRAND_FWD, STRAND_FWD},
	};

	const idx_t N = 2;
	for (idx_t i = 0; i < N; i++) {
		char *p_line = strdup(p_lines[i]);
		if (p_line == NULL) {
			fprintf(stderr, "Memory allocation failed\n");
			exit(EXIT_FAILURE);
		}
		u32 line_len = (u32)strlen(p_line);

		ref *r = parse_ref_line(P_LINE, p_line, line_len);
		ASSERT_NE(r, nullptr);

		ASSERT_EQ(get_line_prefix(r), P_LINE);
		ASSERT_EQ(get_ref_id_type(r), expected_ref_name_types[i]);
		ASSERT_STREQ(get_sample_name(r), expected_sample_names[i]);

		if (get_ref_id_type(r) == REF_ID_PANSN) {
			// In this case we cannot parse hap_id and contig_name
			// from the name
			ASSERT_EQ(get_hap_id(r), 0);
			ASSERT_STREQ(get_contig_name(r), "tig00000001");
		}

		ASSERT_STREQ(get_tag(r), expected_tags[i]);

		for (idx_t j = 0; j < get_step_count(r); j++) {
			ASSERT_EQ(get_walk_v_ids(r)[j], expected_v_ids[i][j]);
			ASSERT_EQ(get_walk_strands(r)[j], exp_strands[i][j]);
		}
		destroy_ref(&r);
		ASSERT_EQ(r, nullptr);
	}
}

TEST(Refs, AllocatesAndFrees)
{
	idx_t ref_count = 1000;
	idx_t step_count = 5;

	struct ref **refs; // the reference sequences;
	refs = (struct ref **)malloc(sizeof(struct ref *) * ref_count);

	ASSERT_NE(refs, nullptr);
	for (idx_t i = 0; i < ref_count; i++) {
		struct ref_walk *rw = alloc_ref_walk(step_count);
		ASSERT_NE(rw, nullptr);
	}

	for (idx_t i = 0; i < ref_count; i++) {
		destroy_ref(&refs[i]);
		ASSERT_EQ(refs[i], nullptr);
	}
}
