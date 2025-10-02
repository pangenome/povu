#include <gtest/gtest.h>
#include <vector>

#include "../../include/refs/refs.hpp"

std::pair<pt::idx_t, std::vector<std::string>> all_undefined_refs()
{
	std::vector<std::string> test_samples = {
		"chm13__LPA__tig00000001",   "HG002__LPA__tig00000001",
		"HG002__LPA__tig00000005",   "HG00733__LPA__tig00000001",
		"HG00733__LPA__tig00000008", "HG01358__LPA__tig00000002",
		"HG01358__LPA__tig00000010", "HG02572__LPA__tig00000005",
		"HG02572__LPA__tig00000001", "NA19239__LPA__tig00000002",
		"NA19239__LPA__tig00000006", "NA19240__LPA__tig00000001",
		"NA19240__LPA__tig00000012"};

	const pt::idx_t P_LINE_COUNT = test_samples.size();
	return {P_LINE_COUNT, test_samples};
}

std::pair<pt::idx_t, std::vector<std::string>> bad_refs1()
{
	std::vector<std::string> test_samples = {"",
						 " ",
						 "##",
						 "###",
						 "sample#1#contig",
						 "sample##contig",
						 "#1#contig",
						 "sample#1#",
						 "sample#-1#contig",
						 "sample#abc#contig",
						 "sample#1#contig#extra"};

	const pt::idx_t P_LINE_COUNT = test_samples.size();
	return {P_LINE_COUNT, test_samples};
}

std::pair<pt::idx_t, std::vector<std::string>> bad_refs2()
{
	std::vector<std::string> test_samples = {"",
						 " ",
						 "##",
						 "###",
						 "sample#1#contig",
						 "sample##contig",
						 "#1#contig",
						 "sample#1#",
						 "sample#-1#contig",
						 "sample#abc#contig",
						 "sample#1#contig#extra"};

	const pt::idx_t P_LINE_COUNT = test_samples.size();
	return {P_LINE_COUNT, test_samples};
}

TEST(FormatTest, UndefinedRefCount)
{
	auto [P_LINE_COUNT, test_samples] = all_undefined_refs();

	// create refs
	pr::Refs refs(P_LINE_COUNT);

	// add each ref to refs
	for (pt::idx_t i{}; i < P_LINE_COUNT; ++i) {
		const std::string &s = test_samples[i];
		refs.add_ref(s, '#');
	}

	EXPECT_EQ(refs.ref_count(), P_LINE_COUNT);
}

TEST(FormatTest, UndefinedSampleCount)
{

	auto [P_LINE_COUNT, test_samples] = all_undefined_refs();

	// create refs
	pr::Refs refs(P_LINE_COUNT);

	// add each ref to refs
	for (pt::idx_t i{}; i < P_LINE_COUNT; ++i) {
		const std::string &s = test_samples[i];
		refs.add_ref(s, '#');
	}

	for (pt::idx_t i{}; i < P_LINE_COUNT; ++i) {
		const std::string &sample_name = test_samples[i];
		std::set<pt::id_t> ref_ids =
			refs.get_refs_in_sample(sample_name);
		EXPECT_EQ(ref_ids.size(), 1);
	}
}
