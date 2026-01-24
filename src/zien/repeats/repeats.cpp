#include "zien/repeats/repeats.hpp" // for foo

#include <liteseq/refs.h> // for ref_walk, ref

#include "ita/align/align.hpp"	   // for align, aln_level_e
#include "ita/genomics/allele.hpp" // for ia::at_itn

#include "povu/common/core.hpp" // for pt
// #include "povu/common/utils.hpp"     // for pu::concat_with
#include "povu/graph/bidirected.hpp" // for VG
#include "povu/graph/types.hpp"

namespace zien::tui::repeats
{
namespace lq = liteseq;

std::map<pt::id_t, std::vector<pt::slice>>
unroll_tangle(const bd::VG &g, const pt::op_t<ptg::id_or_t> &ef)
{
	std::map<pt::id_t, std::vector<pt::slice>> unrolled_steps;

	auto [u, _] = ef.first;
	auto [v, __] = ef.second;

	for (pt::u32 h_idx{}; h_idx < g.get_hap_count(); h_idx++) {
		const std::vector<pt::u32> &positions_u =
			g.get_vertex_ref_idxs(g.v_id_to_idx(u), h_idx);

		const std::vector<pt::u32> &positions_v =
			g.get_vertex_ref_idxs(g.v_id_to_idx(v), h_idx);

		pt::u32 N = std::min(positions_u.size(), positions_v.size());
		for (pt::u32 i{}; i < N; i++) {
			pt::u32 pos_u = positions_u[i];
			pt::u32 pos_v = positions_v[i];

			if (pos_v < pos_u) // swap
				std::swap(pos_u, pos_v);

			pt::slice s{pos_u, pos_v - pos_u + 1};
			unrolled_steps[h_idx].emplace_back(s);
		}
	}

	return unrolled_steps;
}

ptg::or_e lq_strand_to_or(lq::strand s)
{
	return (s == lq::strand::STRAND_FWD) ? ptg::or_e::forward
					     : ptg::or_e::reverse;
}

ptg::walk_t slice_to_walk(const bd::VG &g, pt::u32 h_idx, const pt::slice &s)
{
	ptg::walk_t w;

	const lq::ref_walk *rw = g.get_ref_vec(h_idx)->walk;
	for (pt::u32 i = s.start(); i < s.start() + s.len(); i++) {
		pt::id_t v_id = rw->v_ids[i];
		ptg::or_e o = lq_strand_to_or(rw->strands[i]);
		w.push_back({v_id, o});
	}

	return w;
}

struct bar {
	pt::u32 h_idx;
	std::vector<pt::slice> itn;

	// when h_idx == ref hap then aln is empty or all matches
	std::string aln; // alignment string or edit transcript

	// methods
	[[nodiscard]]
	std::string baz(const bd::VG &g, pt::u32 h_idx,
			const std::vector<pt::slice> &curr_itn,
			pt::u32 slice_idx) const
	{
		if (slice_idx >= curr_itn.size())
			return "";

		ptg::walk_t w = slice_to_walk(g, h_idx, curr_itn.at(slice_idx));
		return ptg::to_string(w);
	}

	void to_tui_aln(const bd::VG &g, pt::u32 ref_h_idx,
			const std::vector<pt::slice> &ref_itn, pt::u32 pos,
			std::string &ref_at_str, std::string &alt_at_str,
			pt::slice &ref_sl, pt::slice &alt_sl) const
	{
		pt::u32 i{}; // aln index
		pt::u32 j{}; // ref itn index
		pt::u32 k{}; // alt itn index

		pt::u32 N = std::max(this->itn.size(), ref_itn.size());

		// rs ref string, as alt string
		auto get_rs = [&]() -> std::string
		{
			return this->baz(g, ref_h_idx, ref_itn, j++);
		};

		auto get_as = [&]() -> std::string
		{
			return this->baz(g, this->h_idx, this->itn, k++);
		};

		// Replace the character at the middle index with '-'
		auto dashen = [](std::string &s)
		{
			size_t middle = s.size() / 2;
			s[middle] = '-';
		};

		auto pos_covered = [&](pt::u32 j) -> bool
		{
			// std::cerr << "Checking pos_covered at j=" << j
			//	  << " ref itn size " << ref_itn.size() << "\n";

			if (j >= ref_itn.size())
				return false;

			pt::slice ref_sl = ref_itn.at(j);
			pt::u32 ref_start = ref_sl.start();
			// std::cerr << "Ref SL: " << ref_sl.start() << " "
			//	  << ref_sl.len() << "\n";
			pt::u32 ref_end = ref_start + (ref_sl.len() - 1);
			// std::cerr << "ref end " << ref_end;
			const lq::ref_walk *rw = g.get_ref_vec(ref_h_idx)->walk;

			pt::u32 s = rw->loci[ref_start];
			pt::u32 e = rw->loci[ref_end];
			pt::id_t v_id = rw->v_ids[ref_end];
			e += g.get_vertex_by_id(v_id).get_length();

			return (pos >= s && pos < e);
		};

		std::string ref_at_step;
		std::string alt_at_step;

		for (; i < N; i++) {

			ref_at_step.clear();
			alt_at_step.clear();
			char c = this->aln[i];

			if (c == 'I') {
				ref_at_step = get_rs();
			}
			else if (c == 'D') {
				j++;
				alt_at_step = get_as();
			}
			else { // M or X
				ref_at_step = get_rs();
				alt_at_step = get_as();
			}

			if (ref_at_step.empty()) {
				std::string s(alt_at_step.size(), ' ');
				dashen(s);
				ref_at_step = s;
			}
			else if (alt_at_step.empty()) {
				std::string s(ref_at_step.size(), ' ');
				dashen(s);
				alt_at_step = s;
			}
			else if (ref_at_step.length() < alt_at_step.length()) {
				// pad ref
				pt::u32 diff = alt_at_step.length() -
					       ref_at_step.length();
				ref_at_step += std::string(diff, ' ');
			}
			else if (alt_at_step.length() < ref_at_step.length()) {
				// pad alt
				pt::u32 diff = ref_at_step.length() -
					       alt_at_step.length();
				alt_at_step += std::string(diff, ' ');
			}

			if (pos_covered(j - 1)) { // highlight the position
				ref_sl = pt::slice(ref_at_str.length(),
						   ref_at_step.length());
				alt_sl = pt::slice(alt_at_str.length(),
						   alt_at_step.length());
			}

			ref_at_str += ref_at_step;
			alt_at_str += alt_at_step;

			if (i + 1 < N) {
				ref_at_str += std::string(2, ' ');
				alt_at_str += std::string(2, ' ');
			}
		}
	}
};

std::tuple<std::vector<pt::u32>, std::vector<pt::slice>,
	   std::vector<std::string>>
foo(const bd::VG &g, pt::u32 ref_h_idx, const pt::op_t<ptg::id_or_t> &ef,
    pt::u32 pos)
{
	std::vector<std::string> lines;
	std::vector<bar> bars;

	std::map<pt::id_t, std::vector<pt::slice>> unrolled =
		unroll_tangle(g, ef);

	pt::u32 longest_loop{};
	for (auto &[_, slices] : unrolled)
		if (slices.size() > longest_loop)
			longest_loop = slices.size();

	// std::cerr << "A" << "\n";

	bars.emplace_back(bar{ref_h_idx, unrolled.at(ref_h_idx), ""});

	auto lvl = ita::align::aln_level_e::at;

	std::map<pt::u32, ia::at_itn> itns;

	std::map<pt::u32, std::string> hap_walk_to_alns;

	for (const auto &[h_idx, slices] : unrolled) {
		std::vector<ptg::walk_t> walks;
		for (auto s : slices)
			walks.emplace_back(slice_to_walk(g, h_idx, s));

		ia::at_itn at(std::move(walks));
		itns[h_idx] = at;
	}

	// std::cerr << "B" << "\n";

	for (pt::u32 i{}; i < g.get_hap_count(); i++) {
		if (i == ref_h_idx || !pv_cmp::contains(itns, i))
			continue;

		const ia::at_itn &ref_itn = itns.at(ref_h_idx);
		const ia::at_itn &alt_itn = itns.at(i);
		std::string et = ita::align::align(ref_itn, alt_itn, lvl);
		hap_walk_to_alns[i] = et;

		// std::cerr << i << " " << et << "\n";

		bars.emplace_back(bar{i, unrolled.at(i), et});
	}

	// std::cerr << "C" << "\n";

	std::vector<pt::u32> header_lens;
	std::vector<pt::slice> highlight_slices; // highlight slices

	lines.reserve(bars.size());
	std::string ref_at_str;
	std::string alt_at_str;

	pt::slice ref_sl(0, 0);
	pt::slice alt_sl(0, 0);

	for (pt::u32 i{}; i < bars.size(); i++) {
		const bar &b = bars.at(i);
		if (b.h_idx == ref_h_idx)
			continue;

		ref_at_str += g.get_tag(ref_h_idx); // append ref tag
		alt_at_str += g.get_tag(b.h_idx);   // append alt tag

		// std::cerr << "Bar h_idx: " << ref_at_str << " " << alt_at_str
		//	  << "\n";

		header_lens.emplace_back(ref_at_str.length());
		header_lens.emplace_back(alt_at_str.length());

		b.to_tui_aln(g, ref_h_idx, unrolled.at(ref_h_idx), pos,
			     ref_at_str, alt_at_str, ref_sl, alt_sl);

		highlight_slices.push_back(ref_sl);
		highlight_slices.push_back(alt_sl);

		lines.push_back(ref_at_str);
		lines.push_back(alt_at_str);

		ref_at_str.clear();
		alt_at_str.clear();
	}

	// std::cerr << "D" << "\n";

	// std::cerr << "lines\n";
	// pt::u32 k{};
	// for (auto &l : lines)
	//	std::cerr << k++ << " " << l << "\n";

	return {header_lens, highlight_slices, lines};
}
} // namespace zien::tui::repeats
