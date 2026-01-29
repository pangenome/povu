#include <algorithm>
#include <string>
#include <vector>

#include <liteseq/refs.h> // for ref_walk, ref

#include "povu/common/constants.hpp"
#include "povu/common/core.hpp"
#include "povu/graph/bidirected.hpp"	  // for VG
#include "zien/common/common.hpp"	  // for to_char
#include "zien/components/components.hpp" // for display_lines

namespace zien::components::paths
{
using cell = std::vector<std::string>;
using hap_row = std::vector<cell>;
using matrix = std::vector<hap_row>;

std::vector<std::string> comp_order(const bd::VG &g, pt::u32 start, pt::u32 end)
{
	std::vector<std::string> h;
	h.reserve(end - start);
	for (pt::u32 i{start}; i < end; i++)
		h.emplace_back(std::to_string(g.v_idx_to_id(i)));

	return h;
}

void update_display_lines(const bd::VG &g, const matrix &view_matrix,
			  const std::vector<pt::u32> &col_width,
			  const std::vector<pt::u32> &hap_row_count,
			  display_lines &pd)
{
	auto cell_to_str = [](const cell &cell_) -> std::string
	{
		std::string s;
		for (const auto &part : cell_)
			s += part;

		return s;
	};

	pt::u32 hap_idx{};
	pt::u32 hap_row_idx{};
	std::string line;

	pt::u32 matrix_row_idx{};
	pt::u32 matrix_row_count = view_matrix.size();

	while (matrix_row_idx < matrix_row_count) {
		const hap_row &hap_row = view_matrix[matrix_row_idx];

		if (hap_row_idx == 0) {
			std::string tag = g.get_tag(hap_idx);
			line += tag;
			pd.meta[matrix_row_idx].ref_name_pos = tag.length();
		}
		else {
			pd.meta[matrix_row_idx].ref_name_pos = 0;
		}

		// c_idx = cell index
		for (pt::u32 c_idx{}; c_idx < hap_row.size(); c_idx++) {
			const cell &cell_ = hap_row[c_idx];
			std::string k = cell_to_str(cell_);
			pt::u32 w = col_width[c_idx];

			if (k.length() < w)
				line += k + std::string(w - k.length(), ' ');
			else
				line += k;
		}

		pd.lines.push_back(line);
		line.clear();

		if (hap_row_idx + 1 >= hap_row_count[hap_idx]) {
			hap_idx++;
			hap_row_idx = 0;
			pd.group_lines.insert(matrix_row_idx + 1);
		}
		else {
			hap_row_idx++;
		}

		matrix_row_idx++;
	}
}

matrix comp_hap_rows(const bd::VG &g, const std::vector<std::string> &order,
		     pt::u32 hap_idx, pt::u32 start, pt::u32 end,
		     const liteseq::ref_walk *rw,
		     std::vector<pt::u32> &col_width,
		     std::vector<pt::u32> &hap_row_count)
{
	auto count_rows = [](const std::vector<pt::idx_t> &positions) -> pt::u32
	{
		pt::u32 N = positions.size();
		if (N < 2)
			return positions.size();

		pt::u32 row_count{};
		for (pt::u32 j{}; j < N - 1; j++)
			if (positions[j + 1] > positions[j])
				row_count++;

		return row_count;
	};

	auto baz = [](const liteseq::ref_walk *rw,
		      const std::vector<pt::idx_t> &positions, pt::u32 pos_idx,
		      pt::u32 order_step_idx,
		      const std::vector<std::string> &order, pt::u32 &prev_pos,
		      pt::u32 &w, cell &cell_)
	{
		pt::u32 pos = positions[pos_idx];
		liteseq::strand s = rw->strands[pos];
		char c = zien::common::to_char(s);
		std::string k = c + order[order_step_idx];
		w += k.length();

		cell_.emplace_back(k);
		prev_pos = pos;
	};

	auto check_prev = [](const std::vector<pt::idx_t> &positions,
			     pt::u32 prev_pos, bool row_has_data,
			     pt::u32 row_idx) -> std::vector<pt::u32>
	{
		std::vector<pt::u32> res{};

		pt::u32 N = positions.size();
		if (prev_pos == pc::INVALID_IDX) {
			res.push_back(0);
		}
		else {
			for (pt::u32 pos_idx{}; pos_idx < N; pos_idx++) {
				if ((positions[pos_idx] - 1) == prev_pos) {
					res.push_back(pos_idx);
					break;
				}
				else if (!row_has_data && pos_idx == row_idx) {
					res.push_back(pos_idx);
					break;
				}
			}

			if (res.empty())
				return {};
		}

		for (pt::u32 i{res.front()}; i < N - 1; i++)
			if (positions[i] + 1 == (positions[i + 1]))
				res.push_back(i + 1);

		return res;
	};

	matrix m; // hap rows

	pt::u32 curr_row{};
	pt::u32 max_row_count{};

	pt::u32 prev_pos = pc::INVALID_IDX;
	do {
		bool row_has_data{false};
		hap_row row; // single hap row
		for (pt::u32 v_idx{start}; v_idx < end; v_idx++) {
			const std::vector<pt::idx_t> &positions =
				g.get_vertex_ref_idxs(v_idx, hap_idx);

			if (positions.empty()) {
				row.emplace_back(); // empty cell
				continue;
			}

			// do this only once to get the max row count
			if (curr_row == 0) {
				pt::u32 row_count = count_rows(positions);
				if (row_count > max_row_count)
					max_row_count = row_count;
			}

			// the col for that specific node
			// i - start to get the index in order
			pt::u32 order_step_idx = v_idx - start;

			std::vector<pt::u32> curr_pos_idxs = check_prev(
				positions, prev_pos, row_has_data, curr_row);

			if (!curr_pos_idxs.empty())
				row_has_data = true;

			cell cell_;
			pt::u32 cell_width{};

			for (pt::u32 curr_pos_idx : curr_pos_idxs)
				baz(rw, positions, curr_pos_idx, order_step_idx,
				    order, prev_pos, cell_width, cell_);

			// if empty it adds an empty cell
			row.emplace_back(cell_);

			if (cell_width > col_width[order_step_idx])
				col_width[order_step_idx] = cell_width;
		}
		m.emplace_back(row);
	} while (++curr_row < max_row_count);

	hap_row_count.emplace_back(max_row_count);

	return m;
}

pt::op_t<pt::u32> comp_path_view_range(const bd::VG &g, const ui_state &state)
{
	const pt::u32 min_idx{};	     // minimum valid vertex index
	const pt::u32 half_window_size{250}; // halfway point vertex idx

	pt::u32 start_idx = (state.paths_view_mid > half_window_size)
				    ? state.paths_view_mid - half_window_size
				    : min_idx;

	pt::u32 end_idx = std::min(state.paths_view_mid + half_window_size,
				   g.vtx_count());

	return {start_idx, end_idx};
}

void update_paths(const bd::VG &g, ui_state &state, display_lines &pd)
{
	auto [start, end] = comp_path_view_range(g, state);
	state.paths_view_range = {start, end};
	std::vector<std::string> order = comp_order(g, start, end);

	pt::u32 col_count = order.size();
	std::vector<pt::u32> col_width(col_count, 0); // default col width is 0

	std::vector<pt::u32> hap_row_count;
	hap_row_count.reserve(g.get_hap_count());

	matrix view_matrix;

	// populate view matrix
	for (pt::u32 hap_idx{}; hap_idx < g.get_hap_count(); hap_idx++) {
		// row header
		std::string tag = g.get_tag(hap_idx);
		if (tag.length() > pd.lh)
			pd.lh = tag.length();

		// row content
		const liteseq::ref_walk *rw = g.get_ref_vec(hap_idx)->walk;
		pt::u32 N = std::min(rw->step_count, end);
		matrix hap_rows = comp_hap_rows(g, order, hap_idx, start, N, rw,
						col_width, hap_row_count);

		for (auto &&r : hap_rows)
			view_matrix.emplace_back(r);
	}

	// update display lines from view matrix
	update_display_lines(g, view_matrix, col_width, hap_row_count, pd);
}
} // namespace zien::components::paths
