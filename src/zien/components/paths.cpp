#include <algorithm>
#include <queue>
#include <string>
#include <vector>

#include <liteseq/refs.h> // for ref_walk, ref

#include "meza/matrix/matrix.hpp"	  // for matrix2d
#include "povu/common/constants.hpp"	  // for pv_cmp
#include "povu/common/core.hpp"		  // for pt, pc
#include "povu/graph/bidirected.hpp"	  // for VG
#include "zien/common/common.hpp"	  // for to_char
#include "zien/components/components.hpp" // for display_lines

namespace zien::components::paths
{
using cell = std::vector<std::string>;
using hap_row = std::vector<cell>;
using matrix = std::vector<hap_row>;

struct range_meta {
	pt::u32 min_pos = pc::MAX_IDX;
	pt::u32 row_count = 1;
	pt::u32 wrap_idx = 0;

	range_meta() = default; // default constructor
};

struct disp_matrix : public meza::matrix::path_matrix {
	std::vector<range_meta> meta;
	std::vector<pt::u32> col_width; // default col width is 0

	// ------------
	// constructors
	// ------------
	disp_matrix() = delete;

	disp_matrix(pt::u32 I, pt::u32 J, std::vector<range_meta> &&m)
	    : meza::matrix::path_matrix(I, J), meta(std::move(m)),
	      col_width(J, 0)
	{}
};

std::vector<std::string> comp_order(const bd::VG &g, pt::u32 start, pt::u32 end)
{
	std::vector<std::string> h;
	h.reserve(end - start);
	for (pt::u32 i{start}; i < end; i++)
		h.emplace_back(std::to_string(g.v_idx_to_id(i)));

	return h;
}

/**
 * Update the display lines based on the provided display matrix and metadata.
 *
 * @param g The bidirected graph containing the haplotypes and vertices.
 * @param meta A vector of range_meta structures containing metadata for each
 * haplotype.
 * @param dm The display matrix containing the cell data to be displayed.
 * @param pd A reference to the display_lines structure to be updated with the
 * formatted lines.
 */
void update_display_lines(const bd::VG &g, const disp_matrix &dm,
			  display_lines &pd)
{
	auto cell_to_str = [](const cell &cell_) -> std::string
	{
		std::string s;
		for (const auto &part : cell_)
			s += part;

		return s;
	};

	// const matrix &view_matrix = dm.data;
	const std::vector<pt::u32> &col_width = dm.col_width;
	const std::vector<range_meta> &meta = dm.meta;

	pt::u32 hap_idx{};
	pt::u32 hap_row_idx{};
	std::string line;

	pt::u32 matrix_row_idx{};
	pt::u32 matrix_row_count = dm.base().rows();
	while (matrix_row_idx < matrix_row_count) {
		const range_meta &rm = meta[hap_idx];
		pt::u32 hap_row_count = rm.row_count;
		// const hap_row &hap_row = view_matrix[matrix_row_idx];

		if (hap_row_idx == 0) {
			std::string tag = g.get_tag(hap_idx);
			line += tag;
			pd.meta[matrix_row_idx].ref_name_pos = tag.length();
		}
		else {
			pd.meta[matrix_row_idx].ref_name_pos = 0;
		}

		// TODO: better handle blank rows
		// blank
		if (hap_row_count == 0) {
			pd.lines.push_back(line);
			line.clear();
			hap_idx++;
			hap_row_idx = 0;
			pd.group_lines.insert(matrix_row_idx + 1);
			matrix_row_idx++;
			continue;
		}

		// c_idx = cell index
		for (pt::u32 c_idx{}; c_idx < dm.base().cols(); c_idx++) {
			const cell &cell_ = dm.base().at(matrix_row_idx, c_idx);
			std::string k = cell_to_str(cell_);
			pt::u32 w = col_width[c_idx];

			if (k.length() < w)
				line += k + std::string(w - k.length(), ' ');
			else
				line += k;
		}

		pd.lines.push_back(line);
		line.clear();

		if (hap_row_idx + 1 >= hap_row_count) {
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

/**
 * Compute the wrap around index for a given haplotype in the specified vertex
 * range.
 *
 * @param g The bidirected graph containing the haplotypes and vertices.
 * @param hap_idx The index of the haplotype to compute the wrap around for.
 * @param start_v_idx The starting vertex index of the range (inclusive).
 * @param end_v_idx The ending vertex index of the range (exclusive).
 * @param rm A reference to a range_meta structure to store the computed wrap
 * index.
 */
void comp_wrap_around(const bd::VG &g, pt::u32 hap_idx, pt::u32 start_v_idx,
		      pt::u32 end_v_idx, range_meta &rm)
{
	for (pt::u32 v_idx{end_v_idx}; v_idx-- > start_v_idx;) {
		const std::vector<pt::idx_t> &positions =
			g.get_vertex_ref_idxs(v_idx, hap_idx);

		if (positions.empty())
			continue;

		pt::u32 min_pos =
			*std::min_element(positions.begin(), positions.end());

		if (min_pos > rm.wrap_idx) {
			rm.wrap_idx = min_pos;
			break;
		}
	}
}

void scan_range(const bd::VG &g, pt::u32 hap_idx, range_meta &rm,
		pt::u32 start_v_idx, pt::u32 end_v_idx)
{
	comp_wrap_around(g, hap_idx, start_v_idx, end_v_idx, rm);

	for (pt::u32 v_idx{start_v_idx}; v_idx < end_v_idx; v_idx++) {
		const std::vector<pt::idx_t> &positions =
			g.get_vertex_ref_idxs(v_idx, hap_idx);

		if (positions.empty())
			continue;

		// compute min pos
		for (pt::u32 pos : positions)
			if (pos < rm.min_pos)
				rm.min_pos = pos;

		// compute row count
		//
		// wrap around detected
		pt::u32 row_count_positions{1};
		for (pt::u32 i{}; i < positions.size() - 1; i++)
			if (positions[i] + 1 < positions[i + 1])
				row_count_positions++;

		// wrap around row count
		pt::u32 row_count_wrap{0};
		pt::u32 row_idx{};
		for (pt::u32 step_idx : positions) {
			if (step_idx > rm.wrap_idx) {
				row_idx = step_idx / rm.wrap_idx;

				if (row_count_wrap < row_idx)
					row_count_wrap = row_idx;
			}
		}

		// final row count for this haplotype in this range
		pt::u32 row_count =
			std::max(row_count_positions, row_count_wrap);
		if (rm.row_count < row_count)
			rm.row_count = row_count;
	}

	if (rm.row_count > 1)
		rm.row_count++; // for 0 based index

	return;
}

void comp_hap_rows(const bd::VG &g, pt::u32 hap_idx, pt::u32 start, pt::u32 end,
		   disp_matrix &dm, pt::u32 row_offset)
{
	auto fill_cell = [](const liteseq::ref_walk *rw,
			    const std::vector<pt::idx_t> &positions,
			    pt::u32 pos_idx, pt::u32 &prev_pos, pt::u32 &w,
			    disp_matrix &dm, pt::u32 row_idx,
			    pt::u32 col_idx) -> void
	{
		pt::u32 pos = positions[pos_idx];

		liteseq::strand s = rw->strands[pos];
		char c = zien::common::to_char(s);
		std::string k = c + std::to_string(rw->v_ids[pos]);

		// std::string pos_str = "[" + std::to_string(pos) + "]";
		// k += pos_str; // to indicate position

		w += k.length();

		cell &cell_ = dm.base_mut().at(row_idx, col_idx);
		cell_.emplace_back(k);
		dm.base_mut().mark_non_blank(row_idx, col_idx);
		prev_pos = pos;
	};

	// matrix &m = dm.data;
	std::vector<pt::u32> &col_width = dm.col_width;
	const liteseq::ref_walk *rw = g.get_ref_vec(hap_idx)->walk;

	std::map<pt::u32, std::queue<pt::u32>> buff;
	std::set<pt::u32> seen;
	const range_meta &rm = dm.meta[hap_idx];
	pt::u32 row_count = rm.row_count;
	// bool row_has_data{false};
	pt::u32 prev_pos = pc::INVALID_IDX;

	pt::u32 conceptual_row_idx{};
	pt::u32 row_idx{}; // actual row index in the matrix
	for (; conceptual_row_idx < row_count; row_idx++) {
		prev_pos = pc::INVALID_IDX;

		for (pt::u32 v_idx{start}; v_idx < end; v_idx++) {

			const std::vector<pt::idx_t> &positions =
				g.get_vertex_ref_idxs(v_idx, hap_idx);

			if (positions.empty())
				continue;

			// the col for that specific node
			// i - start to get the index in order
			pt::u32 order_step_idx = v_idx - start;
			pt::u32 col_idx = order_step_idx;

			std::set<pt::u32> curr_pos_idxs;

			if (pv_cmp::contains(buff, v_idx)) {
				auto &q = buff[v_idx];
				while (!q.empty()) {
					pt::u32 curr_pos_idx = q.front();
					q.pop();
					curr_pos_idxs.insert(curr_pos_idx);

					if (!q.empty() &&
					    positions[curr_pos_idx] + 1 <
						    positions[q.front()]) {
						break;
					}
				}

				if (q.empty())
					buff.erase(v_idx);
			}

			for (pt::u32 i{}; i < positions.size(); i++) {
				if (rm.wrap_idx == 0)
					continue;

				pt::u32 prev_step_idx{pc::INVALID_IDX};
				if (!curr_pos_idxs.empty()) {
					pt::u32 i = *curr_pos_idxs.rbegin();
					prev_step_idx = positions[i];
				};

				pt::u32 curr_step_idx = positions[i];

				pt::u32 expected_row_idx =
					positions[i] / rm.wrap_idx;

				if (expected_row_idx > 0 &&
				    positions[i] % rm.wrap_idx == 0)
					expected_row_idx--;

				if (expected_row_idx < conceptual_row_idx || //
				    pv_cmp::contains(seen, curr_step_idx))   //
					continue;

				if (expected_row_idx > conceptual_row_idx)
					break;

				// going forward
				// expected_row_idx == conceptual_row_idx

				// buffer
				if (prev_pos == pc::INVALID_IDX &&	// new
				    prev_step_idx != pc::INVALID_IDX && // x
				    prev_step_idx + 1 < curr_step_idx	// x
				) {
					buff[v_idx].push(i);
				}
				else if (prev_pos == pc::INVALID_IDX ||
					 prev_step_idx + 1 == curr_step_idx ||
					 (prev_pos != pc::INVALID_IDX &&
					  (prev_pos + 1 == curr_step_idx))) {
					curr_pos_idxs.insert(i);
				}

				seen.insert(curr_step_idx);
			}

			// extend from previous col
			// checks if a position is adjacent to prev_pos
			if (prev_pos != pc::INVALID_IDX) {
				for (pt::u32 i{}; i < positions.size(); i++) {
					if (prev_pos + 1 == positions[i] ||
					    prev_pos - 1 == positions[i]) {
						curr_pos_idxs.insert(i);
						break;
					}
				}
			}

			for (pt::u32 i{1}; i < positions.size(); i++) {
				if (positions[i - 1] + 1 < positions[i] &&
				    prev_pos == positions[i - 1]) {
					curr_pos_idxs.insert(i);
				}
			}

			// if (!row_has_data && !curr_pos_idxs.empty())
			//	row_has_data = true;

			pt::u32 cell_width{};

			for (pt::u32 curr_pos_idx : curr_pos_idxs)
				fill_cell(rw, positions, curr_pos_idx, prev_pos,
					  cell_width, dm, row_idx + row_offset,
					  col_idx);

			if (cell_width > col_width[order_step_idx])
				col_width[order_step_idx] = cell_width;
		}

		if (buff.empty())
			conceptual_row_idx++;
	}
}

pt::op_t<pt::u32> comp_path_view_range(const bd::VG &g, const ui_state &state)
{
	const pt::u32 min_idx{}; // minimum valid vertex index
	// halfway point vertex idx
	const pt::u32 half_window_size{state.half_window_size};

	pt::u32 start_idx = (state.paths_view_mid > half_window_size)
				    ? state.paths_view_mid - half_window_size
				    : min_idx;

	pt::u32 end_idx = std::min(state.paths_view_mid + half_window_size,
				   g.vtx_count());

	return {start_idx, end_idx};
}

/**
  pre comp metadata for each haplotype in the range
*/
std::pair<pt::u32, std::vector<range_meta>>
pre_comp_meta(const bd::VG &g, pt::u32 start, pt::u32 end)
{
	std::vector<range_meta> meta(g.get_hap_count());
	pt::u32 total_row_count{};
	for (pt::u32 hap_idx{}; hap_idx < g.get_hap_count(); hap_idx++) {
		// row content
		// idx in the vertex is the hap idx the value is the
		// range meta for that hap
		pt::u32 N = std::min(g.vtx_count(), end);
		range_meta &rm = meta[hap_idx];
		scan_range(g, hap_idx, rm, start, N);
		total_row_count += rm.row_count;
	}

	return {total_row_count, meta};
}

void update_paths(const bd::VG &g, ui_state &state, display_lines &pd)
{
	auto [start, end] = comp_path_view_range(g, state);
	state.paths_view_range = {start, end};
	std::vector<std::string> order = comp_order(g, start, end);

	pt::u32 col_count = order.size();

	std::vector<pt::u32> hap_row_count;
	hap_row_count.reserve(g.get_hap_count());

	auto [total_row_count, meta] = pre_comp_meta(g, start, end);
	disp_matrix dm(total_row_count, col_count, std::move(meta));

	pt::u32 row_offset{};
	// populate view matrix
	for (pt::u32 hap_idx{}; hap_idx < g.get_hap_count(); hap_idx++) {
		// row header
		std::string tag = g.get_tag(hap_idx);
		if (tag.length() > pd.lh)
			pd.lh = tag.length();

		// row content
		pt::u32 N = std::min(g.vtx_count(), end);
		comp_hap_rows(g, hap_idx, start, N, dm, row_offset);
		row_offset += dm.meta[hap_idx].row_count;
	}

	update_display_lines(g, dm, pd);
}
} // namespace zien::components::paths
