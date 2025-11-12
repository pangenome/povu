#include "povu/variation/rov.hpp"

// #include <cstddef>  // for size_t
// #include <cstdlib>  // for exit, EXIT_FAILURE
// #include <iterator> // for back_insert_iterator, bac...
// #include <optional> // for optional, operator==
// #include <string>   // for basic_string, string

#include <map>	   // for map
#include <set>	   // for set
#include <utility> // for move
#include <vector>  // for vector

// #include "fmt/core.h"		  // for format_to
// #include "indicators/setting.hpp" // for PostfixText
// #include "povu/common/app.hpp"	  // for config
// #include "povu/common/constants.hpp"
// #include "povu/common/log.hpp"
// #include "povu/common/progress.hpp" // for set_progress_bar_com...
// #include "povu/common/utils.hpp"
// #include "povu/genomics/graph.hpp" // for RoV, find_walks, pgt
// #include "povu/genomics/vcf.hpp"
// #include "povu/graph/pvst.hpp" // for Tree, VertexBase
// #include "povu/variation/color.hpp"
// #include "povu/common/utils.hpp" // for print_with_comma

#include "povu/common/compat.hpp" // for pv_cmp::contains
#include "povu/common/core.hpp"	  // for pt
#include "povu/graph/types.hpp"	  // for or_e, id_or_t, walk_t

namespace povu::var::rov
{
// namespace pgg = povu::genomics::graph;
// namespace pvst = povu::pvst;
// using namespace povu::progress;

// constexpr var_type_e ins = var_type_e::ins;
// constexpr var_type_e del = var_type_e::del;
// constexpr var_type_e sub = var_type_e::sub;

struct walk_matrix {
private:
	std::vector<pt::u8> data;
	// std::map<matrix_idx, pt::u32> idx;
	pt::u32 vtx_count;

	pt::u32 walk_count; // also row count
	pt::u32 col_count;

	std::map<pt::u32, pt::u32> col_to_vid;
	std::map<pt::u32, pt::u32> vid_to_col;

	walk_matrix(pt::u32 N, const std::set<pt::u32> &vtxs)
	    : data(N, 0), vtx_count(vtxs.size()), col_count(vtxs.size())
	{
		pt::u32 col_idx{};
		for (pt::u32 v_id : vtxs) {
			vid_to_col[v_id] = col_idx;
			col_to_vid[col_idx] = v_id;
			col_idx++;
		}
	}

	[[nodiscard]]
	bool can_extend(const std::vector<pt::u32> &from,
			const std::vector<pt::u32> &to) const
	{
		for (pt::u32 i{}; i < walk_count; i++)
			if (from[i] == 1 && to[i] == 0)
				return false;

		return true;
	}

	[[nodiscard]]
	bool try_ext_left(pt::u32 u) const
	{
		if (u == 0)
			return false;

		std::vector<pt::u32> u_walks = this->sweep_col(u);
		std::vector<pt::u32> u_walks_prev = this->sweep_col(u - 1);

		// std::cerr << "Trying extend left for col " << u
		//	  << " with col count " << col_count << "\n";

		return can_extend(u_walks, u_walks_prev);
	}

	[[nodiscard]]
	bool try_ext_right(pt::u32 v) const
	{
		if (v >= col_count - 1)
			return false;

		std::vector<pt::u32> v_walks = this->sweep_col(v);
		std::vector<pt::u32> v_walks_nxt = this->sweep_col(v + 1);

		return can_extend(v_walks, v_walks_nxt);
	}

public:
	[[nodiscard]]
	walk_matrix static create_blank(pt::u32 walk_count, pt::u32 vtx_count,
					const std::set<pt::u32> &vtxs)

	{
		walk_matrix wm{walk_count * vtx_count, vtxs};
		wm.walk_count = walk_count;
		return wm;
	}

	void fill(const std::set<pt::u32> &all_vertices,
		  const std::map<pt::u32, std::set<pt::u32>> &v_to_walks)
	{
		for (pt::u32 v_id : all_vertices)
			for (pt::u32 w_idx : v_to_walks.at(v_id))
				this->set_val(w_idx, v_id);
	}

	void print_matrix() const
	{
		for (pt::u32 w_idx{}; w_idx < walk_count; w_idx++) {
			for (pt::u32 col_idx{}; col_idx < col_count;
			     col_idx++) {
				pt::u32 matrix_idx =
					w_idx * vtx_count + col_idx;
				std::cerr << static_cast<pt::u32>(
						     data[matrix_idx])
					  << " ";
			}
			std::cerr << "\n";
		}
	}

	[[nodiscard]]
	pt::u32 get_col(pt::u32 v_id) const
	{
		return this->vid_to_col.at(v_id);
	}

	[[nodiscard]]
	pt::u32 matrix_idx(pt::u32 walk_idx, pt::u32 v_id) const
	{
		pt::u32 col_idx = this->get_col(v_id);
		return walk_idx * vtx_count + col_idx;
	}

	void set_val(pt::u32 walk_idx, pt::u32 v_id)
	{
		pt::u32 col_idx = vid_to_col[v_id];
		pt::u32 matrix_idx = walk_idx * vtx_count + col_idx;
		data[matrix_idx] = 1;
	}

	[[nodiscard]]
	std::vector<pt::u32> sweep_col(pt::u32 col_idx) const
	{
		std::vector<pt::u32> walk_idxs;
		for (pt::u32 w_idx{}; w_idx < walk_count; w_idx++) {
			pt::u32 matrix_idx = w_idx * vtx_count + col_idx;
			walk_idxs.push_back(data[matrix_idx]);
		}
		return walk_idxs;
	}

	[[nodiscard]]
	std::vector<pt::u32> sweep_row(pt::u32 row_idx) const
	{
		std::vector<pt::u32> row;
		for (pt::u32 col_idx{}; col_idx < col_count; col_idx++) {
			pt::u32 matrix_idx = row_idx * vtx_count + col_idx;
			row.push_back(data[matrix_idx]);
		}
		return row;
	}

	[[nodiscard]]
	pt::op_t<pt::u32> comp_u_v(pt::u32 a, pt::u32 b) const
	{
		pt::u32 a_col_idx = this->get_col(a);
		pt::u32 b_col_idx = this->get_col(b);

		pt::u32 u = std::min(a_col_idx, b_col_idx);
		pt::u32 v = std::max(a_col_idx, b_col_idx);

		return {u, v};
	}

	[[nodiscard]]
	bool is_nested(pt::u32 a, pt::u32 b) const
	{
		pt::u32 a_col_idx = this->get_col(a);
		pt::u32 b_col_idx = this->get_col(b);

		pt::u32 u = std::min(a_col_idx, b_col_idx);
		pt::u32 v = std::max(a_col_idx, b_col_idx);

		if (u == 0 && v >= col_count - 1)
			return false;

		return try_ext_left(u) || try_ext_right(v);
	}

	[[nodiscard]]
	std::set<pt::u32> shared_walks(pt::u32 u, pt::u32 v) const
	{
		auto find_col_walks = [&](const std::vector<pt::u32> &col)
			-> std::set<pt::u32>
		{
			std::set<pt::u32> walks;
			for (pt::u32 i{}; i < walk_count; i++)
				if (col[i] == 1)
					walks.insert(i);
			return walks;
		};

		std::set<pt::u32> shared;
		std::vector<pt::u32> u_col = sweep_col(u);
		std::set<pt::u32> u_walks = find_col_walks(u_col);

		std::vector<pt::u32> v_col = sweep_col(v);
		std::set<pt::u32> v_walks = find_col_walks(v_col);

		std::set_intersection(u_walks.begin(), u_walks.end(),
				      v_walks.begin(), v_walks.end(),
				      std::inserter(shared, shared.begin()));
		return shared;
	}

	// maybe check in the actual walks and not just in the walk matrix?
	[[nodiscard]]
	bool adj_cols(pt::u32 a, pt::u32 b)
	{
		pt::u32 a_col_idx = this->get_col(a);
		pt::u32 b_col_idx = this->get_col(b);

		pt::u32 u = std::min(a_col_idx, b_col_idx);
		pt::u32 v = std::max(a_col_idx, b_col_idx);

		std::set<pt::u32> sw = shared_walks(u, v);

		std::vector<pt::u32> walk_idxs;
		for (pt::u32 w_idx{}; w_idx < walk_count; w_idx++) {
			if (!pv_cmp::contains(sw, w_idx))
				continue;

			pt::u32 row = w_idx * vtx_count;
			for (pt::u32 col_idx = u + 1; col_idx < v; col_idx++) {
				pt::u32 matrix_idx = row + col_idx;
				if (data[matrix_idx] == 1)
					return false;
			}
		}

		return true;
	}

	// [[nodiscard]]
	// std::set<pt::u32> shared_walks_deep(pt::u32 a, pt::u32 b) const
	// {
	//	auto [u, v] = comp_u_v(a, b);
	//	std::set<pt::u32> all_vertices;

	//	return all_vertices;
	// }
};

struct walk_guide {

	std::map<pt::u32, std::set<pt::u32>> v_to_walks;

	[[nodiscard]]
	static walk_guide create(const std::vector<pgt::walk_t> &walks)
	{
		std::map<pt::u32, std::set<pt::u32>> v_to_walks;

		for (pt::u32 i{}; i < walks.size(); i++) {
			const pgt::walk_t &w = walks[i];
			for (pt::u32 s{}; s < w.size(); s++) {
				v_to_walks[w[s].v_id].insert(i);
			}
		}

		return walk_guide{std::move(v_to_walks)};
	}

	[[nodiscard]]
	std::set<pt::u32> get_all_vertices() const
	{
		std::set<pt::u32> all_vertices;
		for (const auto &[v_id, _] : v_to_walks)
			all_vertices.insert(v_id);

		return all_vertices;
	}

	[[nodiscard]]
	const std::map<pt::u32, std::set<pt::u32>> &get_data() const
	{
		return v_to_walks;
	}

	[[nodiscard]]
	pt::u32 vtx_count() const
	{
		return v_to_walks.size();
	}
};

std::set<pt::up_t<pt::u32>>
find_non_planar(const std::vector<pgt::walk_t> &walks)
{
	walk_guide wg = walk_guide::create(walks);
	std::set<pt::u32> all_vertices = wg.get_all_vertices();

	std::set<pt::u32> flanks;
	for (const auto &[v_id, walk_set] : wg.get_data())
		if (walk_set.size() > 1)
			flanks.insert(v_id);

	pt::u32 vtx_count = wg.vtx_count();

	// walk idx, v id
	// using matrix_idx = std::pair<pt::u32, pt::u32>;

	walk_matrix wm = walk_matrix::create_blank(walks.size(), vtx_count,
						   all_vertices);

	wm.fill(all_vertices, wg.get_data());

	std::set<pt::up_t<pt::u32>> all_bounds;
	for (pt::u32 u : flanks) {
		for (pt::u32 v : flanks) {
			pt::up_t<pt::u32> k{u, v};
			if (u == v || pv_cmp::contains(all_bounds, k))
				continue;

			if (!wm.is_nested(u, v) && !wm.adj_cols(u, v))
				all_bounds.insert({u, v});
		}
	}

	// std::cerr << "Walk matrix:\n";
	// wm.print_matrix();

	// std::cerr << "------\n";

	if (all_bounds.empty())
		std::exit(EXIT_FAILURE);

	// std::cerr << "Non-planar pairs in Rov:\n";
	// for (auto [u, v] : all_bounds)
	//	std::cerr << "(" << u << ", " << v << ")\n";

	// std::exit(EXIT_SUCCESS);

	return all_bounds;
}

} // namespace povu::var::rov
