#ifndef POVU_GENOMICS_GRAPH_HPP
#define POVU_GENOMICS_GRAPH_HPP

#include <algorithm>
#include <deque>
#include <set>
#include <stack>
#include <string>
#include <sys/types.h>
#include <utility>
#include <vector>

"#include "povu/common/compat.hpp"
"#include "povu/common/core.hpp"
"#include "povu/common/log.hpp"
"#include "povu/graph/bidirected.hpp"
"#include "povu/graph/pvst.hpp"
"#include "povu/graph/types.hpp"

namespace povu::genomics::graph
{
inline constexpr std::string_view MODULE = "povu::genomics::graph";

namespace pvst = povu::pvst;
namespace pgt = povu::types::graph;

// Maximum number of steps to take from flubble start to end
const pt::idx_t MAX_FLUBBLE_STEPS{20};

// direction for traversing a vertex in a bidirected graph
enum class dir_e {
	in,
	out
};
std::string_view to_str(dir_e d);
const dir_e IN = dir_e::in;
const dir_e OUT = dir_e::out;
auto format_as(dir_e d);

typedef pgt::id_or_t idx_or_t; // specifically for idx instead of id

/**
 * a collection of walks within a region of variation from start to end
 */
class RoV
{
	std::vector<pgt::walk_t> walks_;
	const pvst::VertexBase *pvst_vtx;

public:
	// print when RoV is moved
	RoV(RoV &&other) noexcept
	    : walks_(std::move(other.walks_)), pvst_vtx(other.pvst_vtx)
	{
		other.pvst_vtx = nullptr;
		// INFO("Moved RoV {}", this->as_str());
	}

	RoV(const RoV &other) = delete;		   // disable copy constructor
	RoV &operator=(const RoV &other) = delete; // disable copy assignment
	RoV &
	operator=(RoV &&other) noexcept = delete; // disable move assignment

	~RoV()
	{
		if (this->pvst_vtx != nullptr) {
			// INFO("Destroyed RoV {}", this->as_str());
		}
	}

	// --------------
	// constructor(s)
	// --------------
	RoV(const pvst::VertexBase *v) : walks_(), pvst_vtx(v)
	{}

	// ---------
	// getter(s)
	// ---------

	pt::idx_t walk_count() const
	{
		return this->walks_.size();
	}

	const pvst::VertexBase *get_pvst_vtx() const
	{
		return this->pvst_vtx;
	}

	const std::vector<pgt::walk_t> &get_walks() const
	{
		return this->walks_;
	}

	std::vector<pgt::walk_t> &get_walks_mut()
	{
		return this->walks_;
	}

	// ---------
	// setter(s)
	// ---------

	void set_walks(std::vector<pgt::walk_t> &&walks)
	{
		this->walks_ = std::move(walks);
	}

	// --------
	// other(s)
	// --------

	std::string as_str() const
	{
		return this->pvst_vtx->as_str();
	}
};

void find_walks(const bd::VG &g, RoV &rov);
} // namespace povu::genomics::graph

#endif // POVU_GENOMICS_GRAPH_HPP
