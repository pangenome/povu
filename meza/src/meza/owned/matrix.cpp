#include "meza/owned/matrix.hpp"

#include <cstddef>
#include <stdexcept>
#include <string>
#include <string_view>

#include <log/location.hpp> // for LOG_HERE
#include <quilt/shim.hpp>   // for qs::contains, qs::format
#include <quilt/types.hpp>  // for qt::u32, qt::u8, qt::op_t

namespace meza::matrix
{
[[nodiscard]]
std::size_t get_idx(qt::u32 i, qt::u32 j, layout l, const qt::u32 I,
		    const qt::u32 J)
{
	if (i >= I || j >= J) {
		std::string err = qs::format("{} out of bounds access [{},{}]",
					     LOG_HERE, i, j);
		throw std::out_of_range(err);
	}

	switch (l) {
	case layout::DenseRowMajor:
		return static_cast<std::size_t>(i) * J + j;
	case layout::LowerSymmetricSquare:
		if (I != J) {
			std::string err =
				qs::format("{} LowerSymmetricSquare requires I "
					   "== J got {} x {} ",
					   LOG_HERE, I, J);
			throw std::invalid_argument(err);
		}
		if (i < j)
			std::swap(i, j);
		return static_cast<std::size_t>(i) * (i + 1) / 2 + j;
	case layout::RepeatedRow:
		return j;
	}

	std::string err = qs::format("{} Invalid layout: {}", LOG_HERE,
				     static_cast<int>(l));
	throw std::invalid_argument(err);
}
} // namespace meza::matrix
