// #include "meza/view/matrix_view.hpp"

// #include <cstddef>
// #include <stdexcept>
// #include <string>
// #include <string_view>

// #include "quilt/types.hpp" // for tc

// namespace meza::matrix_view
// {
// [[nodiscard]]
// std::size_t get_idx(qt::u32 i, qt::u32 j, layout l, const qt::u32 I,
//		    const qt::u32 J)
// {
//	if (i >= I || j >= J) {
//		std::string err; //  = pv_cmp::format(
//				 // "{} out of bounds access [{},{}]", MODULE,
//				 // i, j);
//		throw std::out_of_range(err);
//	}

//	switch (l) {
//	case layout::DenseRowMajor:
//		return static_cast<std::size_t>(i) * J + j;
//	case layout::LowerSymmetricSquare:
//		if (I != J) {
//			std::string err; // =
//					 // pv_cmp::format("{}
//					 // LowerSymmetricSquare "
//					 // "requires I==J got {}x{}", MODULE,
//					 // I, J);
//			throw std::invalid_argument(err);
//		}
//		if (i < j)
//			std::swap(i, j);
//		return static_cast<std::size_t>(i) * (i + 1) / 2 + j;
//	case layout::RepeatedRow:
//		return j;
//	default:
//		throw std::invalid_argument("");
//		// pv_cmp::format(
//		//"{} invalid layout {}", MODULE, static_cast<int>(l))
//	}
// }
// } // namespace meza::matrix_view
