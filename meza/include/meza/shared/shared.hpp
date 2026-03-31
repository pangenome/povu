#ifndef MZ_SHARED_HPP
#define MZ_SHARED_HPP

#include "quilt/types.hpp" // for qt::u32, qt::u8, qt::op_t

namespace meza::shared
{
inline constexpr std::string_view MODULE = "meza::shared";

enum class layout : qt::u8 {
	DenseRowMajor,
	LowerSymmetricSquare,
	RepeatedRow,
};
}; // namespace meza::shared
#endif // MZ_SHARED_HPP
