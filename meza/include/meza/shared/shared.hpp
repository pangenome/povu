#ifndef MZ_SHARED_HPP
#define MZ_SHARED_HPP

#include "quilt/types.hpp" // for qt::u8

namespace meza::shared
{

enum class layout : qt::u8 {
	DenseRowMajor,
	LowerSymmetricSquare,
	RepeatedRow,
};
}; // namespace meza::shared
#endif // MZ_SHARED_HPP
