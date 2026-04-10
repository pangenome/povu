#ifndef MEZA_MATRIX_POOL_SPLIT_TYPES_HPP
#define MEZA_MATRIX_POOL_SPLIT_TYPES_HPP

#include <cstddef>
#include <ostream>
#include <string>
#include <string_view>

#include "meza/shared/shared.hpp" // for layout
#include "meza/view/view.hpp"
#include "quilt/types.hpp"

namespace meza::pool
{
using layout = meza::shared::layout;
using qt::u32;
using qt::u8;

using ov_mat_t = meza::view::ov_matrix<u8, std::string, std::string>;

enum comparison_op : qt::u8 {
	sum,
	bitwise_xor,
};

enum class pool_region : qt::u8 {
	Reference,
	Filter,
	Xor,
};

inline constexpr std::string_view to_string_view(pool_region region)
{
	switch (region) {
	case pool_region::Reference:
		return "Reference";
	case pool_region::Filter:
		return "Filter";
	case pool_region::Xor:
		return "Xor";
	default:
		return "Unknown";
	}
}

inline std::string to_string(pool_region region)
{
	return std::string(to_string_view(region));
}

inline std::ostream &operator<<(std::ostream &os, pool_region region)
{
	return os << to_string_view(region);
}

struct comparison_matrices {
	ov_mat_t ref;
	ov_mat_t filter;
	ov_mat_t xor_result;

	u32 j_offset = 0; // offset in the pool for the filter matrix

	void dbg_print() const
	{
		std::cerr << "Reference Matrix:\n";
		ref.dbg_print();
		std::cerr << "Filter Matrix:\n";
		filter.dbg_print();
		std::cerr << "XOR Result Matrix:\n";
		xor_result.dbg_print();
	}
};

using rov_mat_set = comparison_matrices;

template <typename T>
void print_ptr(std::ostream &os, T *ptr)
{
	os << "addr=0x" << std::hex << reinterpret_cast<std::uintptr_t>(ptr)
	   << std::dec << "\n";
}

/**
 * the starting offsets for each region of the pool
 */
struct partition_offsets {
	std::size_t ref_start = 0;
	std::size_t filter_start = 0;
	std::size_t xor_start = 0;

	friend std::ostream &operator<<(std::ostream &os,
					const partition_offsets &po)
	{
		os << "PoolOffsets {ref_start= " << po.ref_start
		   << ", filter_start= " << po.filter_start
		   << ", xor_start= " << po.xor_start << "}";
		return os;
	}
};

/**
 * the number of elements for each region of the pool
 */
struct partition_sizes {
	std::size_t ref_len;
	std::size_t filter_len;
	std::size_t xor_len;

	friend std::ostream &operator<<(std::ostream &os,
					const partition_sizes &ps)
	{
		os << "PoolSplit {ref_len= " << ps.ref_len
		   << ", filter_len= " << ps.filter_len
		   << ", xor_len= " << ps.xor_len << "}";
		return os;
	}
};

struct partition_usage {
	std::size_t ref_used;
	std::size_t filter_used;
	std::size_t xor_used;

	friend std::ostream &operator<<(std::ostream &os,
					const partition_usage &u)
	{
		os << "Usage {ref_used= " << u.ref_used
		   << ", filter_used= " << u.filter_used
		   << ", xor_used= " << u.xor_used << "}";
		return os;
	}
};

} // namespace meza::pool

#endif // MEZA_MATRIX_POOL_SPLIT_TYPES_HPP
