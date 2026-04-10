#ifndef MEZA_MATRIX_POOL_SPLIT_HPP
#define MEZA_MATRIX_POOL_SPLIT_HPP

#include <cstddef>
#include <stdexcept>

#include "meza/pool/split_pool_types.hpp"
#include "meza/view/view.hpp"
#include "quilt/types.hpp"

namespace meza::pool
{

template <typename T>
struct matrix_pool {
public:
	/* ================= ================== */
	partition_sizes reset()
	{
		auto [ref_len, filter_len, xor_len] = default_lengths();

		ref_start_ = 0;
		filter_start_ = ref_start_ + ref_len;
		xor_start_ = filter_start_ + filter_len;

		head_ref_ = ref_start_;
		head_filter_ = filter_start_;
		head_xor_ = xor_start_;

		return {ref_len, filter_len, xor_len};
	}

	[[nodiscard]] partition_sizes default_lengths() const
	{
		std::size_t ref_len = capacity_ / 3;
		std::size_t filter_len = capacity_ / 3;
		std::size_t xor_len = capacity_ - ref_len - filter_len;

		return {ref_len, filter_len, xor_len};
	}

	[[nodiscard]] partition_offsets default_offsets() const
	{
		auto [ref_len, filter_len, _] = default_lengths();

		partition_offsets po;
		po.ref_start = 0;
		po.filter_start = po.ref_start + ref_len;
		po.xor_start = po.filter_start + filter_len;

		return po;
	}

	[[nodiscard]] partition_offsets offsets() const
	{
		return {ref_start_, filter_start_, xor_start_};
	}

	[[nodiscard]] partition_usage used() const
	{
		return {head_ref_ - ref_start_, head_filter_ - filter_start_,
			head_xor_ - xor_start_};
	}

	void clear()
	{
		// overwrite with zeros for safety/debugging
		std::fill(host_storage_, host_storage_ + capacity_, T{});
		reset();
	}

	[[nodiscard]] size_t empty() const
	{
		return head_ref_ == ref_start_ &&
		       head_filter_ == filter_start_ && head_xor_ == xor_start_;
	}

	/**
	 * is full at 90% of any of the regions
	 * to trigger an early flush before we run out of memory
	 */
	[[nodiscard]] bool is_full() const
	{
		auto [ref_used, filter_used, xor_used] = this->used();
		auto [ref_capacity, filter_capacity, xor_capacity] =
			this->default_lengths();

		double ref_util = static_cast<double>(ref_used) /
				  static_cast<double>(ref_capacity);
		double filter_util = static_cast<double>(filter_used) /
				     static_cast<double>(filter_capacity);
		double xor_util = static_cast<double>(xor_used) /
				  static_cast<double>(xor_capacity);

		return ref_util >= 0.9 || filter_util >= 0.9 || xor_util >= 0.9;
	}

	/* ================= ================== */
	[[nodiscard]] std::size_t ref_start() const noexcept
	{
		return ref_start_;
	}

	[[nodiscard]] std::size_t filter_start() const noexcept
	{
		return filter_start_;
	}

	[[nodiscard]] std::size_t xor_start() const noexcept
	{
		return xor_start_;
	}

	[[nodiscard]] T *ref_start_ptr() noexcept
	{
		return host_storage_ + ref_start_;
	}

	[[nodiscard]] T *filter_start_ptr() noexcept
	{
		return host_storage_ + filter_start_;
	}

	[[nodiscard]] T *xor_start_ptr() noexcept
	{
		return host_storage_ + xor_start_;
	}

	[[nodiscard]] std::size_t ref_used() const noexcept
	{
		return head_ref_ - ref_start_;
	}

	[[nodiscard]] std::size_t filter_used() const noexcept
	{
		return head_filter_ - filter_start_;
	}

	[[nodiscard]] std::size_t xor_used() const noexcept
	{
		return head_xor_ - xor_start_;
	}

	T *host_data() noexcept
	{
		return host_storage_;
	}

	const T *host_data() const noexcept
	{
		return host_storage_;
	}

	/* ================= ================== */

	template <typename U, typename W>
	[[nodiscard]] meza::view::ov_matrix<T, U, W>
	alloc_ov_matrix(qt::u32 I, qt::u32 J, pool_region region)
	{
		T *ptr = alloc_full(I, J, region);
		return meza::view::ov_matrix<T, U, W>{I, J, ptr};
	}

	/* ========================== constructors ====================== */

	// we delete copy constructor and copy assignment operator to prevent
	// copying of the pool, which could lead to double frees and other
	// issues
	matrix_pool(const matrix_pool &) = delete;
	matrix_pool &operator=(const matrix_pool &) = delete;
	matrix_pool() = delete;

	// move operators
	matrix_pool(matrix_pool &&other) noexcept
	    : capacity_(other.capacity_), ref_start_(other.ref_start_),
	      filter_start_(other.filter_start_), xor_start_(other.xor_start_),
	      head_ref_(other.head_ref_), head_filter_(other.head_filter_),
	      head_xor_(other.head_xor_), host_storage_(other.host_storage_)
	{
		other.host_storage_ = nullptr;
		other.capacity_ = 0;
		other.ref_start_ = other.filter_start_ = other.xor_start_ = 0;
		other.head_ref_ = other.head_filter_ = other.head_xor_ = 0;
	}

	matrix_pool &operator=(matrix_pool &&other) noexcept
	{
		if (this != &other) {
			delete[] host_storage_;
			capacity_ = other.capacity_;
			ref_start_ = other.ref_start_;
			filter_start_ = other.filter_start_;
			xor_start_ = other.xor_start_;
			head_ref_ = other.head_ref_;
			head_filter_ = other.head_filter_;
			head_xor_ = other.head_xor_;
			host_storage_ = other.host_storage_;

			other.host_storage_ = nullptr;
			other.capacity_ = 0;
			other.ref_start_ = other.filter_start_ =
				other.xor_start_ = 0;
			other.head_ref_ = other.head_filter_ = other.head_xor_ =
				0;
		}
		return *this;
	}

	/**
	 * Initialise the pool with a given capacity (number of elements)
	 *
	 * @param N the total number of elements in the pool
	 */
	explicit matrix_pool(size_t N) : capacity_(N)
	{
		host_storage_ = new T[N];
		reset();
	}

	/**
	 * Factory method to create a matrix pool with a specified size in
	 * megabytes.
	 */
	static matrix_pool create_from_megabytes(std::size_t space_mb = 512)
	{
		std::size_t bytes = space_mb * 1024 * 1024;
		std::size_t elements{bytes / sizeof(T)};
		return matrix_pool{elements};
	}

	~matrix_pool()
	{
		delete[] host_storage_;
	}

private:
	/* ============= private data members ======================== */
	// max no. of elements that can be stored in the pool
	std::size_t capacity_ = 0;

	// The starting offsets for each region of the pool.
	// These are set at the beginning and do not change, while the "head"
	// offsets track how much of that region has been allocated.
	std::size_t ref_start_ = 0;
	std::size_t filter_start_ = 0;
	std::size_t xor_start_ = 0;

	// The "head" tracks how much of that region has been allocated.
	// These should always be greater than or equal to the
	// corresponding start offsets, and less than the next region's start
	// offset (or the total capacity for the xor region)
	std::size_t head_ref_ = 0;
	std::size_t head_filter_ = 0;
	std::size_t head_xor_ = 0;

	T *host_storage_ = nullptr;

	/* ============= private methods ======================== */

	[[nodiscard]] static size_t required_elements(std::size_t I,
						      std::size_t J, layout lo)
	{
		switch (lo) {
		case layout::DenseRowMajor:
			return I * J;
		case layout::LowerSymmetricSquare:
			if (I != J)
				throw std::invalid_argument(
					"symmetric needs square");

			return I * (I + 1) / 2;
		case layout::RepeatedRow:
			return J;
		}

		throw std::logic_error("unhandled layout");
	}

	[[nodiscard]] bool can_allocate(std::size_t need,
					pool_region region) const
	{
		switch (region) {
		case pool_region::Reference:
			return this->head_ref_ + need < head_filter_;
		case pool_region::Filter:
			return this->head_filter_ + need < head_xor_;
		case pool_region::Xor:
			return this->head_xor_ + need < capacity_;
		}

		throw std::logic_error("invalid layout");
	}

	[[nodiscard]] T *alloc(std::size_t need, pool_region region)
	{
		if (!can_allocate(need, region))
			throw std::runtime_error("not enough space in region " +
						 to_string(region));

		T *ptr;
		switch (region) {
		case pool_region::Reference:
			ptr = host_storage_ + head_ref_;
			head_ref_ += need;
			return ptr;
		case pool_region::Filter:
			ptr = host_storage_ + head_filter_;
			head_filter_ += need;
			return ptr;
		case pool_region::Xor:
			ptr = host_storage_ + head_xor_;
			head_xor_ += need;
			return ptr;
		}

		throw std::logic_error("invalid pool region");
	}

	[[nodiscard]] T *alloc_full(qt::u32 I, qt::u32 J, pool_region region)
	{
		constexpr layout DRM = layout::DenseRowMajor;
		std::size_t need = required_elements(I, J, DRM);
		return alloc(need, region);
	}
};

} // namespace meza::pool

#endif // MEZA_MATRIX_POOL_SPLIT_HPP
