#ifndef MZ_MATRIX_POOL_JOINT_HPP
#define MZ_MATRIX_POOL_JOINT_HPP

#include <cstddef>
#include <driver_types.h>
#include <vector>
#include <vector_types.h>

#include "quilt/types.hpp"

namespace meza::pool::joint
{

template <typename T>
struct full_view;

/**
 * a simple joint pool implementation that allows for allocating contiguous
 * blocks of memory for matrices. The pool keeps track of the free and used
 * space, and can extend itself when needed. The full_view struct provides a
 * view into the pool for a specific matrix, allowing for easy access to the
 * elements.
 *
 * Note: this is a simple implementation and does not include any bounds
 * checking or deallocation of individual blocks. It is designed for simplicity
 * and may not be the most efficient implementation for all use cases.
 */
template <typename T>
struct joint_pool {
private:
	std::size_t free_;
	std::size_t used_ = 0;
	std::size_t start_ = 0; // position of the start of the free region
	std::size_t capacity_ = 0;

	std::vector<T> data_;

	// -----------
	// constructor
	// -----------

	// N is the number of elements in the pool
	explicit joint_pool(std::size_t N)
	    : data_(N, T{}), free_(N), capacity_(N)
	{}

	void extend_pool()
	{
		std::size_t new_capacity = capacity_ * 2;
		data_.resize(new_capacity, T{});
		free_ += (new_capacity - capacity_);
		capacity_ = new_capacity;
	}

	std::size_t alloc(std::size_t need)
	{
		if (!can_allocate(need))
			extend_pool();

		std::size_t alloc_start = start_;
		start_ += need;
		used_ += need;
		free_ -= need;

		return alloc_start; // return the starting index of the
				    // allocated block
	}

public:
	// -----------
	// constructor
	// -----------

	static joint_pool<T> init(std::size_t N)
	{
		return joint_pool<T>{N};
	}

	full_view<T> alloc_full(qt::u32 I, qt::u32 J)
	{
		std::size_t need = static_cast<std::size_t>(I) * J;
		std::size_t start_idx = alloc(need);
		return full_view<T>{this, start_idx, I, J};
	}

	// -------
	// getters
	// -------

	[[nodiscard]]
	std::size_t capacity() const
	{
		return capacity_;
	}

	[[nodiscard]]
	std::size_t used() const
	{
		return used_;
	}

	[[nodiscard]]
	std::size_t free() const
	{
		return free_;
	}

	[[nodiscard]]
	std::size_t start() const
	{
		return start_;
	}

	[[nodiscard]]
	bool can_allocate(std::size_t need) const
	{
		return need <= free();
	}

	const T &at(std::size_t i) const
	{
		return data_[i];
	}

	T &at_mut(std::size_t i)
	{
		return data_[i];
	}

	// ---------
	// modifiers
	// ---------

	void reset()
	{
		std::fill(data_.begin(), data_.end(), T{});
		free_ = data_.size();
		used_ = 0;
		start_ = 0;
	}
};

/**
 * a matrix view that uses a joint pool for storage. The matrix does not own
 * the data, but rather holds a reference to the pool and an index into the
 * pool where its data starts. This allows us to have multiple matrices that
 * share the same underlying storage, which can be more efficient in terms of
 * memory usage and allocation overhead.
 *
 * Note: this is a simple implementation and does not include any bounds
 * checking
 */
template <typename T>
struct full_view {
	using pool_t = joint_pool<T>;

	pool_t *pool_ptr_ = nullptr;
	qt::u32 I_;
	qt::u32 J_;

	std::size_t start_; // start index into the joint pool
	std::size_t len_;   // end index (+1 of the last element/exclusive) into
			    // the joint pool

	full_view(pool_t *pool_ptr, std::size_t start, qt::u32 I, qt::u32 J)
	    : pool_ptr_(pool_ptr), I_(I), J_(J), start_(start)
	{
		std::size_t need = static_cast<std::size_t>(I) * J;
		this->len_ = need;
	}

	[[nodiscard]]
	qt::u32 rows() const
	{
		return I_;
	}

	[[nodiscard]]
	qt::u32 cols() const
	{
		return J_;
	}

	[[nodiscard]]
	std::size_t get_idx(qt::u32 i, qt::u32 j) const
	{
		return start_ + static_cast<std::size_t>(i) * J_ + j;
	}

	std::pair<const T *, qt::u32> get_slice() const
	{
		const T *ptr = &pool_ptr_->at(start_);
		return {ptr, this->len_};
	}

	[[nodiscard]]
	const T &at(qt::u32 i, qt::u32 j) const
	{
		std::size_t idx = this->get_idx(i, j);
		return this->pool_ptr_->at(start_ + idx);
	}

	[[nodiscard]]
	T &at_mut(qt::u32 i, qt::u32 j)
	{
		std::size_t idx = this->get_idx(i, j);
		return this->pool_ptr_->at_mut(start_ + idx);
	}
};

} // namespace meza::pool::joint
#endif // MZ_MATRIX_POOL_JOINT_HPP
