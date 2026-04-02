#ifndef MZ_MATRIX_POOL_SPLIT_HPP
#define MZ_MATRIX_POOL_SPLIT_HPP

#include <cstddef>
#include <ostream>
// #include <set>
#include <stdexcept>
#include <string_view>
#include <vector>

#if MEZA_USE_CUDA
#include <cuda_runtime.h>
#include <cuda_runtime_api.h>
#include <driver_types.h>
#include <vector_types.h>

#include "meza/ops/ops.cuh"
#else
#include "meza/ops/ops.hpp"
#endif

#include "meza/shared/shared.hpp" // for layout
#include "meza/view/view.hpp"
#include "quilt/types.hpp"

namespace meza::pool::split
{
inline constexpr std::string_view MODULE = "meza::matrix_pool";

// -------
// aliases
// -------

using layout = meza::shared::layout;
using qt::u32;
using qt::u8;

using ov_mat_t = meza::view::ov_matrix<u8, std::string, std::string>;

// --------------
// helper structs
// --------------

enum comparison_op : qt::u8 {
	sum,
	bitwise_xor,
};

enum class pool_region : qt::u8 {
	Reference,
	Filter,
	Xor,
};

struct rov_mat_set {
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

using comparison_matrices = rov_mat_set;

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
	qt::u32 ref_len;
	qt::u32 filter_len;
	qt::u32 xor_len;

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

/**
 *
 * collapse a set of matrices into a single matrix
 *
 * 2^28 (268,435,456 elements)
 * 2^30 (1,073,741,824 elements)
 * 2^32 (4,294,967,296 elements)
 *
 */
template <typename T>
struct matrix_pool {
private:
	// ------------
	// data members
	// ------------

	std::size_t capacity; // no. of elements in the pool

	std::size_t ref_start = 0;    // start of reference matrices
	std::size_t filter_start = 0; // start of filter matrices
	std::size_t xor_start = 0;    // start of xor matrices

	std::size_t head_ref = 0;    // current head for reference matrices
	std::size_t head_filter = 0; // current head for filter matrices
	std::size_t head_xor = 0;    // current head for xor matrices

	T *host_storage = nullptr;

	// optional-used for CUDA:
	// separate device storage for matrices: reference, filter & xor
	T *d_ref_store = nullptr;
	T *d_filter_store = nullptr;
	T *d_xor_store = nullptr;

	// haps xor is a different set of xors
	T *h_haps_xor_store = nullptr;
	T *d_haps_xor_store = nullptr;

	// haps sum is a different set of xors
	T *h_haps_sum_store = nullptr;
	T *d_haps_sum_store = nullptr;

	// ------------------------
	// private helper functions
	// ------------------------

	static size_t required_elements(std::size_t I, std::size_t J, layout lo)
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

		return 0;
	}

	[[nodiscard]] partition_sizes default_lens() const
	{
		std::size_t ref_len = capacity / 3;
		std::size_t filter_len = capacity / 3;
		std::size_t xor_len = capacity - ref_len - filter_len;

		return {(qt::u32)ref_len, (qt::u32)filter_len,
			(qt::u32)xor_len};
	}

	void cuda_setup_haps_xor()
	{
		std::size_t default_bytes = 50; // we want 50MiB of u8s
		std::size_t exp_bytes = default_bytes * 1024 * 1024;
		cudaError_t e;

		// --------------------------------
		// allocate device mem for haps xor
		// --------------------------------
		e = cudaMalloc((void **)&d_haps_xor_store, exp_bytes);
		if (e != cudaSuccess)
			throw std::runtime_error(cudaGetErrorString(e));

		// allocate pinned mem on the host
		e = cudaMallocHost((void **)&h_haps_xor_store, exp_bytes);
		if (e != cudaSuccess)
			throw std::runtime_error(cudaGetErrorString(e));

		// --------------------------------
		// allocate device mem for haps sum
		// --------------------------------
		e = cudaMalloc((void **)&d_haps_sum_store, exp_bytes);
		if (e != cudaSuccess)
			throw std::runtime_error(cudaGetErrorString(e));

		e = cudaMallocHost((void **)&h_haps_sum_store, exp_bytes);
		if (e != cudaSuccess)
			throw std::runtime_error(cudaGetErrorString(e));
	}

public:
	// -----------------------
	// public helper functions
	// -----------------------

	[[nodiscard]] partition_offsets default_offsets() const
	{
		partition_sizes ps = default_lens();
		auto [ref_len, filter_len, xor_len] = ps;

		partition_offsets po;
		po.ref_start = 0;
		po.filter_start = po.ref_start + ref_len;
		po.xor_start = po.filter_start + filter_len;

		return po;
	}

	/**
	 * reset the pool by moving the head back to the beginning
	 * Note: this does not actually free any memory, but allows us
	 * to reuse the existing storage for new matrices. The caller is
	 * responsible for ensuring that any matrices allocated from the
	 * pool are no longer in use before calling reset.
	 */
	partition_sizes reset()
	{
		partition_sizes ps = default_lens();
		auto [ref_len, filter_len, xor_len] = ps;

		ref_start = 0;
		filter_start = ref_start + ref_len;
		xor_start = filter_start + filter_len;

		head_ref = ref_start;
		head_filter = filter_start;
		head_xor = xor_start;

		return ps;
	}

	void clear()
	{
		// overwrite with zeros for safety/debugging
		// auto [ref_used, filter_used, xor_used]this->used();
		std::fill(host_storage, host_storage + capacity, T{});
		std::size_t exp_bytes = 50 * 1024 * 1024; // we want 5MiB of u8s
		std::fill(h_haps_xor_store, h_haps_xor_store + exp_bytes, T{});
		std::fill(h_haps_sum_store, h_haps_sum_store + exp_bytes, T{});

		reset();
	}

	[[nodiscard]] partition_usage usage_split()
	{
		std::size_t ref_len = capacity / 3;
		std::size_t filter_len = capacity / 3;
		std::size_t xor_len = capacity - ref_len - filter_len;

		return {ref_len, filter_len, xor_len};
	}

	/**
	 * returns the number of elements used in each region of the pool
	 */
	[[nodiscard]] partition_usage used() const
	{
		return {
			head_ref - ref_start,	    //
			head_filter - filter_start, //
			head_xor - xor_start	    //
		};
	}

	[[nodiscard]] size_t free() const
	{
		return capacity - used();
	}

	[[nodiscard]] bool empty() const
	{
		auto [exp_head_ref, exp_head_filter, exp_head_xor] =
			default_offsets();

		return head_ref == exp_head_ref &&
		       head_filter == exp_head_filter &&
		       head_xor == exp_head_xor;
	}

	[[nodiscard]] bool is_full()
	{
		auto [ref_used, _, __] = this->used();
		auto [ref_capacity, ___, ____] = this->usage_split();

		// is full at 80 %
		return (static_cast<double>(ref_used) /
			static_cast<double>(ref_capacity)) >= 0.8;
	}

	// ------------
	// constructors
	// ------------

	// delete copy constructor and copy assignment operator to prevent
	// copying of the pool, which could lead to double frees and other
	// issues
	matrix_pool(const matrix_pool &) = delete;
	matrix_pool &operator=(const matrix_pool &) = delete;

	matrix_pool() = default;

	explicit matrix_pool(size_t capacity_elems) : capacity(capacity_elems)
	{
#if MEZA_USE_CUDA
		auto alloc_cuda = [](void **ptr, size_t bytes)
		{
			cudaError_t e = cudaMalloc(ptr, bytes);
			if (e != cudaSuccess) {
				*ptr = nullptr;
				throw std::runtime_error(cudaGetErrorString(e));
			}
		};

		size_t bytes = capacity * sizeof(T);
		cudaError_t err = cudaMallocHost((void **)&host_storage, bytes);
		if (err != cudaSuccess)
			throw std::runtime_error(cudaGetErrorString(err));

		// bytes div by three
		// allocate device storage for the pool
		auto [ref_len, filter_len, xor_len] = reset();

		size_t bytes0 = ref_len * sizeof(T);
		size_t bytes1 = filter_len * sizeof(T);
		size_t bytes2 = xor_len * sizeof(T);

		try {
			alloc_cuda((void **)&d_ref_store, bytes0);
			alloc_cuda((void **)&d_filter_store, bytes1);
			alloc_cuda((void **)&d_xor_store, bytes2);
		}
		catch (...) {
			// Free any allocated resources before rethrowing
			if (d_ref_store)
				cudaFree(d_ref_store);
			if (d_filter_store)
				cudaFree(d_filter_store);
			if (d_xor_store)
				cudaFree(d_xor_store);
			if (host_storage)
				cudaFreeHost(host_storage);
			throw; // Rethrow the original exception
		}

#else
		// fallback if you want non-CUDA builds to work
		host_storage = (T *)std::malloc(capacity * sizeof(T));
		if (!host_storage)
			throw std::bad_alloc{};
#endif
	}

	~matrix_pool()
	{
#if MEZA_USE_CUDA
		auto free_cuda = [](void *ptr)
		{
			if (ptr)
				cudaFree(ptr);
		};

		free_cuda(d_ref_store);
		free_cuda(d_filter_store);
		free_cuda(d_xor_store);

		if (host_storage)
			cudaFreeHost(host_storage);
#else
		std::free(host_storage);
#endif
	}

	// ----------
	// allocators
	// ----------

	T *alloc(std::size_t need, pool_region region)
	{
		if (region == pool_region::Reference) {
			T *ptr = host_storage + head_ref;
			head_ref += need;
			return ptr;
		}
		else if (region == pool_region::Filter) {
			T *ptr = host_storage + head_filter;
			head_filter += need;
			return ptr;
		}
		else if (region == pool_region::Xor) {
			T *ptr = host_storage + head_xor;
			head_xor += need;
			return ptr;
		}
		else {
			std::string err =
				qs::format("{} Invalid pool region: {}", MODULE,
					   static_cast<int>(region));
			throw std::logic_error(err);
		}
	}

	T *alloc_full(qt::u32 I, qt::u32 J, pool_region region)
	{
		std::size_t need =
			required_elements(I, J, layout::DenseRowMajor);
		return alloc(need, region);
	}

	T *alloc_repeated_row(qt::u32 I, qt::u32 J, pool_region region)
	{
		std::size_t need = required_elements(I, J, layout::RepeatedRow);
		return alloc(need, region);
	}

	static matrix_pool &init()
	{
		// 512 MB is a reasonable default size for the pool,
		// also		static matrix_pool pool{1u << 28};
		static constexpr std::size_t space_mb = 512;
		static constexpr size_t BYTES = space_mb * 1024 * 1024;
		std::size_t elements{BYTES / sizeof(T)};
		static matrix_pool pool{elements};

		pool.cuda_setup_haps_xor();

		return pool;
	}

	[[nodiscard]] bool can_allocate(qt::u32 I, qt::u32 J, layout lo)
	{
		std::size_t need = required_elements(I, J, lo);
		if (layout::DenseRowMajor == lo)
			return this->head_ref + need <
			       head_filter; // check against filter head
		else if (layout::LowerSymmetricSquare == lo)
			return this->head_filter + need <
			       head_xor; // check against xor head
		else if (layout::RepeatedRow == lo)
			return this->head_xor + need <
			       capacity; // check against total capacity
		else
			throw std::logic_error("invalid layout");
	}

	template <typename U, typename W>
	meza::view::ov_matrix<T, U, W> alloc_ov_matrix(qt::u32 I, qt::u32 J,
						       pool_region region)
	{
		T *ptr = alloc_full(I, J, region);
		return meza::view::ov_matrix<T, U, W>{I, J, ptr};
	}

	template <typename U, typename W>
	meza::view::ref_matrix<T, U, W> alloc_ref_matrix(qt::u32 I, qt::u32 J,
							 pool_region region)
	{
		T *ptr = alloc_repeated_row(I, J, region);
		return meza::view::ref_matrix<T, U, W>{I, J, ptr};
	}

	// -----------------
	// device operations
	// -----------------

	void copy_to_mem(cudaStream_t stream = 0)
	{
		// how many elements used in each region
		const size_t used_ref = head_ref - ref_start;
		const size_t used_filter = head_filter - filter_start;
		const size_t used_xor = head_xor - xor_start;

		const std::size_t ref_capacity = filter_start - ref_start;
		const std::size_t filter_capacity = xor_start - filter_start;
		const std::size_t xor_capacity = capacity - xor_start;

		const T *h_ref = host_storage + ref_start;
		const T *h_filter = host_storage + filter_start;
		const T *h_xor = host_storage + xor_start;

		std::vector<T> ref_vec{ref_capacity, T{}};
		std::vector<T> filter_vec{filter_capacity, T{}};
		std::vector<T> xor_vec{xor_capacity, T{}};

		for (size_t i = 0; i < used_ref; i++)
			ref_vec[i] = h_ref[i];

		for (size_t i = 0; i < used_filter; i++)
			filter_vec[i] = h_filter[i];

		for (size_t i = 0; i < used_xor; i++)
			xor_vec[i] = h_xor[i];
	}

	void copy_to_device(cudaStream_t stream = 0)
	{
#if MEZA_USE_CUDA
		// how many elements used in each region
		const size_t used_ref = head_ref - ref_start;
		const size_t used_filter = head_filter - filter_start;
		const size_t used_xor = head_xor - xor_start;

		const T *h_ref = host_storage + ref_start;
		const T *h_filter = host_storage + filter_start;
		const T *h_xor = host_storage + xor_start;

		cudaError_t e;

		if (used_ref) {
			e = cudaMemcpyAsync(d_ref_store, h_ref,
					    used_ref * sizeof(T),
					    cudaMemcpyHostToDevice, stream);
			if (e != cudaSuccess)
				throw std::runtime_error(cudaGetErrorString(e));
		}

		if (used_filter) {
			e = cudaMemcpyAsync(d_filter_store, h_filter,
					    used_filter * sizeof(T),
					    cudaMemcpyHostToDevice, stream);
			if (e != cudaSuccess)
				throw std::runtime_error(cudaGetErrorString(e));
		}

		if (used_xor) {
			e = cudaMemcpyAsync(d_xor_store, h_xor,
					    used_xor * sizeof(T),
					    cudaMemcpyHostToDevice, stream);
			if (e != cudaSuccess)
				throw std::runtime_error(cudaGetErrorString(e));
		}
#endif
	}

	//	void run_filter(u32 I, u32 J)
	//	{
	// #if MEZA_USE_CUDA
	//		xor_on_device(I, J);
	// #else
	//		std::cerr << qs::format("{}{}", MODULE, __func__)
	//			  << "operation on device is not "
	//			     "implemented in non-CUDA mode\n";
	// #endif
	//	}

	void run_filter(u32 N)
	{
#if MEZA_USE_CUDA
		meza::cuda_ops::cuda_mat_xor(d_ref_store, d_filter_store,
					     d_xor_store, N);
#else
		meza::cpu_ops::cpu_mat_xor(ref_start, filter_start, xor_start,
					   N);
#endif
	}

	// void xor_on_device(qt::u32 N)
	//	{
	// #if MEZA_USE_CUDA
	//		meza::cuda_ops::cuda_mat_xor(d_ref_store,
	// d_filter_store, d_xor_store, N); #endif
	//	}

	[[nodiscard]]
	T *get_haps_xor()
	{
		return h_haps_xor_store;
	}

	[[nodiscard]]
	T *get_haps_sum()
	{
		return h_haps_sum_store;
	}

	void copy_slice_to_device(u32 J, T *h_xor = nullptr,
				  cudaStream_t stream = 0)
	{
		std::size_t size_bytes = J * sizeof(T);

		cudaError_t e;
		e = cudaMemcpyAsync(h_xor, d_xor_store, size_bytes,
				    cudaMemcpyDeviceToHost, stream);
		if (e != cudaSuccess)
			throw std::runtime_error(cudaGetErrorString(e));
	}

	/**
	 * essentially a slice sum
	 */
	void haps_sum(u32 k, u32 pool_j_offset, u32 mat_j, const u32 I,
		      std::vector<T> &hap_res, cudaStream_t stream = 0)
	{
		const u32 J = mat_j;

		u32 a_row = 0; // row
		u32 b_row = k; // col

		u32 k_len = I - k;
		u32 N = k_len * J;

		std::size_t mat_filter_start =
			filter_start + (pool_j_offset * I);

		std::size_t a_off = mat_filter_start + (a_row * J);
		std::size_t b_off = mat_filter_start + (b_row * J);

		for (std::size_t j = 0; j < N; j++) {
			T val_a = host_storage[a_off + j];
			T val_b = host_storage[b_off + j];

			hap_res.emplace_back(val_a + val_b);
		}
	}

	void copy_haps_sum_to_host(std::size_t bytes, cudaStream_t stream = 0)
	{
#if MEZA_USE_CUDA
		cudaMemcpyAsync(h_haps_sum_store, d_haps_sum_store, bytes,
				cudaMemcpyDeviceToHost, stream);
		cudaStreamSynchronize(stream);
#endif
	}

	void copy_haps_xor_to_host(std::size_t bytes, cudaStream_t stream = 0)
	{
#if MEZA_USE_CUDA
		cudaMemcpyAsync(h_haps_xor_store, d_haps_xor_store, bytes,
				cudaMemcpyDeviceToHost, stream);
		cudaStreamSynchronize(stream);
#endif
	}

	void cuda_free_haps_xor()
	{
		if (d_haps_xor_store)
			cudaFree(d_haps_xor_store);
		d_haps_xor_store = nullptr;

		if (h_haps_xor_store)
			cudaFreeHost(h_haps_xor_store);
		h_haps_xor_store = nullptr;
	}

	void hap_sum_gpu(u32 mat_off, u32 len, u32 col_shift, u32 res_shift,
			 cudaStream_t stream) const
	{
		cuda_ops::cuda_haps_sum(d_filter_store, d_haps_sum_store,
					mat_off, col_shift, res_shift, len,
					stream);
	}

	void hap_xor_gpu(u32 mat_off, u32 len, u32 col_shift, u32 res_shift,
			 cudaStream_t stream) const
	{
		cuda_ops::cuda_haps_xor(d_filter_store, d_haps_xor_store,
					mat_off, col_shift, res_shift, len,
					stream);
	}

	void hap_xor_cpu(u32 pool_off, u32 len, u32 col_shift,
			 u32 res_shift) const
	{
		u32 mat_filter_start = filter_start + pool_off;
		u32 mat_filter_end = mat_filter_start + len;

		qt::u32 j{res_shift};
		qt::u32 i{mat_filter_start};

		for (; i < mat_filter_end; i++, j++) {
			T val_a = host_storage[i];
			T val_b = host_storage[i + col_shift];

			h_haps_xor_store[j] = val_a ^ val_b;
		}
	}

	/**
	 * essentially a slice xor
	 */
	void haps_xor(u32 pool_off, u32 len, u32 col_shift, u32 res_shift,
		      cudaStream_t stream) const
	{
#if MEZA_USE_CUDA
		hap_xor_gpu(pool_off, len, col_shift, res_shift, stream);
#else
		hap_xor_cpu(pool_off, len, col_shift, res_shift);
#endif
	}

	/**
	 * essentially a slice sum
	 */
	void hap_sum(u32 pool_off, u32 len, u32 col_shift, u32 res_shift,
		     cudaStream_t stream) const
	{
		hap_sum_gpu(pool_off, len, col_shift, res_shift, stream);
		// #if CONVO_USE_CUDA
		//		hap_sum_gpu(pool_off, len, col_shift,
		// res_shift); #else		hap_xor_cpu(pool_off,
		// len, col_shift, res_shift); #endif
	}

	void copy_to_host_thirds(pool_region region, cudaStream_t stream = 0)
	{
#if MEZA_USE_CUDA
		const size_t used_ref = head_ref - ref_start;
		const size_t used_filter = head_filter - filter_start;
		const size_t used_xor = head_xor - xor_start;

		T *h_ref = host_storage + ref_start;
		T *h_filter = host_storage + filter_start;
		T *h_xor = host_storage + xor_start;

		cudaError_t e;

		if (region == pool_region::Reference && used_ref) {
			e = cudaMemcpyAsync(h_ref, d_ref_store,
					    used_ref * sizeof(T),
					    cudaMemcpyDeviceToHost, stream);
			if (e != cudaSuccess)
				throw std::runtime_error(cudaGetErrorString(e));
		}

		if (region == pool_region::Filter && used_filter) {
			e = cudaMemcpyAsync(h_filter, d_filter_store,
					    used_filter * sizeof(T),
					    cudaMemcpyDeviceToHost, stream);
			if (e != cudaSuccess)
				throw std::runtime_error(cudaGetErrorString(e));
		}

		if (region == pool_region::Xor && used_xor) {
			std::size_t bytes = used_xor * sizeof(T);
			// std::cerr << "copying " << used_xor << " " <<
			// bytes
			//	  << " bytes from device to host for Xor
			//"	     "region\n";
			e = cudaMemcpyAsync(h_xor, d_xor_store, bytes,
					    cudaMemcpyDeviceToHost, stream);
			if (e != cudaSuccess)
				throw std::runtime_error(cudaGetErrorString(e));
		}
#endif
	}

	void sync_device(cudaStream_t stream = 0)
	{
#if MEZA_USE_CUDA
		cudaError_t e = cudaStreamSynchronize(stream);
		if (e != cudaSuccess)
			throw std::runtime_error(cudaGetErrorString(e));
#endif
	}

	// -------
	// getters
	// -------

	void dbg_print(std::ostream &os) const
	{

		auto to_u32 = [](T v) -> u32
		{
			return static_cast<u32>(v);
		};

		auto [ref_used, filter_used, xor_used] = used();

		os << "ref: ";
		for (size_t i = 0; i < ref_used; i++)
			os << to_u32(host_storage[ref_start + i]) << ", ";
		std::cerr << "\n";

		os << "\nfilter: ";
		for (size_t i = 0; i < filter_used; i++)
			os << to_u32(host_storage[filter_start + i]) << ", ";
		std::cerr << "\n";

		os << "\nxor: ";
		for (size_t i = 0; i < xor_used; i++)
			os << to_u32(host_storage[xor_start + i]) << ", ";
		std::cerr << "\n";
	}
};

// /**
//  * haplotype comparison matrix
//  *
//  *
//  * a square matrix
//  *
//  * stores from diagonal 1 upwards
//  */
// struct hap_comp_matrix {
// private:
//	// a reference to the filter matrix
//	const ov_mat_t &filter;

//	// the offset in the pool where the haplotype comparison matrix
//	// starts.
//	// This is used to calculate the correct offsets for
//	// accessing the data in the pool when performing comparisons.
//	qt::u32 pool_offset;

//	// number of haplotypes
//	// also no of rows in the filter matrix
//	// also the number of rows and cols in the haplotype comparison matrix
//	//
//	// matrix the number of rows and cols in the hap comp matrix
//	qt::u32 H;

//	// the total number of elements in the haplotype comparison matrix
//	// diagonal 1 upwards
//	qt::u32 data_size_;

//	// pointers to the data in the pool for the xor and sum results
//	// these are used to store the results of the comparisons for each
//	// pair of haplotypes. The data is stored in a flattened format,
//	// where the comparisons for each pair of haplotypes are stored
//	// contiguously in memory. The offsets for accessing the correct
//	// data for each pair of haplotypes are calculated based on the
//	// k value and the number of comparisons that come before it in the
//	// upper triangle of the matrix.
//	qt::u8 *xor_data_;
//	qt::u8 *sum_data_;

//	std::vector<qt::u32> xor_ps; // prefix sum

//	// -----------------
//	// helpers (private)
//	// -----------------
//	/**
//	 * calculates the number of elements in the upper triangle of an
//	 * n x n matrix excluding the diagonal
//	 */
//	[[nodiscard]]
//	static qt::u32 triangle_number(qt::u32 n)
//	{
//		return n * (n + 1) / 2;
//	}

//	/**
//	 * total filter size
//	 */
//	[[nodiscard]]
//	qt::u32 filter_size() const
//	{
//		return filter.rows() * filter.rows();
//	}

//	[[nodiscard]]
//	qt::u32 k_len(qt::u32 k) const
//	{
//		return this->H - k;
//	}

//	/**
//	 * Ensure H is set before calling this fn
//	 *
//	 * calculates the expected size of the haplotype comparison
//	 * matrix based on the number of haplotypes (H) and the number
//	 * of elements per row in the filter matrix
//	 *
//	 */
//	[[nodiscard]] qt::u32 expected_size(const ov_mat_t &f)
//	{
//		qt::u32 exp_comparisons =
//			(k_len(1) * H) - triangle_number(H - 1);
//		qt::u32 elements_per_row = f.cols();
//		qt::u32 exp_size = exp_comparisons * elements_per_row;

//		// std::cerr << __func__ << " l_k(1) " << k_len(1)
//		//	  << " exp comparisons " << exp_comparisons
//		//	  << " elements per row " << elements_per_row
//		//	  << " exp size " << exp_size << "\n";

//		return exp_size;
//	}

//	/**
//	 * @brief Given k & k-offset, calculate the row & col pair in
//	 * the matrix
//	 */
//	static qt::up_t<qt::u32> comp_hap_pair(qt::u32 k, qt::u32 k_offset)
//	{
//		return {k_offset, k + k_offset}; // {h1, h2}
//	}

// public:
//	// ------------
//	// constructors
//	// ------------

//	// delete default constructor
//	hap_comp_matrix() = delete;

//	hap_comp_matrix(qt::u8 *pool_xor_ptr, qt::u8 *pool_sum_ptr,
//			const ov_mat_t &filter, u32 pool_offset)
//	    : filter(filter), pool_offset(pool_offset), H(filter.rows()),
//	      data_size_(expected_size(filter)), xor_data_(pool_xor_ptr),
//	      sum_data_(pool_sum_ptr)
//	{}

//	// ----------------
//	// helpers (public)
//	// ----------------

//	/**
//	 * calculates the offset in the haplotype comparison matrix for
//	 * a given k value. This is based on the number of comparisons
//	 * that come before the k-th diagonal in the upper triangle of
//	 * the matrix.
//	 */
//	[[nodiscard]]
//	qt::u32 k_offset(qt::u32 k) const
//	{
//		return (k * k_len(1)) - triangle_number(k - 1) - (k_len(k));
//	}

//	void run_in_haps(matrix_pool<qt::u8> &ov_pool, comparison_op op,
//			 cudaStream_t stream = 0)
//	{
//		qt::u32 J = filter.cols();
//		qt::u32 K = H;

//		// std::cerr << __func__ << " data size: " << data_size_
//		// <<
//		// "\n";

//		for (qt::u32 k{1}; k < K; k++) {
//			qt::u32 col_shift = k * J;
//			qt::u32 xor_shift = k_offset(k) * J;
//			qt::u32 len = k_len(k) * J;

//			// std::cerr << "k " << k << " pool offset " <<
//			// pool_offset
//			//	  << " k offset " << k_offset(k)
//			//	  << " col_shift " << col_shift << "
//			// xor_shift "
//			//	  << xor_shift << " len " << len <<
//			//"\n";

//			// ov_pool.hap_xor_new(pool_offset, len,
//			// col_shift,
//			//		    xor_shift);

//			if (op == comparison_op::bitwise_xor)
//				ov_pool.haps_xor(pool_offset, len, col_shift,
//						 xor_shift, stream);
//			else if (op == comparison_op::sum)
//				ov_pool.hap_sum(pool_offset, len, col_shift,
//						xor_shift, stream);
//		}

//		if (op == comparison_op::bitwise_xor)
//			ov_pool.copy_haps_xor_to_host(data_size_);
//		else if (op == comparison_op::sum)
//			ov_pool.copy_haps_sum_to_host(data_size_);

//		ov_pool.sync_device(stream); // TODO: don't do this in prod

//		if (op == comparison_op::sum)
//			return;

//		xor_ps.assign(xor_data_, xor_data_ + data_size_);

//		meza::cuda_ops::prefix_sum(xor_ps.data(), data_size_);
//	}

//	qt::op_t<std::set<qt::up_t<qt::u32>>> explore_pairs()
//	{
//		std::set<qt::up_t<qt::u32>> matches;
//		std::set<qt::up_t<qt::u32>> mismatches;

//		auto no_inc = [&](qt::u32 start, qt::u32 end) -> bool
//		{
//			if (start == 0)
//				return xor_ps[end] == 0;

//			// start > 0
//			start--;
//			return (xor_ps[end] - xor_ps[start]) == 0;
//		};

//		qt::u32 J = filter.cols();
//		qt::u32 K = H;
//		for (qt::u32 k{1}; k < K; k++) {
//			qt::u32 k_len = this->k_len(k);
//			for (qt::u32 k_off{}; k_off < k_len; k_off++) {
//				auto [h1, h2] = comp_hap_pair(k, k_off);

//				if (filter.base().is_row_blank(h1) ||
//				    filter.base().is_row_blank(h2)) {
//					continue;
//				}

//				qt::u32 start = (k_offset(k) + k_off) * J;
//				qt::u32 end = (start + J) - 1;

//				if (no_inc(start, end))
//					matches.insert({h1, h2});
//				else
//					mismatches.insert({h1, h2});
//			}
//		}

//		return {matches, mismatches};
//	}

//	std::set<qt::up_t<qt::u32>> find_reversals()
//	{
//		std::set<qt::up_t<qt::u32>> reversals;

//		// std::cerr << __func__ << " sum data:\n";
//		// for (qt::u32 i{}; i < data_size_; i++)
//		//	std::cerr << static_cast<u32>(sum_data_[i]) <<
//		//", ";
//		// std::cerr << "\n";

//		std::vector<qt::u32> mask(data_size_, 0);
//		for (qt::u32 i{}; i < data_size_; i++)
//			if (sum_data_[i] == 1 || sum_data_[i] == 2)
//				mask[i] = sum_data_[i];

//		meza::cuda_ops::prefix_sum(mask.data(), data_size_);

//		// Precompute the `presence` array
//		// - tracks whether `1` or `2` exist up to each index
//		std::vector<qt::u32> presence(data_size_, 0);
//		for (qt::u32 i{}; i < data_size_; i++)
//			presence[i] = (i > 0 ? presence[i - 1] : 0) +
//				      (sum_data_[i] == 1 || sum_data_[i] == 2);

//		// Lambda for the range check
//		auto no_inc = [&mask, &presence](qt::u32 start,
//						 qt::u32 end) -> bool
//		{
//			if (start > 0) {
//				// Check if prefix sum is constant in
//				// range and `1`/`2` exists
//				return (mask[end] - mask[start - 1]) == 0 &&
//				       (presence[end] - presence[start - 1]) >
//					       0;
//			}
//			else {
//				// Special case when start is 0
//				return (mask[end] == 0) && presence[end] > 0;
//			}
//		};

//		// std::cerr << __func__ << " sum data summed:\n";
//		// for (qt::u32 i{}; i < data_size_; i++)
//		//	std::cerr << static_cast<u32>(mask[i]) << ", ";
//		// std::cerr << "\n";

//		qt::u32 J = filter.cols();
//		qt::u32 K = H;
//		for (qt::u32 k{1}; k < K; k++) {
//			qt::u32 k_len = this->k_len(k);
//			for (qt::u32 k_off{}; k_off < k_len; k_off++) {
//				auto [h1, h2] = comp_hap_pair(k, k_off);

//				if (filter.base().is_row_blank(h1) ||
//				    filter.base().is_row_blank(h2)) {
//					continue;
//				}

//				qt::u32 start = (k_offset(k) + k_off) * J;
//				qt::u32 end = (start + J) - 1;

//				if (no_inc(start, end))
//					reversals.insert({h1, h2});
//			}
//		}

//		return reversals;
//	}
// };

// struct haps_comp_set {
//	std::set<qt::up_t<qt::u32>> reversals;
//	std::set<qt::up_t<qt::u32>> matches;
//	std::set<qt::up_t<qt::u32>> mismatches;
// };

} // namespace meza::pool::split
#endif // MZ_MATRIX_POOL_SPLIT_HPP
