#ifdef MEZA_USE_CUDA
#ifndef MZ_INFO_CUH
#define MZ_INFO_CUH

#include <iostream>

#include <cuda_runtime.h>

namespace meza::info
{

void print_cuda_device_info();

} // namespace meza::info

#endif // MZ_INFO_CUH
#endif // MEZA_USE_CUDA
