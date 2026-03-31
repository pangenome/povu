#include <iostream>

#include <cuda_runtime.h>

namespace meza::info
{
void print_cuda_device_info()
{
	int count = 0;
	cudaError_t err = cudaGetDeviceCount(&count);
	if (err != cudaSuccess) {
		std::cerr << "cudaGetDeviceCount failed: "
			  << cudaGetErrorString(err) << "\n";
		return;
	}

	std::cerr << "CUDA devices found: " << count << "\n";

	int current = 0;
	err = cudaGetDevice(&current);
	if (err == cudaSuccess) {
		std::cerr << "Current device: " << current << "\n";
	}

	for (int d = 0; d < count; ++d) {
		cudaDeviceProp p{};
		err = cudaGetDeviceProperties(&p, d);
		if (err != cudaSuccess) {
			std::cerr << "cudaGetDeviceProperties(" << d
				  << ") failed: " << cudaGetErrorString(err)
				  << "\n";
			continue;
		}

		std::cerr << "\n=== Device " << d << " ===\n";
		std::cerr << "Name: " << p.name << "\n";
		std::cerr << "Compute capability: " << p.major << "." << p.minor
			  << "\n";
		std::cerr << "Total global memory: "
			  << (p.totalGlobalMem / (1024.0 * 1024.0)) << " MiB\n";
		std::cerr << "SMs (multiProcessorCount): "
			  << p.multiProcessorCount << "\n";
		std::cerr << "Max threads per block: " << p.maxThreadsPerBlock
			  << "\n";
		std::cerr << "Max threads dim: (" << p.maxThreadsDim[0] << ", "
			  << p.maxThreadsDim[1] << ", " << p.maxThreadsDim[2]
			  << ")\n";
		std::cerr << "Max grid size: (" << p.maxGridSize[0] << ", "
			  << p.maxGridSize[1] << ", " << p.maxGridSize[2]
			  << ")\n";
		std::cerr << "Warp size: " << p.warpSize << "\n";
		std::cerr << "Shared mem per block: "
			  << (p.sharedMemPerBlock / 1024.0) << " KiB\n";
		std::cerr << "Regs per block: " << p.regsPerBlock << "\n";

		// std::cerr << "Memory clock rate: "
		//	  << (p.memoryClockRate / 1000.0)
		//	  << " MHz\n"; // kHz->MHz
		std::cerr << "Memory bus width: " << p.memoryBusWidth
			  << " bits\n";

		std::cerr << "Async engines (copy engines): "
			  << p.asyncEngineCount << "\n";
		// std::cerr << "Device overlap supported: "
		//	  << (p.deviceOverlap ? "yes" : "no") << "\n";

		std::cerr << "Unified addressing: "
			  << (p.unifiedAddressing ? "yes" : "no") << "\n";
		std::cerr << "Managed memory: "
			  << (p.managedMemory ? "yes" : "no") << "\n";
		std::cerr << "Concurrent managed access: "
			  << (p.concurrentManagedAccess ? "yes" : "no") << "\n";

		std::cerr << "Can map host memory: "
			  << (p.canMapHostMemory ? "yes" : "no") << "\n";
		std::cerr << "ECC enabled: " << (p.ECCEnabled ? "yes" : "no")
			  << "\n";
	}
}
} // namespace meza::info
