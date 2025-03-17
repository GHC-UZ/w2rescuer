#include <iostream>
#include <cuda_runtime.h>

__global__ void hello() {
    printf("\n************");
    printf("\nHi from GPU!");
    printf("\n************\n\n");
}

int main() {

    // Step 1: Get the number of GPUs available
    int numGPUs = 0;
    cudaGetDeviceCount(&numGPUs);
    std::cout << "Number of GPUs available: " << numGPUs << std::endl;

    if (numGPUs == 0) {
        std::cerr << "No GPUs found!" << std::endl;
        return -1;
    }



    // Step 2: List properties of each GPU
    for (int i = 0; i < numGPUs; i++) {
        cudaDeviceProp props;
        cudaGetDeviceProperties(&props, i);

        std::cout << "GPU " << i << ": " << props.name << std::endl;
    } 



    // Step 3: Select a GPU (e.g., GPU 0)
    int selectedGPU = 0;
    cudaSetDevice(selectedGPU);  




    // Step 4: Verify the selected GPU
    cudaDeviceProp selectedProps;
    cudaGetDeviceProperties(&selectedProps, selectedGPU);

    std::cout << "\nSelected GPU: " << selectedProps.name << std::endl; 
    std::cout << "  Compute Capability: " << selectedProps.major << "." << selectedProps.minor << std::endl; 
    std::cout << "  Multiprocessors (SM): " << selectedProps.multiProcessorCount << std::endl;
    std::cout << "  Maximun threads per block: " << selectedProps.maxThreadsPerBlock << std::endl;
    std::cout << "  Maximun number of blocks (in x): " << selectedProps.maxGridSize[0] << std::endl;
    std::cout << "  Maximun threads in each SM: " << selectedProps.maxThreadsPerMultiProcessor << std::endl;
    



    // Step 5: Get total and free memory
    size_t freeMemory, totalMemory;
    cudaMemGetInfo(&freeMemory, &totalMemory);

    std::cout << "  Total Global Memory: " << totalMemory / (1024 * 1024) << " MB" << std::endl;
    std::cout << "  Free Memory: " << freeMemory / (1024 * 1024) << " MB" << std::endl;



    // Launch kernel
    hello <<<1,1>>> ();

    // Check for launch errors
    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess) {
        std::cerr << "Error launching kernel: " << cudaGetErrorString(err) << std::endl;
        return -1;
    }




    // Synchronize and check for execution errors
    cudaError_t sync_err = cudaDeviceSynchronize();
    if (sync_err != cudaSuccess) {
        std::cerr << "Error synchronizing: " << cudaGetErrorString(sync_err) << std::endl;
        return -1;
    }




    // Flush stdout to ensure printf output is displayed
    fflush(stdout);

    return 0;
}
