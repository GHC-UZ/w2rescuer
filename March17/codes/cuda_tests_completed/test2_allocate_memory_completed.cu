#include <iostream>
#include <time.h>
#include <cuda_runtime.h>

// Kernel functions for memory allocation
__global__ void checkMemoryAllocate(int N, size_t size, double *d_A) {

    int ithread = threadIdx.x;
    int iblock = blockIdx.x; 
    int idx = ( blockIdx.x * blockDim.x ) + threadIdx.x;

    if (idx < N) {

        printf("<<gpu>> Thread %d in block %d - d_A[%d] = %.0lf \n", ithread, iblock, idx, d_A[idx]);

    }
       
}


__global__ void registerMemoryAllocate(int N, size_t size, double *d_A) {

    int ithread = threadIdx.x;
    int iblock = blockIdx.x; 
    int idx = ( blockIdx.x * blockDim.x ) + threadIdx.x;

    //Declare registers
    double rA;

    if (idx < N) {

        // Use the register memory
        rA = (double)(idx);

        // Write the result back to global memory
        d_A[idx] = rA;

        printf("<<gpu>> Thread %d in block %d - d_A[%d] = %.0lf \n", ithread, iblock, idx, d_A[idx]);

    }
       
}


__global__ void sharedMemoryAllocate(int N, size_t size, double *d_A) {

    int ithread = threadIdx.x;
    int iblock = blockIdx.x; 
    int idx = ( blockIdx.x * blockDim.x ) + threadIdx.x;

    // Declare dynamic shared memory
    extern __shared__ double shared_A[];

    // Allocate shared memory
    if (idx < N) {
        shared_A[ithread] = (double)(idx);
    }

    // Synchronize to ensure all threads have written to shared memory
    __syncthreads();        

    if (idx < N) {
        // Write the result back to global memory
        d_A[idx] = shared_A[ithread];

        printf("<<gpu>> Thread %d in block %d - d_A[%d] = %.0lf \n", ithread, iblock, idx, d_A[idx]);

    }
    
}
                                       


int main() {

    cudaError_t err; // CUDA error
    double exec_time=0.0; //timer 
    clock_t stime1, stime2;

    // Size of the vector
    int N = 1000;

    // Allocate memory on the host
    size_t size = N*sizeof(double);
    double *h_A = (double*) malloc(size);

    // Initialize host vectors
    for (int i = 0; i < N; i++) {
        h_A[i] = (double)i;
    }



    // Select a GPU
    int selectedGPU = 0; // Change this to select your asigned GPU
    cudaSetDevice(selectedGPU);  

    // Get the device propoerties
    cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop, selectedGPU);    

    // Get total and free memory
    size_t freeMemory, totalMemory;
    cudaMemGetInfo(&freeMemory, &totalMemory);
    std::cout << "  Total Global Memory: " << totalMemory / (1024 * 1024) << " MB" << std::endl;
    std::cout << "  Free Memory: " << freeMemory / (1024 * 1024) << " MB" << std::endl;
           
    getchar(); //Pause point




    // Allocate memory on the device
    double *d_A;
    err = cudaMalloc((void**) &d_A, size);

    // Check if the allocation was successful
    if (err != cudaSuccess) {
        std::cerr << "cudaMalloc failed: " << cudaGetErrorString(err) << std::endl;
        return 1;
    }    

    // Check the updated free memory
    cudaMemGetInfo(&freeMemory, &totalMemory);
    std::cout << "  Updated free Memory: " << freeMemory / (1024 * 1024) << " MB" << std::endl;

    getchar(); //Pause point





    //Start IO time .....................................
    stime1=clock();

    // Launch free-load kernel
    // int threadsPerBlock = 256;
    // int blocksPerGrid = N/threadsPerBlock+1;

    // checkMemoryAllocate <<<blocksPerGrid, threadsPerBlock>>> (N, size, d_A);




    //Launch register allocation kernel
    // int threadsPerBlock = 256;
    // int blocksPerGrid = N/threadsPerBlock+1;

    // registerMemoryAllocate <<<blocksPerGrid, threadsPerBlock>>> (N, size, d_A);



    //Launch shared allocation kernel
    int threadsPerBlock = 256;
    int blocksPerGrid = N/threadsPerBlock+1;    
    size_t requiredSharedMemory = threadsPerBlock*sizeof(double);
    std::cout << "  Maximum shared memory per block: " << prop.sharedMemPerBlock << " bytes" << std::endl; 
    std::cout << "  Required shared memory per block: " << requiredSharedMemory << " bytes" << std::endl;
    
    sharedMemoryAllocate <<<blocksPerGrid, threadsPerBlock, requiredSharedMemory>>> (N, size, d_A);



    // Synchronize device
    err = cudaDeviceSynchronize();
    if (err != cudaSuccess) {
        std::cerr << "Error en cudaDeviceSynchronize: " << cudaGetErrorString(err) << std::endl;
    }

    stime2=clock();
    exec_time += double(stime2-stime1)/CLOCKS_PER_SEC;    
    std::cerr << "  Execution time: " << exec_time << std::endl; 
    //End IO time .....................................   

    getchar(); //Pause point





    // Free device memory
    cudaFree(d_A);

    // Check the updated free memory
    cudaDeviceSynchronize();
    cudaMemGetInfo(&freeMemory, &totalMemory);
    std::cout << "  Final free Memory: " << freeMemory / (1024 * 1024) << " MB" << std::endl;

    // Free host memory
    free(h_A);

    return 0;
}

