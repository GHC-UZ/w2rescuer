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
    //******* complete code here *******/

    //**********************************/

    if (numGPUs == 0) {
        std::cerr << "No GPUs found!" << std::endl;
        return -1;
    }



    // Step 2: List properties of each GPU
    for (int i = 0; i < numGPUs; i++) {
        //******* complete code here *******/

        //**********************************/
    } 



    // Step 3: Select a GPU (e.g., GPU 0)
    //******* complete code here *******/

    //**********************************/




    // Step 4: Verify the selected GPU
    //******* complete code here *******/

    //**********************************/
    



    // Step 5: Get total and free memory
    size_t freeMemory, totalMemory;
    //******* complete code here *******/

    //**********************************/




    // Launch kernel
    hello <<<1,1>>> ();

    // Check for launch errors
    //******* complete code here *******/

    //**********************************/



    // Synchronize and check for execution errors
    //******* complete code here *******/

    //**********************************/



    // Flush stdout to ensure printf output is displayed
    fflush(stdout);

    return 0;
}
