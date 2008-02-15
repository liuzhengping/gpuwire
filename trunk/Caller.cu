#include <stdio.h>
__global__ void
emptyKernel(){
}



extern "C" void 
call(){
    dim3  grid( 1, 1, 1);
    dim3  threads( 256, 1, 1);
	emptyKernel    <<<grid, threads >>> ();
	cudaThreadSynchronize();
	printf("Called\n");
}

