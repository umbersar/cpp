    // computed by the thread// MP 3: Due Sunday, Dec 30, 2012 at 11:59 p.m. PST
#include    <wb.h>

#define wbCheck(stmt) do {                                 \
        cudaError_t err = stmt;                            \
        if (err != cudaSuccess) {                          \
            wbLog(ERROR, "Failed to run stmt ", #stmt);    \
            return -1;                                     \
        }                                                  \
    } while(0)

// Compute C = A * B
__global__ void matrixMultiplyShared(float * A, float * B, float * C,
			             int numARows, int numAColumns,
			             int numBRows, int numBColumns,
			             int numCRows, int numCColumns) {
    //@@ Insert code to implement matrix multiplication here
    //@@ You have to use shared memory for this MP
    const int TILE_WIDTH = 16;
    __shared__ float ds_A[TILE_WIDTH][TILE_WIDTH];
   	__shared__ float ds_B[TILE_WIDTH][TILE_WIDTH];
  
	//int bx = blockIdx.x; int by = blockIdx.y;
  	int tx = threadIdx.x; 
    int ty = threadIdx.y;
  
   	int x = blockIdx.x*blockDim.x+threadIdx.x;//col
    int y = blockIdx.y*blockDim.y+threadIdx.y;//row
 
	//if ((y < numARows) || (x < numBColumns)) {
    // value stores the element that is 

    float Pvalue = 0;
    //for (int i = 0; i < numAColumns; ++i)
    for (int m = 0; m < (numAColumns-1)/TILE_WIDTH + 1 ; ++m)
    {
      if((m*TILE_WIDTH+tx < numARows) && y<numAColumns)
      {
      	ds_A[tx][ty] = A[(m*TILE_WIDTH+tx)*numAColumns + y];
      }
      else
        ds_A[tx][ty] = 0;
      
      if(x<numBRows && m*TILE_WIDTH+ty < numBColumns )
      {
	  	ds_B[tx][ty] = B[x*numBColumns + m*TILE_WIDTH+ty];
      }
      else
        ds_B[tx][ty]=0;
        
	  __syncthreads();
      
      //float elementA = A[tx * numAColumns + k];
      //float elementB = B[k * numBColumns + ty];
      //Pvalue += elementA * elementB;
      for (int k = 0; k < TILE_WIDTH; ++k)
		Pvalue += ds_A[tx][k] * ds_B[k][ty];
      __syncthreads();
    }
 
    // Write the matrix to device memory each 
    // thread writes one element
    if ((x < numCRows) && (y < numCColumns)) 
    C[x * numCColumns + y] = Pvalue;
   //}   
   
}

int main(int argc, char ** argv) {
    wbArg_t args;
    float * hostA; // The A matrix
    float * hostB; // The B matrix
    float * hostC; // The output C matrix
    float * deviceA;
    float * deviceB;
    float * deviceC;
    int numARows; // number of rows in the matrix A
    int numAColumns; // number of columns in the matrix A
    int numBRows; // number of rows in the matrix B
    int numBColumns; // number of columns in the matrix B
    int numCRows; // number of rows in the matrix C (you have to set this)
    int numCColumns; // number of columns in the matrix C (you have to set this)

    args = wbArg_read(argc, argv);

    wbTime_start(Generic, "Importing data and creating memory on host");
    hostA = (float *) wbImport(wbArg_getInputFile(args, 0), &numARows, &numAColumns);
    hostB = (float *) wbImport(wbArg_getInputFile(args, 1), &numBRows, &numBColumns);
    //@@ Set numCRows and numCColumns
    numCRows = numARows;
    numCColumns = numBColumns;
    //@@ Allocate the hostC matrix
    unsigned int size_C = numARows * numBColumns;
    unsigned int mem_size_C = sizeof(float) * size_C;
    hostC = (float*) malloc(mem_size_C);
    wbTime_stop(Generic, "Importing data and creating memory on host");

    wbLog(TRACE, "The dimensions of A are ", numARows, " x ", numAColumns);
    wbLog(TRACE, "The dimensions of B are ", numBRows, " x ", numBColumns);

    wbTime_start(GPU, "Allocating GPU memory.");
    //@@ Allocate GPU memory here
    unsigned int mem_size_A = sizeof(float) * numARows * numAColumns;
    cudaMalloc((void**) &deviceA, mem_size_A);
  
  	unsigned int mem_size_B = sizeof(float) * numBRows * numBColumns;  
  	cudaMalloc((void**) &deviceB, mem_size_B);
        
    cudaMalloc((void**) &deviceC, mem_size_C);
  
    wbTime_stop(GPU, "Allocating GPU memory.");

    wbTime_start(GPU, "Copying input memory to the GPU.");
    //@@ Copy memory to the GPU here
	cudaMemcpy(deviceA, hostA, mem_size_A, 
               cudaMemcpyHostToDevice);
    cudaMemcpy(deviceB, hostB, mem_size_B, 
               cudaMemcpyHostToDevice);
  
    wbTime_stop(GPU, "Copying input memory to the GPU.");
    
    //@@ Initialize the grid and block dimensions here
    dim3 threads(16, 16);
    dim3 grid(ceil(numCRows / (float)threads.x), ceil(numCColumns / (float)threads.y));
  
    wbTime_start(Compute, "Performing CUDA computation");
    //@@ Launch the GPU Kernel here
    matrixMultiplyShared<<< grid, threads >>>(deviceA, deviceB, 
                                   deviceC, numARows, numAColumns, numBRows, numBColumns, numCRows, numCColumns);
      
    cudaThreadSynchronize();
     
    wbTime_stop(Compute, "Performing CUDA computation");
    
    wbTime_start(Copy, "Copying output memory to the CPU");
    //@@ Copy the GPU memory back to the CPU here
    cudaMemcpy(hostC, deviceC, mem_size_C, 
               cudaMemcpyDeviceToHost);

    wbTime_stop(Copy, "Copying output memory to the CPU");

    wbTime_start(GPU, "Freeing GPU Memory");
    //@@ Free the GPU memory here
    cudaFree(deviceA);
    cudaFree(deviceB);
    cudaFree(deviceC);
  
    wbTime_stop(GPU, "Freeing GPU Memory");

    wbSolution(args, hostC, numCRows, numCColumns);

    free(hostA);
    free(hostB);
    free(hostC);

    return 0;
}

