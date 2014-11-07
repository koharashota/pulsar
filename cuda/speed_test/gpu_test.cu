/*
 * Copyright 1993-2010 NVIDIA Corporation.  All rights reserved.
 *
 * NVIDIA Corporation and its licensors retain all intellectual property and 
 * proprietary rights in and to this software and related documentation. 
 * Any use, reproduction, disclosure, or distribution of this software 
 * and related documentation without an express license agreement from
 * NVIDIA Corporation is strictly prohibited.
 *
 * Please refer to the applicable NVIDIA end user license agreement (EULA) 
 * associated with this source code for terms and conditions that govern 
 * your use of this NVIDIA software.
 * 
 */


#include "../cuda_by_example/common/book.h"
#include <time.h>
#define imin(a,b) (a<b?a:b)

#define N 50000
#define M 10000
const int threadsPerBlock = 500;
const int blocksPerGrid =
            imin( 60000, (N+threadsPerBlock-1) / threadsPerBlock );

__global__ void add( float *a, float *b, float *c ) {
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    if (tid < N)
      c[tid] = (a[tid] / b[tid] + 120)*a[tid]/b[tid] + (a[tid] / b[tid] + 120)*a[tid]/b[tid] + (a[tid] / b[tid] + 120)*a[tid]/b[tid] +  (a[tid] / b[tid] + 120)*a[tid]/b[tid] + (a[tid] / b[tid] + 120)*a[tid]/b[tid];
      c[tid] += (a[tid] / b[tid] + 120)*a[tid]/b[tid] + (a[tid] / b[tid] + 120)*a[tid]/b[tid] + (a[tid] / b[tid] + 120)*a[tid]/b[tid] +  (a[tid] / b[tid] + 120)*a[tid]/b[tid] + (a[tid] / b[tid] + 120)*a[tid]/b[tid];
      c[tid] += (a[tid] / b[tid] + 120)*a[tid]/b[tid] + (a[tid] / b[tid] + 120)*a[tid]/b[tid] + (a[tid] / b[tid] + 120)*a[tid]/b[tid] +  (a[tid] / b[tid] + 120)*a[tid]/b[tid] + (a[tid] / b[tid] + 120)*a[tid]/b[tid];
      c[tid] = c[tid]/N/N;
      tid += blockDim.x * gridDim.x;
}

int main( void ) {
    clock_t start, end;
    float a[N], b[N], c[N];
    int i,j;
    float *dev_a, *dev_b, *dev_c;
    start = clock();

    for(j=0; j<M;j++) {
      // allocate the memory on the GPU
      HANDLE_ERROR( cudaMalloc( (void**)&dev_a, N * sizeof(float)) );
      HANDLE_ERROR( cudaMalloc( (void**)&dev_b, N * sizeof(float)) );
      HANDLE_ERROR( cudaMalloc( (void**)&dev_c, N * sizeof(float)) );
      
      // fill the arrays 'a' and 'b' on the CPU
      for (i=0; i<N; i++) {
        a[i] =  (i);
        b[i] =  (i+j);
        if (j==230){
          //printf( "%f / %f = \n", a[i], b[i] );
          //printf( "%f / %f = %f\n", a[i], b[i], c[i] );
        }
      }
      
      // copy the arrays 'a' and 'b' to the GPU
      HANDLE_ERROR( cudaMemcpy( dev_a, a, N * sizeof(float),
                                cudaMemcpyHostToDevice ) );
      HANDLE_ERROR( cudaMemcpy( dev_b, b, N * sizeof(float),
                                cudaMemcpyHostToDevice ) );

      //add<<<N,1>>>( dev_a, dev_b, dev_c );
      add<<<blocksPerGrid,threadsPerBlock>>>( dev_a, dev_b, dev_c );

      // copy the array 'c' back from the GPU to the CPU
      HANDLE_ERROR( cudaMemcpy( c, dev_c,   N * sizeof(float),
                                cudaMemcpyDeviceToHost ) );
      //if (j==230){
        for (int i=0; i<N; i++) {
          //printf( "%f + %f = %f\n", a[i], b[i], c[i] );
          //printf( "%f \n", c[i] );
        }
      //}
      // free the memory allocated on the GPU
      HANDLE_ERROR( cudaFree( dev_a ) );
      HANDLE_ERROR( cudaFree( dev_b ) );
      HANDLE_ERROR( cudaFree( dev_c ) );
    }
    
    end = clock();
    printf("%.2f秒かかりました\n",(double)(end-start)/CLOCKS_PER_SEC);
    return 0;
}
