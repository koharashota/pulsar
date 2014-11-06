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

#define N 50000
#define M 10000
/**
void add( int *a, int *b, int *c ) {
    int tid = 0;    // this is CPU zero, so we start at zero
    while (tid < N) {
        c[tid] = a[tid] * b[tid];
        tid += 1;   // we have one CPU, so we increment by one
    }
}
**/

int main( void ) {
    clock_t start, end;
    float a[N], b[N], c[N];
    int i,j;
    start = clock();
    
    for(j=0; j<M;j++) {
      for (i=0; i<N; i++) {
        a[i] =  (i);
        b[i] =  (i+j);
        c[i] =  (a[i] / b[i] + 120)*a[i]/b[i] + (a[i] / b[i] + 120)*a[i]/b[i] + (a[i] / b[i] + 120)*a[i]/b[i] + (a[i] / b[i] + 120)*a[i]/b[i] + (a[i] / b[i] + 120)*a[i]/b[i] ;
        c[i] +=  (a[i] / b[i] + 120)*a[i]/b[i] + (a[i] / b[i] + 120)*a[i]/b[i] + (a[i] / b[i] + 120)*a[i]/b[i] + (a[i] / b[i] + 120)*a[i]/b[i] + (a[i] / b[i] + 120)*a[i]/b[i] ;
        c[i] +=  (a[i] / b[i] + 120)*a[i]/b[i] + (a[i] / b[i] + 120)*a[i]/b[i] + (a[i] / b[i] + 120)*a[i]/b[i] + (a[i] / b[i] + 120)*a[i]/b[i] + (a[i] / b[i] + 120)*a[i]/b[i] ;
        c[i] = c[i]/N/N;
        //c[i] =  a[i] / b[i];
        if (i==230){
          //printf( "%f / %f = %lf\n", a[i], b[i], c[i] );
        }
      }
    }
    end = clock();
    printf("%.2f秒かかりました\n",(double)(end-start)/CLOCKS_PER_SEC);
    return 0;
}
