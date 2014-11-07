#include "../cuda_by_example/common/book.h"
#define imin(a,b) (a<b?a:b)

//総データ数
const int N = 100 * 1024;
//一分割(スペクトルのグラフ)毎のデータ数(=総データ数/時間毎の分割数)
const int M =  100;

//一度に使うスレッド数は512を超えないようにする
const int threadsPerBlock = 500;
//const int threadsPerBlock = 256;
// 一度につかうブロック数は6,5535を超えないようにする
const int blocksPerGrid =
            imin( 60000, (N+threadsPerBlock-1) / threadsPerBlock );


__global__ void PowerSpectral( float *t, float *f, float *X, float *S ) {
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    
    //ここにそれぞれのパラメーターに対しての処理を書く
    while (tid < N) {
      S[tid] = X[tid] / 5;
      tid += blockDim.x * gridDim.x;
    }

    
    //同期的にスレッド間のやり取りが必要な場合、共有メモリを使う。非同期なら不要
    /** 
    __shared__ float cache[threadsPerBlock];
    int cacheIndex = threadIdx.x;
    
    float   temp = 0;
    while (tid < N) {
        temp += a[tid] * b[tid];
        tid += blockDim.x * gridDim.x;
    }
    
    // set the cache values
    cache[cacheIndex] = temp;
    
    // synchronize threads in this block
    __syncthreads();

    // for reductions, threadsPerBlock must be a power of 2
    // because of the following code
    int i = blockDim.x/2;
    while (i != 0) {
        if (cacheIndex < i)
            cache[cacheIndex] += cache[cacheIndex + i];
        __syncthreads();
        i /= 2;
    }

    if (cacheIndex == 0)
        c[blockIdx.x] = cache[0];
    **/
}


int main( void ) {
    
    //t 時間, f 周波数, X 強度, S パワースペクトラム
    float   *t, *f, *X, *S;
    float   *dev_t, *dev_f, *dev_X, *dev_S;

    // allocate memory on the cpu side
    t = (float*)malloc( N*sizeof(float) );
    f = (float*)malloc( N*sizeof(float) );
    X = (float*)malloc( N*sizeof(float) );
    S = (float*)malloc( N*sizeof(float) );

    // allocate the memory on the GPU
    HANDLE_ERROR( cudaMalloc( (void**)&dev_t,
                              N*sizeof(float) ) );
    HANDLE_ERROR( cudaMalloc( (void**)&dev_f,
                              N*sizeof(float) ) );
    HANDLE_ERROR( cudaMalloc( (void**)&dev_X,
                              N*sizeof(float) ) );
    HANDLE_ERROR( cudaMalloc( (void**)&dev_S,
                              blocksPerGrid*sizeof(float) ) );

    // fill in the host memory with data
    for (int j=0; j<N/M; j++) {
      for (int i=0; i<M; i++) {
        //0.2秒毎のデータと仮定 
        t[i + M*j] = 0.2*j;
        //1400MHzから2MHz毎にデータを取ってると仮定
        f[i + M*j] = 2*i + 1400;
        //強度を適当に仮定
        X[i + M*j] = i*2 + 100;
      }
    }

    // copy the arrays 'a' and 'b' to the GPU
    HANDLE_ERROR( cudaMemcpy( dev_t, t, N*sizeof(float),
                              cudaMemcpyHostToDevice ) );
    HANDLE_ERROR( cudaMemcpy( dev_f, f, N*sizeof(float),
                              cudaMemcpyHostToDevice ) ); 
    HANDLE_ERROR( cudaMemcpy( dev_X, X, N*sizeof(float),
                              cudaMemcpyHostToDevice ) ); 

    PowerSpectral<<<blocksPerGrid,threadsPerBlock>>>( dev_t, dev_f, dev_X,
                                                      dev_S );

    // copy the array 'c' back from the GPU to the CPU
    HANDLE_ERROR( cudaMemcpy( S, dev_S,
                              N*sizeof(float),
                              cudaMemcpyDeviceToHost ) );

    // finish up on the CPU side
    //for (int i=0; i<blocksPerGrid; i++) {
    for (int i=0; i<N; i++) {
      printf( "時間%f, 強度%f, パワースペクトラム%f, 周波数%f\n", t[i], X[i], S[i], f[i]);
    }


    // free memory on the gpu side
    HANDLE_ERROR( cudaFree( dev_t ) );
    HANDLE_ERROR( cudaFree( dev_f ) );
    HANDLE_ERROR( cudaFree( dev_X ) );
    HANDLE_ERROR( cudaFree( dev_S ) );

    // free memory on the cpu side
    free( t );
    free( f );
    free( X );
    free( S );
}
