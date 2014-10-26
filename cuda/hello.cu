#include <stdio.h>
/* GPU用strcpy */
__device__ void
dev_strcpy(char *dst, const char *src)
{
  while (*dst++ = *src++);
}
/* GPU側エントリ */
__global__ void gen_hello(char *A)
{
  dev_strcpy(A, "hello");
}

int main()
{
  char *d_hello;

  /* ホストのメモリを確保 */
  char hello[128];

  /* GPU側のメモリ(デバイスメモリ)確保 */
  cudaMalloc((void**)&d_hello, 128);

  /* gen_hello 呼び出し */
  gen_hello<<<1,1>>>(d_hello);

  /* GPU側のデータを取得 */
  cudaMemcpy(hello, d_hello, 128, cudaMemcpyDeviceToHost);

  /* 確保したメモリを解放 */
  cudaFree(d_hello);

  /* 出力 */
  puts(hello);
}
