#include <stdio.h>

/* 1-98番目までのテキスト。__ の部分を今の数字に、 ## の部分をひとつ減らした数字に置き換える */
__device__ char text[] = "__ bottles of beer on the wall, __ bottles of beer!\n"
    "Take one down, and pass it around, ## bottles of beer on the wall!\n\n";

/* 99番目のテキスト。そのまま表示する */
__device__ char end[] = 
    "01 bottle of beer on the wall, 01 bottle of beer.\n"
    "Take one down and pass it around, no more bottles of beer on the wall.\n"
    "\n"
    "No more bottles of beer on the wall, no more bottles of beer.\n"
    "Go to the store and buy some more, 99 bottles of beer on the wall.";

#define SIZE_TEXT (sizeof(text)-1) /* 98番目までのテキストの長さ */
#define SIZE_END (sizeof(end)-1) /* 99番目のテキストの長さ */

__global__
void bottle99(char *addr)
{
    /* スレッド ID(threadIdx.x) の取得 */
    int x = threadIdx.x;

    /* 結果の格納アドレスを求める */
    addr += x * SIZE_TEXT;

    /* 残りのボトル数を求める */
    int bottle = 99 - x;

    if (bottle == 1) {
  /* 99番目, ボトルが無くなったのでendを表示 */
        for (int i=0; i<SIZE_END; i++) {
            addr[i] = end[i];
        }
        addr[SIZE_END] = '\0';
    } else {
  /* 1-98番目 */

  /* 残りのボトル数のASCII表現 */
        char c1 = (bottle/10) + '0';
        char c2 = (bottle%10) + '0';

  /* ひとつ減ったあとのボトル数のASCII表現 */
        char d1 = ((bottle-1)/10) + '0';
        char d2 = ((bottle-1)%10) + '0';

        for (int i=0; i<SIZE_TEXT; i++) {
            int c = text[i];

            if (c == '_') {
    /* '__' の置き換え */
                addr[i] = c1;
                addr[i+1] = c2;
                i++;
            } else if (c == '#') {
    /* '##' の置き換え */
                addr[i] = d1;
                addr[i+1] = d2;
                i++;
            } else {
    /* 文字列のコピー */
                addr[i] = text[i];
            }
        }
    }
}

int main()
{
    char *buffer;
    char *d_buffer;

    /* 必要な領域を計算する */
    int size = SIZE_TEXT * 98 + SIZE_END + 1;

    buffer = new char[size];
    cudaMalloc((void**)&d_buffer, size);

    /* 99 bottles of beer を求めるために99個のスレッドを起動する */
    bottle99<<<1, 99>>>(d_buffer);

    /* 結果を取得 */
    cudaMemcpy(buffer, d_buffer, size, cudaMemcpyDeviceToHost);
    cudaFree(d_buffer);

    /* 表示 */
    puts(buffer);
    free(buffer);
}

