#include <stdio.h>
#include <math.h>

int main(void)
{
  int i, j, ret;
  double tmp;
  FILE *fp;
  FILE *fp1;
  char *inputfname = "data.csv";
  char *outputfname = "integrate.csv";
  
  //合計値をいれる
  double s;
  
  //何秒積分か入力
  double sec;
  //1秒あたりの観測ポイント数を定義(1ポイント毎の秒数としても利用)
  double n;
  //観測合計時間
  double len;
  printf("何秒積分したいですか？(double): ");
  scanf("%lf",&sec);
  printf("毎秒何ポイント観測してますか？(double): ");
  scanf("%lf",&n);
  printf("何秒分の観測データですか？(double): ");
  scanf("%lf",&len);
  
  //積分の配列に入れるポイント数を定義
  int p = ceil(sec*n);
  //合計ポイント数を定義
  int N = ceil(len*n);
  //ループ処理回数
  int loop = floor(N/p);
  //配列を定義。p個分観測データを入れる
  double y[p+1];
  double data[N+1];

  //読み込みファイルを開く
  fp1 = fopen( inputfname, "r" );
  if( fp1 == NULL ){
    printf( "%sファイルが開けません\n", inputfname );
    return -1;
  }
  for(i=0;i<N;i++) {
    fscanf(fp1, "%lf", &data[i]);
  }
  //読み込み終了
  fclose( fp1 );

  //書き込みファイルを開く
  fp = fopen( outputfname, "w" );
  if( fp == NULL ){
    printf( "%sファイルが開けません\n", outputfname );
    return -1;
  }
  
  for(j=0; j < loop; j++){
    //初期化
    s=0.0;
    
    //printf("%d\n",p*j);

    //この辺はシンプソン公式通り
    for(i=0; i < p; i++)
    {
      s+=data[i+p*j];
    }
    printf("Intergral : %lf(%lfs)\n", s/n, (p)*(j+1)/n);
    
    //csvファイルに書き込み
    fprintf( fp, "%lf\n", s/n);
  }
  //編集終了
  fclose( fp );
  printf( "%sファイル書き込みが終わりました\n", outputfname );
  return 0;
}

