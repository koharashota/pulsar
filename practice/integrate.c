//シンプソン公式利用
//http://www.opt.utsunomiya-u.ac.jp/~yatagai/Lectures/jouhou/Chap6.htm

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
  
  //シンプソン公式のためのsを定義
  double s,s1,s2;
  
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
  /**
  while(( ret = fscanf( fp1, "%lf", &tmp )) != EOF ){
    printf("%lf\n",tmp);
    printf("%d\n",i);
    //data[i] = tmp;
    i = i + 1;
  }
  **/
  //読み込み終了
  fclose( fp1 );

  //書き込みファイルを開く
  fp = fopen( outputfname, "w" );
  if( fp == NULL ){
    printf( "%sファイルが開けません\n", outputfname );
    return -1;
  }
  
  for(j=0; j<=loop; j++){
    //初期化
    s=0.0; s1=0.0; s2=0.0;
    
    //この辺はシンプソン公式通り
    for(i=0; i<=p; i++)
    {
      y[i]=data[i+(p+1)*j];
    }
    for(i=1; i<=p-1;i+=2)
    {
      s1+=y[i];
    }
    for(i=2;i<=p-2;i+=2)
    {
      s2+=y[i];
    }
    s=(y[0]+4.0*s1+2.0*s2+y[p])/n/3.0;
    
    printf("Intergral : %lf(%lfs)\n", s, (p+1)*j/n);
    
    //csvファイルに書き込み
    fprintf( fp, "%lf\n", s);
  }
  //編集終了
  fclose( fp );
  printf( "%sファイル書き込みが終わりました\n", outputfname );
  return 0;
}

