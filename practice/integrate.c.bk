#include <stdio.h>
//シンプソン公式を利用
//観測データを先ず読み込む
//一秒分のポイント数(N)を定義してその分だけの要素数をもった観測強度の配列を作成
#define N 100
double func(double x)
{
  return(1.0/x);
}
//作成された配列の分だけループをまわす
//[["21","231","12", ..Nの数だけ.. , "1224"], ... , ["", ... , ""] ]
int main(void)
{
  int i;
  double a,b,s,s1,s2;
  double x,y[N+1],d;
  //a, bは必要ないね。xも必要ない
  //dに(観測秒/観測ポイント) = (1/N) を入れたら解決？ 
  a=1,0; b=2.0;
  s=0.0; s1=0.0; s2=0.0;
  d=(b-a)/(double)N;
  for(i=0; i<=N; i++)
  {
    x=(double)i*d+a;
    //yに観測データの配列からとってきた観測強度を入れる
    y[i]=func(x);
  }
  for(i=1; i<=N-1;i+=2)
  {
    s1+=y[i];
  }
  for(i=2;i<=N-2;i+=2)
  {
    s2+=y[i];
  }
  s=(y[0]+4.0*s1+2.0*s2+y[N])*d/3.0;

  printf("Intergral : %f\n",s);
}

//作成した1秒積分の配列を一行ずつprintfして、それをグラフ化する
