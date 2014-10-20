#include <iostream>  // for cout
#include <math.h>    // for sin(), cos()
#include <stdio.h>   // for printf()

#define N 100                // 分割数
#define CSV_DFT  "DFT.csv"   // 出力ファイル (DFT)

using namespace std;

/*
 *  * 計算クラス
 *   */
class Calc
{
  double SRC_re[N];   //元データの実部
  double SRC_im[N];   //元データの虚部
  double DFT_re[N];   //DFTの実部
  double DFT_im[N];   //DFTの虚部
  double Power[N];  //PowerSpectun
  double Phase_re[N];  //Phaser
  double Phase_im[N];  //Phaser

  public:
  void makeSourceData();  // 元データ作成
  void executeDFT();      // 離散フーリエ変換

  private:
  double calcTerm(int n, double x);  //各項計算
};

/*
 *  * 元データ作成
 *   */
void Calc::makeSourceData()
{
  int i;

  for (i = 0; i < N; i++) {
    SRC_re[i] = 2 * sin(4 * (2 * M_PI / N) * i)
      + 3 * cos(2 * (2 * M_PI / N) * i);
    SRC_im[i] = 0.0;
  }
}

void Calc::executeDFT()
{
    int k, n;  // LOOPインデックス
    FILE *pf;  // ファイルポインタ

    // 出力ファイルOPEN
    pf = fopen(CSV_DFT, "w");

    // ヘッダ出力 ( k,PowerSpectum,Phaser_re,Phaser_im,  DFT(実部), DFT(虚部) )
    fprintf(pf, "k,PowerSpectum,Phaser_re,Phaser_im,X_re,X_im\n");

    // 計算・結果出力
    for (k = 0; k < N; k++) {
        DFT_re[k] = 0.0;
        DFT_im[k] = 0.0;
        for (n = 0; n < N; n++) {
            DFT_re[k] += SRC_re[n] * ( cos((2 * M_PI / N) * k * n))
                       + SRC_im[n] * ( sin((2 * M_PI / N) * k * n));
            DFT_im[k] += SRC_re[n] * (-sin((2 * M_PI / N) * k * n))
                       + SRC_im[n] * ( cos((2 * M_PI / N) * k * n));
        }
        Power[k] = DFT_re[k]*DFT_re[k]+DFT_im[k]*DFT_im[k];
        Phase_re[k] = sqrt(Power[k])*(cos(atan(DFT_re[k]/DFT_re[k])));
        Phase_im[k] = sqrt(Power[k])*(sin(atan(DFT_re[k]/DFT_re[k])));
        
        fprintf(pf, "%d,%lf,%lf,%lf,%lf,%lf\n",
            k, Power[k] ,Phase_re[k],Phase_im[k] ,DFT_re[k], DFT_im[k]);
    }

    // 出力ファイルCLOSE
    fclose(pf);
}
/*
 * メイン処理
 */
int main()
{
    try
    {
        // 計算クラスインスタンス化
        Calc objCalc;
        // 元データ作成
        objCalc.makeSourceData();
        // 離散フーリエ変換
        objCalc.executeDFT();
    }
    catch (...) {
        cout << "例外発生！" << endl;
        return -1;
    }

    // 正常終了
    return 0;
}
