
// 1データファイル分の積分データを配列に格納し、main 最後あたりでファイル出力
// 不要な文は全て消す
// ファイル末尾に中途半端に作った積分プログラムがあるので適当に改良せよ
// void time_series_integration_complex (const float integration_time_s, ta::FileInformation &fileinfo, ta::BinData &data, ta::Config &config, const path dest_dir)

//argv[0] 実行ファイル名 argv[1] 観測データ argv[2] 設定ファイル config.txt?

#include <cstdio>
#include <iostream>
#include <fstream>
#include <string>
#include <ctime>
#include "Functions.h"
//c++のライブラリを使ってディレクトリとかファイルパスの操作 fs
#include "FileInformation.h"
//データの項とかの説明
#include "Config.h"
#include "BinData.h"
#include "Stats.h"

#include <cuda.h>
#include <cufft.h>
#include <cuda_runtime.h>

//c++のファイルやディレクトリを操作ライブラリ
#include <boost/filesystem.hpp>

#pragma comment(lib, "cudart")
#pragma comment(lib, "cufft")

#define STR_LEN_MAX 2048
#define CL

#ifdef CL
//#include <Windows.h>
#define PLOT "pgnuplot"
#elif defined(GCC)
#include <unistd.h>
#define PLOT "gnuplot"
#endif

using std::cout;
using std::endl;
using std::string;
using boost::filesystem::path;


void time_series_integration_complex (int integration_point, ta::BinData &data, int num_entire_data_pt, int datum_bit_length_B, ta::FileInformation &fileinfo, ta::Config &config, const path dest_dir);

int main(int argc, char *argv[])
{
    
    bool flag_create_spectrogram_image = false;
    try{
        // Input File の引数に誤りがあったらエラーをはく
        if(argc != 3){
            printf("Usage: %s <file> <confiig>¥n", argv[0]);
            throw "Invalid Arguments";
        }
        
        //fs ファイル情報宣言
        ta::FileInformation fileinfo (argv[1], '_');
        printf ("Data:   %s¥n", fileinfo.getFilePath().string().c_str());
        printf ("Config: %s¥n", path(argv[2]).string().c_str());
        
        // Observation Parameters
        ta::Config config (argv[2]);
        const bool   flag_complex_data   = true; //config.isComplexData(); // Data format: <Re><Re><Re>... or <Re><Im><Re><Im>...
        const int    datum_bit_length_B  = config.getDatumBitLength_B(); // byte
        const double sampling_freq_Hz    = config.getSamplingFrequency_Hz(); // Hz, sample/sec
        const double sampling_interval_s = config.getSamplingInterval_s(); // sec
        
        const float  pulsar_p0_s         = config.getPulsarP0_s();
        const float  pulsar_w50_s        = config.getPulsarW50_s();
        
        
        int tmp_num_entire_data_pt = fileinfo.getFilesize_B() / datum_bit_length_B;
        if (flag_complex_data == true){
            tmp_num_entire_data_pt /= 2;
        }
        
        const int num_entire_data_pt = tmp_num_entire_data_pt;
        
        // Display Parameters
        printf ("¥nPulsar Information¥n");
        printf ("- Period      P0      = %f ms¥n", pulsar_p0_s  * 1000);
        printf ("- Pulse Width W50     = %f ms¥n", pulsar_w50_s * 1000);
        printf ("¥nData Information¥n");
        printf ("- Sampling frequency  = %f MHz¥n", sampling_freq_Hz / 1E+6);
        printf ("- Sampling interval   = %f ns¥n",  sampling_interval_s * 1E+9);
        
        printf ("¥nAnalysis Information¥n");
        printf ("- Entire  Data Size   = %f ms = %d pt¥n",                  sampling_interval_s * num_entire_data_pt * 1000, num_entire_data_pt);
        //printf ("- Number of Analysis Partitions = %d¥n", num_partition);
        
        printf ("¥nInitiate the process?");
        if (!ta::stdin_yes_or_no()) {
            throw "Task Terminated";
        }
        printf ("¥nAnalyzing...¥n");
        
        
        // Data
        ta::BinData data;
        if (flag_complex_data == true) {
            clock_t t0 = clock();
            data.load_binary_as_double (fileinfo.getFilePath().string(), num_entire_data_pt, datum_bit_length_B);
            clock_t t1 = clock();
            ta::print_elapse ("Data Loading", t0, t1);
            
        } else {
            throw "not complex data";
        }
        
        // Create Directory
        //config.txtに設定されたパスを設定
        const path dest_dir = path (config.getOutputDirectory()) / fileinfo.getFileStemSubstr(); // = /opt/data/summation
        //c++のライブラリを使ってパスの作成
        boost::filesystem::create_directories (dest_dir);
        boost::filesystem::current_path (dest_dir);
        
        
        //1s積分の読み出し
        time_series_integration_complex (20000000, data, num_entire_data_pt, datum_bit_length_B, fileinfo, config, dest_dir);
        
        // Delete
        return 0;
    }
    
    catch (const char *err) {
        fprintf(stderr, "%s¥n", err);
        system ("pause");
        return -1;
    }
    catch (const string err) {
        fprintf(stderr, "%s¥n", err.c_str());
        system ("pause");
        return -1;
    }
}


/**
 ・・くみこめばうごくかなー
 integration_point 何ポイントの合計値を取るか(一秒積分なら1/50^[-9]なので20000000ポイント)
 data バイナリから読み取ったデータ配列をもらう。スペクトルの大きさの配列を受け取りたい(ta::BinData data.load_binary_as_doubleでよんだデータ)
 num_entire_data_pt データの全ポイント数
 tmp_num_partial_data 分割データのポイント数(いらないかも、てかいらなそう)
 fileinfo fileのパス情報とか受け取る
 config configの読み取り
 dest_dir 保存先のディレクトリかな？
 complex_flagはオンと仮定
 **/

void time_series_integration_complex (int integration_point, ta::BinData &data, int num_entire_data_pt, int datum_bit_length_B, ta::FileInformation &fileinfo, ta::Config &config, const path dest_dir)
{
    
    int tmp_num_partial_data = 100000000; //適当にセット(必要ない。寧ろintegrate_window_width軽くした方がよさげ)
    const int   integrate_window_width  = integration_point;
    const int   integrate_batch = std::floor(static_cast<float>(tmp_num_partial_data) / integrate_window_width);
    const int   num_partial_data       = integrate_window_width * integrate_batch;
    const int   num_partition          = std::floor( static_cast<float>(num_entire_data_pt) / num_partial_data);
    //これ入れないと最後のデータ一部とばしちゃう
    //const int   skipped_data_size_pts  = num_entire_data_pt - num_partial_data * num_partition;
    
    //1秒積分のために用意
    float *integrate_re = new float[integrate_batch * num_partition];
    float *integrate_im = new float[integrate_batch * num_partition];
    for(int f=0; f<integrate_batch * num_partition; f++){
        integrate_re[f] = 0;
        integrate_im[f] = 0;
    }
    
    //float *h_idata   = new float[num_partial_data];
    cufftComplex *h_idata   = new cufftComplex[num_partial_data];
    for(int i=0; i<num_partial_data; i++){
        h_idata[i].x = 0.0;
        h_idata[i].y = 0.0;
    }
    
    //
    // MAIN
    //
    
    for (int pos = 0; pos < num_partition; pos++) {
        const int partial_data_begin_pt = pos * num_partial_data;
        //hi_dataの作成        
        data.load_binary_as_double (fileinfo.getFilePath().string(), num_entire_data_pt, datum_bit_length_B);
        // ここの処理でreとimのそれぞれの強度をポイント毎に入れられれば良い
        
        // Confirm Data For Debug
        if (pos == 0) {
            std::ofstream fout_test ( (dest_dir / (fileinfo.getFileStem().string() + ".txt")).string() );
            for (int i = 0; i < 100; i++) {
                fout_test << h_idata[i].x << "¥t" << h_idata[i].y << "¥n";
            }
            fout_test.close();
        }
        
        //一秒積分
        for(int j=0; j < integrate_batch ; j++){
            //初期化
            double re=0.0;
            double im=0.0;
            int integrate_first_point = j * integrate_window_width + partial_data_begin_pt;
            
            //for(int i=0; i < integrate_window_width; i++)
            for(int i=0; i < integrate_window_width; i++)
            {
                re+=h_idata[i+integrate_first_point].x;
                im+=h_idata[i+integrate_first_point].y;
            }
            integrate_re[j + pos * integrate_batch] = re;
            integrate_im[j + pos * integrate_batch] = im;
            //デバッグ
            printf("Intergral re : %lf, im : %lf¥n", re , im);
        }
        
    } // Next pos
    
    // Save the Integrate as a text file
    char str[STR_LEN_MAX + 1];
    
    sprintf (str, "%s.integrate.im", fileinfo.getFileStem().string().c_str());
    const path output_im_filepath   = dest_dir / str; // Specify the destination directory and avoid depending on the current directory for Gnuplot configuration.
    FILE *ofp_im = fopen(output_im_filepath.string().c_str(), "w"); if(ofp_im == NULL){throw ta::messageFileOpenError(output_im_filepath.string());}
    
    sprintf (str, "%s.integrate.re", fileinfo.getFileStem().string().c_str());
    const path output_re_filepath   = dest_dir / str; // Specify the destination directory and avoid depending on the current directory for Gnuplot configuration.
    FILE *ofp_re = fopen(output_re_filepath.string().c_str(), "w"); if(ofp_re == NULL){throw ta::messageFileOpenError(output_re_filepath.string());}
    
    for(int t=0; t<integrate_batch * num_partition; t++){
        fprintf(ofp_im,  "%f¥n", integrate_im[t]);
        fprintf(ofp_re,  "%f¥n", integrate_re[t]);
    }
    
    fclose(ofp_im);
    fclose(ofp_re);
    
    delete [] h_idata;
    delete [] integrate_re;
    delete [] integrate_im;
}
