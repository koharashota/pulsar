
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
#include <cuda_runtime.h>

//c++のファイルやディレクトリを操作ライブラリ
#include <boost/filesystem.hpp>


#define STR_LEN_MAX 2048
#define CL

#ifdef CL
	#include <Windows.h>
	#define PLOT "pgnuplot"
#elif defined(GCC)
	#include <unistd.h>
	#define PLOT "gnuplot"
#endif

using std::cout;
using std::endl;
using std::string;
using boost::filesystem::path;

typedef struct {
	string xmin, xmax;
	string ymin, ymax;
	string zmin, zmax;
} plot_range;


int gnuplot(const char *data_filepath_char, const plot_range range, ta::Config &config);
void time_series_integration_complex (const float integration_point, float *data, int num_entire_data_pt,int tmp_num_partial_data, ta::FileInformation &fileinfo, ta::Config &config, const path dest_dir)

int main(int argc, char *argv[])
{

	bool flag_create_spectrogram_image = false;
	try{
		// Input File の引数に誤りがあったらエラーをはく
		if(argc != 3){
			printf("Usage: %s <file> <confiig>\n", argv[0]);
			throw "Invalid Arguments";
		}

		//fs ファイル情報宣言
    ta::FileInformation fileinfo (argv[1], '_');
		printf ("Data:   %s\n", fileinfo.getFilePath().string().c_str());
		printf ("Config: %s\n", path(argv[2]).string().c_str());	

		// Observation Parameters
    // config.txtで設定した情報を読み取って、定数を宣言
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
    
    //1/50*10^-9 = 20000000
    //tmp_num_partial_data>20000000
    int tmp_num_partial_data = static_cast<int>(config.getPulsarP0_pt() * 5); // 分割データ内に少なくとも1つのパルスが入るように、分割データの点数を設定。
    if (tmp_num_partial_data < 20000000){
      tmp_num_partial_data = 20000000
    }
		if (tmp_num_partial_data > num_entire_data_pt){ // ただし全データ数がそれに満たない場合、分割データ数を全データ数に設定。
			tmp_num_partial_data = num_entire_data_pt;
		}

		//const int   integrate_window_width = 1024;
		//const int   integrate_window_width       = config.getFFTWindowWidth();
		//何ポイントで一秒になるか
    const int   integrate_window_width       = 20000000;
    const int   integrate_batch        = std::floor(static_cast<float>(tmp_num_partial_data) / integrate_window_width);
    const int   num_partial_data       = integrate_window_width * integrate_batch;
		const int   num_partition          = std::floor( static_cast<float>(num_entire_data_pt) / num_partial_data);
		const int   skipped_data_size_pts  = num_entire_data_pt - num_partial_data * num_partition;
		
    
    
    // Display Parameters
		printf ("\nPulsar Information\n");
		printf ("- Period      P0      = %f ms\n", pulsar_p0_s  * 1000);
		printf ("- Pulse Width W50     = %f ms\n", pulsar_w50_s * 1000);
		printf ("\nData Information\n");
		printf ("- Sampling frequency  = %f MHz\n", sampling_freq_Hz / 1E+6);
		printf ("- Sampling interval   = %f ns\n",  sampling_interval_s * 1E+9);
		
    printf ("\nAnalysis Information\n");
		printf ("- Interate Window Width    = %d pt = %f ms\n", integrate_window_width, integrate_window_width * sampling_interval_s * 1000);
		printf ("- Partial Data Size   = %f ms = %d pt = %d MB (on RAM)\n", sampling_interval_s * num_partial_data * 1000, num_partial_data, sizeof(float)*num_partial_data/static_cast<int>(1E+6));
		printf ("- Entire  Data Size   = %f ms = %d pt\n",                  sampling_interval_s * num_entire_data_pt * 1000, num_entire_data_pt);
		printf ("- Skipped Data Size   = %f ms = %d pt\n",                  sampling_interval_s * skipped_data_size_pts * 1000, skipped_data_size_pts);
		printf ("- Number of Analysis Partitions = %d\n", num_partition);

		printf ("\nInitiate the process?");
		if (!ta::stdin_yes_or_no()) {
			throw "Task Terminated";
		}

		printf ("\nAnalyzing...\n");

		
    
    
    
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

    float *h_idata   = new float[num_partial_data];
		for(int i=0; i<num_partial_data; i++){
			h_idata[i].x = 0.0;
			h_idata[i].y = 0.0;
		}

    //1秒積分のために用意
    float *integrate_re = new float[integrate_batch * num_partition];
    float *integrate_im = new float[integrate_batch * num_partition];
		for(int f=0; f<integrate_batch * num_partition; f++){
			integrate_re[f] = 0;
			integrate_im[f] = 0;
		}

		ta::Stats stats; // Statistic parameters including mean & variance values

		//
		// MAIN
		//
    
    /**
    //1/50*10^-9 = 20000000
		int tmp_num_partial_data = static_cast<int>(config.getPulsarP0_pt() * 5); // 分割データ内に少なくとも1つのパルスが入るように、分割データの点数を設定。
    const int num_entire_data_pt = tmp_num_entire_data_pt;
		const int   integrate_window_width = 1024;
    const int   integrate_batch        = std::floor(static_cast<float>(tmp_num_partial_data) / integrate_window_width);
    const int   num_partial_data       = fft_window_width * fft_batch;
		const int   num_partition          = std::floor( static_cast<float>(num_entire_data_pt) / num_partial_data);
		const int   skipped_data_size_pts  = num_entire_data_pt - num_partial_data * num_partition;
		**/

    for (int pos = 0; pos < num_partition; pos++) {
			const int partial_data_begin_pt = pos * num_partial_data;
      
      //hi_dataの作成
      if (flag_complex_data == true) {
        const{ 
          float *cudata = h_idata
          const int extraction_first_point = partial_data_begin_pt
          const int extraction_width = num_partial_data
            
          if( !load_d_executed ){return -1;}
          
          if (extraction_first_point < 0 || extraction_width < 0 || extraction_width > num_data){
            throw "Exception: 関数 Data::extractData() の引数不正";
          }
          if (extraction_first_point + extraction_width> num_data){
            for(int i=0; i<extraction_width; i++){
              cudata[i].x = 0;
              cudata[i].y = 0;
            }
            for(int i=0; i<num_data - extraction_first_point; i++){
              //ここでデータの時間と強度がよめれば良い多分re, imで出来てる?
              cudata[i].x = static_cast<float>(static_cast<unsigned>(data[extraction_first_point + i].x));
              cudata[i].y = static_cast<float>(static_cast<unsigned>(data[extraction_first_point + i].y));
            }
          }else{
            for(int i=0; i<extraction_width; i++){
              //ここでデータの時間と強度がよめれば良い多分re, imで出来てる?
              cudata[i].x = static_cast<float>(static_cast<unsigned>(data[extraction_first_point + i].x)); // start から start+width-1 ポイントのデータを抽出
              cudata[i].y = static_cast<float>(static_cast<unsigned>(data[extraction_first_point + i].y)); // start から start+width-1 ポイントのデータを抽出
            }
          }
          return 0;
        }
      } else {
        throw "not complex data";
      }
     


			// Confirm Data For Debug
			if (pos == 0) {
				std::ofstream fout_test ( (dest_dir / (fileinfo.getFileStem().string() + ".txt")).string() );
				for (int i = 0; i < 100; i++) {
					fout_test << h_idata[i].x << "\t" << h_idata[i].y << "\n";
				}
				fout_test.close();
			}
      
      //一秒積分
      for(j=0; j < integrate_batch ; j++){
        //初期化
        double re=0.0;
        double im=0.0;
        int integrate_first_point = j * integrate_window_width + pos * num_partial_data;
        
        for(i=0; i < integrate_window_width; i++)
        {
          re+=h_idata[i+integrate_first_point];
          im+=h_idata[i+integrate_first_point];
        }
        integrate_re[j + pos * integrate_batch] = re;
        integrate_im[j + pos * integrate_batch] = im;
        //デバッグ
        printf("Intergral re : %lf, im : %lf\n", re , im);
      }

      
      //integrate_re,integrate_imを保存
      if (flag_complex_data == true) {
        path fout = dest_dir / (fileinfo.getFileStem().string() + ".intergrate.2d");
        //saveData2D (const char *filename, const float *data, const int data_length, const float x_resolution)
        ta::saveData2D(fout.string().c_str(),re , fft_window_width);	
      }
      else {
        path fout = dest_dir / (fileinfo.getFileStem().string() + ".integrate.1d");
        ta::saveData1D(fout.string().c_str(), band_amp, fft_window_width);
      }

		} // Next pos		
		
		// Delete
		delete [] h_idata;
		delete [] integrate_re;
		delete [] integrate_im;

		return 0;
	}
	catch (const char *err) {
		fprintf(stderr, "%s\n", err);
		system ("pause");
		return -1;
	}
	catch (const string err) {
		fprintf(stderr, "%s\n", err.c_str());
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

void time_series_integration_complex (const float integration_point, ta::BinData &data, int num_entire_data_pt, ta::FileInformation &fileinfo, ta::Config &config, const path dest_dir)
{
  
  int tmp_num_partial_data = 100000000 //適当にセット(多分必要ない。寧ろintegrate_window_width軽くした方がよさげ)
  const int   integrate_window_width  = integration_point;
  const int   integrate_batch = std::floor(static_cast<float>(tmp_num_partial_data) / integrate_window_width);
  const int   num_partial_data       = integrate_window_width * integrate_batch;
  const int   num_partition          = std::floor( static_cast<float>(num_entire_data_pt) / num_partial_data);
  const int   skipped_data_size_pts  = num_entire_data_pt - num_partial_data * num_partition;

  //1秒積分のために用意
  float *integrate_re = new float[integrate_batch * num_partition];
  float *integrate_im = new float[integrate_batch * num_partition];
  for(int f=0; f<integrate_batch * num_partition; f++){
    integrate_re[f] = 0;
    integrate_im[f] = 0;
  }
  
  float *h_idata   = new float[num_partial_data];
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
    // ここの処理でreとimのそれぞれの強度をポイント毎に入れられれば良い
    if (flag_complex_data == true) {
      const{ 
        float *cudata = h_idata
        const int extraction_first_point = partial_data_begin_pt
        const int extraction_width = num_partial_data
          
        if( !load_d_executed ){return -1;}
        
        if (extraction_first_point < 0 || extraction_width < 0 || extraction_width > num_data){
          throw "Exception: 関数 Data::extractData() の引数不正";
        }
        if (extraction_first_point + extraction_width> num_data){
          for(int i=0; i<extraction_width; i++){
            cudata[i].x = 0;
            cudata[i].y = 0;
          }
          for(int i=0; i<num_data - extraction_first_point; i++){
            //ここでデータの時間と強度がよめれば良い多分re, imで出来てる?
            cudata[i].x = static_cast<float>(static_cast<unsigned>(data[extraction_first_point + i].x));
            cudata[i].y = static_cast<float>(static_cast<unsigned>(data[extraction_first_point + i].y));
          }
        }else{
          for(int i=0; i<extraction_width; i++){
            //ここでデータの時間と強度がよめれば良い多分re, imで出来てる?
            cudata[i].x = static_cast<float>(static_cast<unsigned>(data[extraction_first_point + i].x)); // start から start+width-1 ポイントのデータを抽出
            cudata[i].y = static_cast<float>(static_cast<unsigned>(data[extraction_first_point + i].y)); // start から start+width-1 ポイントのデータを抽出
          }
        }
        return 0;
      }
    } else {
      throw "not complex data";
    }

    // Confirm Data For Debug
    if (pos == 0) {
      std::ofstream fout_test ( (dest_dir / (fileinfo.getFileStem().string() + ".txt")).string() );
      for (int i = 0; i < 100; i++) {
        fout_test << h_idata[i].x << "\t" << h_idata[i].y << "\n";
      }
      fout_test.close();
    }
    
    //一秒積分
    for(j=0; j < integrate_batch ; j++){
      //初期化
      double re=0.0;
      double im=0.0;
      int integrate_first_point = j * integrate_window_width + pos * num_partial_data;
      
      for(i=0; i < integrate_window_width; i++)
      {
        re+=h_idata[i+integrate_first_point];
        im+=h_idata[i+integrate_first_point];
      }
      integrate_re[j + pos * integrate_batch] = re;
      integrate_im[j + pos * integrate_batch] = im;
      //デバッグ
      printf("Intergral re : %lf, im : %lf\n", re , im);
    }

    //データセーブ integrate_re integrate_imを保存すればいいね
    path fout = dest_dir / (fileinfo.getFileStem().string() + ".intergrate.2d");
    //saveData2D (const char *filename, const float *data, const int data_length, const float x_resolution)
    ta::saveData2D(fout.string().c_str(),re , fft_window_width);	

  } // Next pos		
	
  delete [] h_idata;
  delete [] integrate_re;
  delete [] integrate_im;
}
