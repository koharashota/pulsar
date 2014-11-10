// 1�f�[�^�t�@�C�����̐ϕ��f�[�^��z��Ɋi�[���Amain �Ō゠����Ńt�@�C���o��
// �s�v�ȕ��͑S�ď���
// �t�@�C�������ɒ��r���[�ɍ�����ϕ��v���O����������̂œK���ɉ��ǂ���
// void time_series_integration_complex (const float integration_time_s, ta::FileInformation &fileinfo, ta::BinData &data, ta::Config &config, const path dest_dir)

//argv[0] ���s�t�@�C���� argv[1] �ϑ��f�[�^ argv[2] �ݒ�t�@�C�� config.txt?

#include <cstdio>
#include <iostream>
#include <fstream>
#include <string>
#include <ctime>
#include "Functions.h"
//c++�̃��C�u�������g���ăf�B���N�g���Ƃ��t�@�C���p�X�̑��� fs
#include "FileInformation.h"
//�f�[�^�̍��Ƃ��̐���
#include "Config.h"
#include "BinData.h"


#include <cuda.h>
#include <cuda_runtime.h>
//c++�̃t�@�C����f�B���N�g���𑀍색�C�u����
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


void time_series_integration_complex (int integration_point, ta::BinData &data, int num_entire_data_pt, ta::FileInformation &fileinfo, ta::Config &config, const path dest_dir)

int main(int argc, char *argv[])
{

	bool flag_create_spectrogram_image = false;
	try{
		// Input File �̈����Ɍ�肪��������G���[���͂�
		if(argc != 3){
			printf("Usage: %s <file> <confiig>\n", argv[0]);
			throw "Invalid Arguments";
		}

		//fs �t�@�C�����錾
    ta::FileInformation fileinfo (argv[1], '_');
		printf ("Data:   %s\n", fileinfo.getFilePath().string().c_str());
		printf ("Config: %s\n", path(argv[2]).string().c_str());	

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
		printf ("\nPulsar Information\n");
		printf ("- Period      P0      = %f ms\n", pulsar_p0_s  * 1000);
		printf ("- Pulse Width W50     = %f ms\n", pulsar_w50_s * 1000);
		printf ("\nData Information\n");
		printf ("- Sampling frequency  = %f MHz\n", sampling_freq_Hz / 1E+6);
		printf ("- Sampling interval   = %f ns\n",  sampling_interval_s * 1E+9);
		
    printf ("\nAnalysis Information\n");
		printf ("- Entire  Data Size   = %f ms = %d pt\n",                  sampling_interval_s * num_entire_data_pt * 1000, num_entire_data_pt);
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
		//config.txt�ɐݒ肳�ꂽ�p�X��ݒ�
    const path dest_dir = path (config.getOutputDirectory()) / fileinfo.getFileStemSubstr(); // = /opt/data/summation
		//c++�̃��C�u�������g���ăp�X�̍쐬
		boost::filesystem::create_directories (dest_dir);
		boost::filesystem::current_path (dest_dir);

    
    //1s�ϕ��̓ǂݏo�� 
    time_series_integration_complex (20000000, data, num_entire_data_pt, fileinfo, config, dest_dir) 

    // Delete
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
  �E�E���݂��߂΂��������ȁ[
  integration_point ���|�C���g�̍��v�l����邩(��b�ϕ��Ȃ�1/50^[-9]�Ȃ̂�20000000�|�C���g)
  data �o�C�i������ǂݎ�����f�[�^�z������炤�B�X�y�N�g���̑傫���̔z����󂯎�肽��(ta::BinData data.load_binary_as_double�ł�񂾃f�[�^)
  num_entire_data_pt �f�[�^�̑S�|�C���g��
  tmp_num_partial_data �����f�[�^�̃|�C���g��(����Ȃ������A�Ă�����Ȃ���)
  fileinfo file�̃p�X���Ƃ��󂯎��
  config config�̓ǂݎ��
  dest_dir �ۑ���̃f�B���N�g�����ȁH
  complex_flag�̓I���Ɖ���
**/

void time_series_integration_complex (int integration_point, ta::BinData &data, int num_entire_data_pt, ta::FileInformation &fileinfo, ta::Config &config, const path dest_dir)
{
  
  int tmp_num_partial_data = 100000000 //�K���ɃZ�b�g(�K�v�Ȃ��B�J��integrate_window_width�y�����������悳��)
  const int   integrate_window_width  = integration_point;
  const int   integrate_batch = std::floor(static_cast<float>(tmp_num_partial_data) / integrate_window_width);
  const int   num_partial_data       = integrate_window_width * integrate_batch;
  const int   num_partition          = std::floor( static_cast<float>(num_entire_data_pt) / num_partial_data);
  const int   skipped_data_size_pts  = num_entire_data_pt - num_partial_data * num_partition;

  //1�b�ϕ��̂��߂ɗp��
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
    //hi_data�̍쐬
    // �����̏�����re��im�̂��ꂼ��̋��x���|�C���g���ɓ������Ηǂ�
    if (flag_complex_data == true) {
      const{

        float *cudata = h_idata
        const int extraction_first_point = partial_data_begin_pt
        const int extraction_width = num_partial_data
          
        if( !load_d_executed ){return -1;}
        
        if (extraction_first_point < 0 || extraction_width < 0 || extraction_width > num_data){
          throw "Exception: �֐� Data::extractData() �̈����s��";
        }
        if (extraction_first_point + extraction_width> num_data){
          for(int i=0; i<extraction_width; i++){
            cudata[i].x = 0;
            cudata[i].y = 0;
          }
          for(int i=0; i<num_data - extraction_first_point; i++){
            //�����Ńf�[�^�̎��ԂƋ��x����߂�Ηǂ�����re, im�ŏo���Ă�?
            cudata[i].x = static_cast<float>(static_cast<unsigned>(data[extraction_first_point + i].x));
            cudata[i].y = static_cast<float>(static_cast<unsigned>(data[extraction_first_point + i].y));
          }
        }else{
          for(int i=0; i<extraction_width; i++){
            //�����Ńf�[�^�̎��ԂƋ��x����߂�Ηǂ�����re, im�ŏo���Ă�?
            cudata[i].x = static_cast<float>(static_cast<unsigned>(data[extraction_first_point + i].x)); // start ���� start+width-1 �|�C���g�̃f�[�^�𒊏o
            cudata[i].y = static_cast<float>(static_cast<unsigned>(data[extraction_first_point + i].y)); // start ���� start+width-1 �|�C���g�̃f�[�^�𒊏o
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
    
    //��b�ϕ�
    for(j=0; j < integrate_batch ; j++){
      //������
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
      //�f�o�b�O
      printf("Intergral re : %lf, im : %lf\n", re , im);
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
    fprintf(ofp_im,  "%f\n", integrate_im[t]);
    fprintf(ofp_re,  "%f\n", integrate_re[t]);
	}
	
  fclose(ofp_im);
	fclose(ofp_re);
  
  delete [] h_idata;
  delete [] integrate_re;
  delete [] integrate_im;
}
