
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
#include "Stats.h"

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

typedef struct {
	string xmin, xmax;
	string ymin, ymax;
	string zmin, zmax;
} plot_range;


int gnuplot(const char *data_filepath_char, const plot_range range, ta::Config &config);
void time_series_integration_complex (const float integration_time_s, ta::FileInformation &fileinfo, ta::BinData &data, ta::Config &config, const path dest_dir)

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
    // config.txt�Őݒ肵������ǂݎ���āA�萔��錾
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
		int tmp_num_partial_data = static_cast<int>(config.getPulsarP0_pt() * 5); // �����f�[�^���ɏ��Ȃ��Ƃ�1�̃p���X������悤�ɁA�����f�[�^�̓_����ݒ�B
		if (tmp_num_partial_data > num_entire_data_pt){ // �������S�f�[�^��������ɖ����Ȃ��ꍇ�A�����f�[�^����S�f�[�^���ɐݒ�B
			tmp_num_partial_data = num_entire_data_pt;
		}

		//const int   integrate_window_width = 1024;
		const int   integrate_window_width       = config.getFFTWindowWidth();
    const int   integrate_batch        = std::floor(static_cast<float>(tmp_num_partial_data) / integrate_window_width);
    const int   num_partial_data       = fft_window_width * fft_batch;
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
		//config.txt�ɐݒ肳�ꂽ�p�X��ݒ�
    const path dest_dir = path (config.getOutputDirectory()) / fileinfo.getFileStemSubstr(); // = /opt/data/summation
		//c++�̃��C�u�������g���ăp�X�̍쐬
		boost::filesystem::create_directories (dest_dir);
		boost::filesystem::current_path (dest_dir);

		// Initialize Parameters
		//num_partial_data�̓o�b�`�����œǂݍ��ރf�[�^��
    //h_idata�͉��Ɏg��?
    //spectral_amp[t], spectral_phase[t]�̔����������Bfft_batch������t������Ă�
    /**
    cufftComplex *h_idata   = new cufftComplex[num_partial_data];
		float **spectral_amp    = new float*[fft_batch]; // Power Spectral Density, psd[time][freq]
		float **spectral_phase  = new float*[fft_batch];
		for(int t=0; t<fft_batch; t++){
			spectral_amp[t]   = new float[fft_window_width];
			spectral_phase[t] = new float[fft_window_width];
		}
		for(int i=0; i<num_partial_data; i++){
			h_idata[i].x = 0.0;
			h_idata[i].y = 0.0;
		}
    **/

		//amp��phase�ɂ���鐔�l�̏�����
    float *h_idata   = new float[num_partial_data];
    float *gain_amp   = new float[integrate_window_width];
		float *gain_phase = new float[integrate_window_width];
		float **spectral_amp    = new float*[integrate_batch]; // Power Spectral Density, psd[time][freq]
		float **spectral_phase  = new float*[integrate_batch];
		for(int f=0; f<integrate_window_width; f++){
			gain_amp[f] = 0;
			gain_phase[f] = 0;
		}
		for(int i=0; i<num_partial_data; i++){
			h_idata[i].x = 0.0;
			h_idata[i].y = 0.0;
		}

		ta::Stats stats; // Statistic parameters including mean & variance values

		//
		// MAIN
		//
    //num_partition�͕����t�@�C����
    //�����Ńo�C�i���f�[�^�̓ǂݎ������Ă��܂���
    
    /**
    //1/50*10^-9 = 20000000
		int tmp_num_partial_data = static_cast<int>(config.getPulsarP0_pt() * 5); // �����f�[�^���ɏ��Ȃ��Ƃ�1�̃p���X������悤�ɁA�����f�[�^�̓_����ݒ�B
    const int num_entire_data_pt = tmp_num_entire_data_pt;
		const int   integrate_window_width = 1024;
    const int   integrate_batch        = std::floor(static_cast<float>(tmp_num_partial_data) / integrate_window_width);
    const int   num_partial_data       = fft_window_width * fft_batch;
		const int   num_partition          = std::floor( static_cast<float>(num_entire_data_pt) / num_partial_data);
		const int   skipped_data_size_pts  = num_entire_data_pt - num_partial_data * num_partition;
		**/

    for (int pos = 0; pos < num_partition; pos++) {
			const int partial_data_begin_pt = pos * num_partial_data;
      if (flag_complex_data == true) {
        //data.extract_binary_data_xy (h_idata, partial_data_begin_pt, num_partial_data);
        //int BinData::extract_binary_data_xy(cufftComplex *cudata, const int extraction_first_point, const int extraction_width) const
        const{ 
          float *cudata = h_idata
          const int extraction_first_point = partial_data_begin_pt
          const int extraction_width num_partial_data
            
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
              //�����Ńf�[�^�̎��ԂƋ��x����߂�Ηǂ�
              cudata[i].x = static_cast<float>(static_cast<unsigned>(data[extraction_first_point + i].x));
              cudata[i].y = static_cast<float>(static_cast<unsigned>(data[extraction_first_point + i].y));
            }
          }else{
            for(int i=0; i<extraction_width; i++){
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

      //�����悤�Ȋ����ň�b�ϕ�
      //stft (h_idata, fft_window_width, fft_batch, false);

		  /**
        �����悤�Ȋ�����gnuplot
      clock_t spec0 = clock();
			spectrogram (spectral_amp, spectral_phase, fft_batch, fft_window_width, h_idata);
			clock_t spec1 = clock();
			ta::print_elapse ("Spectrogram", spec0, spec1);
      **/

		} // Next pos		
		
		// Delete
		delete [] h_idata;
		for(int t=0; t<fft_batch; t++){
			delete [] spectral_amp[t];
			delete [] spectral_phase[t];
		}
		delete [] spectral_amp;
		delete [] spectral_phase;

		delete [] gain_amp;
		delete [] gain_phase;

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
void time_series_integration_complex (const float integration_time_s, ta::FileInformation &fileinfo, ta::BinData &data, ta::Config &config, const path dest_dir)
{
	const int num_entire_data_pt = data.getNumData();
	const int integration_time_pt = integration_time_s / config.getSamplingInterval_s();
	const int num_integ_data_pt = num_entire_data_pt / integration_time_pt + 1; // Important +1
	const int num_residue_pt = num_entire_data_pt % integration_time_pt;

	cufftComplex *data_xy = new cufftComplex[integration_time_pt];
	float *tmp_data_x = new float[integration_time_pt];
	float *tmp_data_y = new float[integration_time_pt];
	float *integ_data_x = new float[num_integ_data_pt];
	float *integ_data_y = new float[num_integ_data_pt];

	for (int pos = 0; pos < num_integ_data_pt; pos++) {
		data.extract_binary_data_xy (data_xy, pos * integration_time_pt, integration_time_pt);
		for (int i = 0; i < integration_time_pt; i++) {
			tmp_data_x[i] = data_xy[i].x;
			tmp_data_y[i] = data_xy[i].y;
		}

		if (pos < num_integ_data_pt - 1) {
			integ_data_x[pos] = ta::mean (tmp_data_x, integration_time_pt);
			integ_data_y[pos] = ta::mean (tmp_data_y, integration_time_pt);
		} else {
			integ_data_x[pos] = ta::mean (tmp_data_x, num_residue_pt);
			integ_data_y[pos] = ta::mean (tmp_data_y, num_residue_pt);
		}
	}

//	path filepath_re = 
//	path filepath_im = 
//	ta::saveData1D (filepath_re.string().c_str(), integ_data_x, num_integ_data_pt);
//	ta::saveData1D (filepath_re.string().c_str(), integ_data_x, num_integ_data_pt);
	//ta::saveData2D (filepath_re.string().c_str(), integ_data_x, num_integ_data_pt, integration_time_pt);
	//ta::saveData2D (filepath_im.string().c_str(), integ_data_y, num_integ_data_pt, integration_time_pt);


	
	delete [] data_xy;
	delete [] tmp_data_x;
	delete [] tmp_data_y;
	delete [] integ_data_x;
	delete [] integ_data_y;
}
**/
