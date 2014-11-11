/**
*
* Notes:
*	- Use 'float' instead of 'double' because 'cufftComplex' is 'float'.
*/

#include <cstdio>
#include <fstream>
#include <string>
#include <cuda.h>
#include <cufft.h>
#include <cuda_runtime.h>
#include <boost/filesystem.hpp>
#include "Functions.h"
#include "FileInformation.h"
#include "Config.h"

#include <ctime>

#pragma comment(lib, "cudart")
#pragma comment(lib, "cufft")

#define STR_LEN_MAX 2048
#define CL

#ifdef CL
	#include <Windows.h>
	#define PLOT "pgnuplot"
#elif defined(GCC)
	#include <unistd.h>
	#define PLOT "gnuplot"
#endif

using std::string;
using boost::filesystem::path;
typedef struct {
	int xmin, xmax;
	int ymin, ymax;
	int zmin, zmax;
	//string xmin, xmax;
	//string ymin, ymax;
	//string zmin, zmax;
} plot_range;

int loadBinaryDataComplex(float *data_re, float *data_im, const int data_length, const char *filename, const int offset_B, const int datum_bit_length_B);
int loadBinaryData(float *data, const int data_length, const char *filename, const int offset_B, const int datum_bit_length_B);
int STFT(cufftComplex *h_idata, const int fft_window_width, const int fft_batch, bool inverse);
int gnuplot(const char *data_filepath_char, const plot_range range, const float sampling_interval_s, const float fft_time_resolution_s, const float fft_freq_resolution_Hz);


int main(int argc, char *argv[])
{
	//bool create_spectrogram_image = true;
	bool create_spectrogram_image = false;
	bool complex_data   = true; // True if data format is <Re 1 byte><Im 1 byte><Re 1 byte><Im 1 byte>...

	try{
		char buf1[STR_LEN_MAX + 1];	
		char buf2[STR_LEN_MAX + 1];	
		
		// Input File
		if(argc != 3){
			fprintf(stderr, "Usage: %s <file> <confiig>\n", argv[0]);
			throw "Invalid Arguments";
		}
		Aoki::FileInformation fileinfo(argv[1], '_');
		printf("Data:   %s\n", fileinfo.getFilePath().string().c_str());
		printf("Config: %s\n", path(argv[2]).string().c_str());	

		// Observation Parameters
		Aoki::Config config(argv[2]);
		const int    datum_bit_length_B  = config.getDatumBitLength_B(); // byte
		const float  sampling_freq_Hz    = config.getSamplingFrequency_Hz(); // Hz, sample/sec
		const float  sampling_interval_s = config.getSamplingInterval_s(); // sec
		//const float  beam_width_deg      = config.getBeamWidth_deg(); // degree
		const int    entire_data_point   = fileinfo.getFilesize_B() / datum_bit_length_B;

		const float  pulsar_p0_s         = config.getPulsarP0_s(); // sec, Barycentric period of the pulsar
		const float  pulsar_w50_s        = config.getPulsarW50_s(); // sec, Width of pulse at 50% of peak.
		const int    pulsar_p0_point     = static_cast<int>(pulsar_p0_s / sampling_interval_s);
		const int    pulsar_w50_point    = static_cast<int>(pulsar_w50_s / sampling_interval_s);


		// Analysis Parameters
		// データファイルは大容量なので、適当なデータサイズに分割して解析する。
		// 規定の分割データサイズに満たないデータは読み捨てられるので、1ファイルの最後の方は解析されない。
		// 無限に分割し続けるわけにはいかないので、最大分割数を設定する。
		//int        tmp_partition_num = 0; // Division number of data points
		//const int  max_partition_num = 1000; 
		//const float allowed_pointing_error_s   = 240 * beam_width_deg / 40; // Pointing error must be less than twentieth or fortieth part of the beam width. 1 deg = 24/360 hour = 24/360 * 3600 sec = 240 sec
		//const int   allowed_pointing_error_pts = std::floor (allowed_pointing_error_s / sampling_interval_s);

		const float partition_width_approx_s = config.getDataPartitioningWdithApprox_s(); // sec
		const int   fft_window_width       = config.getFFTWindowWidth(); // point
		const float fft_time_resolution_s  = sampling_interval_s * fft_window_width; // sec
		const float fft_freq_resolution_Hz = 1.0 / fft_time_resolution_s; // Hz
		const int   fft_batch              = std::floor (partition_width_approx_s / fft_time_resolution_s);
		const int   partial_data_width_pts = fft_window_width * fft_batch;	

		
		//int tmp_partial_data_width_pts = pow( 2.0, ceil(log(5.0 * pulsar_w50_s / sampling_interval_s)/log(2.0)) ); // points; データ数を2^mになるよう設定。パルス全幅の約5倍のデータをとる。
		//int tmp_partial_data_width_pts = static_cast<int>(pulsar_p0_point * 2); // 分割データ内に少なくとも1つのパルスが入るように、分割データの点数を設定。
		if (partial_data_width_pts > entire_data_point){ 
			fprintf(stderr, "The data partitioning width exceeds the entire data width of the input file. Edit the config file: %s\n", argv[2]);
			throw "Invalid Configuration";
		}
		
		const int   num_partition          = std::floor( static_cast<float>(entire_data_point) / partial_data_width_pts);
		const int   skipped_data_size_pts  = entire_data_point - partial_data_width_pts * num_partition;

		// Display Parameters
		printf("\nPulsar Information\n");
		printf("- Period      P0      = %f ms\n", pulsar_p0_s  * 1000);
		printf("- Pulse Width W50     = %f ms\n", pulsar_w50_s * 1000);
		printf("\nData Information\n");
		printf("- Sampling frequency  = %f MHz\n", sampling_freq_Hz / 1E+6);
		printf("- Sampling interval   = %f ns\n",  sampling_interval_s * 1E+9);
		printf("\nAnalysis Information\n");
		printf("- FFT Window Width    = %d pts = %f ms\n", fft_window_width, fft_window_width * sampling_interval_s * 1000);
		printf("- FFT Time resolution = %f ms\n",  fft_time_resolution_s  * 1000);
 		printf("- FFT Freq resolution = %f kHz\n", fft_freq_resolution_Hz / 1000);
		if(fft_time_resolution_s > pulsar_w50_s){
			printf("\n### CAUTION ###\nFFT Time Resolution > Pulse Width.\nThis analysis can not resolve the pulse profile.\n");
		}
		printf("- Partial Data Size   = %f ms = %d pts = %d MB (on RAM)\n", sampling_interval_s * partial_data_width_pts * 1000, partial_data_width_pts, sizeof(float) * partial_data_width_pts / static_cast<int>(1E+6));
		printf("- Entire  Data Size   = %f ms = %d pts\n",                  sampling_interval_s * entire_data_point * 1000,entire_data_point);
		printf("- Skipped Data Size   = %f ms = %d pts\n",                  sampling_interval_s * skipped_data_size_pts * 1000, skipped_data_size_pts);
		printf("- Number of Data Partitions = %d\n", num_partition);
		printf("\nAnalyzing...\n");

		// Create Directory
		const path dest_dir = path(config.getOutputDirectory()) / fileinfo.getFileStemSubstr(); // = /opt/data/summation
		boost::filesystem::create_directories(dest_dir);
		boost::filesystem::current_path(dest_dir);

		// Initialize Parameters
		cufftComplex *h_idata   = new cufftComplex[partial_data_width_pts];
		float *partial_data_re  = new float[partial_data_width_pts];
		float *partial_data_im  = new float[partial_data_width_pts];
		float *fft_re           = new float[fft_window_width];
		float *fft_im           = new float[fft_window_width];
		float **spectral_amp    = new float*[fft_batch]; // Power Spectral Density, psd[time][freq]
		float **spectral_phase  = new float*[fft_batch];
		for(int t=0; t<fft_batch; t++){
			spectral_amp[t]   = new float[fft_window_width];
			spectral_phase[t] = new float[fft_window_width];
		}
		for(int i=0; i<partial_data_width_pts; i++){
			partial_data_re[i] = 0.0;
			partial_data_im[i] = 0.0;
			h_idata[i].x       = 0.0;
			h_idata[i].y       = 0.0;
		}
		float mean_re, stdev_re, div_by_stdev_re;
		float mean_im, stdev_im, div_by_stdev_im;

		float *spectral_amp_at_t   = new float[fft_window_width];
		float *spectral_phase_at_t = new float[fft_window_width];

		//
		// MAIN
		//
		// DEBUG
		

		int loop_break_flag = 1;
		for (int pos=0; pos<num_partition; pos++) {

			clock_t start_time = clock();

			// Load Partial Data
			for (int i=0; i<partial_data_width_pts; i++) {
				partial_data_re[i] = 0;
				partial_data_im[i] = 0;
			}
			if (complex_data == true) {
				loop_break_flag = loadBinaryDataComplex(partial_data_re, partial_data_im, partial_data_width_pts, fileinfo.getFilePath().string().c_str(),  pos * datum_bit_length_B * partial_data_width_pts, datum_bit_length_B);
			} else {
				loop_break_flag = loadBinaryData(partial_data_re, partial_data_width_pts, fileinfo.getFilePath().string().c_str(),  pos * datum_bit_length_B * partial_data_width_pts, datum_bit_length_B);
			}

			clock_t end1 = clock();
			printf("Data Loading: %.2f s\n", (double)(end1 - start_time)/CLOCKS_PER_SEC);
			
			/* If skipping the last residue part of the data file, the following description is needed.
			if (loop_break_flag != 0) { // Need to break also here.
				break;
			}
			*/


			// For debug
			if (pos == 0) {
				std::ofstream fout ("test.txt");
				for (int i=0; i<100; i++){
					fout << partial_data_re[i] << "\t" << partial_data_im[i] << "\n";
				}
				fout.close();
			}
			
			// Normalize data using the mean & variance of the first data block (pos = 0)
			if (pos == 0) {
				mean_re  = Aoki::mean (partial_data_re, partial_data_width_pts);
				mean_im  = Aoki::mean (partial_data_im, partial_data_width_pts);
				stdev_re = std::sqrt  (Aoki::unbiased_variance(partial_data_re, partial_data_width_pts, mean_re));
				stdev_im = std::sqrt  (Aoki::unbiased_variance(partial_data_im, partial_data_width_pts, mean_im));
				div_by_stdev_re = 1.0 / stdev_re;
				div_by_stdev_im = 1.0 / stdev_im;
			}
			for (int i=0; i<partial_data_width_pts; i++) {
				partial_data_re[i] = (partial_data_re[i] - mean_re) * div_by_stdev_re;
				partial_data_im[i] = (partial_data_im[i] - mean_im) * div_by_stdev_im;
			}
			
			// Create input data on CPU
			// h_idata = input data on host CPU
			for(int i=0; i<partial_data_width_pts; i++){
				h_idata[i].x = partial_data_re[i];
				h_idata[i].y = partial_data_im[i];
			}

			// Forward FFT
			clock_t s = clock();
			STFT(h_idata, fft_window_width, fft_batch, false);
			clock_t e = clock();
			printf("STFT: %.2f s\n", (double)(e-s)/CLOCKS_PER_SEC);


			


			// Create a spectrogram, psd[t][f]
			//h_idata[0].x = 0; h_idata[0].y = 0; // Remove the direct current (DC) component
			int index = 0; 
			const float div_by_root_N = 1.0 / std::sqrt (fft_window_width);
			for (int t=0; t<fft_batch; t++) {
				for (int f=0; f<fft_window_width; f++) {
					fft_re[f]            = h_idata[index].x * div_by_root_N; 
					fft_im[f]            = h_idata[index].y * div_by_root_N; // This assignment is only for code readability.
					spectral_amp[t][f]   = std::hypot (fft_re[f], fft_im[f]);
					spectral_phase[t][f] = std::atan2 (fft_im[f], fft_re[f]);
					index++;
				}
			}


			// Spectrum File
			sprintf(buf1, "%s.%d-%dS.averaged_amplitude_spectrum.2d", fileinfo.getFileStem().string().c_str(), pos * partial_data_width_pts, (pos+1) * partial_data_width_pts);
			sprintf(buf2, "%s.%d-%dS.averaged_phase_spectrum.2d", fileinfo.getFileStem().string().c_str(), pos * partial_data_width_pts, (pos+1) * partial_data_width_pts);
			const path output_ampspe_filepath = dest_dir / buf1;
			const path output_phaspe_filepath = dest_dir / buf2;
			const float div_by_fft_batch = 1.0 / fft_batch;
			spectral_amp_at_t[0]   = 0;
			spectral_phase_at_t[0] = 0;
			for(int f=1; f<fft_window_width; f++){
				spectral_amp_at_t[f]   = 0;
				spectral_phase_at_t[f] = 0;
				for(int t=0; t<fft_batch; t++){
					spectral_amp_at_t[f]   += spectral_amp[t][f];
					spectral_phase_at_t[f] += spectral_phase[t][f];
				}
				spectral_amp_at_t[f]   *= div_by_fft_batch;
				spectral_phase_at_t[f] *= div_by_fft_batch;
			}
			Aoki::saveSpectrumOfComplexData(output_ampspe_filepath.string().c_str(), spectral_amp_at_t,   fft_window_width);
			//Aoki::saveSpectrumOfComplexData(output_phaspe_filepath.string().c_str(), spectral_phase_at_t, fft_window_width);
		
		



			// Create a spectrogram image
			if(create_spectrogram_image == true){

				// Save the spectrogram as a text file
				sprintf(buf1, "%s.%d-%dPT.spectrogram.3d", fileinfo.getFileStem().string().c_str(), pos * partial_data_width_pts, (pos+1) * partial_data_width_pts);
				sprintf(buf2, "%s.%d-%dPT.phase.3d", fileinfo.getFileStem().string().c_str(), pos * partial_data_width_pts, (pos+1) * partial_data_width_pts);
				const path output_spectrogram_filepath = dest_dir / buf1; // Specify the destination directory and avoid depending on the current directory for Gnuplot configuration.
				const path output_phase_filepath       = dest_dir / buf2; // Specify the destination directory and avoid depending on the current directory for Gnuplot configuration.
				FILE *ofp_spec  = fopen(output_spectrogram_filepath.string().c_str(), "w"); if(ofp_spec == NULL) {throw Aoki::messageFileOpenError(output_spectrogram_filepath.string());}
				FILE *ofp_phase = fopen(output_phase_filepath.string().c_str(), "w");       if(ofp_phase == NULL){throw Aoki::messageFileOpenError(output_phase_filepath.string());}

				for(int t=0; t<fft_batch; t++){
					for(int f=0; f<fft_window_width; f++){
						fprintf(ofp_spec,  "%d\t%d\t%.3f\n", pos * partial_data_width_pts + t * fft_window_width, f, spectral_amp[t][f]*spectral_amp[t][f]);
						fprintf(ofp_phase, "%d\t%d\t%.2f\n", pos * partial_data_width_pts + t * fft_window_width, f, spectral_phase[t][f]);
					}
					fprintf(ofp_spec, "\n");
					fprintf(ofp_phase, "\n");
				}
				fclose(ofp_spec);
				fclose(ofp_phase);


				// Make an PNG image of the spectrogram
				plot_range range;
				range.xmin = pos       * partial_data_width_pts;
				range.xmax = (pos + 1) * partial_data_width_pts;
				range.ymin = 0;
				range.ymax = fft_window_width;
				range.zmin = 0;
				range.zmax = 30;
				gnuplot(output_spectrogram_filepath.string().c_str(), range, sampling_interval_s, fft_time_resolution_s, fft_freq_resolution_Hz);

				range.zmin = - 4;
				range.zmax = 4;
				gnuplot(output_phase_filepath.string().c_str(), range, sampling_interval_s, fft_time_resolution_s, fft_freq_resolution_Hz);
			}

			clock_t end_time = clock();
			printf("Total Elapse: %.2f s\n", (double)(end_time - start_time)/CLOCKS_PER_SEC);
			//printf("Number of analyzed data partition: %d of %d\r", pos, num_partition); fflush(stdout);
			printf("Number of analyzed data partition: %d of %d\n", pos + 1, num_partition); fflush(stdout);
		} // Next pos		
		printf("\n");
				


		// Delete
		delete [] partial_data_re;
		delete [] partial_data_im;
		delete [] h_idata;
		for(int t=0; t<fft_batch; t++){
			delete [] spectral_amp[t];
			delete [] spectral_phase[t];
		}
		delete [] spectral_amp;
		delete [] spectral_phase;
		return 0;
	
	}catch(const char *err){
		fprintf(stderr, "%s\n", err);
		return -1;
	}catch(const string err){
		fprintf(stderr, "%s\n", err.c_str());
		return -1;
	}
}

int loadBinaryDataComplex(float *data_re, float *data_im, const int data_length, const char *filename, const int offset_B, const int datum_bit_length_B)
{
	FILE *ifp = fopen(filename, "rb");
	if(ifp == NULL){
		throw Aoki::messageFileOpenError(filename);
	}

	int error_flag = 1;
	unsigned char unsigned_re, unsigned_im;
	fseek(ifp, offset_B, SEEK_SET);
	for(int i=0; i<data_length; i++){
		error_flag = fread(&unsigned_re, datum_bit_length_B, 1, ifp); if(error_flag == 0) break;
		error_flag = fread(&unsigned_im, datum_bit_length_B, 1, ifp); if(error_flag == 0) break;
		data_re[i] = static_cast<float>(unsigned_re);
		data_im[i] = static_cast<float>(unsigned_im);
		//printf("(%f, %f), ", data_re[i], data_im[i]);
	}
	fclose(ifp);

	if(error_flag == 0){
		return -1;
	}
	return 0;
}

int loadBinaryData(float *data, const int data_length, const char *filename, const int offset_B, const int datum_bit_length_B)
{
	FILE *ifp = fopen(filename, "rb");
	if(ifp == NULL){
		throw Aoki::messageFileOpenError(filename);
	}

	int error_flag = 1;
	unsigned char unsigned_re;
	fseek(ifp, offset_B, SEEK_SET);
	for(int i=0; i<data_length; i++){
		error_flag = fread(&unsigned_re, datum_bit_length_B, 1, ifp); if(error_flag == 0) break;
		data[i] = static_cast<float>(unsigned_re);
	}
	fclose(ifp);

	if(error_flag == 0){
		return -1;
	}
	return 0;
}

int STFT(cufftComplex *h_idata, const int fft_window_width, const int fft_batch, bool inverse)
{
	const int NX = fft_window_width * fft_batch;

	// Transfer the data from CPU to GPU
	// d_idata = input data on device GPU
	cufftComplex *d_idata;
	cudaMalloc( (void**)&d_idata, sizeof(cufftComplex) * NX);
	cudaMemcpy(d_idata, h_idata,  sizeof(cufftComplex) * NX, cudaMemcpyHostToDevice);
			
	// Create an FFT plan
	cufftHandle plan;
	cufftPlan1d(&plan, fft_window_width, CUFFT_C2C, fft_batch);

	// Execute the FFT
	if(inverse == false){
		cufftExecC2C(plan, d_idata, d_idata, CUFFT_FORWARD);
	}else{
		cufftExecC2C(plan, d_idata, d_idata, CUFFT_INVERSE);
	}

	// Transfer data from GPU to CPU 
	cudaMemcpy(h_idata, d_idata, sizeof(cufftComplex) * NX, cudaMemcpyDeviceToHost);

	// Erase the data on GPU
	cudaFree(d_idata);
	cufftDestroy(plan);
	return 0;
}

int gnuplot(const char *data_filepath_char, const plot_range range, const float sampling_interval_s, const float fft_time_resolution_s, const float fft_freq_resolution_Hz)
{
	// Gnuplot cannot recognize '\' in a file path, and therefore '\' must be repalced to '/'.
	string data_filepath   = data_filepath_char; Aoki::strReplace(data_filepath, "\\", "/");
	string config_filepath = data_filepath + ".plt";

	// Create a config file
	FILE *ofp = fopen(config_filepath.c_str(), "w");
	fprintf(ofp, "reset\n");
	fprintf(ofp, "set pm3d map corners2color max\n");
	fprintf(ofp, "set palette gray\n");
	fprintf(ofp, "set lmargin 0\n");
	fprintf(ofp, "set rmargin 0\n");
	fprintf(ofp, "set tmargin 0\n");
	fprintf(ofp, "set bmargin 0\n");
	
	fprintf(ofp, "set title \"STFT from #%d to #%d; Time resolution %f ms; Freq resolution %f kHz\"\n", range.xmin, range.xmax, fft_time_resolution_s * 1E+3, fft_freq_resolution_Hz / 1000);
	fprintf(ofp, "set xlabel \"Time [%f ns]\"\n", sampling_interval_s * 1E+9);
	fprintf(ofp, "set ylabel \"Frequency [%f kHz]\"\n", fft_freq_resolution_Hz / 1000);

	//fprintf(ofp, "set xrange [%d:%d]\n", 0, xend - xstart);
	fprintf(ofp, "set xrange [%d:%d]\n", range.xmin, range.xmax);
	fprintf(ofp, "set yrange [%d:%d]\n", range.ymin, range.ymax);
	fprintf(ofp, "set zrange [%d:%d]\n", range.zmin, range.zmax); // When data ~ N(0,1), spectral power ~ Exp(1/2). p-value = 28 for a significance level of 1E-6.
	fprintf(ofp, "set ytics 128\n");

	fprintf(ofp, "set terminal png font \"Verdana\" 10 size 4800,3000\n"); // font = Verdana or Serif; Consolas looks dirty.
	fprintf(ofp, "set output \"%s.png\"\n", data_filepath.c_str());
	fprintf(ofp, "splot \"%s\" using 1:2:3 notitle\n", data_filepath.c_str());
	fprintf(ofp, "# If skipping some data, use an splot option 'every 5:5'.\n");
	//fprintf(ofp, "splot \"%s\" using ($1 - %d):2:3 notitle\n", input_data_file_path, xstart);

	fprintf(ofp, "\n"); // Important \n. If it does not exist, gnuplot may not progress.
	fclose(ofp);

	// Load the config file and plot the data
	FILE *gp = _popen(PLOT, "w");
	fprintf(gp, "load \"%s\"\n", config_filepath.c_str());
	fprintf(gp, "exit\n");
	_pclose(gp);
	
	return 0;
}
