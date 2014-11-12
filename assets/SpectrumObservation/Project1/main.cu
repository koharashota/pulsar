/*
* Usage:
*	<this program> <data file> <config file>
*
* Binary Data Format must be real or complex number data such as
*	- <Re><Re><Re>...
*	- <Re><Im><Re><Im>...
*
*/

#include <cstdio>
#include <iostream>
#include <fstream>
#include <string>
#include <ctime>
#include "Functions.h"
#include "FileInformation.h"
#include "Config.h"
#include "BinData.h"
#include "Stats.h"

#include <cuda.h>
#include <cufft.h>
#include <cuda_runtime.h>
#include <boost/filesystem.hpp>

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

using std::cout;
using std::endl;
using std::string;
using boost::filesystem::path;

typedef struct {
	string xmin, xmax;
	string ymin, ymax;
	string zmin, zmax;
} plot_range;

void stft (cufftComplex *h_idata, const int fft_window_width, const int fft_batch, bool inverse);
void normalize_data (cufftComplex *partial_data, const int num_partial_data, const ta::Stats &stats);
void spectrogram (float **spectral_amp, float **spectral_phase, const int num_time_index, const int num_freq_index, const cufftComplex *partial_data);
void mid_band_characteristics (float *band_amp, const int fft_batch, const int fft_window_width, float **spectral_amp);
void band_characteristics (float *band_amp, const int num_partition, const int fft_window_width);
void save_band_characteristics (const float *band_amp, const int fft_window_width, const bool flag_complex_data, ta::FileInformation &fileinfo, const path dest_dir);

void create_spectrogram_image (float **spectral_amp, float **spectral_phase, const int fft_batch, const int fft_window_width,  const int data_range_begin_pt, const int data_range_end_pt, const path dest_dir, ta::FileInformation &fileinfo, ta::Config &config);
int gnuplot(const char *data_filepath_char, const plot_range range, ta::Config &config);
//void time_series_integration_complex (int integration_point, ta::BinData &data, int num_entire_data_pt, ta::FileInformation &fileinfo, ta::Config &config, const path dest_dir)


int main(int argc, char *argv[])
{
	bool flag_create_spectrogram_image = false;

	try{
		// Input File
		if(argc != 3){
			printf("Usage:\n\t%s <data file> <confiig file>\n", argv[0]);
			throw "Invalid Arguments";
		}
		ta::FileInformation fileinfo (argv[1], '_');

		// Observation Parameters
		ta::Config config (argv[2]);
		const bool   flag_complex_data   = true; //config.isComplexData(); 
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

		// Data Handling Parameters
		// データファイルは大容量なので、適当なデータサイズに分割して解析する。
		int tmp_num_partial_data = static_cast<int>(config.getPulsarP0_pt() * 5); // 分割データ内に少なくとも1つのパルスが入るように、分割データの点数を設定。
		if (tmp_num_partial_data > num_entire_data_pt){ // ただし全データ数がそれに満たない場合、分割データ数を全データ数に設定。
			tmp_num_partial_data = num_entire_data_pt;
		}
		const int   fft_window_width       = config.getFFTWindowWidth();
		const int   fft_batch              = std::floor(static_cast<float>(tmp_num_partial_data) / fft_window_width);
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
		printf ("- FFT Window Width    = %d pt = %f ms\n", fft_window_width, fft_window_width * sampling_interval_s * 1000);
		printf ("- FFT Time resolution = %f ms\n", config.getFFTTimeResolution_s() * 1000);
		printf ("- FFT Freq resolution = %f kHz\n", config.getFFTFreqResolution_Hz() / 1000);
		if (config.getFFTTimeResolution_s() > pulsar_w50_s){
			printf("\n### CAUTION ###\nFFT Time Resolution > Pulse Width.\nThis analysis can not resolve the pulse profile.\n");
		}
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
		const path dest_dir = path (config.getOutputDirectory()) / fileinfo.getFileStemSubstr();
		boost::filesystem::create_directories (dest_dir);
		boost::filesystem::current_path (dest_dir);
    
    //1s積分の読み出し 
    //time_series_integration_complex (20000000, data, num_entire_data_pt, fileinfo, config, dest_dir) 

		// Initialize Parameters
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

		float *gain_amp   = new float[fft_window_width];
		float *gain_phase = new float[fft_window_width];
		for(int f=0; f<fft_window_width; f++){
			gain_amp[f] = 0;
			gain_phase[f] = 0;
		}
		ta::Stats stats; // Statistic parameters including mean & variance values

		//
		// MAIN
		//
		for (int pos = 0; pos < num_partition; pos++) {
			const int partial_data_begin_pt = pos * num_partial_data;

			// Extract patial data
			if (flag_complex_data == true) {
				data.extract_binary_data_xy (h_idata, partial_data_begin_pt, num_partial_data);
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

			// Normalize the data using the mean & variance of the first partial data block (pos = 0)
			if (pos == 0) {
				stats.calcParams (h_idata, num_partial_data); // Create a Stats instance, calculating the mean & variance of the data.
			} 
			normalize_data (h_idata, num_partial_data, stats);

			// Forward FFT
			stft (h_idata, fft_window_width, fft_batch, false);

			// Create a spectrogram: spectral_amp[t][f], spectral_phase[t][f]
			clock_t spec0 = clock();
			spectrogram (spectral_amp, spectral_phase, fft_batch, fft_window_width, h_idata);
			clock_t spec1 = clock();
			ta::print_elapse ("Spectrogram Creating", spec0, spec1);

			// Derive a band characteristics gain_amp[f] from spectral_amp[t][f]
			mid_band_characteristics (gain_amp, fft_batch, fft_window_width, spectral_amp);
			
			if(pos == 0){
				//ta::saveSpectrumOfComplexData("gain amplitude at first data block.2d", gain_amp, fft_window_width);
				//ta::saveSpectrumOfComplexData("gain phase at first data block.2d", gain_phase, fft_window_width);
			}
			
			
			// Create a spectrogram image
			if(flag_create_spectrogram_image == true){
				create_spectrogram_image (spectral_amp, spectral_phase, fft_batch, fft_window_width, partial_data_begin_pt, partial_data_begin_pt + num_partial_data, dest_dir, fileinfo, config);
			}
		} // Next pos		
		
		// Band Characteristics
		band_characteristics (gain_amp, num_partition, fft_window_width);
		save_band_characteristics (gain_amp, fft_window_width, flag_complex_data, fileinfo, dest_dir);
		
/*
		//
		// Calibrate Data
		//

		float *calib_gain_amp = new float [fft_window_width];
		float *calib_gain_phase = new float [fft_window_width];
		for(int f=0; f<fft_window_width; f++){
			calib_gain_amp[f] = 0;
			calib_gain_phase[f] = 0;
		}

		for(int pos=0; pos<num_partition; pos++){

			// Load Data
			if(complex_data == true){
				loop_break_flag = loadBinaryDataComplex(partial_data_re, partial_data_im, num_partial_data, fileinfo.getFilePath().string().c_str(),  pos * datum_bit_length_B * num_partial_data, datum_bit_length_B);
			}else{
				loop_break_flag = loadBinaryData(partial_data_re, num_partial_data, fileinfo.getFilePath().string().c_str(),  pos * datum_bit_length_B * num_partial_data, datum_bit_length_B);
			}
			if(loop_break_flag != 0){ // Need to break also here.
				break;
			}

			

			// Normalize data using the mean & variance of the first data block (pos = 0)
			if(pos == 0){
				mean_re  = ta::mean(partial_data_re, num_partial_data);
				mean_im  = ta::mean(partial_data_im, num_partial_data);
				stdev_re = std::sqrt(ta::unbiased_variance(partial_data_re, num_partial_data, mean_re));
				stdev_im = std::sqrt(ta::unbiased_variance(partial_data_im, num_partial_data, mean_im));
			}	
			for(int i=0; i<num_partial_data; i++){
				partial_data_re[i] = (partial_data_re[i] - mean_re) / stdev_re;
				partial_data_im[i] = (partial_data_im[i] - mean_im) / stdev_im;
			}

			
			// Create input data on CPU
			// h_idata = input data on host CPU
			for(int i=0; i<num_partial_data; i++){
				h_idata[i].x = partial_data_re[i];
				h_idata[i].y = partial_data_im[i];
			}

			// Forward FFT
			STFT(h_idata, fft_window_width, fft_batch, false);

			// Create a spectrogram
			//h_idata[0].x = 0; h_idata[0].y = 0; // Remove the direct current (DC) component
			int index = 0; 
			const float root_N = std::sqrt(fft_window_width);
			for(int t=0; t<fft_batch; t++){
				for(int f=0; f<fft_window_width; f++){
					fft_re[f] = h_idata[index].x / root_N; 
					fft_im[f] = h_idata[index].y / root_N; // This assignment is only for readability.

					// Remove gain
					fft_re[f] /= gain_amp[f];
					fft_im[f] /= gain_amp[f];
					spectral_amp[t][f]   = std::hypot(fft_re[f], fft_im[f]);
					spectral_phase[t][f] = std::atan2(fft_im[f], fft_re[f]);

					h_idata[index].x = fft_re[f];
					h_idata[index].y = fft_im[f];
					index++;
				}
			}

			// Inverse FFT
			STFT(h_idata, fft_window_width, fft_batch, true);

			// Save Data
			
			if(pos < 5){
			unsigned int *test = new unsigned int [num_partial_data];
			for (int i=0; i<num_partial_data; i++){
				//test[i] = static_cast<unsigned int>(h_idata[i].x);
				//printf("(float %f, uint %u) ", h_idata[i].x, test[i]);
				fprintf(ofp_calib, "%.3f\t%.3f\n", h_idata[i].x / root_N, h_idata[i].y / root_N);
			}
			//fwrite(test, sizeof(float), num_partial_data, ofp_calib);
			delete [] test;
			}
			

			// Derive a band characteristics, gain_psd[f]
			for(int f=0; f<fft_window_width; f++){
				for(int t=0; t<fft_batch; t++){
					calib_gain_amp[f]   += spectral_amp[t][f] / fft_batch;
					calib_gain_phase[f] += spectral_phase[t][f] / fft_batch;
				}
				// Do not plase the division by fft_batch here.
			}
		} // Next pos	
		
		// Band Characteristics
		for(int f=0; f<fft_window_width; f++){
			calib_gain_amp[f]   /= num_partition;
			calib_gain_phase[f] /= num_partition;
		}
		calib_gain_amp[0] = 0;
		ta::smoothUsingBPF(calib_gain_amp, fft_window_width, 0, std::floor(static_cast<float>(fft_window_width)/2.0));
		if(complex_data == true){
			ta::saveSpectrumOfComplexData("calib_gain_amp.2d", calib_gain_amp, fft_window_width);		
			ta::saveSpectrumOfComplexData("calib_gain_phase.2d", calib_gain_phase, fft_window_width);
		}else{
			ta::saveData1D("calib_gain_amp.1d", calib_gain_amp, fft_window_width);
			ta::saveData1D("calib_gain_phase.1d", calib_gain_phase, fft_window_width);
		}
*/

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
		//delete [] calib_gain_amp;
		//delete [] calib_gain_phase;

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

void stft (cufftComplex *h_idata, const int fft_window_width, const int fft_batch, bool inverse)
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
}



int gnuplot(const char *data_filepath_char, const plot_range range, ta::Config &config)
{
	// Gnuplot cannot recognize '\' in a file path, and therefore '\' must be repalced to '/'.
	string data_filepath   = data_filepath_char; ta::strReplace(data_filepath, "\\", "/");
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
	
	fprintf(ofp, "set title \"STFT from #%s to #%s; Time resolution %f ms; Freq resolution %f kHz\"\n", range.xmin.c_str(), range.xmax.c_str(), config.getFFTTimeResolution_s() * 1E+3, config.getFFTFreqResolution_Hz() / 1000);
	fprintf(ofp, "set xlabel \"Time [%f ns]\"\n", config.getSamplingInterval_s() * 1E+9);
	fprintf(ofp, "set ylabel \"Frequency [%f kHz]\"\n", config.getFFTFreqResolution_Hz() / 1000);

	//fprintf(ofp, "set xrange [%d:%d]\n", 0, xend - xstart);
	fprintf(ofp, "set xrange [%s:%s]\n", range.xmin.c_str(), range.xmax.c_str());
	fprintf(ofp, "set yrange [%s:%s]\n", range.ymin.c_str(), range.ymax.c_str());
	fprintf(ofp, "set zrange [%s:%s]\n", range.zmin.c_str(), range.zmax.c_str()); // When data ~ N(0,1), spectral power ~ Exp(1/2). p-value = 28 for a significance level of 1E-6.
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

void normalize_data (cufftComplex *partial_data, const int num_partial_data, const ta::Stats &stats)
{
	const float mean_x = stats.getMean_x ();
	const float mean_y = stats.getMean_y ();
	const float inverse_stdev_x = 1.0 / stats.getStdev_x ();
	const float inverse_stdev_y = 1.0 / stats.getStdev_y ();
	for (int i=0; i<num_partial_data; i++) {
		partial_data[i].x = (partial_data[i].x - mean_x) * inverse_stdev_x;
		partial_data[i].y = (partial_data[i].y - mean_y) * inverse_stdev_y;
	}
}

void spectrogram (float **spectral_amp, float **spectral_phase, const int num_time_index, const int num_freq_index, const cufftComplex *partial_data)
{
	int index = 0; 
	for (int t = 0; t < num_time_index; t++) {
		for (int f = 0; f < num_freq_index; f++) {
			spectral_amp[t][f]   = static_cast<float>(std::hypot (partial_data[index].x, partial_data[index].y));
			spectral_phase[t][f] = static_cast<float>(std::atan2 (partial_data[index].y, partial_data[index].x));
			index++;
		}
	}
}

void mid_band_characteristics (float *band_amp, const int fft_batch, const int fft_window_width, float **spectral_amp)
{
	// Do NOT initialize band_amp[f], becuuse the parameter has been and will be integrated.
	const float inverse_fft_batch = 1.0 / fft_batch;
	for (int f = 0; f < fft_window_width; f++) {
		for (int t = 0; t < fft_batch; t++) {
			band_amp[f] += spectral_amp[t][f] * inverse_fft_batch;
			//band_phase[f] += spectral_phase[t][f] * inverse_fft_batch; No practical meaning. The spectral phase is random noise.
		}
	}
}

void band_characteristics (float *band_amp, const int num_partition, const int fft_window_width)
{
	// Derive a band characteristics
	const float inverse_num_partition = 1.0 / static_cast<float>(num_partition);
	for (int f = 0; f < fft_window_width; f++){
		band_amp[f] *= inverse_num_partition;
	}

	// Smoothe the band characteristics
	ta::smoothUsingBPF (band_amp, fft_window_width, 0, std::floor(static_cast<float>(fft_window_width)/2.0));

	// Normalize data using the minimum value
	const float min_band_amp = ta::minimum (band_amp, fft_window_width);
	const float inverse_min_gain_amp = 1.0 / min_band_amp;
	for(int f=0; f<fft_window_width; f++){
		band_amp[f] *= inverse_min_gain_amp;
	}
	band_amp[0] = 0;
}

void save_band_characteristics (const float *band_amp, const int fft_window_width, const bool flag_complex_data, ta::FileInformation &fileinfo, const path dest_dir)
{
	if (flag_complex_data == true) {
		path fout = dest_dir / (fileinfo.getFileStem().string() + ".spectral_amplitude.2d");
		ta::saveSpectrumOfComplexData(fout.string().c_str(), band_amp, fft_window_width);	
	}
	else {
		path fout = dest_dir / (fileinfo.getFileStem().string() + ".spectral_amplitude.1d");
		ta::saveData1D(fout.string().c_str(), band_amp, fft_window_width);
	}
}

void create_spectrogram_image (float **spectral_amp, float **spectral_phase, const int fft_batch, const int fft_window_width,  const int data_range_begin_pt, const int data_range_end_pt, const path dest_dir, ta::FileInformation &fileinfo, ta::Config &config)
{
	char str[STR_LEN_MAX + 1];

	// Save the spectrogram as a text file
	const int data_range_begin_B = config.getDatumBitLength_B() * data_range_begin_pt;
	const int data_range_end_B   = config.getDatumBitLength_B() * data_range_end_pt;

	sprintf (str, "%s.%d-%dB.spectrogram.3d", fileinfo.getFileStem().string().c_str(), data_range_begin_B, data_range_end_B);
	const path output_specamp_filepath   = dest_dir / str; // Specify the destination directory and avoid depending on the current directory for Gnuplot configuration.
	FILE *ofp_spec = fopen(output_specamp_filepath.string().c_str(), "w"); if(ofp_spec == NULL){throw ta::messageFileOpenError(output_specamp_filepath.string());}

	sprintf (str, "%s.%d-%dB.phase.3d", fileinfo.getFileStem().string().c_str(), data_range_begin_B, data_range_end_B);
	const path output_specphase_filepath = dest_dir / str; // Specify the destination directory and avoid depending on the current directory for Gnuplot configuration.
	FILE *ofp_phase = fopen(output_specphase_filepath.string().c_str(), "w"); if(ofp_phase == NULL){throw ta::messageFileOpenError(output_specphase_filepath.string());}

	for(int t=0; t<fft_batch; t++){
		for(int f=0; f<fft_window_width; f++){
			fprintf(ofp_spec,  "%d\t%d\t%.3f\n", data_range_begin_pt + t * fft_window_width, f, spectral_amp[t][f] * spectral_amp[t][f]);
			fprintf(ofp_phase, "%d\t%d\t%.2f\n", data_range_begin_pt + t * fft_window_width, f, spectral_phase[t][f]);
		}
		fprintf(ofp_spec, "\n");
		fprintf(ofp_phase, "\n");
	}
	fclose(ofp_spec);
	fclose(ofp_phase);


	// Make an PNG image of the spectrogram
	plot_range range;
	range.xmin = data_range_begin_pt;
	range.xmax = data_range_end_pt;
	range.ymin = "0";
	range.ymax = fft_window_width;
	range.zmin = "0";
	range.zmax = "";
	gnuplot (output_specamp_filepath.string().c_str(), range, config);

	range.zmin = "- 4";
	range.zmax = "4";
	gnuplot (output_specphase_filepath.string().c_str(), range, config);
}

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

/**
void time_series_integration_complex (int integration_point, ta::BinData &data, int num_entire_data_pt, ta::FileInformation &fileinfo, ta::Config &config, const path dest_dir)
{
  
  int tmp_num_partial_data = 100000000 //適当にセット(必要ない。寧ろintegrate_window_width軽くした方がよさげ)
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
**/
