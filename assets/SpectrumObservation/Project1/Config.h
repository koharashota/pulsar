#ifndef TA_CONFIG
#define TA_CONFIG

#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <string>

namespace ta
{
class Config
{
private:

	int   datum_bit_length_B;
	float sampling_freq_Hz;
	float sampling_interval_s;
	bool  complex_data;

	// Integration Config
	float integration_time_s;

	// FFT Config
	float pulsar_p0_s;   // second, Barycentric period of the pulsar
	float pulsar_w50_s;  // second, Width of pulse at 50% of peak.
	int   pulsar_p0_pt;  // data point, depending on the sampling interval
	int   pulsar_w50_pt; // data point
	//float pulsar_DM_pccm3;

	// FFT Config
	int   fft_window_width;
	float fft_time_resolution_s;
	float fft_freq_resolution_Hz;

	// Output Config
	std::string output_dir;

	void initialize();

public:
	Config();
	Config(const std::string filepath);
	~Config();
	int    loadConfig(const std::string filepath);

	int   getDatumBitLength_B();
	float getSamplingFrequency_Hz();
	float getSamplingInterval_s();
	bool  isComplexData()             {return complex_data;};

	float getIntegrationTime_s();
	float getPulsarP0_s();
	float getPulsarW50_s();
	int   getPulsarP0_pt();
	int   getPulsarW50_pt();

	//float getPulsarDM
	int   getFFTWindowWidth();
	float getFFTTimeResolution_s();
	float getFFTFreqResolution_Hz();
	std::string getOutputDirectory();

};
}
#endif
