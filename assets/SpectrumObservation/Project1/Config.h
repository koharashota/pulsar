#ifndef AOKI_CONFIG
#define AOKI_CONFIG

#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <string>

namespace Aoki
{
class Config
{
private:

	int   datum_bit_length_B;
	float sampling_freq_Hz;
	float sampling_interval_s;
	float beam_width_deg;
	bool  complex_data;

	float data_partitioning_width_approx_s;

	// Integration Config
	float integration_time_s;

	// FFT Config
	float pulsar_p0_s;
	float pulsar_w50_s;
	//double pulsar_DM_pccm3;
	int    fft_window_width;

	// Output Config
	std::string output_dir;

public:
	Config();
	Config(const std::string filepath);
	~Config();
	int    loadConfig(const std::string filepath);

	int   getDatumBitLength_B()       {return datum_bit_length_B;};
	float getSamplingFrequency_Hz()   {return sampling_freq_Hz;};
	float getSamplingInterval_s()     {return sampling_interval_s;};
	float getBeamWidth_deg()          {return beam_width_deg;};
	bool  isComplexData()             {return complex_data;};

	float getDataPartitioningWdithApprox_s() {return data_partitioning_width_approx_s;};

	float getIntegrationTime_s()      {return integration_time_s;};
	float getPulsarP0_s()             {return pulsar_p0_s;};
	float getPulsarW50_s()            {return pulsar_w50_s;};
	//double getPulsarDM
	int    getFFTWindowWidth()         {return fft_window_width;};
	std::string getOutputDirectory()   {return output_dir;};
	
};
}
#endif
