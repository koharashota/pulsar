#include "Config.h"

namespace ta
{

	void Config::initialize()
	{
		datum_bit_length_B = 0;
		sampling_freq_Hz = 0;
		sampling_interval_s = 0;
		complex_data = false;

		integration_time_s = 0;

		pulsar_p0_s = 0;
		pulsar_w50_s = 0;
		pulsar_p0_pt = 0; 
		pulsar_w50_pt = 0;
		//pulsar_DM_pccm3;

		fft_window_width = 0;
		fft_time_resolution_s = 0;
		fft_freq_resolution_Hz = 0;

		std::string output_dir = "";	
	}

	Config::Config()
	{
		//initialize();
	}

	Config::Config(const std::string filepath)
	{
		//initialize();
		if (loadConfig(filepath) != 0){
			throw "Configuration Error";
		}

		pulsar_p0_pt  = static_cast<int>(pulsar_p0_s  / sampling_interval_s);
		pulsar_w50_pt = static_cast<int>(pulsar_w50_s / sampling_interval_s);
		fft_time_resolution_s  = sampling_interval_s * fft_window_width;
		fft_freq_resolution_Hz = static_cast<float>(1.0 / fft_time_resolution_s);
	}

	Config::~Config(){}



	int Config::loadConfig(const std::string filepath)
	{
		std::ifstream fin(filepath.c_str());
		if(!fin){
			fprintf(stderr, "File Open Error: %s\n", filepath.c_str());
			return 1;
		}

		std::string buf;
		while(std::getline(fin, buf)){

			// コメント文は読み飛ばし
			if(buf.substr(0,1) == "#") continue;
		
			// 空白文字で変数と値を分割
			const size_t divpos = buf.find_first_of(" \t");
			if(divpos == 0 || divpos == std::string::npos){
				fprintf(stderr, "Invalid Config Format in %s:\n  %s", filepath.c_str(), buf.c_str());
				return 1;
			}

			std::string param = buf.substr(0, divpos);
			std::string value = buf.substr(divpos + 1, buf.length());

			// 設定項目はここに追加していく
			// 追加方法: if (param == "設定項目名")
			if(param == "datum_bit_length_B"){
				datum_bit_length_B = atoi(value.c_str());
				if(datum_bit_length_B <= 0){
					fprintf(stderr, "Invalid Config: datum_bit_length_B = %d\n", datum_bit_length_B);
					return 1;
				}
			}
			else if(param == "sampling_frequency_MHz"){
				sampling_freq_Hz = static_cast<float>(atof(value.c_str()) * 1E+6);
				if(sampling_freq_Hz <= 0){
					fprintf(stderr, "Invalid Config: sampling_frequency_MHz = %f\n", sampling_freq_Hz / 1E+6);
					return 1;
				}
				sampling_interval_s = 1 / sampling_freq_Hz;
			}
			else if(param == "comeplex_data"){
				if(value == "true"){
					complex_data = true;
				}else{
					complex_data = false;
				}
			}
			else if(param == "integration_time_s"){
				integration_time_s = static_cast<float>(atof(value.c_str()));
				if(integration_time_s < 0){
					fprintf(stderr, "Invalid Config: integration_time_s = %f\n", integration_time_s);
					return 1;
				}

			}
			else if(param == "pulsar_p0_s"){
				pulsar_p0_s = static_cast<float>(atof(value.c_str()));
				if(pulsar_p0_s <= 0){
					fprintf(stderr, "Invalid Config: pulsar_p0_s = %f\n", pulsar_p0_s);
					return 1;
				}
			}
			else if(param == "pulsar_w50_ms"){
				pulsar_w50_s = static_cast<float>(atof(value.c_str()) / 1000);
				if(pulsar_w50_s <= 0){
					fprintf(stderr, "Invalid Config: pulsar_w50_ms = %f\n", pulsar_w50_s * 1000);
					return 1;
				}
			}
			else if(param == "fft_window_width"){
				fft_window_width = atoi(value.c_str());
				if(fft_window_width < 4){
					fprintf(stderr, "Invalid Config: fft_window_width = %d\n", fft_window_width);
					return 1;
				}
			}
			else if(param == "output_directory"){
				for(size_t i=0; i<value.length(); i++){
					const size_t erasepos = value.find_first_of("\t\"");
					if(erasepos == std::string::npos) break;
					value.erase(erasepos, 1);
				}
				output_dir = value;
			}
			else{
				fprintf(stderr, "Skipping Invalid Config Entry: %s\n", param.c_str());
			}
		}
		return 0;
	}

	int Config::getDatumBitLength_B()
	{
		if (datum_bit_length_B < 1) {
			throw "Config Error";
		}
		return datum_bit_length_B;
	}


	float Config::getSamplingFrequency_Hz()
	{
		if (sampling_freq_Hz <= 0) {
			throw "config error";
		}
		return sampling_freq_Hz;
	}

	float Config::getSamplingInterval_s()
	{
		if (sampling_interval_s <= 0) {
			throw "config error";
		}
		return sampling_interval_s;
	}

	float Config::getIntegrationTime_s()
	{
		if (integration_time_s < 0) {
			throw "config error";
		}
		return integration_time_s;
	}

	float Config::getPulsarP0_s()
	{
		if (pulsar_p0_s <= 0) {
			throw "config error";
		}
		return pulsar_p0_s;
	}

	float Config::getPulsarW50_s()
	{
		if (pulsar_w50_s <= 0) {
			throw "config error";
		}
		return pulsar_w50_s;
	}

	int Config::getPulsarP0_pt()
	{
		if (pulsar_p0_pt <= 0) {
			throw "config error";
		}
		return pulsar_p0_pt;
	}

	int Config::getPulsarW50_pt()
	{
		if (pulsar_w50_pt <= 0) {
			throw "config error";
		}
		return pulsar_w50_pt;
	}

	//float getPulsarDM

	int Config::getFFTWindowWidth()
	{
		if (fft_window_width <= 0) {
			throw "config error";
		}
		return fft_window_width;
	}

	float Config::getFFTTimeResolution_s()
	{
		if (fft_time_resolution_s <= 0) {
			throw "config error";
		}
		return fft_time_resolution_s;
	}

	float Config::getFFTFreqResolution_Hz()
	{
		if (fft_freq_resolution_Hz <= 0) {
			throw "CONFIG ERROR";
		}
		return fft_freq_resolution_Hz;
	}

	std::string Config::getOutputDirectory()
	{
		if (output_dir == "") {
			throw "config error";
		}
		return output_dir;
	}
}
