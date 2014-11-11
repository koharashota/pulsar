#include "Config.h"

namespace Aoki
{

	Config::Config()
	{
		// Obs Param
		datum_bit_length_B  = -1;
		sampling_freq_Hz    = -1;
		sampling_interval_s = -1;
		beam_width_deg      = -1;
		complex_data        = false;

		data_partitioning_width_approx_s = -1;

		// Target Param
		pulsar_p0_s  = -1;
		pulsar_w50_s = -1;

		// Analysis Param
		integration_time_s = -1;
		fft_window_width   = -1;

		output_dir = "";
	}

	Config::Config(const std::string filepath)
	{
		if(loadConfig(filepath) != 0){
			throw "Configuration Error";
		}
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
			const int divpos = buf.find_first_of(" \t");
			if(divpos <= 0){
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
				sampling_freq_Hz = atof(value.c_str()) * 1E+6;
				if(sampling_freq_Hz <= 0){
					fprintf(stderr, "Invalid Config: sampling_frequency_MHz = %f\n", sampling_freq_Hz / 1E+6);
					return 1;
				}
				sampling_interval_s = 1 / sampling_freq_Hz;
			}
			else if (param == "beam_width_arcmin") {
				beam_width_deg = atof (value.c_str()) / 60;
				if (beam_width_deg <= 0) {
					fprintf (stderr, "Invalid Config: beam_width_arcmin = %f\n", beam_width_deg * 60);
					return 1;
				}
			}
			else if (param == "data_partitioning_width_approx_s") {
				data_partitioning_width_approx_s = atof (value.c_str());
				if (data_partitioning_width_approx_s <= 0) {
					fprintf (stderr, "Invalid Config: data_partitioning_width_approx_s = %f\n", data_partitioning_width_approx_s);
					return 1;
				}
			}
			
			else if (param == "comeplex_data") {
				if(value == "true"){
					complex_data = true;
				}else{
					complex_data = false;
				}
			}
			else if(param == "integration_time_s"){
				integration_time_s = atof(value.c_str());
				if(integration_time_s < 0){
					fprintf(stderr, "Invalid Config: integration_time_s = %f\n", integration_time_s);
					return 1;
				}

			}
			else if(param == "pulsar_p0_s"){
				pulsar_p0_s = atof(value.c_str());
				if(pulsar_p0_s <= 0){
					fprintf(stderr, "Invalid Config: pulsar_p0_s = %f\n", pulsar_p0_s);
					return 1;
				}
			}
			else if(param == "pulsar_w50_ms"){
				pulsar_w50_s = atof(value.c_str()) / 1000;
				if(pulsar_w50_s <= 0){
					fprintf(stderr, "Invalid Config: pulsar_w50_ms = %f\n", pulsar_w50_s * 1000);
					return 1;
				}
			}
			else if(param == "fft_window_width"){
				fft_window_width = atoi(value.c_str());
				if(fft_window_width < 0){
					fprintf(stderr, "Invalid Config: fft_window_width = %d\n", fft_window_width);
					return 1;
				}
			}
			else if(param == "output_directory"){
				for(size_t i=0; i<value.length(); i++){
					const int erasepos = value.find_first_of("\t\"");
					if(erasepos < 0) break;
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
}
