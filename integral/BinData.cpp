#include "BinData.h"

namespace ta
{

	BinData::BinData()
	{
		load_s_executed = false;
		load_d_executed = false;
	}

	BinData::~BinData()
	{
	}

	int BinData::load_binary_as_single (const string file, const int data_length_pt, const size_t datum_size_B)
	{
		num_data = data_length_pt;
		bin_data_x = new unsigned char[num_data];
		
		for (int i = 0; i < num_data; i++){
			bin_data_x[i] = 0.0;
		}

		FILE *ifp = fopen(file.c_str(), "rb");
		if(ifp == NULL){
			std::string err = "File Open Error: " + file;
			throw err;
		}

		int error_flag = 1;
		unsigned char unsigned_x;
		fseek (ifp, 0, SEEK_SET);
		for (int i = 0; i < num_data; i++) {
			error_flag = fread(&unsigned_x, datum_size_B, 1, ifp); if(error_flag == 0) break;
			bin_data_x[i] = unsigned_x;
		}
		fclose(ifp);

		if(error_flag == 0){
			return -1;
		}
		load_s_executed =  true;
		return 0;
	}

	int BinData::load_binary_as_double_2 (const string file, const int data_length_pt, const size_t datum_size_B)
	{
		num_data = data_length_pt;
		FILE *ifp = fopen(file.c_str(), "rb");
		if(ifp == NULL){
			std::string err = "File Open Error: " + file;
			throw err;
		}

		
		bin_data_x = new unsigned char[num_data];
		bin_data_y = new unsigned char[num_data];
		for (int i = 0; i < num_data; i++){
			bin_data_x[i] = '0';
			bin_data_y[i] = '0';
		}


		int error_flag = 1;
		unsigned char unsigned_x, unsigned_y;
		fseek (ifp, 0, SEEK_SET);
		for (int i = 0; i < num_data; i++) {
			error_flag = fread(&unsigned_x, datum_size_B, 1, ifp); if(error_flag == 0) break;
			error_flag = fread(&unsigned_y, datum_size_B, 1, ifp); if(error_flag == 0) break;
			bin_data_x[i] = unsigned_x;
			bin_data_y[i] = unsigned_y;
		}
		fclose(ifp);

		if(error_flag == 0){
			return -1;
		}
		load_d_executed =  true;
		return 0;
	}

	int BinData::load_binary_as_double (const string file, const int data_length_pt, const size_t datum_size_B)
	{
		num_data = data_length_pt;
		data = new data_xy [num_data];

		FILE *ifp = fopen(file.c_str(), "rb");
		if(ifp == NULL){
			std::string err = "File Open Error: " + file;
			throw err;
		}

		for (int i = 0; i < num_data; i++){
			data[i].x = '0';
			data[i].y = '0';
		}

		fseek (ifp, 0, SEEK_SET);
		num_data = fread (data, sizeof(unsigned char) * 2, num_data, ifp);
		if (num_data == 0) {
			std::string err = "No data in " + file;
			throw err;
		}
		else if (num_data != data_length_pt){
			fprintf (stderr, "Number of data points = %d, which is less than the specified data length.\n", num_data);
		}
		fclose(ifp);

		load_d_executed =  true;
		return 0;
	}


	int BinData::extract_binary_data_xy_2 (cufftComplex *cudata, const int extraction_first_point, const int extraction_width) const
	{
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
				cudata[i].x = static_cast<float>(static_cast<unsigned>(bin_data_x[extraction_first_point + i]));
				cudata[i].y = static_cast<float>(static_cast<unsigned>(bin_data_y[extraction_first_point + i]));
			}

		}else{
			for(int i=0; i<extraction_width; i++){
				cudata[i].x = static_cast<float>(static_cast<unsigned>(bin_data_x[extraction_first_point + i])); // start から start+width-1 ポイントのデータを抽出
				cudata[i].y = static_cast<float>(static_cast<unsigned>(bin_data_y[extraction_first_point + i])); // start から start+width-1 ポイントのデータを抽出
			}
		}

		return 0;
	}

	int BinData::extract_binary_data_xy(cufftComplex *cudata, const int extraction_first_point, const int extraction_width) const
	{
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
				cudata[i].x = static_cast<float>(static_cast<unsigned>(data[extraction_first_point + i].x));
				cudata[i].y = static_cast<float>(static_cast<unsigned>(data[extraction_first_point + i].y));
			}

		}else{
			for(int i=0; i<extraction_width; i++){
				cudata[i].x = static_cast<float>(static_cast<unsigned>(data[extraction_first_point + i].x)); // start から start+width-1 ポイントのデータを抽出
				cudata[i].y = static_cast<float>(static_cast<unsigned>(data[extraction_first_point + i].y)); // start から start+width-1 ポイントのデータを抽出
			}
		}

		return 0;
	}
//------------------------------------------------------------------------------
}
