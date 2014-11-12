#include <string>
#include <fstream>
#include <cctype>
#include <cuda.h>
#include <cufft.h>

#ifndef TA_BINDATA
#define TA_BINDATA

using std::string;

namespace ta
{
	class BinData
	{
		private:
			int  num_data;   // Number of data samples, which output by getNumData().

			unsigned char *bin_data_x; // 1 series binary data
			unsigned char *bin_data_y; // 1 series binary data

			typedef struct {
				unsigned char x;
				unsigned char y;
			} data_xy;
			data_xy *data;
			
			bool load_s_executed;
			bool load_d_executed;

		public:
			BinData();
			~BinData();

			// データファイルを読み込む
			int load_binary_as_single   (const string file, const int data_length_pt, const size_t datum_size_B);
			int load_binary_as_double   (const string file, const int data_length_pt, const size_t datum_size_B);
			int load_binary_as_double_2 (const string file, const int data_length_pt, const size_t datum_size_B); // メモリ上のデータ配置が違うバージョン

			// データを返す
			int extract_binary_data_xy   (cufftComplex *cudata, const int extraction_first_point, const int extraction_width) const;
			int extract_binary_data_xy_2 (cufftComplex *cudata, const int extraction_first_point, const int extraction_width) const;
		
			// データ数を返す
			int getNumData() const {return num_data;};

	};
}
#endif
