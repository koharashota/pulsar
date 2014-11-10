/*******************************************************************************
�g����

	#include "Data.h"

	// �f�[�^�C���X�^���X�����Bdata�Ƀf�[�^���i�[�����B
	Aoki::Data data;
	data.load_s("./single.dat"); // �f�[�^��1��̏ꍇ
�܂���
	data.load_d("./double.dat"); // �f�[�^��2��̏ꍇ
�܂���
	data.load_binary_s("./single.bin", 8); // �f�[�^�� <x0><x1><x2>... �� 8 bit ������ł���`���̏ꍇ
�܂���
	data.load_binary_d("./double.bin", 8); // �f�[�^�� <x0><y0><x1><y1><x2><y2>... �� 8 bit ������ł���`���̏ꍇ

	// �f�[�^���𒲂ׂ�
	int num_data = data.getNumData(); // �Ⴆ��86400�̂悤�Ȓl���Ԃ�

	// �f�[�^�𒊏o
	double x[500];
	data.extract_data_s( x, 1000, 500 ); // �f�[�^�t�@�C���� 1000 �Ԃ��� 1499 �Ԃ܂ł̃f�[�^�𒊏o�� data[] �Ɋi�[
�܂���
	double t[500];
	double x[500];
	data.extract_data_d( t, x, 1000, 500 );

*******************************************************************************/
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

			// �f�[�^�t�@�C����ǂݍ���
			int load_binary_as_single   (const string file, const int data_length_pt, const size_t datum_size_B);
			int load_binary_as_double   (const string file, const int data_length_pt, const size_t datum_size_B);
			int load_binary_as_double_2 (const string file, const int data_length_pt, const size_t datum_size_B); // ��������̃f�[�^�z�u���Ⴄ�o�[�W����

			// �f�[�^��Ԃ�
			int extract_binary_data_xy   (cufftComplex *cudata, const int extraction_first_point, const int extraction_width) const;
			int extract_binary_data_xy_2 (cufftComplex *cudata, const int extraction_first_point, const int extraction_width) const;
		
			// �f�[�^����Ԃ�
			int getNumData() const {return num_data;};

	};
}
#endif
