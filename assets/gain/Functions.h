#ifndef TA_FUNCTIONS
#define TA_FUNCTIONS

#include <cmath>
#include <ctime>
#include <string>
#include <iostream>
#include <fstream>
//#include <gsl/gsl_matrix.h>

namespace ta
{
	/* �~���� */
	static const double pi     = 6.0 * std::asin(0.5);
	static const double twopi  = 2.0 * pi;
	static const double halfpi = pi / 2.0;

	/* ��{�֐� */
	int bitReversal(int N, int L);
	void swap(double *x, double *y);

	/* �f�[�^���o */
	void extract_data(double *extracted_data, const double *original_data, const int original_data_num, const int extraction_start_point, const int extraction_width);

	/* �e��֐� */
	double sinc(double x);
	double rayleigh(const double x, const double sigma);
	double theta_time_converter(const double time, const double declination);
	double antenna_radiation_pattern(const double theta, const double wave_length, const double aperture_radius);
	double fringe(const double theta, const double wave_length, const double aperture_radius, const double baseline_length);

	/* �M������ */
	void windowFunction(double *x, int N, int id);
	void autocorrelation(double *R, double *x, int N);
	void movingAverage(double *ma, double *x, int datNum, int m);
	void fft(double *xr, double *xi, int L, int id);
	int  dft(double *xr, double *xi, const int N, const bool inverse);
	int  dft(float  *xr, float  *xi, const int N, const bool inverse);
	int  bpf(double *xr, double *xi, const int width, const int inf_freq_channel, const int sup_freq_channel);
	int  bpf(float  *xr, float  *xi, const int width, const int inf_freq_channel, const int sup_freq_channel);
	int  smoothUsingBPF(double *data, const int width, const int inf_freq_channel, const int sup_freq_channel);
	int  smoothUsingBPF(float  *data, const int width, const int inf_freq_channel, const int sup_freq_channel);
	void normalize_data(float  *data, const unsigned int data_size);
	void normailze_data(double *data, const unsigned int data_size);
	//void create_spectrogram(const float *dft_re, const float *dft_im, const int 

	/* ���v */
	float  summation (const float *data, const int num_data);
	float  minimum(const float *data, const int num_data);
	double mean(const double *data, const int data_num);
	float  mean(const float  *data, const int data_num);
	double unbiased_variance(const double *data, const int data_num, const double mean); 
	float  unbiased_variance(const float  *data, const int data_num, const float  mean);
	double sample_variance(const double *data, const int data_num, const double mean);
	double sample_covariance(const double *x, const double *y, const int data_num, const double x_average, const double y_average);

	double gradient_of_regression_line(const double *x, const double *y, const int data_num);
	double intercept_of_regression_line(const double *x, const double *y, const int data_num);

	/* ���� ************************************************************************
	#include <time>
	time_t timer = atot("19991231000000", 0);
	time_t timer2 = ttot(1999, 12, 31, 0, 0, 0, 0);	
	*******************************************************************************/
	time_t atot(const char *local_YYYYMMDDhhmmss, const int isdst);
	time_t ttot(const int local_year, const int local_month, const int local_day, const int local_hour, const int local_minute, const int local_second, const int isdst);
	void ttoa(char* str, const time_t timer);

	/* �s�� ************************************************************************
	3�����̉�]�s��𐶐�����
	#include <gsl/gsl_matrix.h>
	double theta = 3.14;
	gsl_matrix * Rx = gsl_matrix_alloc(3, 3);
	gsl_matrix * Ry = gsl_matrix_alloc(3, 3);
	gsl_matrix * Rz = gsl_matrix_alloc(3, 3);
	Aoki::createRotationMatrix_3D(Rx, theta, 'x');
	Aoki::createRotationMatrix_3D(Ry, theta, 'y');
	Aoki::createRotationMatrix_3D(Rz, theta, 'z');
	gsl_matrix_free(Rx);
	gsl_matrix_free(Ry);
	gsl_matrix_free(Rz);
	*******************************************************************************/
//	void createRotationMatrix_3D(gsl_matrix * matrix, const double angle, const char axis);

	/* �p�x �� �� 0 �� �� �� 2�� �ɔ[�߂� */
	double rangeAngle(const double angle_radian);

	/* �p�x�̒P�ʂ� radian ���� hms/dms �ɕϊ����� */
	void rad_to_hms(double* hms, const double radian); /* radian ���� hour-minute-second �ɕϊ� */
	void rad_to_dms(double* dms, const double radian); /* radian ���� deg-arcmin-arcsec �ɕϊ� */
	double hms_to_rad(const double hour, const double minute, const double second);         /* hour-minute-second ���� radian �ɕϊ� */
	double dms_to_rad(const double degree, const double arcminute, const double arcsecond); /* deg-arcmin-arcsec ���� radian �ɕϊ� */

	/* �����񑀍� */
	void strReplace (std::string& str, const std::string& deleted_str, const std::string& added_str);

	/* ���o�� */
	size_t read_stdin(char * str, const size_t STRLENMAX);
	bool stdin_yes_or_no();
	int saveSpectrumOfComplexData (const char *filename, const float *data, const int data_length);
	int saveData1D (const char *filename, const float *data, const int data_length);
	int saveData2D (const char *filename, const float *data, const int data_length, const float x_resolution);
	int saveDataForGnuplot3D (const char *filename, float **data, const int x_length, const float x_unit, const int y_length, const float y_unit);

	/* ���b�Z�[�W���� */
	std::string messageFileOpenError(const std::string filename);
	void print_elapse (std::string str, clock_t start_time, clock_t end_time);

}
#endif
