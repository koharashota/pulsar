#include "Functions.h"

namespace ta
{
	//------------------------------------------------------------------------------
	// ���� N �̃r�b�g���t�]���������� brN �����߂�v���O����
	//------------------------------------------------------------------------------
	int bitReversal(int N, int L) // �r�b�g���t�]������������ N �ƁA���̃r�b�g�� L
	{
		int n1, n2;
		int brN; // ���� N �̃r�b�g���t�]���������� bit reversal N

		n1 = N;
		brN = 0;

		// �P�� L ��J��Ԃ�
		for(int i=0; i<L; i++){
			n2 = int(n1 / 2);
			brN = 2*brN + n1 - 2*n2;
			n1 = n2;
		}

		return brN;
	}

	//------------------------------------------------------------------------------
	// ���l x, y ����������v���O����
	//------------------------------------------------------------------------------
	void swap(double *x, double *y)
	{
		double temp;
		temp = *x;
		*x = *y;
		*y = temp;
	}

	//------------------------------------------------------------------------------
	// Window Function
	//------------------------------------------------------------------------------
	void windowFunction(double *x, int N, int id)
	{
		// x[i]   : ���̎����f�[�^�i�����f�[�^������ꍇ�̓v���O�����������K�v����j
		//          ���̌��f�[�^�͎�����
		// N      : �f�[�^��

		// id = 1 : �n�~���O��
		// id = 2 : �n�j���O��
		// id = 3 : �u���b�N�}����

		double omega = 2 * 4*std::atan(1) / (N-1); // 2��/(N-1) �͌v�Z���Ă���

		// Hamming
		if(id == 1){
			for(int i=0; i<N; i++){
				x[i] = x[i] * ( 0.54 - 0.46 * std::cos(omega*i) );
			}
		}
		// Hann
		else if(id == 2){
			for(int i=0; i<N; i++){
				x[i] = x[i] * ( 0.5 - 0.5 * std::cos(omega*i) );
			}
		}
		// Blackman
		else if(id == 3){
			for(int i=0; i<N; i++){
				x[i] = x[i] * ( 0.42 - 0.5*std::cos(omega*i) + 0.08*cos(2*omega*i) );
			}
		}
		// Gaussian
		else if(id == 4){
			double mu    = N / 2.0;
			double sigma = N / 4.0; // N = 4�� �ƃЂ��`����
			double index;

			for(int i=0; i<N; i++){
				index = (x[i] - mu) / sigma;
				index = index * index;
				index = - index / 2.0;
				x[i] = x[i] * std::exp(index);
			}
		}

	}

	//------------------------------------------------------------------------------
	// SINC Function
	//------------------------------------------------------------------------------
	double sinc(double x)
	{
		double value;

		if(x == 0.0){
			value = 1.0;
		}else{
			value = std::sin(x) / x;
		}

		return value;
	}

	//------------------------------------------------------------------------------
	// Rayleigh Distribution Function
	//------------------------------------------------------------------------------
	double rayleigh(const double x, const double sigma)
	{
		double sigma_squared = sigma * sigma;
		double f = x * std::exp( - x * x / 2.0 / sigma_squared) / sigma_squared;
		return f;
	}

	//------------------------------------------------------------------------------
	// theta-time converter
	//------------------------------------------------------------------------------
	double theta_time_converter(const double time, const double declination)
	{
		const double omega = 2 * pi / 86164.091;
		return omega * time * std::cos(declination);
	}

	//------------------------------------------------------------------------------
	// single antenna radiation pattern
	//------------------------------------------------------------------------------
	double antenna_radiation_pattern(const double theta, const double wave_length, const double aperture_radius)
	{
		const double K = 2 * pi / wave_length;
		const double A = aperture_radius;

		double kasin = K * A * std::sin(theta);
		double beam = sinc(kasin);

		return beam * beam;
	}

	//------------------------------------------------------------------------------
	// Fringe Function
	//------------------------------------------------------------------------------

	double fringe(const double theta, const double wave_length, const double aperture_radius, const double baseline_length)
	{
		const double K = 2 * pi / wave_length;
		const double A = aperture_radius;
		const double D = baseline_length;

		double kasin = K * A * std::sin(theta);
		double kDsin = K * D * std::sin(theta);
		double beam = sinc(kasin);
		beam = beam * beam;

		return beam * std::cos(kDsin);
	}

	//------------------------------------------------------------------------------
	// Moving Average
	//------------------------------------------------------------------------------
	void movingAverage(double *ma, double *x, int datNum, int m)
	{
		/* 2���E3���������K���ړ����� (Savitzky-Golay�@) ���v�Z
		 ma[i] = 1/W sigma_{j=-m}^{m} w[j] x[i+j]
		 m      : �������_�� 2m+1
		 datNum : ���f�[�^��
		 �������\�͈� i = m+1, ... , datNum-m-1 �̌v datNum-2m ��
		 */

		// �d�݌W�� w ��ݒ�
		double *w = new double[2*m+1];
		double W = (4*m*m - 1) * (2*m + 3) / 3;
		for(int jj=0; jj <= 2*m; jj++){ // j = -m ,., m �̑���� jj = j + m = 0 ,., 2m �Ƃ���
			w[jj] = 3*m*(m+1) - 1 - 5*(jj-m)*(jj-m);
			w[jj] = w[jj] / W;
		}

		// i = 0 ,., m-1 ��m�̃f�[�^�͕������ł��Ȃ��̂Ő��f�[�^���R�s�[
		for(int i=0; i<m; i++){
			ma[i] = x[i];
		}

		// i = datNum-m ,., datNum-1 ��m�̃f�[�^�͕������ł��Ȃ��̂Ő��f�[�^���R�s�[
		for(int i=datNum-m; i<datNum; i++){
			ma[i] = x[i];
		}

		// i = m ,., datNum-m-1 ��datNum-2m�̃f�[�^�𕽊���
		for(int i=m; i<datNum-m; i++){
			// �����ł����� ma[i] �͏��������Ȃ�����ˁI
			ma[i] = 0;
			for(int jj=0; jj<= 2*m; jj++){ // j = -m ,., m �̑���� jj = j + m = 0 ,., 2m �Ƃ���
				ma[i] = ma[i] + ( w[jj] * x[i + jj - m] );
			}
		}

		// w[j] �͂����g����
		delete [] w;

	}

	//------------------------------------------------------------------------------
	// Discrete Fourier Transform
	//------------------------------------------------------------------------------
	int dft(double *xr, double *xi, const int N, const bool inverse)
	{
		// xr : ���̃f�[�^�̎����i�����Ɍ��ʂ�����B���̃f�[�^�͎�����j
		// xi : ���̃f�[�^�̋����i�����Ɍ��ʂ�����B���̃f�[�^�͎�����j
		// N  : �f�[�^��

		// ���� Xr Xi �͕K�v�B�A���S���Y���� xr xi �ɒ��ڌ��ʂ������邱�Ƃ͂ł��Ȃ��B
		double *Xr = new double[N]; // �t�[���G�ϊ��̎���
		double *Xi = new double[N]; // �t�[���G�ϊ��̋���
		if(!Xr || !Xi){
			fprintf(stderr, "Memory Allocation Error\n");
			return -1;
		}
		
		double omega = 2 * 6*asin(0.5) / N; // 2��/N �� omega �Ƃ���
		double root_N = std::sqrt(N);
		double cos, sin;
		
		if(inverse == false){
			// DFT
			for(int k=0; k<N; k++){
				Xr[k] = 0;
				Xi[k] = 0;
				for(int j=0; j<N; j++){
					cos = std::cos(omega*j*k);
					sin = std::sin(omega*j*k);
					Xr[k] += xr[j] * cos + xi[j] * sin;
					Xi[k] += - xr[j] * sin + xi[j] * cos;
				}
				Xr[k] /= root_N;
				Xi[k] /= root_N;
			}
		}else{
			// Inverse DFT
			for(int k=0; k<N; k++){
				Xr[k] = 0;
				Xi[k] = 0;
				for(int j=0; j<N; j++){
					cos = std::cos(omega*j*k);
					sin = std::sin(omega*j*k);
					Xr[k] += xr[j] * cos - xi[j] * sin;
					Xi[k] += xr[j] * sin + xi[j] * cos;
				}
				Xr[k] /= root_N;
				Xi[k] /= root_N;
			}

		}

		// x �� X ������ X ������
		for(int k=0; k<N; k++){
			xr[k] = Xr[k];
			xi[k] = Xi[k];
		}
		delete [] Xr;
		delete [] Xi;
		return 0;
	}

	int dft(float *xr, float *xi, const int N, const bool inverse)
	{
		float *Xr = new float[N]; // �t�[���G�ϊ��̎���
		float *Xi = new float[N]; // �t�[���G�ϊ��̋���
		if(!Xr || !Xi){
			fprintf(stderr, "Memory Allocation Error\n");
			return -1;
		}
		
		float omega = 2 * 6*asin(0.5) / N; // 2��/N �� omega �Ƃ���
		float root_N = std::sqrt(N);
		float cos, sin;
		
		if(inverse == false){
			// DFT
			for(int k=0; k<N; k++){
				Xr[k] = 0;
				Xi[k] = 0;
				for(int j=0; j<N; j++){
					cos = std::cos(omega*j*k);
					sin = std::sin(omega*j*k);
					Xr[k] += xr[j] * cos + xi[j] * sin;
					Xi[k] += - xr[j] * sin + xi[j] * cos;
				}
				Xr[k] /= root_N;
				Xi[k] /= root_N;
			}
		}else{
			// Inverse DFT
			for(int k=0; k<N; k++){
				Xr[k] = 0;
				Xi[k] = 0;
				for(int j=0; j<N; j++){
					cos = std::cos(omega*j*k);
					sin = std::sin(omega*j*k);
					Xr[k] += xr[j] * cos - xi[j] * sin;
					Xi[k] += xr[j] * sin + xi[j] * cos;
				}
				Xr[k] /= root_N;
				Xi[k] /= root_N;
			}

		}

		for(int k=0; k<N; k++){
			xr[k] = Xr[k];
			xi[k] = Xi[k];
		}
		delete [] Xr;
		delete [] Xi;
		return 0;
	}

	//------------------------------------------------------------------------------
	// Autocorrelation
	//------------------------------------------------------------------------------
	void autocorrelation(double *R, double *x, int N)
	{
		// �̗p���Ă��鐔����
		// R[dt] = 1/(N-dt) * sum_{i=0}^{N-1-dt} x[i] * x[i+dt]
		// ���̎��̏ꍇ�Adt ���傫���Ȃ�Ƒ��ւ̃o���c�L���傫���Ȃ��Ă��܂��B

		// x  : ���̃f�[�^
		// N  : �f�[�^ x �̌�
		// R  : ���߂�ׂ� x �̎��ȑ���
		// dt : delay time (dt = 0, 1, 2, ..., N-1)

		int dt;

		// ���ȑ��։��Z�J�n
		for(dt=0; dt<N; dt++){
			// �܂��͏�����
			R[dt] = 0.0;

			// �e�x������ dt �ɂ����鎩�ȑ��ւ̌v�Z
			for(int i=0; i< N-dt; i++){
				R[dt] += x[i] * x[i+dt];
			}
			R[dt] = R[dt] / (N-dt);
		}
		// ���ȑ��։��Z�I��

		// R[0] = 1 �ƂȂ�悤�ɐ��K���������Ȃ̂ł��B
		// ���[�v�̍Ō�� R[0] = 1 �Ƃ��Ȃ���΂Ȃ�Ȃ����Ƃɒ��ӁB
		for(int dt=N-1; dt>=0; dt--){
			R[dt] = R[dt] / R[0];
		}
	}

	//------------------------------------------------------------------------------
	// Fast Fourier Transform
	//------------------------------------------------------------------------------
	void fft(double *xr, double *xi, int L, int id)
	{
		const int N = std::pow(2, L);
		// xr[k] : �����l�f�[�^�i�����Ɍ��ʂ�����B���̃f�[�^�͎�����j
		// xi[k] : �����l�f�[�^�i�����Ɍ��ʂ�����B���̃f�[�^�͎�����j
		// N     : �f�[�^�̌��i2�̗ݏ�B256, 512, 1024,...�j
		// L     : N = 2^L �� L
		// id    : 1 or -1 �i1�̂Ƃ��t�[���G�ϊ��A-1�̂Ƃ��t�t�[���G�ϊ��j

		int k = 0; // �f�[�^�̓Y���� x[k], x[kk]
		int kk;
		int l;     // �K�� l = 0,1,2,...,L �� L+1 ��

		double omega = 2 * pi / N;      // exp(-2��/N) �� 2��/N
		double p;                         // W^p �� p

		int distance = N;                 // �o�^�t���C���Z����2�f�[�^�̊Ԋu
		int shiftnum = L-1;               // p ���߂邽�߂̉E�V�t�g��

		// �v�Z�ʌ��炷���߂̂ǂ��ł�������Ɨp�ϐ�
		double arg;
		double c, s;
		double alpha, beta;    

		// FFT
		for(l=1; l<=L; l++){

			// �o�^�t���C���Z����2�f�[�^�Ԃ̋�����ݒ�
			// 1�K�ڂł� N/2
			// 2�K�ڂł� N/4
			// 3�K�ڂł� N/8
			distance = distance / 2;

			// ����K l �ɂ����鉉�Z�J�n
			while(k < N){
				for(int i=0; i<distance; i++){

					// p �����߂�
					p = k >> shiftnum;
					p = bitReversal(p, L);

					// x[k] �ɑ΂���o�^�t���C���Z�̑��� x[k+N/2^l] �� k+N/2^l
					kk = k + distance;

					// �o�^�t���C���Z�J�n
					arg = omega * p;
					c = cos(arg);
					s = sin(id * arg);
					alpha = xr[kk] * c - xi[kk] * s;
					beta  = xr[kk] * s + xi[kk] * c;

					xr[kk] = xr[k] - alpha;
					xi[kk] = xi[k] - beta;
					xr[k]  = xr[k] + alpha;
					xi[k]  = xi[k] + beta;
					// �o�^�t���C���Z�I��

					k++;
				}
				k += distance;

			} // ����K l �ɂ����鉉�Z�I��

			// ���̊K�̉��Z�ɓ���O�ɁAk ��������
			k = 0;

			// p �����߂�ۂ� k �̉E�V�t�g����ݒ�
			shiftnum--;
			// ������Ȃ��Ď��̊K�̃o�^�t���C���Z�Ɉȍ~

		} // FFT �I��

		// FFT �����ɂ���ċt�]�����r�b�g�����ɖ߂�
		int brk; // bit reversal of k
		for(k=0; k<N; k++){
			brk = bitReversal(k, L);
			if(brk > k){
				swap(&xr[k], &xr[brk]);
				swap(&xi[k], &xi[brk]);
			}
		}

		// �t�t�[���G�ϊ��ł� N �Ŋ���
		if(id == -1){
			for(k=0; k<N; k++){
				xr[k] = xr[k] / N;
				xi[k] = xi[k] / N;
			}
		}
	}

	//------------------------------------------------------------------------------
	// Band Pass Filter
	//------------------------------------------------------------------------------
	int bpf(double *xr, double *xi, const int width, const int inf_channel, const int sup_channel)
	{
		// ftr             ... real part of Fourier transform (���� ftr �z����t�B���^�����O��̎����f�[�^�Ƃ��ĕԂ��B���� ftr �z��f�[�^�͏�����B)
		// fti             ... imaginary part of Fourier transform
		// WIDTH           ... ��������f�[�^��
		// INF_PASSED_FREQ ... infimum frequency channel of the passed band
		// SUP_PASSED_FREQ ... supremum frequency channel of the passed band

		if( inf_channel < 0 || inf_channel > sup_channel || sup_channel > width){
			fprintf(stderr, "Band Pass Filtering Configuration Error\n");
			return -1;
		}

		// �t�B���^�����O (INF �` SUP �� WID-SUP-1 �` WID-INF-1 �̓X���[)
		for(int i=0; i<inf_channel; i++){
			xr[i] = 0;
			xi[i] = 0;
		}
		for(int i = sup_channel + 1; i < width - sup_channel - 2; i++){
			xr[i] = 0;
			xi[i] = 0;
		}
		for(int i = width - inf_channel; i < width; i++){
			xr[i] = 0;
			xi[i] = 0;
		}

		// �t�[���G�t�ϊ�
		//fft(ftr, fti, std::floor( std::log(WIDTH) / std::log(2.0) ), -1);
		return 0;
	}

	int bpf(float *xr, float *xi, const int width, const int inf_channel, const int sup_channel)
	{
		// ftr             ... real part of Fourier transform (���� ftr �z����t�B���^�����O��̎����f�[�^�Ƃ��ĕԂ��B���� ftr �z��f�[�^�͏�����B)
		// fti             ... imaginary part of Fourier transform
		// WIDTH           ... ��������f�[�^��
		// INF_PASSED_FREQ ... infimum frequency channel of the passed band
		// SUP_PASSED_FREQ ... supremum frequency channel of the passed band

		if( inf_channel < 0 || inf_channel > sup_channel || sup_channel > width){
			fprintf(stderr, "Band Pass Filtering Configuration Error\n");
			return -1;
		}

		// �t�B���^�����O (INF �` SUP �� WID-SUP-1 �` WID-INF-1 �̓X���[)
		for(int i=0; i<inf_channel; i++){
			xr[i] = 0;
			xi[i] = 0;
		}
		for(int i = sup_channel + 1; i < width - sup_channel - 2; i++){
			xr[i] = 0;
			xi[i] = 0;
		}
		for(int i = width - inf_channel; i < width; i++){
			xr[i] = 0;
			xi[i] = 0;
		}

		// �t�[���G�t�ϊ�
		//fft(ftr, fti, std::floor( std::log(WIDTH) / std::log(2.0) ), -1);
		return 0;
	}

	int smoothUsingBPF(double *data, const int width, const int inf_freq_channel, const int sup_freq_channel)
	{
		double *xr = new double[width];
		double *xi = new double[width];
		if(!xr || !xi){
			fprintf(stderr, "Memory Allocation Error\n");
			return -1;
		}
		for(int i=0; i<width; i++){
			xr[i] = data[i];
			xi[i] = 0;
		}

		dft(xr, xi, width, false);
		bpf(xr, xi, width, inf_freq_channel, sup_freq_channel);
		dft(xr, xi, width, true);

		for(int i=0; i<width; i++){
			data[i] = xr[i];
		}

		delete [] xr;
		delete [] xi;
		return 0;
	}

	int smoothUsingBPF(float *data, const int width, const int inf_freq_channel, const int sup_freq_channel)
	{
		float *xr = new float[width];
		float *xi = new float[width];
		if(!xr || !xi){
			fprintf(stderr, "Memory Allocation Error\n");
			return -1;
		}
		for(int i=0; i<width; i++){
			xr[i] = data[i];
			xi[i] = 0;
		}

		dft(xr, xi, width, false);
		bpf(xr, xi, width, inf_freq_channel, sup_freq_channel);
		dft(xr, xi, width, true);

		for(int i=0; i<width; i++){
			data[i] = xr[i];
		}

		delete [] xr;
		delete [] xi;
		return 0;
	}

	/**
	* Normalize data such that y = (x - mean)/stdev
	*/
	void normalize_data(double *data, const unsigned int data_size)
	{
		const double mu    = mean(data, data_size);
		const double sigma = std::sqrt(unbiased_variance(data, data_size, mu));
		for(unsigned int i=0; i<data_size; i++){
			data[i] = (data[i] - mu) / sigma;
		}
	}
	void normalize_data(float *data, const unsigned int data_size)
	{
		const float mu    = mean(data, data_size);
		const float sigma = std::sqrt(unbiased_variance(data, data_size, mu));
		for(unsigned int i=0; i<data_size; i++){
			data[i] = (data[i] - mu) / sigma;
		}
	}

	/*
	*  Summation
	*/
	float summation (const float *data, const int num_data)
	{
		float sum = 0;
		for (int i=0; i<num_data; i++) {
			sum += data[i];
		}
		return sum;
	}

	/**
	* Derive the minimum of data.
	*/
	float  minimum(const float *data, const int num_data)
	{
		float min = data[0];
		for (int i=1; i<num_data; i++) {
			if (min > data[i]){
				min = data[i];
			}
		}
		return min;	
	}


	/**
	* Calculates a mean value.
	* @param data     the target data array
	* @param data_num the size of the data
	*/
	double mean(const double *data, const int data_num)
	{
		double ave = 0;
		for(int i=0; i<data_num; i++){
			ave += data[i];
		}
		ave /= data_num;
		return ave;
	}

	float mean(const float *data, const int data_num)
	{
		float ave = 0;
		for(int i=0; i<data_num; i++){
			ave += data[i];
		}
		ave /= data_num;
		return ave;
	}

	/**
	* Calculates an unbiased variance.
	* @param data     the target data array
	* @param data_num the size of the data
	* @param mean     the mean value of the data
	*/
	double unbiased_variance(const double *data, const int data_num, const double mean)
	{
		double var = 0;
		double err;
		for(int i=0; i<data_num; i++){
			err = data[i] - mean;
			var += (err * err);			
		}
		var /= (data_num - 1); // �s�Ε��U���߂�̂� n �łȂ� n-1 �Ŋ���

		return var;
	}

	float unbiased_variance(const float *data, const int data_num, const float mean)
	{
		float var = 0;
		float err;
		for(int i=0; i<data_num; i++){
			err = data[i] - mean;
			var += (err * err);			
		}
		var /= (data_num - 1); // �s�Ε��U���߂�̂� n �łȂ� n-1 �Ŋ���

		return var;
	}
	//------------------------------------------------------------------------------
	// �W�{���U (n �Ŋ���)
	//------------------------------------------------------------------------------
	double sample_variance(const double *data, const int data_num, const double mean)
	{
		double svar = 0;
		double err;
		for(int i=0; i<data_num; i++){
			err = data[i] - mean;
			svar += (err * err);			
		}
		svar /= data_num; // �W�{���U�Ȃ̂� n �Ŋ���

		return svar;
	}

	//------------------------------------------------------------------------------
	// �W�{�����U
	//------------------------------------------------------------------------------
	double sample_covariance(const double *x, const double *y, const int data_num, const double x_average, const double y_average)
	{
		double cov = 0;
		for(int i=0; i<data_num; i++){
			cov += (x[i] - x_average) * (y[i] - y_average);
		}
		cov /= data_num;

		return cov;
	}

	//------------------------------------------------------------------------------
	// ��A�����̌X����y�ؕ� �i�X���ƐؕЂ̌v�Z��ʂ̊֐��ɕ������Ă���̂ŁA�������߂�ꍇ���Ȃ薳�ʂȌv�Z���s���j
	//------------------------------------------------------------------------------
	double gradient_of_regression_line(const double *x, const double *y, const int data_num)
	{
		double ave_x = mean(x, data_num);
		double ave_y = mean(y, data_num);
		double var_x = sample_variance(x, data_num, ave_x);
		double cov_xy = sample_covariance(x, y, data_num, ave_x, ave_y);

		return cov_xy / var_x ;
	}

	double intercept_of_regression_line(const double *x, const double *y, const int data_num)
	{
		double ave_x = mean(x, data_num);
		double ave_y = mean(y, data_num);
		double var_x = sample_variance(x, data_num, ave_x);
		double cov_xy = sample_covariance(x, y, data_num, ave_x, ave_y);
		double a = cov_xy / var_x ;

		return ave_y - a * ave_x ;
	}


	//------------------------------------------------------------------------------
	// ������\����������J�����_�[���� time_t �ɕϊ�
	//------------------------------------------------------------------------------
	time_t atot(const char *local_YYYYMMDDhhmmss, const int isdst)
	{
		// �e������擾
		const std::string local_time = local_YYYYMMDDhhmmss;
		const std::string yyyy = local_time.substr(0, 4);
		const std::string mm   = local_time.substr(4, 2);
		const std::string dd   = local_time.substr(6, 2);
		const std::string t_hh = local_time.substr(8, 2);
		const std::string t_mm = local_time.substr(10, 2);
		const std::string t_ss = local_time.substr(12, 2);

		// �擾������� tm �\���̂ɓ����
		// mktime �֐��̎d�l��A���[�J�����Ԃ̂܂܂ɂ��Ȃ���΂Ȃ�Ȃ�
		struct tm tm_local_time;
		tm_local_time.tm_year = std::atoi( yyyy.c_str() ) - 1900;
		tm_local_time.tm_mon  = std::atoi( mm.c_str() ) - 1;
		tm_local_time.tm_mday = std::atoi( dd.c_str() );
		tm_local_time.tm_hour = std::atoi( t_hh.c_str() );
		tm_local_time.tm_min  = std::atoi( t_mm.c_str() );
		tm_local_time.tm_sec  = std::atoi( t_ss.c_str() );
		tm_local_time.tm_isdst = isdst;

		// ���[�J������ tm ���J�����_�[���� time_t �ɕϊ�
		time_t time = mktime( &tm_local_time );

		return time;
	}

	//------------------------------------------------------------------------------
	// �������J�����_�[���� time_t �ɕϊ�
	//------------------------------------------------------------------------------
	time_t ttot(const int local_year, const int local_month, const int local_day, const int local_hour, const int local_minute, const int local_second, const int isdst)
	{
		// tm �\���̂ɓ����
		// mktime �֐��̎d�l��A���[�J�����Ԃ̂܂܂ɂ��Ȃ���΂Ȃ�Ȃ�
		struct tm tm_local_time;
		tm_local_time.tm_year = local_year - 1900;
		tm_local_time.tm_mon  = local_month - 1;
		tm_local_time.tm_mday = local_day;
		tm_local_time.tm_hour = local_hour;
		tm_local_time.tm_min  = local_minute;
		tm_local_time.tm_sec  = local_second;
		tm_local_time.tm_isdst = isdst;

		// ���[�J������ tm ���J�����_�[���� time_t �ɕϊ�
		time_t time = mktime( &tm_local_time );

		return time;
	}

	//------------------------------------------------------------------------------
	// �J�����_�[���� time_t �𕶎���ɕϊ�
	//------------------------------------------------------------------------------
	void ttoa(char* str, const time_t timer)
	{
		struct tm tm_local_time = *localtime(&timer);
		int year = tm_local_time.tm_year + 1900;
		int mon  = tm_local_time.tm_mon + 1;
		int mday = tm_local_time.tm_mday;
		int hour = tm_local_time.tm_hour;
		int min  = tm_local_time.tm_min;
		int sec  = tm_local_time.tm_sec;

		sprintf(str, "%4d%02d%02d%02d%02d%02d", year, mon, mday, hour, min, sec);
	}
	//------------------------------------------------------------------------------
	// 3�����̉�]�s��𐶐�����
	//------------------------------------------------------------------------------
/*	void createRotationMatrix_3D(gsl_matrix * matrix, const double angle, const char axis)
	{
		// matrix : �߂�l�ƂȂ�3�����̐����s��
		// angle  : ��]�p�x [rad]
		// axis   : ��]��

		// axis = 'x' : x ���܂��� angle ��]�������]�s��
		// axis = 'y' : y ���܂��ɁE�E�E
		// axis = 'z' : z ���܂��ɁE�E�E

		// matrix ��3�������m�F
		if(matrix->size1 != 3 || matrix->size2 != 3){
			throw "�֐� createRotationMatrix_3D() �̑� 1 �����̍s��T�C�Y���s��";
		}

		// ��]�s��̗v�f���`
		if(axis == 'x'){

			// x ���܂��� angle ��]������s��
			double m[] = { 1.0 ,   0.0             , 0.0             ,
				0.0 ,   std::cos(angle) , std::sin(angle) ,
				0.0 , - std::sin(angle) , std::cos(angle) };

			// gsl_matrix_view �s�������Ă���A����� matrix �ɃR�s�[����B
			// gsl_matrix_set() �֐���p����΁A view_matrix �Ȃ�č�炸���� matrix �ɒl��ݒ�ł��邪�A
			// �R�[�h�̌��₷���̂��ߔz�� m[] �����A�����N�b�V�����u���B
			gsl_matrix_view view_matrix = gsl_matrix_view_array(m, 3, 3);

			// ���������s�� view_matrix ��߂�l matrix �ɃR�s�[
			gsl_matrix_memcpy(matrix, &view_matrix.matrix);


		}else if(axis == 'y'){

			// y ���܂��� angle ��]������s��
			double m[] = { std::cos(angle) , 0.0 , - std::sin(angle) ,
				0.0             , 1.0 ,   0.0             ,
				std::sin(angle) , 0.0 ,   std::cos(angle) };

			gsl_matrix_view view_matrix = gsl_matrix_view_array(m, 3, 3);
			gsl_matrix_memcpy(matrix, &view_matrix.matrix);


		}else if(axis == 'z'){

			// z ���܂��� angle ��]������s��
			double m[] = {   std::cos(angle) , std::sin(angle) , 0.0 ,
				- std::sin(angle) , std::cos(angle) , 0.0 ,
				0.0               , 0.0             , 1.0 };

			gsl_matrix_view view_matrix = gsl_matrix_view_array(m, 3, 3);
			gsl_matrix_memcpy(matrix, &view_matrix.matrix);

		}else{
			throw "�֐� createRotationMatrix_3D() �̑� 3 �����s��";
		}
	} // end of function
*/

	//------------------------------------------------------------------------------
	// �p�x �� �� 0 �� �� �� 2�� �ɔ[�߂�
	//------------------------------------------------------------------------------
	double rangeAngle(const double angle_radian)
	{
		// �� = �� + 2n��, 0 �� �� �� 2��
		double theta = angle_radian;	
		int n = std::floor( theta / twopi ); // ��/2�� = ��/2�� + n �� ��/2�� �͕K������
		double phi = theta - twopi * n;

		return phi;
	}

	//------------------------------------------------------------------------------
	// radian �� hour-minute-second �P�ʂɕϊ�
	//------------------------------------------------------------------------------
	void rad_to_hms(double* hms, const double radian)
	{
		double rad;
		double temp;

		// radian��(h, m, s)�ɕϊ�
		rad = 12.0 * radian / pi; // rad �� h �ɒ����B?.?? h �Ƃ����l�ɂȂ��Ă���B
		hms[0] = floor(rad);      // �������������o���B? h �̂ł�������B
		temp = rad - hms[0];      // �������������o���Btemp = 0.?? h �Ƃ����l�ɂȂ��Ă���B
		temp = temp * 60.0;       // h �� m �ɒ����Btemp = 0.?? h = 0.??*60 m = #.## m �Ƃ����l�ɂȂ��Ă���B
		hms[1] = floor(temp);     // �������������o���B# m �̂ł�������B 
		temp = temp - hms[1];     // �������������o���Btemp = 0.## m �Ƃ����l�ɂȂ��Ă���B
		hms[2] = temp * 60.0;     // m �� s �ɒ����B0.## m = 0.## * 60 s = $.$$ s �Ƃ����l�ɂȂ��Ă���B����łł�������B
	}

	//------------------------------------------------------------------------------
	// radian �� degree-arcminute-arcsecond �P�ʂɕϊ�
	//------------------------------------------------------------------------------
	void rad_to_dms(double* dms, const double radian)
	{
		double rad;
		double temp;

		// radian��(deg, arcmin, arcsec)�ɕϊ�
		rad = 180.0 * radian / pi; // rad �� deg �ɒ����B       ?.?? deg �Ƃ����l�ɂȂ��Ă���B
		dms[0] = floor(rad);       // �������������o���B      ? deg �̂ł�������B
		temp = rad - dms[0];       // �������������o���B      temp = 0.?? deg �Ƃ����l�ɂȂ��Ă���B
		temp = temp * 60.0;        // deg �� arcmin �ɒ����B    0.?? deg = 0.?? * 60 arcmin = #.## arcmin �Ƃ����l�ɂȂ��Ă���B
		dms[1] = floor(temp);      // �������������o���B      # arcmin �̂ł�������B
		temp = temp - dms[1];      // �������������o���B      temp = 0.## arcmin �Ƃ����l�ɂȂ��Ă���B
		dms[2] = temp * 60.0;      // arcmin �� arcsec �ɒ����B 0.## arcmin = 0.## * 60 arcsec = $.$$ arcsec �Ƃ����l�ɂȂ��Ă���B
	}

	//------------------------------------------------------------------------------
	// hour-minute-second �� radian �P�ʂɕϊ�
	//------------------------------------------------------------------------------
	double hms_to_rad(const double hour, const double minute, const double second)
	{
		double h = hour + minute / 60.0 + second / 3600.0; // �Ⴆ�� 5h30m00s �� 5.5h �ɂȂ�
		double rad = h * twopi / 24.0;
		return rad;
	}

	//------------------------------------------------------------------------------
	// degree-arcminute-arcsecond �� radian �P�ʂɕϊ�
	//------------------------------------------------------------------------------
	double dms_to_rad(const double degree, const double arcminute, const double arcsecond)
	{
		double d = degree + arcminute / 60.0 + arcsecond / 3600.0; // �Ⴆ�� 30d30'00'' �� 30.5d �ɂȂ�
		double rad = d * twopi / 360.0;
		return rad;	
	}


	/**
	 * In String_A, searches String_B and replaces it with String_C.
	 * This code is from the Internet.
	 * @param str         String_A
	 * @param deleted_str String_B
	 * @param added_str   String_C
	 * 
	 */
	void strReplace (std::string& str, const std::string& deleted_str, const std::string& added_str) {
		std::string::size_type pos = 0;
		while(pos = str.find(deleted_str, pos), pos != std::string::npos) {
			str.replace(pos, deleted_str.length(), added_str);
			pos += added_str.length();
		}
	}


	/**
	 * Read STANDARD INPUT.
	 * @param str       the input string
	 * @param STRLENMAX the maximum readable length of the input string
	 * @return          the length of the input string
	 *
	 * Usage:
	 * 	printf("Please input: ");
	 *	const size_t STRLENMAX = 3;
	 *	char   str[STRLENMAX + 1];
	 *	size_t strlength = ta_read_stdin(str, STRLENMAX);
	 *	printf("String: %s\n", str);
	 *	printf("Length: %zd\n", strlength);
	 */
	size_t read_stdin(char * str, const size_t STRLENMAX)
	{
		/* getchar() ��int�^�Œl���󂯎��*/
		int c = 0;
		size_t strlength = 0;
		while((c = getchar()) != EOF){
			if(c == '\n'){
				break;
			}
			if(strlength < STRLENMAX){
				str[strlength] = c;
			}else{
				break;
			}
			strlength++;
		}
		str[strlength] = '\0';
		return strlength;
	}	

	bool stdin_yes_or_no()
	{
		printf(" [y/n]: ");
		int c = getchar();
		bool confirm = false;
		if (c == 'y') {
			confirm = true;
		} else {
			confirm = false;
		}
		/*
		if(c == 'y'){
			confirm = 1;
		}else if(c == 'n'){
			confirm = 0;
		}else{
			confirm = -1;
		}
		*/
		// Throw residue on stdin
		while((c = getchar()) != EOF){
			if(c == '\n'){
				break;
			}
		}
		return confirm;
	}

	int saveSpectrumOfComplexData (const char *filename, const float *data, const int data_length)
	{
		std::ofstream fout (filename);
		const int positive_region = std::floor (data_length / 2.0);
		for (int f = 1 + positive_region; f < data_length; f++) {
			fout << f - data_length << "\t" << data[f] << "\n";
		}
		for (int f = 0; f <= positive_region; f++) {
			fout << f << "\t" << data[f] << "\n";
		}
		fout.close();
		std::cout << "Output: " << filename << std::endl;
		return 0;
	}

	int saveData1D (const char *filename, const float *data, const int data_length)
	{
		std::ofstream fout (filename);
		for(int i = 0; i < data_length; i++){
			fout << data[i] << "\n";
		}
		fout.close();
		std::cout << "Output: " << filename << std::endl;
		return 0;
	}

	int saveData2D (const char *filename, const float *data, const int data_length, const float x_resolution)
	{
		std::ofstream fout (filename);
		for (int i = 0; i < data_length; i++) {
			fout << i * x_resolution << "\t" << data[i] << "\n";
		}
		fout.close();
		std::cout << "Output: " << filename << std::endl;
		return 0;
	}
    /**
	int saveDataForGnuplot3D(const char *filename, float **data, const int x_length, const float x_unit, const int y_length, const float y_unit)
	{
		FILE *ofp = fopen(filename, "w");
		if(ofp == NULL){
			fprintf(stderr, "%s\n", messageFileOpenError(filename));
			return -1;
		}
		for(int x=0; x<x_length; x++){
			for(int y=0; y<y_length; y++){
				fprintf(ofp, "%f\t%f\t%f\n", x * x_unit, y * y_unit, data[x][y]);
			}
			fprintf(ofp, "\n");
		}
		fclose(ofp);
		return 0;
	}
     **/


	/* ���b�Z�[�W���� */
	std::string messageFileOpenError(const std::string filename)
	{
		return "File Open Error: " + filename;
	}

	void print_elapse (std::string str, clock_t start_time, const clock_t end_time)
	{
		printf("%s: %.2f s\n", str.c_str(), (double)(end_time - start_time)/CLOCKS_PER_SEC);
		return;
	}



} // end of namespace
