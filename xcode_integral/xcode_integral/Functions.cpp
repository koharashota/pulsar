#include "Functions.h"

namespace ta
{
	//------------------------------------------------------------------------------
	// 整数 N のビットを逆転させた整数 brN を求めるプログラム
	//------------------------------------------------------------------------------
	int bitReversal(int N, int L) // ビットを逆転させたい整数 N と、そのビット数 L
	{
		int n1, n2;
		int brN; // 整数 N のビットを逆転させた整数 bit reversal N

		n1 = N;
		brN = 0;

		// 単に L 回繰り返す
		for(int i=0; i<L; i++){
			n2 = int(n1 / 2);
			brN = 2*brN + n1 - 2*n2;
			n1 = n2;
		}

		return brN;
	}

	//------------------------------------------------------------------------------
	// 数値 x, y を交換するプログラム
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
		// x[i]   : 元の実数データ（虚数データもある場合はプログラムし直す必要あり）
		//          この元データは失われる
		// N      : データ数

		// id = 1 : ハミング窓
		// id = 2 : ハニング窓
		// id = 3 : ブラックマン窓

		double omega = 2 * 4*std::atan(1) / (N-1); // 2π/(N-1) は計算しておく

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
			double sigma = N / 4.0; // N = 4σ とσを定義する
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
		/* 2次・3次多項式適合移動平均 (Savitzky-Golay法) を計算
		 ma[i] = 1/W sigma_{j=-m}^{m} w[j] x[i+j]
		 m      : 平滑化点数 2m+1
		 datNum : 元データ数
		 平滑化可能範囲 i = m+1, ... , datNum-m-1 の計 datNum-2m 個
		 */

		// 重み係数 w を設定
		double *w = new double[2*m+1];
		double W = (4*m*m - 1) * (2*m + 3) / 3;
		for(int jj=0; jj <= 2*m; jj++){ // j = -m ,., m の代わりに jj = j + m = 0 ,., 2m とする
			w[jj] = 3*m*(m+1) - 1 - 5*(jj-m)*(jj-m);
			w[jj] = w[jj] / W;
		}

		// i = 0 ,., m-1 のm個のデータは平滑化できないので生データをコピー
		for(int i=0; i<m; i++){
			ma[i] = x[i];
		}

		// i = datNum-m ,., datNum-1 のm個のデータは平滑化できないので生データをコピー
		for(int i=datNum-m; i<datNum; i++){
			ma[i] = x[i];
		}

		// i = m ,., datNum-m-1 のdatNum-2m個のデータを平滑化
		for(int i=m; i<datNum-m; i++){
			// ここでちゃんと ma[i] は初期化しなくちゃね！
			ma[i] = 0;
			for(int jj=0; jj<= 2*m; jj++){ // j = -m ,., m の代わりに jj = j + m = 0 ,., 2m とする
				ma[i] = ma[i] + ( w[jj] * x[i + jj - m] );
			}
		}

		// w[j] はもう使わんよ
		delete [] w;

	}

	//------------------------------------------------------------------------------
	// Discrete Fourier Transform
	//------------------------------------------------------------------------------
	int dft(double *xr, double *xi, const int N, const bool inverse)
	{
		// xr : 元のデータの実部（ここに結果が入る。元のデータは失われる）
		// xi : 元のデータの虚部（ここに結果が入る。元のデータは失われる）
		// N  : データ数

		// この Xr Xi は必要。アルゴリズム上 xr xi に直接結果を代入することはできない。
		double *Xr = new double[N]; // フーリエ変換の実部
		double *Xi = new double[N]; // フーリエ変換の虚部
		if(!Xr || !Xi){
			fprintf(stderr, "Memory Allocation Error\n");
			return -1;
		}
		
		double omega = 2 * 6*asin(0.5) / N; // 2π/N を omega とする
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

		// x に X を代入し X を消去
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
		float *Xr = new float[N]; // フーリエ変換の実部
		float *Xi = new float[N]; // フーリエ変換の虚部
		if(!Xr || !Xi){
			fprintf(stderr, "Memory Allocation Error\n");
			return -1;
		}
		
		float omega = 2 * 6*asin(0.5) / N; // 2π/N を omega とする
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
		// 採用している数式は
		// R[dt] = 1/(N-dt) * sum_{i=0}^{N-1-dt} x[i] * x[i+dt]
		// この式の場合、dt が大きくなると相関のバラツキが大きくなってしまう。

		// x  : 元のデータ
		// N  : データ x の個数
		// R  : 求めるべき x の自己相関
		// dt : delay time (dt = 0, 1, 2, ..., N-1)

		int dt;

		// 自己相関演算開始
		for(dt=0; dt<N; dt++){
			// まずは初期化
			R[dt] = 0.0;

			// 各遅延時間 dt における自己相関の計算
			for(int i=0; i< N-dt; i++){
				R[dt] += x[i] * x[i+dt];
			}
			R[dt] = R[dt] / (N-dt);
		}
		// 自己相関演算終了

		// R[0] = 1 となるように正規化するもんなのです。
		// ループの最後で R[0] = 1 としなければならないことに注意。
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
		// xr[k] : 実数値データ（ここに結果が入る。元のデータは失われる）
		// xi[k] : 虚数値データ（ここに結果が入る。元のデータは失われる）
		// N     : データの個数（2の累乗。256, 512, 1024,...）
		// L     : N = 2^L の L
		// id    : 1 or -1 （1のときフーリエ変換、-1のとき逆フーリエ変換）

		int k = 0; // データの添え字 x[k], x[kk]
		int kk;
		int l;     // 階数 l = 0,1,2,...,L の L+1 個

		double omega = 2 * pi / N;      // exp(-2π/N) の 2π/N
		double p;                         // W^p の p

		int distance = N;                 // バタフライ演算する2データの間隔
		int shiftnum = L-1;               // p 求めるための右シフト回数

		// 計算量減らすためのどうでもいい作業用変数
		double arg;
		double c, s;
		double alpha, beta;    

		// FFT
		for(l=1; l<=L; l++){

			// バタフライ演算する2データ間の距離を設定
			// 1階目では N/2
			// 2階目では N/4
			// 3階目では N/8
			distance = distance / 2;

			// 或る階 l における演算開始
			while(k < N){
				for(int i=0; i<distance; i++){

					// p を求める
					p = k >> shiftnum;
					p = bitReversal(p, L);

					// x[k] に対するバタフライ演算の相方 x[k+N/2^l] の k+N/2^l
					kk = k + distance;

					// バタフライ演算開始
					arg = omega * p;
					c = cos(arg);
					s = sin(id * arg);
					alpha = xr[kk] * c - xi[kk] * s;
					beta  = xr[kk] * s + xi[kk] * c;

					xr[kk] = xr[k] - alpha;
					xi[kk] = xi[k] - beta;
					xr[k]  = xr[k] + alpha;
					xi[k]  = xi[k] + beta;
					// バタフライ演算終了

					k++;
				}
				k += distance;

			} // 或る階 l における演算終了

			// 次の階の演算に入る前に、k を初期化
			k = 0;

			// p を求める際の k の右シフト数を設定
			shiftnum--;
			// これを以って次の階のバタフライ演算に以降

		} // FFT 終了

		// FFT 処理によって逆転したビットを元に戻す
		int brk; // bit reversal of k
		for(k=0; k<N; k++){
			brk = bitReversal(k, L);
			if(brk > k){
				swap(&xr[k], &xr[brk]);
				swap(&xi[k], &xi[brk]);
			}
		}

		// 逆フーリエ変換では N で割る
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
		// ftr             ... real part of Fourier transform (この ftr 配列をフィルタリング後の実数データとして返す。元の ftr 配列データは消える。)
		// fti             ... imaginary part of Fourier transform
		// WIDTH           ... 処理するデータ数
		// INF_PASSED_FREQ ... infimum frequency channel of the passed band
		// SUP_PASSED_FREQ ... supremum frequency channel of the passed band

		if( inf_channel < 0 || inf_channel > sup_channel || sup_channel > width){
			fprintf(stderr, "Band Pass Filtering Configuration Error\n");
			return -1;
		}

		// フィルタリング (INF ～ SUP と WID-SUP-1 ～ WID-INF-1 はスルー)
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

		// フーリエ逆変換
		//fft(ftr, fti, std::floor( std::log(WIDTH) / std::log(2.0) ), -1);
		return 0;
	}

	int bpf(float *xr, float *xi, const int width, const int inf_channel, const int sup_channel)
	{
		// ftr             ... real part of Fourier transform (この ftr 配列をフィルタリング後の実数データとして返す。元の ftr 配列データは消える。)
		// fti             ... imaginary part of Fourier transform
		// WIDTH           ... 処理するデータ数
		// INF_PASSED_FREQ ... infimum frequency channel of the passed band
		// SUP_PASSED_FREQ ... supremum frequency channel of the passed band

		if( inf_channel < 0 || inf_channel > sup_channel || sup_channel > width){
			fprintf(stderr, "Band Pass Filtering Configuration Error\n");
			return -1;
		}

		// フィルタリング (INF ～ SUP と WID-SUP-1 ～ WID-INF-1 はスルー)
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

		// フーリエ逆変換
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
		var /= (data_num - 1); // 不偏分散求めるので n でなく n-1 で割る

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
		var /= (data_num - 1); // 不偏分散求めるので n でなく n-1 で割る

		return var;
	}
	//------------------------------------------------------------------------------
	// 標本分散 (n で割る)
	//------------------------------------------------------------------------------
	double sample_variance(const double *data, const int data_num, const double mean)
	{
		double svar = 0;
		double err;
		for(int i=0; i<data_num; i++){
			err = data[i] - mean;
			svar += (err * err);			
		}
		svar /= data_num; // 標本分散なので n で割る

		return svar;
	}

	//------------------------------------------------------------------------------
	// 標本共分散
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
	// 回帰直線の傾き＆y切片 （傾きと切片の計算を別の関数に分離しているので、両方求める場合かなり無駄な計算を行う）
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
	// 時刻を表す文字列をカレンダー時間 time_t に変換
	//------------------------------------------------------------------------------
	time_t atot(const char *local_YYYYMMDDhhmmss, const int isdst)
	{
		// 各文字列取得
		const std::string local_time = local_YYYYMMDDhhmmss;
		const std::string yyyy = local_time.substr(0, 4);
		const std::string mm   = local_time.substr(4, 2);
		const std::string dd   = local_time.substr(6, 2);
		const std::string t_hh = local_time.substr(8, 2);
		const std::string t_mm = local_time.substr(10, 2);
		const std::string t_ss = local_time.substr(12, 2);

		// 取得文字列を tm 構造体に入れる
		// mktime 関数の仕様上、ローカル時間のままにしなければならない
		struct tm tm_local_time;
		tm_local_time.tm_year = std::atoi( yyyy.c_str() ) - 1900;
		tm_local_time.tm_mon  = std::atoi( mm.c_str() ) - 1;
		tm_local_time.tm_mday = std::atoi( dd.c_str() );
		tm_local_time.tm_hour = std::atoi( t_hh.c_str() );
		tm_local_time.tm_min  = std::atoi( t_mm.c_str() );
		tm_local_time.tm_sec  = std::atoi( t_ss.c_str() );
		tm_local_time.tm_isdst = isdst;

		// ローカル時間 tm をカレンダー時間 time_t に変換
		time_t time = mktime( &tm_local_time );

		return time;
	}

	//------------------------------------------------------------------------------
	// 時刻をカレンダー時間 time_t に変換
	//------------------------------------------------------------------------------
	time_t ttot(const int local_year, const int local_month, const int local_day, const int local_hour, const int local_minute, const int local_second, const int isdst)
	{
		// tm 構造体に入れる
		// mktime 関数の仕様上、ローカル時間のままにしなければならない
		struct tm tm_local_time;
		tm_local_time.tm_year = local_year - 1900;
		tm_local_time.tm_mon  = local_month - 1;
		tm_local_time.tm_mday = local_day;
		tm_local_time.tm_hour = local_hour;
		tm_local_time.tm_min  = local_minute;
		tm_local_time.tm_sec  = local_second;
		tm_local_time.tm_isdst = isdst;

		// ローカル時間 tm をカレンダー時間 time_t に変換
		time_t time = mktime( &tm_local_time );

		return time;
	}

	//------------------------------------------------------------------------------
	// カレンダー時間 time_t を文字列に変換
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
	// 3次元の回転行列を生成する
	//------------------------------------------------------------------------------
/*	void createRotationMatrix_3D(gsl_matrix * matrix, const double angle, const char axis)
	{
		// matrix : 戻り値となる3次元の正方行列
		// angle  : 回転角度 [rad]
		// axis   : 回転軸

		// axis = 'x' : x 軸まわりに angle 回転させる回転行列
		// axis = 'y' : y 軸まわりに・・・
		// axis = 'z' : z 軸まわりに・・・

		// matrix が3次元か確認
		if(matrix->size1 != 3 || matrix->size2 != 3){
			throw "関数 createRotationMatrix_3D() の第 1 引数の行列サイズが不正";
		}

		// 回転行列の要素を定義
		if(axis == 'x'){

			// x 軸まわりに angle 回転させる行列
			double m[] = { 1.0 ,   0.0             , 0.0             ,
				0.0 ,   std::cos(angle) , std::sin(angle) ,
				0.0 , - std::sin(angle) , std::cos(angle) };

			// gsl_matrix_view 行列を作ってから、それを matrix にコピーする。
			// gsl_matrix_set() 関数を用いれば、 view_matrix なんて作らず直接 matrix に値を設定できるが、
			// コードの見やすさのため配列 m[] を作り、ワンクッション置く。
			gsl_matrix_view view_matrix = gsl_matrix_view_array(m, 3, 3);

			// 生成した行列 view_matrix を戻り値 matrix にコピー
			gsl_matrix_memcpy(matrix, &view_matrix.matrix);


		}else if(axis == 'y'){

			// y 軸まわりに angle 回転させる行列
			double m[] = { std::cos(angle) , 0.0 , - std::sin(angle) ,
				0.0             , 1.0 ,   0.0             ,
				std::sin(angle) , 0.0 ,   std::cos(angle) };

			gsl_matrix_view view_matrix = gsl_matrix_view_array(m, 3, 3);
			gsl_matrix_memcpy(matrix, &view_matrix.matrix);


		}else if(axis == 'z'){

			// z 軸まわりに angle 回転させる行列
			double m[] = {   std::cos(angle) , std::sin(angle) , 0.0 ,
				- std::sin(angle) , std::cos(angle) , 0.0 ,
				0.0               , 0.0             , 1.0 };

			gsl_matrix_view view_matrix = gsl_matrix_view_array(m, 3, 3);
			gsl_matrix_memcpy(matrix, &view_matrix.matrix);

		}else{
			throw "関数 createRotationMatrix_3D() の第 3 引数不正";
		}
	} // end of function
*/

	//------------------------------------------------------------------------------
	// 角度 θ を 0 ≦ θ ＜ 2π に納める
	//------------------------------------------------------------------------------
	double rangeAngle(const double angle_radian)
	{
		// θ = φ + 2nπ, 0 ≦ φ ＜ 2π
		double theta = angle_radian;	
		int n = std::floor( theta / twopi ); // θ/2π = φ/2π + n で φ/2π は必ず小数
		double phi = theta - twopi * n;

		return phi;
	}

	//------------------------------------------------------------------------------
	// radian を hour-minute-second 単位に変換
	//------------------------------------------------------------------------------
	void rad_to_hms(double* hms, const double radian)
	{
		double rad;
		double temp;

		// radianを(h, m, s)に変換
		rad = 12.0 * radian / pi; // rad を h に直す。?.?? h という値になっている。
		hms[0] = floor(rad);      // 整数部分を取り出す。? h のできあがり。
		temp = rad - hms[0];      // 小数部分を取り出す。temp = 0.?? h という値になっている。
		temp = temp * 60.0;       // h を m に直す。temp = 0.?? h = 0.??*60 m = #.## m という値になっている。
		hms[1] = floor(temp);     // 整数部分を取り出す。# m のできあがり。 
		temp = temp - hms[1];     // 小数部分を取り出す。temp = 0.## m という値になっている。
		hms[2] = temp * 60.0;     // m を s に直す。0.## m = 0.## * 60 s = $.$$ s という値になっている。これでできあがり。
	}

	//------------------------------------------------------------------------------
	// radian を degree-arcminute-arcsecond 単位に変換
	//------------------------------------------------------------------------------
	void rad_to_dms(double* dms, const double radian)
	{
		double rad;
		double temp;

		// radianを(deg, arcmin, arcsec)に変換
		rad = 180.0 * radian / pi; // rad を deg に直す。       ?.?? deg という値になっている。
		dms[0] = floor(rad);       // 整数部分を取り出す。      ? deg のできあがり。
		temp = rad - dms[0];       // 小数部分を取り出す。      temp = 0.?? deg という値になっている。
		temp = temp * 60.0;        // deg を arcmin に直す。    0.?? deg = 0.?? * 60 arcmin = #.## arcmin という値になっている。
		dms[1] = floor(temp);      // 整数部分を取り出す。      # arcmin のできあがり。
		temp = temp - dms[1];      // 小数部分を取り出す。      temp = 0.## arcmin という値になっている。
		dms[2] = temp * 60.0;      // arcmin を arcsec に直す。 0.## arcmin = 0.## * 60 arcsec = $.$$ arcsec という値になっている。
	}

	//------------------------------------------------------------------------------
	// hour-minute-second を radian 単位に変換
	//------------------------------------------------------------------------------
	double hms_to_rad(const double hour, const double minute, const double second)
	{
		double h = hour + minute / 60.0 + second / 3600.0; // 例えば 5h30m00s は 5.5h になる
		double rad = h * twopi / 24.0;
		return rad;
	}

	//------------------------------------------------------------------------------
	// degree-arcminute-arcsecond を radian 単位に変換
	//------------------------------------------------------------------------------
	double dms_to_rad(const double degree, const double arcminute, const double arcsecond)
	{
		double d = degree + arcminute / 60.0 + arcsecond / 3600.0; // 例えば 30d30'00'' は 30.5d になる
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
		/* getchar() はint型で値を受け取る*/
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


	/* メッセージ生成 */
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
