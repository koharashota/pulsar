#include "Stats.h"

namespace ta
{
	Stats::Stats () 
	{
		num_data = 0;
		mean_x  = 0;
		mean_y  = 0; 
		var_x   = 0; 
		var_y   = 0; 
		stdev_x = 0; 
		stdev_y = 0;
	}

	Stats::~Stats () {}

	void Stats::calcParams (const cufftComplex *data, const int num_data)
	{
		float *re = new float [num_data];
		float *im = new float [num_data];
		for (int i=0; i<num_data; i++) {
			re[i] = data[i].x;
			im[i] = data[i].y;
		}

		mean_x  = ta::mean (re, num_data);
		mean_y  = ta::mean (im, num_data);
		var_x   = ta::unbiased_variance (re, num_data, mean_x);
		var_y   = ta::unbiased_variance (im, num_data, mean_y);
		stdev_x = std::sqrt (var_x);
		stdev_y = std::sqrt (var_y);
	
		delete [] re;
		delete [] im;
	}


}