#include <cuda.h>
#include <cufft.h>
#include <cmath>
#include "Functions.h"

#ifndef TA_STATS
#define TA_STATS

namespace ta
{
	class Stats
	{
		private:
			int  num_data;   // Number of data samples, which output by getNumData().
			float mean_x, mean_y;
			float var_x, var_y;
			float stdev_x, stdev_y;

		public:
			Stats ();
			~Stats ();
			void calcParams (const cufftComplex *data, const int num_data);
			float getMean_x () const {return mean_x;};
			float getMean_y () const {return mean_y;};
			float getStdev_x () const {return stdev_x;};
			float getStdev_y () const {return stdev_y;};
			int getNumData () const {return num_data;};
	};
}
#endif
