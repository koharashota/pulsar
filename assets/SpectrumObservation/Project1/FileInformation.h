/*
	Usage:

	Aoki::FileInformation fileinfo(argv[1], '_');
	const path dest_dir = fileinfo.getOutputDestinationDirectory() / "hogehoge"; // = /opt/data/hogehoge
	printf("Input: %s\n", fileinfo.getFilepath().string().c_str());
*/

#ifndef TA_FILEINFORMATION
#define TA_FILEINFORMATION

#include <ctime>
#include <string>
#include <boost/filesystem.hpp>

namespace fs = boost::filesystem;

namespace ta
{
	class FileInformation
	{
		private:
			fs::path filepath;                     // = /otp/data_0001.bin
			fs::path filestem;                     // = data_0001
			fs::path filestem_substr;              // = data
			fs::path output_destination_directory; // = /opt/data
			size_t   filesize_B;
			
		public:
			FileInformation();
			FileInformation(const char *file, const char delimiter);
			~FileInformation();
			
			fs::path getFilePath()       {return filepath;}
			fs::path getOutputDestinationDirectory(){return output_destination_directory;}
			fs::path getFileStem()       {return filestem;}
			fs::path getFileStemSubstr() {return filestem_substr;}
			size_t   getFilesize_B()     {return filesize_B;}
			

		
	};
}
#endif