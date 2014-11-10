#include "FileInformation.h"

#define BUF_LEN_MAX 2048

namespace ta
{
	// Constructor
	FileInformation::FileInformation(){
	}
	
	FileInformation::FileInformation(const char *file, const char delimiter){
		filepath = fs::path(file);
		
		if(!fs::exists(filepath)){
			fprintf(stderr, "File Not Exist: %s\n", filepath.string().c_str());
			throw "File Open Error";
		}
		
		filestem = filepath.stem().string(); // = data_0001
		const int filestem_division_pos = filestem.string().rfind(delimiter); // = 4
        
		if (filestem_division_pos <= 0){
			fprintf(stderr, "Invalid File Name: An input filename must be <string 1>%c<string 2>.<extension>. The <string 1> will become the name of the ouput directory.\n", delimiter);
			throw "Invalid Input";
		}
		
		filestem_substr = filestem.string().substr(0, filestem_division_pos); // = data
		output_destination_directory = filepath.parent_path() / filestem_substr; // = /opt/data

		filesize_B = fs::file_size(filepath);
	}
	
	// Destructor
	FileInformation::~FileInformation(){
	}



}