#include <fstream>  //input/output from files
#include <stdlib.h> //exit
#include "binary_read_write.h"  //Wrappers to read and write to binary file
#include "generators.h"

using namespace std;

int main(int argc, char** argv){

	
	if(argc < 2){
		printf("You must give the program a configuration file. Check the readme for how the config file should be formatted");
		exit(1); // we cannot continue without the config
	}



	ifstream config_file(argv[1]);

	
	double delta_space; //The spatial step length (evenly sized in x and y)
	int grid_width; //The number of x steps
	int grid_height; //The number of y steps
	string birth_rate_filename;
	string t_filename;
	string i_filename;
	string v_filename;
	string tissue_types_filename;

	config_file >> delta_space;
	config_file >> grid_width;
	config_file >> grid_height;
	config_file >> birth_rate_filename;
	config_file >> t_filename;
	config_file >> i_filename;
	config_file >> v_filename;
	config_file >> tissue_types_filename;
	
	config_file.close();


	ofstream generated_t_file(t_filename.c_str(), ios::out |ios::binary);
	ofstream generated_i_file(i_filename.c_str(), ios::out |ios::binary);
	ofstream generated_v_file(v_filename.c_str(), ios::out |ios::binary);
	ofstream generated_birth_file(birth_rate_filename.c_str(), ios::out |ios::binary);
	ofstream generated_tissue_types_file(tissue_types_filename.c_str(), ios::out |ios::binary);


	for(int i = 0; i < grid_height; ++i){
		for(int j = 0; j < grid_width; ++j){
			binary_write(generated_t_file , generate_T(j, i, grid_width, grid_height, delta_space));
			binary_write(generated_i_file , generate_I(j, i, grid_width, grid_height, delta_space));
			binary_write(generated_v_file , generate_V(j, i, grid_width, grid_height, delta_space));
			binary_write(generated_birth_file, generate_birth(j, i, grid_width, grid_height, delta_space));
			binary_write(generated_tissue_types_file, generate_tissues(j, i, grid_width, grid_height, delta_space));
		}
	}


	generated_t_file.close();
	generated_i_file.close();
	generated_v_file.close();
	generated_birth_file.close();
	generated_tissue_types_file.close();


	return 0;
}

