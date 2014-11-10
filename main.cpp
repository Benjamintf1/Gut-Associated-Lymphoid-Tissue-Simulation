#include <fstream>  //input/output from files
#include <cmath>    //pow
#include <stdlib.h> //atof

using namespace std;

//struct tiv:
//Stores three double values representing the
//population of T-cells, Infected cells, and Viroids
//per volume at a point.
struct tiv{
	double t;
	double i;
	double v;

};
int main (int argc, char** argv){
	//INPUT VARIABLES
	double diffusion_tcells; //D_T: The rate at which T-cells diffuse spatially
	double diffusion_infected; //D_I: The rate at which infected cells diffuse spatially
	double diffusion_virus; //D_V: The rate at which viroids diffuse spatially (generally, > D_I or D_T)
	double death_tcells; //d_T: The death rate of T-cells (from natural deaths)
	double death_infected; //d_I: The death rate of infected cells
	double death_virus; //d_V: The death rate of viroids
	double burst_rate; //N: The amount of viroids created when an infected cell bursts
	double transmission_vt; //k1: The rate of infection when T-cells and Viroids are near
	double transmission_it; //k2: The rate of (cell-to-cell) infection when T-cells and Infected T-cells are near
	double delta_t; //The time step length
	double delta_space; //The spatial step length (evenly sized in x and y)
	int grid_width; //The number of x steps
	int grid_height; //The number of y steps
	int number_of_timesteps; //The number of time steps
	
	tiv** TIV; //A matrix containing the population/volume unit of T-cells, Infected T-cells, and Viroids at each point
	tiv** TIV_next;	//Used to temporarily store the next time step's data.

	double** tcell_birth_rate; //Lambda: A given matrix of where T-cells are produced/generated into the system.
	
	//Variables for reading/writing materices from files:
	string birth_rate_filename;
	string t_filename;
	string i_filename;
	string v_filename;
	
	string result_t_filename;
	string result_i_filename;
	string result_v_filename;
	
	ifstream birth_rate_file;
	
	
	

	if(argc < 2){
		printf("You must give the program a configuration file. Check the readme for how the config file should be formatted");
		exit(1); // we cannot continue without the config
	}



	ifstream config_file(argv[1]);
	
	config_file >> diffusion_tcells;
	config_file >> diffusion_infected;
	config_file >> diffusion_virus;
	config_file >> death_tcells;
	config_file >> death_infected;
	config_file >> death_virus;
	config_file >> transmission_vt;
	config_file >> transmission_it;
	config_file >> delta_t;
	config_file >> delta_space;
	config_file >> grid_width;
	config_file >> grid_height;
	config_file >> number_of_timesteps;
	config_file >> birth_rate_filename;
	config_file >> t_filename;
	config_file >> i_filename;
	config_file >> v_filename;
	config_file >> result_t_filename;
	config_file >> result_i_filename;
	config_file >> result_v_filename;
	
	config_file.close();



	
	ifstream t_file(t_filename);
	ifstream i_file(i_filename);
	ifstream v_file(v_filename);
	for(int i = 0; i < grid_height; ++i){
                for(int j = 0; j < grid_width; ++j){
                        t_file.read(reinterpret_cast<char *>(&TIV[i][j].t), sizeof(TIV[i][j].t));
			i_file.read(reinterpret_cast<char *>(&TIV[i][j].i), sizeof(TIV[i][j].i));
			v_file.read(reinterpret_cast<char *>(&TIV[i][j].v), sizeof(TIV[i][j].v));
                }
        }
	t_file.close();
	i_file.close();
	v_file.close();
	
	birth_rate_file.open(birth_rate_filename.c_str(), std::ios::in | std::ios::binary);
	for(int i = 1; i < grid_height-1; ++i){
		for(int j = 1; j < grid_width-1; ++j){
			birth_rate_file.read(reinterpret_cast<char *>(&tcell_birth_rate[i][j]), sizeof(tcell_birth_rate[i][j]));
			tcell_birth_rate[i][j] *= delta_t;
		}
	}
	birth_rate_file.close();
	
	//VARIABLES: constants in the equations (to avoid repeated multiplications/additions of our input variables together)
	double a4 = delta_t * diffusion_tcells / pow(delta_space, 2);
	double a3 = delta_t * transmission_it;
	double a2 = delta_t * transmission_vt;
	double a1 = 1 - 4 * a4 - delta_t * death_tcells;

	double b4 = delta_t * diffusion_infected / pow(delta_space, 2);
	double b3 = a2;
	double b2 = a3;
	double b1 = 1 - 4 * b4 - delta_t * death_infected;
	
	double c3 = delta_t * diffusion_virus / pow(delta_space, 2);
	double c2 = delta_t * burst_rate * death_infected; 
	double c1 = 1 - 4 * c3 - delta_t * death_virus; 
	
	//Initializing TIV_next
	TIV_next = new tiv*[grid_height];
	for(int i = 0; i < grid_height; ++i ){
		TIV_next[i] = new tiv[grid_width];
	}
	
	// The brunt of the code (TIV_next from TIV)
	for(int n = 0; n < number_of_timesteps; ++n){

		for(int i = 1; i < grid_height-1; ++i){
                	for(int j = 1; j < grid_width-1; ++j){
				TIV_next[i][j].t = TIV[i][j].t * (a1 - a2 * TIV[i][j].v - a3 * TIV[i][j].i) + tcell_birth_rate[i][j] 
						 + a4 * (TIV[i+1][j].t + TIV[i-1][j].t + TIV[i][j+1].t + TIV[i][j-1].t );
				TIV_next[i][j].i = TIV[i][j].i * (b1 + b2*TIV[i][j].t) + b3 * TIV[i][j].t * TIV[i][j].v 
						 + b4 * (TIV[i+1][j].i + TIV[i-1][j].i + TIV[i][j+1].i + TIV[i][j-1].i ); 			
				TIV_next[i][j].v = TIV[i][j].v * c1 + c2 * TIV[i][j].i + c3 * (TIV[i+1][j].v + TIV[i-1][j].v + TIV[i][j+1].v + TIV[i][j-1].v );
	                       
                	}
        	}
	}
	ofstream result_t_file(result_t_filename.c_str(), ios::out |ios::binary);
	ofstream result_i_file(result_i_filename.c_str(), ios::out |ios::binary);
	ofstream result_v_file(result_v_filename.c_str(), ios::out |ios::binary);
	
	for(int i = 0; i < grid_height; ++i){
                for(int j = 0; j < grid_width; ++j){
                        result_t_file.write(reinterpret_cast<char *>(&TIV[i][j].t), sizeof(TIV[i][j].t));
			result_i_file.write(reinterpret_cast<char *>(&TIV[i][j].i), sizeof(TIV[i][j].i));
			result_v_file.write(reinterpret_cast<char *>(&TIV[i][j].v), sizeof(TIV[i][j].v));
                }
        }

	result_t_file.close();
	result_i_file.close();
	result_v_file.close();


	
        birth_rate_file.close();
	
	return 0;
}
