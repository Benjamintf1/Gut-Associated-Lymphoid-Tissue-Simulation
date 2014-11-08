#include <fstream>


using namespace std;

int main (int argc, char** argv){
	double diffusion_tcells;
	double diffusion_infected;
	double diffusion_virus;
	double death_tcells;
	double death_infected;
	double death_virus;
	double burst_rate;
	double transmission_vt;
	double transmission_it;
	double delta_t;
	double delta_space;
	double grid_width;
	double grid_height;
	double number_of_time_steps;
	double** T;
	double** T_next;
	double** V;
	double** V_next;
	double** I;
	double** I_next;
	double** tcell_birth_rate;
	


	string birth_rate_filename;
	ifstream birth_rate_file;
	
	//todo: setup the config 
	birth_rate_file.open(birth_rate_filename.c_str(), std::ios::in | std::ios::binary);
	for(int i = 1; i < grid_height-1; ++i){
		for(int j = 1; j < grid_width-1; ++j){
			birth_rate_file.read(reinterpret_cast<char *>(&tcell_birth_rate[i][j]), sizeof(tcell_birth_rate[i][j]));
			tcell_birth_rate *= delta_t;
		}
	}
	birth_rate_file.close();
	c4 = delta_t * diffusion_tcells / pow(delta_space, 2);
	c3 = delta_t * transmission_it;
	c2 = delta_t * transmission_vt;
	c1 = 1 - 4 * c4 - delta_t * death_tcells;
	

	
	for(int i = 1; i < grid_height-1; ++i){
                for(int j = 1; j < grid_width-1; ++j){
                        
                }
        }

	return 0;
}
