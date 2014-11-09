#include <fstream>


using namespace std;
struct tiv{
	double t;
	double i;
	double v;

}
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
	
	tiv** TIV;
	tiv** TIV_next;	

	double** tcell_birth_rate;
	


	string birth_rate_filename;
	ifstream birth_rate_file;
	
	//todo: setup the config 
	birth_rate_file.open(birth_rate_filename.c_str(), std::ios::in | std::ios::binary);
	for(int i = 1; i < grid_height-1; ++i){
		for(int j = 1; j < grid_width-1; ++j){
			birth_rate_file.read(reinterpret_cast<char *>(&tcell_birth_rate[i][j]), sizeof(tcell_birth_rate[i][j]));
			tcell_birth_rate[i][j] *= delta_t;
		}
	}
	birth_rate_file.close();
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
	
	TIV_next = new tiv*[grid_height];
	for(int i = 0; i < grid_height; ++i ){
		TIV_next[i] = new tiv[grid_width];
	}
	for(int n = 0; n < number_of_time_steps; ++n){

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


	return 0;
}
