/*

	This is the serial version of the code to generat the Virus, Infected and Viral population grid over a number of timesteps.
	
*/
#define master 0
#include <iostream>
#include <mpi.h>
#include <fstream>  //input/output from files
#include <cmath>    //pow
#include <stdlib.h> //atof and exit
#include <algorithm> //max
#include "binary_read_write.h" //Wrappers to read and write to binary file

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
	//START OF PARALLEL
	MPI::Init(argc, argv);
	

	if(argc < 4){
		printf("Not enough command line arguments, check the readme to see how to configure.");
		exit(1); // we cannot continue without the config or divisions
	}




	MPI_Datatype mpi_tiv;


	MPI_Type_contiguous(3, MPI::DOUBLE, &mpi_tiv);  //we need to create a mpi version of this data type so we can send and recieve it
	MPI_Type_commit(&mpi_tiv);

	int rank = MPI::COMM_WORLD.Get_rank();
	int nprocs = MPI::COMM_WORLD.Get_size();

	int nproc_x = atoi(argv[2]);
	int nproc_y = atoi(argv[3]);

	int nprocs_used = nproc_x * nproc_y ;
	if(nprocs_used > nprocs){
		printf("you need to give the program more processors, the input says it requires %d processors, and you only gave it %d",nprocs_used, nprocs );
		exit(2); //We dont have enough processors
	} else if ( nprocs_used != nprocs){
		printf("Warning: not all processors are being utilized");
	}
	
	//INPUT VARIABLES
	double delta_space; //The spatial step length (evenly sized in x and y)
	int grid_width; //The number of x steps
	int grid_height; //The number of y steps
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
	int number_of_timesteps; //The number of time steps
	
	tiv** TIV; //A matrix containing the population/volume unit of T-cells, Infected T-cells, and Viroids at each point

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
	
	
	double broadcast_array[11];
	if(rank == master){
		



		ifstream config_file(argv[1]);

		config_file >> delta_space;
		config_file >> grid_width;
		config_file >> grid_height;
		config_file >> birth_rate_filename;
		config_file >> t_filename;
		config_file >> i_filename;
		config_file >> v_filename;
		config_file >> diffusion_tcells;
		config_file >> diffusion_infected;
		config_file >> diffusion_virus;
		config_file >> death_tcells;
		config_file >> death_infected;
		config_file >> death_virus;
		config_file >> burst_rate;
		config_file >> transmission_vt;
		config_file >> transmission_it;
		config_file >> delta_t;
		config_file >> number_of_timesteps;
		config_file >> result_t_filename;
		config_file >> result_i_filename;
		config_file >> result_v_filename;
	
		config_file.close();

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
		

		broadcast_array[0] = a1;
		broadcast_array[1] = a2;
		broadcast_array[2] = a3;
		broadcast_array[3] = a4;
		broadcast_array[4] = b1;
		broadcast_array[5] = b2;
		broadcast_array[6] = b3;
		broadcast_array[7] = b4;			
		broadcast_array[8] = c1;
		broadcast_array[9] = c2;	
		broadcast_array[10] = c3;
	}
	

	MPI::COMM_WORLD.Bcast(( broadcast_array, 11, MPI::DOUBLE, master );
	//TODO: broadcast 

	a1 = broadcast_array[0];
	a2 = broadcast_array[1];
	a3 = broadcast_array[2];
	a4 = broadcast_array[3];
	b1 = broadcast_array[4];
	b2 = broadcast_array[5];
	b3 = broadcast_array[6];
	b4 = broadcast_array[7];
	c1 = broadcast_array[8];
	c2 = broadcast_array[9];
	c3 = broadcast_array[10];

	
	if(rank == master){
		//Initializing TIV
		TIV = new tiv*[grid_height];
		for(int i = 0; i < grid_height; ++i ){
			TIV[i] = new tiv[grid_width];
		}
		//READING the T, I, and V files into the TIV matrix (opening files, reading for each i/j, closing files)
		ifstream t_file(t_filename.c_str());
		ifstream i_file(i_filename.c_str());
		ifstream v_file(v_filename.c_str());
		for(int i = 0; i < grid_height; ++i){
		        for(int j = 0; j < grid_width; ++j){
		                binary_read(t_file , TIV[i][j].t);
		                binary_read(i_file , TIV[i][j].i);
		                binary_read(v_file , TIV[i][j].v);
		        }
		}
		t_file.close();
		i_file.close();
		v_file.close();

		//Initializing birth rate matrix
	
		tcell_birth_rate = new double*[grid_width];
		for(int i = 0; i < grid_height; ++i ){
			tcell_birth_rate[i] = new double[grid_width];
		} 
	
		//READING the birth rate file (lambda) into a matrix
		birth_rate_file.open(birth_rate_filename.c_str(), std::ios::in | std::ios::binary);
		for(int i = 0; i < grid_height; ++i){
			for(int j = 0; j < grid_width; ++j){
				birth_rate_file.read(reinterpret_cast<char *>(&tcell_birth_rate[i][j]), sizeof(tcell_birth_rate[i][j]));
				tcell_birth_rate[i][j] *= delta_t;
			}
		}
		birth_rate_file.close();
	}
	
	
	

	//TODO: make parallel
	//Initializing TIV_next
	
	// The brunt of the code (TIV_next from TIV)
	for(int n = 0; n < number_of_timesteps; ++n){ //for each time step from 0 to n-1

		//TODO: needs to do what it says in green notes
	}
	
	//TODO: recombine matrix
	
	if(rank == master){
		//WRITING the TIV matrix elements into "result" files t,i, and v.
		ofstream result_t_file(result_t_filename.c_str(), ios::out |ios::binary);
		ofstream result_i_file(result_i_filename.c_str(), ios::out |ios::binary);
		ofstream result_v_file(result_v_filename.c_str(), ios::out |ios::binary);
		for(int i = 0; i < grid_height; ++i){
		        for(int j = 0; j < grid_width; ++j){
		        	binary_write(result_t_file , TIV[i][j].t);
		                binary_write(result_i_file , TIV[i][j].i);
		                binary_write(result_v_file , TIV[i][j].v);
		        }
		}
		result_t_file.close();
		result_i_file.close();
		result_v_file.close();
	}

	

	MPI::Finalize();
	
	return 0;
}