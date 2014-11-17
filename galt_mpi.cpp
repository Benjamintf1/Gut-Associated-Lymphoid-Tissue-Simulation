/*

	This is the pure MPI version of the code to generate the T-Cell, Infected and Viroid population grid over a number of timesteps.
	
*/
//TODO: make parallel
#define master 0

#include <mpi.h>  //This is an mpi program...
#include <fstream>  //input/output from files
#include <cmath>    //pow
#include <stdlib.h> //atof and exit
#include <algorithm> //max
#include "binary_read_write.h" //Wrappers to read and write to binary file
#include <errno.h> //stderr
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
	
	int rank = MPI::COMM_WORLD.Get_rank();
	int nprocs = MPI::COMM_WORLD.Get_size();

	if(rank == master){
		if(argc < 4){
			fprintf(stderr, "Not enough command line arguments, check the readme to see how to configure.\n");
			MPI_Abort(MPI_COMM_WORLD, 1); // we cannot continue without the config or divisions
		}
	}



	MPI_Datatype mpi_tiv;


	MPI_Type_contiguous(3, MPI::DOUBLE, &mpi_tiv);  //we need to create a mpi version of this data type so we can send and recieve it
	MPI_Type_commit(&mpi_tiv);

	

	int nprocs_x = atoi(argv[2]);
	int nprocs_y = atoi(argv[3]);

	int nprocs_used = nprocs_x * nprocs_y ;

	if(rank == master){
		if(nprocs_used > nprocs){
			fprintf(stderr, "you need to give the program more processors, the input says it requires %d processors, and you only gave it %d \n",nprocs_used, nprocs );
			MPI_Abort(MPI_COMM_WORLD, 2); //We dont have enough processors
		} else if ( nprocs_used != nprocs){
			fprintf(stderr, "Warning: not all processors are being utilized \n");
		}

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
	int local_grid_width;
	int local_grid_height;
	
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


	//derived variables
	double a1, a2, a3, a4, b1, b2, b3, b4, c1, c2, c3;
	
	
	double broadcast_array[13];
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



		if( ( grid_width -2 ) % nprocs_x != 0){
			fprintf(stderr, "The width %d(without boundries) is not evenly divisible by %d \n", grid_width-2, nprocs_x);
			MPI_Abort(MPI_COMM_WORLD, 3);
		}


		if( ( grid_height -2 ) % nprocs_y != 0){
			fprintf(stderr, "The height %d(without boundries) is not evenly divisible by %d \n", grid_height-2, nprocs_y);
			MPI_Abort(MPI_COMM_WORLD, 3);
		}


		//VARIABLES: constants in the equations (to avoid repeated multiplications/additions of our input variables together)
		a4 = delta_t * diffusion_tcells / pow(delta_space, 2);
		a3 = delta_t * transmission_it;
		a2 = delta_t * transmission_vt;
		a1 = 1 - 4 * a4 - delta_t * death_tcells;

		b4 = delta_t * diffusion_infected / pow(delta_space, 2);
		b3 = a2;
		b2 = a3;
		b1 = 1 - 4 * b4 - delta_t * death_infected;
	
		c3 = delta_t * diffusion_virus / pow(delta_space, 2);
		c2 = delta_t * burst_rate * death_infected; 
		c1 = 1 - 4 * c3 - delta_t * death_virus; 
		

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

		

		local_grid_width = ((grid_width - 2)/nprocs_x)+2;
		local_grid_height = ((grid_height - 2)/nprocs_y)+2;


		broadcast_array[11] = local_grid_width;
		broadcast_array[12] = local_grid_height;
	}
	

	MPI::COMM_WORLD.Bcast( broadcast_array, 11, MPI::DOUBLE, master );

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

	local_grid_width = broadcast_array[11];
	local_grid_height = broadcast_array[12];
	
	int local_grid_size = local_grid_height * local_grid_width;
	


	
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

		//Setting up Buffers to split the matrices
		tiv** TIV_buffers = new tiv* [nprocs_used];
		for(int i = 0; i < nprocs_used; ++i ){
			TIV_buffers[i] = new tiv[local_grid_size];
		}
		double** birth_rate_buffers = new double* [nprocs_used];
		for(int i = 0; i < nprocs_used; ++i ){
			birth_rate_buffers[i] = new double[local_grid_size];
		}
		//Splitting the Matrices
		int proc_x, proc_y; //The processor's x and y coordinates
		for(int k = 0; k < nprocs_used; ++k) { //for every processor (that we want to use)
			proc_x = k % nprocs_x;
			proc_y = k / nprocs_x;
			
			int n = 0; //location in buffer
			
			for(int i = (proc_y * local_grid_height) - proc_y; i <= (proc_y + 1) * local_grid_height - proc_y; ++i) {
				for(int j = (proc_x * local_grid_width) - proc_x; j <= (proc_x + 1) * local_grid_width - proc_x; ++j) {
					TIV_buffers[k][n] = TIV[i][j];
					birth_rate_buffers[k][n] = tcell_birth_rate[i][j];
					++n;
				}
			}
			
			MPI::COMM_WORLD.Isend(TIV_buffers[k], local_grid_size, mpi_tiv, k, 1); //TAG: 1
			MPI::COMM_WORLD.Isend(birth_rate_buffers[k], local_grid_size, MPI::DOUBLE, k, 2); //TAG: 2
		}
		
		//Deallocating the no longer needed TIV and t_cell_birth_rate matrices:
		for(int i = 0; i < grid_height; ++i ){
			delete[] TIV[i];
			delete[] tcell_birth_rate[i];
		}
		delete[] TIV;
		delete[] tcell_birth_rate;
		
	}
	
	tiv* local_TIV = new tiv[local_grid_size];
	double* local_birth = new double[local_grid_size];

	MPI_Request master_recieve[2];
	master_recieve[0] = MPI::COMM_WORLD.Irecv(local_TIV, local_grid_size, mpi_tiv, master, 1);
	master_recieve[1] = MPI::COMM_WORLD.Irecv(local_birth, local_grid_size, MPI::DOUBLE, master, 2);
	
	
	MPI_Waitall(2, master_recieve, MPI_STATUSES_IGNORE);
	//Initializing TIV_next
	tiv* local_TIV_next = new tiv[local_grid_size];
	
	// The brunt of the code (TIV_next from TIV)

	MPI_Datatype colomn_tiv;
	MPI_Type_vector(local_grid_height-2, 1, local_grid_width, mpi_tiv, &colomn_tiv);
	MPI_Type_commit(&colomn_tiv);

	MPI_Datatype row_tiv;
	MPI_Type_vector(local_grid_width-2, 1, 1, mpi_tiv, &row_tiv);
	MPI_Type_commit(&row_tiv);
	
	bool up, down, left, right;
	int up_proc, down_proc, left_proc, right_proc;
	up = down = left = right = false;

	int proc_x = rank % nprocs_x;
	int proc_y = rank / nprocs_x;

	if(proc_x > 0){
		left = true;
		left_proc = rank - 1;
	} 
	if(proc_x < nprocs_x){
		right = true;
		right_proc = rank + 1;
	} 

	if(proc_y > 0){
		down = true;
		down_proc = rank - nprocs_x;
	} 
	if(proc_y < nprocs_y){
		up = true;
		up_proc = rank + nprocs_x;
	} 
	//Stores our requests so we can wait appropriately
	MPI_request* sends[4];
	MPI_request* receives[4];
	int neighbors = 0;
	for(int n = 0; n < number_of_timesteps; ++n){ //for each time step from 0 to n-1
		
		for(int i = 1; i < local_grid_height-1; ++i){
			//The Brunt of the Math (calculating TIV_next from TIV)
                	for(int j = 1; j < local_grid_width-1; ++j){
                		local_TIV_next[j + i*local_grid_width].t = max(local_TIV[j+i*local_grid_width].t * (a1 - a2 * local_TIV[j+i*local_grid_width].v 
									- a3 * local_TIV[j+i*local_grid_width].i) + local_birth[j+i*local_grid_width]
									+ a4 * (local_TIV[j+(i+1)*local_grid_width].t + local_TIV[j+(i-1)*local_grid_width].t
									+ local_TIV[(j+1)+ i*local_grid_width].t + local_TIV[(j-1)+i*local_grid_width].t), 0.0);

          			local_TIV_next[j + i*local_grid_width].i = max(local_TIV[j+i*local_grid_width].i * (b1 +b2 * local_TIV[j+i*local_grid_width].t) 
									+ b3 * local_TIV[j+i*local_grid_width].t *local_TIV[j+i*local_grid_width].v
                							+ b4 * (local_TIV[j+(i+1)*local_grid_width].i + local_TIV[j+(i-1)*local_grid_width].i 
									+ local_TIV[(j+1)+i*local_grid_width].i + local_TIV[(j-1)+i*local_grid_width].i), 0.0);

                		local_TIV_next[j + i*local_grid_width].v = max(local_TIV[j + i*local_grid_width].v * c1 + c2 * local_TIV[j + i*local_grid_width].i
                							+c3 * (local_TIV[j+(i+1)*local_grid_width].v + local_TIV[j+(i-1)*local_grid_width].v  
									+local_TIV[(j+1)+i*local_grid_width].v + local_TIV[(j-1)+i*local_grid_width].v), 0.0);
                	}
		}
		//Wait for (last time's [n]) Sends:
		if(n > 0) { //Not for the first step (no previous steps)
			MPI_Waitall(neighbors, sends, MPI_STATUSES_IGNORE);
		}
		//Asynchronous Receives and Sends:
		neighbors = 0;
		if(up){
			receives[neighbors] = MPI::COMM_WORLD.Irecv(&local_TIV_next[1], 1, row_tiv, up_proc, 10); 
			sends[neighbors] = MPI::COMM_WORLD.Isend(&local_TIV_next[1+local_grid_width], 1, row_tiv, up_proc, 10); 
			++neighbors;
		}
		if(down){
			receives[neighbors] = MPI::COMM_WORLD.Irecv(&local_TIV_next[1 + (local_grid_height-1)*local_grid_width], 1, row_tiv, down_proc, 10); 
			sends[neighbors] = MPI::COMM_WORLD.Isend(&local_TIV_next[1 + (local_grid_height-2)*local_grid_width], 1, row_tiv, down_proc, 10); 	
			++neighbors;
		}
		if(left){
			receives[neighbors] = MPI::COMM_WORLD.Irecv(&local_TIV_next[local_grid_width], 1, col_tiv, left_proc, 10); 
			sends[neighbors] = MPI::COMM_WORLD.Isend(&local_TIV_next[1+local_grid_width], 1, col_tiv, left_proc, 10); 
			++neighbors;
		}
		if(right){
			receives[neighbors] = MPI::COMM_WORLD.Irecv(&local_TIV_next[2*local_grid_width - 1], 1, col_tiv, right_proc, 10); 
			sends[neighbors] = MPI::COMM_WORLD.Isend(&local_TIV_next[2*local_grid_width-2], 1, col_tiv, right_proc, 10); 
			++neighbors;
		}
		//Swapping TIV and TIV_next pointers
		tiv* temp = local_TIV;
		local_TIV = local_TIV_next;
		local_TIV_next = temp;
		//Wait for (THIS time's [n+1]) Receives:
		MPI_Waitall(neighbors, receives, MPI_STATUSES_IGNORE);
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
