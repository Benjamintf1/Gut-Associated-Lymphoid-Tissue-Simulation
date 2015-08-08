Gut-Associated-Lymphoid-Tissue-Simulation
=========================================

This is a simulation of HIV or other viruses in gut associated lymphoid tissue with a focus on parallelizing for large data sets.

To Compile: Simply run the make file. This should compile the main program as well as a simple generator. If you want to create your own generator, you have to add that to the
makefile if you want to be able to automate the compilation of it, which will be described in the generator section. 

To Run: Before you run the program you need to generate the data files, use ./gen_simple_files HIV.conf for sayers, ./gen_simple_files Scaling.conf for mio

You can either use slurm, or mpiexec to run this program. 

To run it using slurm, you can use the sbatch on the predefined slurm file mio_sim. 

To run it on sayers or another mpi compatible machine, mpiexec -n (number of cores to run on) ./sim_hybrid HIV.conf (number of y nodes) (number of x nodes) for example
mpiexec -n 4 ./sim_hybrid HIV.conf 2 2

If you want to use openmp, export the OMP_NUM_THREADS to equal the number of threads to use. 

The mandatory arguments of the program are (in order), a configuration file, the number of nodes to divide the grid in the x direction, and the number of nodes to divide the grid in the y direction. *Keep in mind the grid given(minus the boundry conditions) must be divisable evenly by the x and y divisors, or the program will output an error to stderr and stop running.*

Configuration Details:

The grid width and height include the boundary spaces.

A configuration file should be a file with one double, integer, or string value per line in the following order:
* (double) delta_space
* (integer) grid_width
* (integer) grid_height
* (string) birth_rate_filename:
* (string) t_filename:
* (string) i_filename:
* (string) v_filename:
* (string) tissue_type_filename
* (double) delta_t:
* (integer) number_of_timesteps:
* (string) result_t_filename:
* (string) result_i_filename:
* (string) result_v_filename:
* (integer) number_of_tissue_types
* (string) sub_config_for_tissue_type_0
* ...
* (string) sub_config_for_tissue_type_n

And for each tissue type, a sub-configuration file (with a name matching what is listed in the main configuration file) should be a file with one double or integer per line in the following order (specialized for the tissue type):
* (double) diffusion_tcells:
* (double) diffusion_infected:
* (double) diffusion_virus:
* (double) death_tcells:
* (double) death_infected:
* (double) death_virus:
* (integer) burst_rate:
* (double) transmission_vt:
* (double) transmission_it:

Visualization: An example way to visualize the data has been included. Use graph_all.m and give the function the configuartion file, for example HIV.conf, or Scaling.conf. If you want to make your own, look at that, and the data files should just be a binary series of the double population values.

Generator: To create your own generator function, make a cpp that implements the funtions described in generators.h(which simply ask for a poplation at a given point, or which tissue type it's using), then compile it to a object file and compile generate_tivb_files.cpp using that object file. 





