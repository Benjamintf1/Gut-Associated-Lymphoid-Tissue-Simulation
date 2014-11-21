all: sim_s gen_files_simple sim_mpi

sim_mpi: galt_mpi.cpp
	mpiCC galt_mpi.cpp -o sim_mpi

sim_s: galt_serial.cpp
	g++ galt_serial.cpp -o sim_s

gen_files_simple: generator1.o generate_tivb_files.cpp
	g++ generate_tivb_files.cpp generator1.o -o gen_files_simple

generator1.o: generator1.cpp
	g++ generator1.cpp -c -o generator1.o

clean:
	rm *.o
	rm *.dat
	rm sim_s
	rm gen_files_simple
