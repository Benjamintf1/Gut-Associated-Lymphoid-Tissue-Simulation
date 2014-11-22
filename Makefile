all: sim_serial gen_files_simple sim_mpi sim_hybrid

sim_mpi: src/galt_mpi.cpp
	mpiCC src/galt_mpi.cpp -o sim_mpi

sim_hybrid: src/galt_hybrid.cpp
	mpiCC -openmp src/galt_hybrid.cpp -o sim_hybrid

sim_serial: src/galt_serial.cpp
	g++ src/galt_serial.cpp -o sim_serial

obj/generator1.o: src/generator1.cpp
	g++ src/generator1.cpp -c -o obj/generator1.o

gen_files_simple: obj/generator1.o src/generate_tivb_files.cpp
	g++ src/generate_tivb_files.cpp obj/generator1.o -o gen_files_simple




