all: sim_s gen_files_simple

sim_s: serial.cpp
	g++ serial.cpp -o sim_s

gen_files_simple: generator1.o generate_tivb_files.cpp
	g++ generate_tivb_files.cpp generator1.o -o gen_files_simple

generator1.o: generator1.cpp
	g++ generator1.cpp -c -o generator1.o

clean:
	rm *.o
	rm *.dat
	rm sim_s
	rm gen_files_simple
