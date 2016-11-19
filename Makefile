CC=gcc-4.8
CCMPI = mpic++
CXXFLAGS = -lstdc++


all: heat_serial heat_omp heat_mpi


heat_serial: heat_serial.cc
	$(CC) -o $@ $^ $(CXXFLAGS)

heat_omp: heat_omp.cc
	$(CC) -fopenmp -o $@ $^ $(CXXFLAGS)

heat_mpi: heat_mpi.cc
	$(CCMPI) -o $@ $<
