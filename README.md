# mpi_omp

Code for the serial, mpi and omp execution of the heat equation using the finite difference approach. 

• To run serial version: ./heat_serial n , where n is the size of the grid (n^2)

• To run OpenMP version: ./heat_omp n nthreads , where n is the size of the gris (n^2) and nthreads is the number of the threads

• To run MPI version with domain decomposition: ./heat_mpi n , where n is the size of the grid (n^2)
