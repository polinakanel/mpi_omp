#include<iostream>
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<mpi.h>
#include<time.h>

double sq(double a){
  return a*a;
}


int main(int argc, char *argv[])
{

  if(argc!=2){
	  printf("input parameters %s <grid size>\n", argv[0]);
	  exit(1);
	}

  const int siz=atof(argv[1]);
  const double kappa=1;

  // time measurement
  clock_t tclock;
  tclock=clock();
  double average_temp = 0;

  // set up MPI
  int nproc;
  int rank;
  MPI_Status stats[4];
  MPI_Request request[4];
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);


  // decompose domain into chunks
  int dim = siz/nproc;

  double** fd_grid = new double*[dim];
  for (int i=0;i<siz;i++){
	  fd_grid[i]=new double[siz];
  }

    for(int i=0;i<dim;i++){
	    for(int j=0;j<siz;j++){
	    	fd_grid[i][j]=0;
	    }
    }
    // boundary conditions
  for(int i=0;i<dim;i++){
	  fd_grid[i][siz-1]=sq(sin((i + rank*dim)*M_PI/siz));
	  fd_grid[i][0]=sq(cos((i + rank*dim)*M_PI/siz));
  }

  double* buff_send_next=new double[siz];
  double* buff_send_prev=new double[siz];
  double* buff_reveive_next=new double[siz];
  double* buff_reveive_prev=new double[siz];
  // buffer grid
  double** fd_grid_t = new double*[dim];
  for (int i=0;i<dim;i++){
	  fd_grid_t[i]=new double[siz];
  }

  // FD parameters
  double dx=M_PI/siz;
  const double dt=sq(dx)/(8*kappa);
  const double time=0.5*sq(M_PI)/kappa;
  const double nsteps=time/dt;


    for(int i=0;i<dim;i++){
	    for(int j=0;j<siz;j++){
	    	fd_grid_t[i][j]=fd_grid[i][j];
	    					  } 	}

    // actual computation
  for(int tt=0;tt<nsteps;tt++){

	  for(int i=1;i<dim-1;i++){
		  for(int j=1;j<siz-1;j++){
			  fd_grid_t[i][j]=fd_grid[i][j]+ kappa*dt*(fd_grid[i-1][j]+fd_grid[i+1][j]+fd_grid[i][j-1]+fd_grid[i][j+1]-4*fd_grid[i][j])/sq(dx);
		  	  	  	  	  	  	  }		}

    for(int i=0;i< dim;i++){
	    for(int j=0;j<siz;j++){
		    fd_grid[i][j]=fd_grid_t[i][j];
	    					  } 	}

    for(int i=0;i<siz;i++){
    	buff_send_prev[i]=fd_grid[1][i];
    	buff_send_next[i]=fd_grid[dim - 1][i];
    					}
    // who send who
    int prev = rank-1;
    int next = rank+1;
    if (rank == 0) prev = nproc - 1;
    if (rank == (nproc - 1)) next = 0;

	int tag1=1, tag2=2;
	// exchange the rows for the next iteration
	MPI_Isend(buff_send_next, siz, MPI_DOUBLE, next , tag1, MPI_COMM_WORLD, &request[0]);
	MPI_Isend(buff_send_prev, siz, MPI_DOUBLE, prev , tag2, MPI_COMM_WORLD, &request[1]);
	MPI_Irecv(buff_reveive_next, siz, MPI_DOUBLE, next , tag2, MPI_COMM_WORLD, &request[2]);
	MPI_Irecv(buff_reveive_prev, siz, MPI_DOUBLE, prev , tag1, MPI_COMM_WORLD, &request[3]);
	MPI_Waitall(4, request, stats);

	for(int i=0;i<siz;i++){
		fd_grid[0][i]= buff_reveive_prev[i];
		fd_grid[siz/nproc+1][i]=buff_reveive_prev[i];
							}
  }

  printf("Average temperature is %15.8f \n", average_temp/(dim*siz));


  MPI_Finalize();

  tclock=clock()-tclock;
  printf("Total time taken is %15.8f \n", float(tclock)/CLOCKS_PER_SEC);

}
