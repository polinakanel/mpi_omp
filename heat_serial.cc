#include<fstream>
#include<stdio.h>
#include<string>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include <iostream>

using namespace std;

double sq(double a){
  return a*a;
}

int main(int argc, char *argv[])
{
  clock_t tclock;
  tclock=clock();
  if(argc!=2){
  printf("input parameters %s <grid size>\n", argv[0]);
  exit(1);
}

  const int siz=atof(argv[1]);
  const double kappa=1;
  double average_temp=0;

  double** fd_grid = new double*[siz];

  for (int i=0;i<siz;i++){
	  fd_grid[i]=new double[siz];
  }

  double** fd_grid_t = new double*[siz];
  for (int i=0;i<siz;i++){
	  fd_grid_t[i]=new double[siz];
  }

  double dx=M_PI/siz; // dx = dy
  const double dt=sq(dx)/(8*kappa);
  const double time=0.5*sq(M_PI)/kappa;
  const double steps=time/dt;


    for(int i=0;i<siz;i++){
	    for(int j=0;j<siz;j++){
	    	fd_grid[i][j]=0;
	    }
    }

  // setting boundaries
  for (int i=0;i<siz;i++){
	  fd_grid[i][0]		=sq(cos(i*dx));
	  fd_grid[i][siz-1]	=sq(sin(i*dx));

  }

  //printf("Size is %d \n", siz);

    for(int i=0;i<siz;i++){
	    for(int j=0;j<siz;j++){
	    	fd_grid_t[i][j]=fd_grid[i][j];
	    }
    }


  for( int timestep=0; timestep<steps; timestep++ ){

	  for(int i=1;i<siz-1;i++){
		  for(int j=1;j<siz-1;j++){
			  fd_grid_t[i][j]=fd_grid[i][j]+ kappa*dt*(fd_grid[i-1][j]+fd_grid[i+1][j]+fd_grid[i][j-1]+fd_grid[i][j+1]-4*fd_grid[i][j])/sq(dx); } }

    for(int i=1;i<siz-1;i++){
    	fd_grid_t[0][i]=fd_grid[0][i]+ kappa*dt*(fd_grid[siz-1][i]+fd_grid[1][i]+fd_grid[0][i-1]+fd_grid[0][i+1]-4*fd_grid[0][i])/sq(dx); }

    for(int i=1;i<siz-1;i++){
    	fd_grid_t[siz-1][i]=fd_grid[siz-1][i]+ kappa*dt*(fd_grid[siz-2][i]+fd_grid[0][i]+fd_grid[siz-1][i-1]+fd_grid[siz-1][i+1]-4*fd_grid[siz-1][i])/sq(dx); }

    // update the grid with the next step approximation
    for(int i=0;i<siz;i++){
	    for(int j=0;j<siz;j++){
	    	fd_grid[i][j]=fd_grid_t[i][j]; } }

  	  	  	  	  	  	  	  }


   for(int i=0;i<siz;i++){
	    for(int j=0;j<siz;j++){
	    	average_temp+=fd_grid[i][j];
	    }
    }

	   tclock=clock()-tclock;
	   printf("Total time taken is %15.8f \n", float(tclock)/CLOCKS_PER_SEC);
	   printf("Average temperature is %15.8f \n", average_temp/(siz*siz));

  // collect the final grid data
  char filename[50];
  sprintf(filename,"map_serial_%d.txt",siz);
  ofstream fout(filename);
    fout<<"dx = dy = "<< dx<< endl;
    for(int i=0;i<siz;i++){
	    for(int j=0;j<siz;j++){
	      fout<< fd_grid[i][j]<<" ";
	    }fout<<endl;
    }

    fout.close();
}
