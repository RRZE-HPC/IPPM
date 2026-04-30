#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "div.h"

#include <mpi.h>

int main(int argc,char **argv) {
    double time_s, time_e, x, y, pi, one_over_rand_max, count, val, count_g;
    unsigned int seed=123;
    int nrank, irank, i, nn, n_local;
    MPI_Status status;
    nn=pow(10,9);
  
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD,&nrank);
    MPI_Comm_rank(MPI_COMM_WORLD,&irank);
  
    n_local=nn/nrank;
    if(irank<nn%nrank) n_local+=1;
    // printf("%d : n_local = %d\n",irank,n_local);
     
    time_s=getTimeStamp();
  
    one_over_rand_max=1.0/(double)RAND_MAX;
    seed=2+irank;
  
    count=0.0;
    for(i=0; i<n_local; ++i) {
      x=rand_r(&seed)*one_over_rand_max;
      y=rand_r(&seed)*one_over_rand_max;
      if(x*x+y*y <1.0) ++count;
    }  
  
    time_e=getTimeStamp();
  
    MPI_FIXME(FIXME,&count_g,FIXME,MPI_DOUBLE,FIXME,FIXME,MPI_COMM_WORLD);
    if(0==irank) {
      pi = 4.0*FIXME/(double)nn; 
      printf("Time: %.3lf sec, accuracy: %14.5E\n",time_e-time_s,fabs(M_PI-pi)/M_PI);
    }
  
    MPI_Finalize();
  
    return 0;
}
