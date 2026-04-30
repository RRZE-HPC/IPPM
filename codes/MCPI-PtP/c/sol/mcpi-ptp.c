#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "div.h"

#include <mpi.h>

int main(int argc,char **argv) {
  double time_s, time_e, x, y, pi, one_over_rand_max, count, val;
  unsigned int seed=123;
  int nrank, irank, i, nn, n_local;
  MPI_Status status;
  nn=pow(10,9);

  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&nrank);
  MPI_Comm_rank(MPI_COMM_WORLD,&irank);

  n_local=nn/nrank;
  //if nrank is not a divisor (factor) of nn, then some ranks should do one more point
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

  if(0==irank) {
    //This is rank=0 and it receives `count' from all other ranks
    for(i=1; i<nrank; ++i) {
      MPI_Recv(&val,1,MPI_DOUBLE,i,0,MPI_COMM_WORLD,&status);
      count+=val; 
    }
    pi=4.0*count/(double)nn; 
    printf("Time: %.3lf sec, accuracy: %14.5E\n",time_e-time_s,fabs(M_PI-pi)/M_PI);
  }
  else {
    //Every rank except rank=0 should send its `count' to rank=0
    MPI_Send(&count,1,MPI_DOUBLE,0,0,MPI_COMM_WORLD);
  }
 
  MPI_Finalize();

  return 0;
}
