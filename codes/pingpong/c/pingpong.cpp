#include <mpi.h>
#include <stdio.h>
int main(int argc, char **argv) {
    int ierr, irank, nrank;
    MPI_Status status;
    double d=0.0;
    ierr=MPI_Init(&argc,&argv);
    ierr=MPI_Comm_rank(MPI_COMM_WORLD,&irank);
    ierr=MPI_Comm_size(MPI_COMM_WORLD,&nrank);
    if(irank==0) d=100.0;
    if(irank==1) d=200.0;
    printf("BEFORE: nrank,irank,d = %5d%5d%8.1f\n",nrank,irank,d);
    if(irank==0) {
        MPI_Send(&d,1,MPI_DOUBLE,1,11,MPI_COMM_WORLD);
        MPI_Recv(&d,1,MPI_DOUBLE,1,22,MPI_COMM_WORLD,&status);
    }
    else if(irank==1) {
        MPI_Recv(&d,1,MPI_DOUBLE,0,11,MPI_COMM_WORLD,&status);
        MPI_Send(&d,1,MPI_DOUBLE,0,22,MPI_COMM_WORLD);
    }
    printf("AFTER:  nrank,irank,d = %5d%5d%8.1f\n",nrank,irank,d);
    ierr=MPI_Finalize();
}

