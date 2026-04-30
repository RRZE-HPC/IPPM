#include <mpi.h>
#include <stdio.h>
int main(int argc, char **argv) {
    int irank, nrank, ierr;

    ierr=MPI_Init(&argc,&argv);
    ierr=MPI_Comm_size(MPI_COMM_WORLD,&nrank);
    ierr=MPI_Comm_rank(MPI_COMM_WORLD,&irank);

    printf("Hello World! I am %d of %d\n",irank,nrank);

    ierr=MPI_Finalize();
}
