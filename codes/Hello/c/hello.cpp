#include <mpi.h>
#include <stdio.h>
int main(int argc, char **argv) {
    int irank, nrank, ierr;

    ierr=MPI_FIXME(FIXME,FIXME);
    ierr=MPI_Comm_FIXME(MPI_COMM_WORLD,&nrank);
    ierr=MPI_Comm_rank(FIXME,FIXME);

    printf("Hello World! I am %d of %d\n",irank,nrank);

    ierr=MPI_FIXME();
}
