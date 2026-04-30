#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include <stdio.h>
int main(int argc, char **argv) {
    int irank_w, nrank_w, icolor, ikey, irank_s, nrank_s, ngroup;
    int irank_w_of_irank_s_equal_zero;
    MPI_Comm subcomm;
    ngroup=3;
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD,&nrank_w);
    MPI_Comm_rank(MPI_COMM_WORLD,&irank_w);

    icolor=10+irank_w%ngroup;
    //key can be used to control of rank assignment within the subcommunicator,
    //setting the keys of a given color the same way, as done in this example, then
    //all the processes in that given color will have the relative rank order as they
    //did in their parent group.
    ikey=1;
    MPI_Comm_split(MPI_COMM_WORLD,FIXME,FIXME,FIXME);
    MPI_Comm_rank(subcomm,&irank_s);
    MPI_Comm_size(subcomm,FIXME);

    //Let's call irank_s=0 of a subcomm the master process.
    //We want all processes of a subcomm to know
    //the irank_w of its master process.
    //The following two lines do it:
    if(irank_s==0) irank_w_of_irank_s_equal_zero=irank_w;
    MPI_Bcast(&irank_w_of_irank_s_equal_zero,1,MPI_INT,0,FIXME);
    printf("%5d     rank in WORLD: %4d%4d         rank in SUBCOMM: %4d%4d\n",
            irank_w_of_irank_s_equal_zero,irank_w,nrank_w,irank_s,nrank_s);

    MPI_Comm_free(FIXME);
    MPI_Finalize();
}
