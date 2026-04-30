#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include <stdio.h>
int main(int argc, char **argv) {
    int nrank, irank;
    int array[5][4];
    int dim_array[2], dim_subarray[2], ind_subarray[2];
    MPI_Datatype type_subarray;
 
    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&irank);
    MPI_Comm_size(MPI_COMM_WORLD,&nrank);
    if(nrank!=2) {
        printf("ERROR: this program must be run with 2 processes!\n");
        MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
    }

    //Creating a new datatype on both processes
    dim_array[0]=5 ; dim_array[1]=4;
    dim_subarray[0]=3 ; dim_subarray[1]=2;
    ind_subarray[0]=1 ; ind_subarray[1]=1;
    MPI_Type_create_subarray(2,dim_array,dim_subarray,ind_subarray,
        MPI_ORDER_C,MPI_INT,&type_subarray);
    MPI_Type_commit(&type_subarray);
 
    if(irank==0) {
        array[0][0]= 1 ; array[0][1]= 6 ; array[0][2]=11 ; array[0][3]=16 ;
        array[1][0]= 2 ; array[1][1]= 7 ; array[1][2]=12 ; array[1][3]=17 ;
        array[2][0]= 3 ; array[2][1]= 8 ; array[2][2]=13 ; array[2][3]=18 ;
        array[3][0]= 4 ; array[3][1]= 9 ; array[3][2]=14 ; array[3][3]=19 ;
        array[4][0]= 5 ; array[4][1]=10 ; array[4][2]=15 ; array[4][3]=20 ;
        printf("%3d%3d%3d%3d\n",array[0][0],array[0][1],array[0][2],array[0][3]);
        printf("%3d%3d%3d%3d\n",array[1][0],array[1][1],array[1][2],array[1][3]);
        printf("%3d%3d%3d%3d\n",array[2][0],array[2][1],array[2][2],array[2][3]);
        printf("%3d%3d%3d%3d\n",array[3][0],array[3][1],array[3][2],array[3][3]);
        printf("%3d%3d%3d%3d\n",array[4][0],array[4][1],array[4][2],array[4][3]);
        MPI_Send(array,1,type_subarray,1,22,MPI_COMM_WORLD);
    }
    else {
        for(int i=0;i<5;i++) for(int j=0;j<4;j++) array[i][j]=0;
        MPI_Recv(array,1,type_subarray,0,22,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        printf("The array on rank=1:\n");
        printf("%3d%3d%3d%3d\n",array[0][0],array[0][1],array[0][2],array[0][3]);
        printf("%3d%3d%3d%3d\n",array[1][0],array[1][1],array[1][2],array[1][3]);
        printf("%3d%3d%3d%3d\n",array[2][0],array[2][1],array[2][2],array[2][3]);
        printf("%3d%3d%3d%3d\n",array[3][0],array[3][1],array[3][2],array[3][3]);
        printf("%3d%3d%3d%3d\n",array[4][0],array[4][1],array[4][2],array[4][3]);
    }
    MPI_Type_free(&type_subarray);
    MPI_Finalize();
}
