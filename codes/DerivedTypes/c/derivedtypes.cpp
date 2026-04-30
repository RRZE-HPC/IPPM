#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include <stdio.h>

struct typ_atom {
    //char useless1[1];
    double xat;
    //char useless2[1];
    char sat[3];
    char useless3[5];
};

int main(int argc, char **argv) {
    int rank, size, ierr, lengths[2], size_notresized, size_resized;
    struct typ_atom atom[2];
    MPI_Datatype types_list[2], typ_atom_mpi_notresized, typ_atom_mpi;
    MPI_Status status;
    MPI_Aint displacements[2], tmp_lb, extent_notresized, extent_true_notresized;
    MPI_Aint extent_resized, extent_true_resized, address_1, address_2;
    bool do_it_unsafe=false;

    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD,&size);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    if(size!=2) {
        printf("This exercise is supposed to be run with 2 MPI processes!\n");
        MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
    }
    lengths[0]=1;
    lengths[1]=3;
    MPI_Get_address(&atom[0].xat,&address_1);
    MPI_Get_address(&atom[0].sat,&address_2);
    displacements[0]=(MPI_Aint) 0;
    displacements[1]=MPI_Aint_diff(address_2,address_1);
    types_list[0]=FIXME;
    types_list[1]=FIXME;
    MPI_Type_create_struct(2,lengths,displacements,types_list,&typ_atom_mpi_notresized);
    if(do_it_unsafe) {
        typ_atom_mpi=typ_atom_mpi_notresized;
    }
    else {
        MPI_Type_create_resized(typ_atom_mpi_notresized,FIXME,
                FIXME,&typ_atom_mpi);
    }
    MPI_Type_commit(&typ_atom_mpi);
    if(rank==0) {
        printf("------------------------------------------------\n");
        MPI_Type_size(typ_atom_mpi_notresized,&size_notresized);
        MPI_Type_get_extent(typ_atom_mpi_notresized,&tmp_lb,&extent_notresized);
        MPI_Type_get_true_extent(typ_atom_mpi_notresized,&tmp_lb,&extent_true_notresized);
        MPI_Type_size(typ_atom_mpi,&size_resized);
        MPI_Type_get_extent(typ_atom_mpi,&tmp_lb,&extent_resized);
        MPI_Type_get_true_extent(typ_atom_mpi,&tmp_lb,&extent_true_resized);
        printf("notresized: size_notresized=%5d\n",size_notresized);
        printf("notresized: extent_notresized=%5d\n",(int) extent_notresized);
        printf("notresized: extent_true_notresized=%5d\n",(int) extent_true_notresized);
        printf("------------------------------------------------\n");
        printf("resized: size_resized=%5d\n",size_resized);
        printf("resized: extent_resized=%5d\n",(int) extent_resized);
        printf("resized: extent_true_resized=%5d\n",(int) extent_true_resized);
        printf("------------------------------------------------\n");
        MPI_Get_address(&atom[0],&address_1);
        MPI_Get_address(&atom[1],&address_2);
        printf("size_from_diff=%5d\n",(int)(MPI_Aint_diff(address_2,address_1)));
        printf("------------------------------------------------\n");
    }
    if(rank==0) {
        atom[0].xat=10.0;
        strcpy(atom[0].sat,"Si");
        atom[1].xat=20.0;
        strcpy(atom[1].sat,"Ti");
        MPI_Send(FIXME,2,typ_atom_mpi,1,11,MPI_COMM_WORLD);
    }
    else {
        MPI_Recv(FIXME,2,typ_atom_mpi,0,11,MPI_COMM_WORLD,&status);
        printf("atom[0]: xat=%lf\n",atom[0].xat);
        printf("atom[0]: sat=%s\n",atom[0].sat);
        printf("atom[1]: xat=%lf\n",atom[1].xat);
        printf("atom[1]: sat=%s\n",atom[1].sat);
    }
    MPI_Finalize();
}
