#include <mpi.h>
#include <cstdio>
#include <cmath>
#include <iostream>
#include <cstring>
#include "utils.h"

void get_input_param(int,int,char **,bool *,char *);

int main(int argc,char *argv[]) {
    int i, n, istep, nstep, irank, nrank, is, ie, nl, irank_l, irank_r;
    double h, errmax, tt, gnrm;
    double *pot, *wa, *pot_global, *errmax_all, *gnrm_all;
    int n_loc_ext, n_tot_left, jrank, n_loc, ios;
    int *n_loc_all;
    char filename[256];
    bool restart;
    FILE *fptr;
    MPI_Status status;
    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&irank);
    MPI_Comm_size(MPI_COMM_WORLD,&nrank);
    get_input_param(irank,argc,argv,&restart,filename);
    if(restart) {
        if(irank==0) {
            fptr=fopen(filename,"r");
            if(fptr==NULL) {
                printf("ERROR: cannot open file %s\n",filename);
                exit(1);
            }
            fscanf(fptr,"%d",&n);
        }
        MPI_Bcast(&n,1,MPI_INT,0,MPI_COMM_WORLD);
    }
    else {
        n=10001;
    }
    nstep=10000000;
    h=1.0/(double) (n-1);
    get_share_per_proc(n-1,irank,nrank,&is,&ie);
    nl=ie-is+3;
    //printf("%8d %8d %8d\n",is,ie,nl);
    irank_l=irank-1;
    irank_r=irank+1;
    if(irank_l<0) irank_l=MPI_PROC_NULL;
    if(irank_r>=nrank) irank_r=MPI_PROC_NULL;
    pot=new double[nl];
    wa=new double[nl];
    n_loc=ie-is+1;
    n_loc_all=new int[nrank];
    for(jrank=0;jrank<nrank;jrank++) n_loc_all[jrank]=0;
    MPI_Allgather(&n_loc,1,MPI_INT,n_loc_all,1,MPI_INT,MPI_COMM_WORLD);
    //n_loc_ext is the extended value of n_loc that it is the largest
    //value of n_loc across all ranks
    n_loc_ext=0;
    for(jrank=0;jrank<nrank;jrank++) n_loc_ext=std::max<int>(n_loc_ext,n_loc_all[jrank]);
    //pot_global stores the potentials on all grid points
    pot_global=new double[n_loc_ext*nrank];
    if(restart) {
        if(irank==0) {
            for(jrank=0;jrank<nrank;jrank++) {
                for(i=0;i<n_loc_all[jrank];i++) {
                    fscanf(fptr,"%lf%lf",&tt,&pot_global[jrank*n_loc_ext+i]);
                }
            }
            fclose(fptr);
        }
        //rank=0 has read the potential on all grid points from the
        //restart file and now it should scatter across all ranks
        MPI_FIXME(pot_global,n_loc_ext,MPI_DOUBLE,
            &pot[1],FIXME,MPI_DOUBLE,FIXME,MPI_COMM_WORLD);
        if(nrank>1) {
            MPI_Send(&pot[1   ],1,MPI_DOUBLE,irank_l,11,MPI_COMM_WORLD);
            MPI_Recv(&pot[nl-1],1,MPI_DOUBLE,irank_r,11,MPI_COMM_WORLD,&status);
            MPI_Send(&pot[nl-2],1,MPI_DOUBLE,irank_r,22,MPI_COMM_WORLD);
            MPI_Recv(&pot[0   ],1,MPI_DOUBLE,irank_l,22,MPI_COMM_WORLD,&status);
        }
    }
    else {
        for(i=0;i<nl;i++) pot[i]=pow((double)(is-1+i)*h,2);
    }
    if(irank==0) pot[0]=0.0;
    if(irank==nrank-1) pot[nl-1]=1.0;
    if(irank==0) wa[0]=0.0;
    if(irank==nrank-1) wa[nl-1]=1.0;
    for(istep=0;istep<=nstep;istep++) {
        for(i=1;i<nl-1;i++) wa[i]=0.5*(pot[i+1]+pot[i-1]);
        if(nrank>1) {
            MPI_Send(&wa[1   ],1,MPI_DOUBLE,irank_l,11,MPI_COMM_WORLD);
            MPI_Recv(&wa[nl-1],1,MPI_DOUBLE,irank_r,11,MPI_COMM_WORLD,&status);
            MPI_Send(&wa[nl-2],1,MPI_DOUBLE,irank_r,22,MPI_COMM_WORLD);
            MPI_Recv(&wa[0   ],1,MPI_DOUBLE,irank_l,22,MPI_COMM_WORLD,&status);
        }
        std::swap(pot,wa);
    }
    //rank=0 will write the entire potential array on all grid points
    //into a restart file and thus it should gather chunks from all ranks
    MPI_FIXME(&pot[1],n_loc_ext,MPI_DOUBLE,
        pot_global,FIXME,MPI_DOUBLE,FIXME,MPI_COMM_WORLD);
    if(irank==0) {
        if(!restart) strcpy(filename,"starting_potential.dat");
        fptr=fopen(filename,"w");
        if(fptr==NULL) {
            printf("ERROR: cannot open file %s\n",filename);
            exit(1);
        }
        fprintf(fptr,"%8d\n",n-1);
        n_tot_left=0;
        for(jrank=0;jrank<nrank;jrank++) {
            for(i=0;i<n_loc_all[jrank];i++) {
                fprintf(fptr,"%26.17E%26.17E\n",(n_tot_left+i+1)*h,pot_global[jrank*n_loc_ext+i]);
            }
        n_tot_left+=n_loc_all[jrank];
        }
        fclose(fptr);
    }
    delete pot_global;
    delete n_loc_all;
    dump_domain(irank,nrank,nstep,nl,pot,n-2);
    errmax=0.0;
    gnrm=0.0;
    for(i=1;i<nl-1;i++) {
        //printf("%14.5E%24.15E\n",(double)i*h,pot[i]);
        errmax=std::max(errmax,std::abs(pot[i]-(double)(is-1+i)*h));
        gnrm=gnrm+pow(pot[i+1]-2.0*pot[i]+pot[i-1],2);
    }
    errmax_all=new double[nrank];
    gnrm_all=new double[nrank];
    MPI_Gather(&errmax,1,MPI_DOUBLE,errmax_all,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Gather(&gnrm,1,MPI_DOUBLE,gnrm_all,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    gnrm=0.0;
    errmax=0.0;
    for(jrank=0;jrank<nrank;jrank++) {
        errmax=std::max(errmax,errmax_all[jrank]);
        gnrm+=gnrm_all[jrank];
    }
    gnrm=pow(gnrm,0.5);
    if(irank==0) {
        printf("Convergence= %14.4E%14.4E\n",errmax,gnrm);
    }
    delete errmax_all;
    delete gnrm_all;
    delete pot;
    delete wa;
    MPI_Finalize();
    return 1;
}

void get_input_param(int irank,int argc, char **argv,bool *restart,char *filename) {
    //int ierr;
    //char argv[256];
    //printf("argc= %d\n",argc);
    //printf("argv0= %s\n",argv[0]);
    //printf("argv1= %s\n",argv[1]);
    //printf("argv2= %s\n",argv[2]);
    if(argc==1) {
        restart[0]=false;
    }
    else if(argc==2) {
        if(irank==0) {
            if(strcmp(argv[1],"-r")==0) {
                printf("ERROR: missing restart filename in the argument list!\n");
            }
            else {
                printf("ERROR: unknown 1st argument!\n");
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
    }
    else if(argc==3) {
        if(strcmp(argv[1],"-r")==0) {
            strcpy(filename,argv[2]);
            restart[0]=true;
        }
        else {
            if(irank==0) printf("ERROR: unknown 1st argument!\n");
            MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
        }
    }
    else {
        if(irank==0) printf("ERROR: too many arguments!\n");
        MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
    }
}
