#include <mpi.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include "utils.h"

int main(int argc,char *argv[]) {
    int i, n, istep, nstep, irank, nrank, ierr, is, ie, nl, irank_l, irank_r;
    double h, errmax;
    double *pot, *wa;
    MPI_Status status;
    ierr=MPI_Init(&argc,&argv);
    ierr=MPI_Comm_rank(MPI_COMM_WORLD,&irank);
    ierr=MPI_Comm_size(MPI_COMM_WORLD,&nrank);
    nstep=100000;
    n=100001;
    h=1.0/(double) (n-1);
    get_share_per_proc(n-1,irank,nrank,&is,&ie);
    //is=1;
    //ie=n-2;
    nl=ie-is+3;
    printf("%8d %8d %8d\n",is,ie,nl);
    irank_l=irank-1;
    irank_r=irank+1;
    if(irank_l<0) irank_l=MPI_PROC_NULL;
    if(irank_r>=nrank) irank_r=MPI_PROC_NULL;
    pot=new double[nl];
    wa=new double[nl];
    for(i=0;i<nl;i++) pot[i]=pow((double)(is-1+i)*h,2);
    if(irank==0) pot[0]=0.0;
    if(irank==nrank-1) pot[nl-1]=1.0;
    if(irank==0) wa[0]=0.0;
    if(irank==nrank-1) wa[nl-1]=1.0;
    for(istep=0;istep<=nstep;istep++) {
        for(i=1;i<nl-1;i++) wa[i]=0.5*(pot[i+1]+pot[i-1]);
        if(nrank>1) {
            ierr=MPI_Sendrecv(&wa[1],1,MPI_DOUBLE_PRECISION,irank_l,11,&wa[nl-1],1,MPI_DOUBLE_PRECISION,irank_r,11,MPI_COMM_WORLD,&status);
            ierr=MPI_Sendrecv(&wa[nl-2],1,MPI_DOUBLE_PRECISION,irank_r,22,wa,1,MPI_DOUBLE_PRECISION,irank_l,22,MPI_COMM_WORLD,&status);
        }
        std::swap(pot,wa);
    }
    dump_domain(irank,nrank,nstep,nl,pot,n-2);
    errmax=0.0;
    for(i=1;i<nl-1;i++) {
        //printf("%14.5E%24.15E\n",(double)i*h,pot[i]);
        errmax=std::max(errmax,abs(pot[i]-(double)(is-1+i)*h));
      }
    printf("Maximum difference relative to analytical solution=%14.5E%6d\n",errmax,irank);
    delete pot;
    delete wa;
    ierr=MPI_Finalize();
    return 1;
}
