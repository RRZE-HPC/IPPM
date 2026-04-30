#include <mpi.h>
#include <cstdio>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <cstring>
#include "utils.h"

int main(int argc,char *argv[]) {
    int i, n, istep, nstep, irank, nrank, is, ie, nl, irank_l, irank_r, istep_print;
    double h, errmax, gmax, gmax_global;
    double *pot, *wa;
    int flag;
    MPI_Status statuses[4], status_conv;
    MPI_Request req[4], req_conv;
    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&irank);
    MPI_Comm_size(MPI_COMM_WORLD,&nrank);
    nstep=pow(10,7);
    n=1001;
    h=1.0/(double) (n-1);
    get_share_per_proc(n-1,irank,nrank,&is,&ie);
    nl=ie-is+3;
    //printf("%8d %8d %8d\n",is,ie,nl);
    //rank of the process working on the left chunk
    irank_l=irank-1;
    //rank of the process working on the right chunk
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
    istep_print=0;
    flag=true;
    for(istep=0;istep<=nstep;istep++) {
        //first computation of start and end elements of the chunk on this rank
        wa[1   ]=0.5*(pot[1+1]+pot[1-1]);
        wa[nl-2]=0.5*(pot[nl-2+1]+pot[nl-2-1]);
        if(nrank>1) {
            MPI_Isend(&wa[1   ],1,MPI_DOUBLE,irank_l,11,MPI_COMM_WORLD,&req[0]);
            MPI_Irecv(&wa[nl-1],1,MPI_DOUBLE,irank_r,11,MPI_COMM_WORLD,&req[1]);
            MPI_Isend(&wa[nl-2],1,MPI_DOUBLE,irank_r,22,MPI_COMM_WORLD,&req[2]);
            MPI_Irecv(&wa[0   ],1,MPI_DOUBLE,irank_l,22,MPI_COMM_WORLD,&req[3]);
        }
        for(i=2;i<nl-2;i++) wa[i]=0.5*(pot[i+1]+pot[i-1]);
        if(nrank>1) {
            MPI_Waitall(4,req,statuses);
        }
        if(istep%10==0) {
            if(flag) {
                gmax=0.0;
                for(i=1;i<nl-1;i++) {
                    gmax=std::max<double>(gmax,std::abs(pot[i+1]-2.0*pot[i]+pot[i-1]));
                }
                MPI_Iallreduce(FIXME,FIXME,FIXME,MPI_DOUBLE,
                    FIXME,MPI_COMM_WORLD,FIXME);
                flag=false;
            }
            else {
                MPI_Wait(FIXME,&status_conv);
                gmax_global/=pow(h,2);
                flag=true;
            }
            if(irank==0 && flag && (istep-istep_print>pow(10,5) || gmax_global<1.E-5)) {
                printf("Convergence= %14.5E%10d\n",gmax_global,istep);
                istep_print=istep;
            }
            if(flag && gmax_global<1.E-5) break;
        }
        std::swap(pot,wa);
    }
    if(!flag) {
        MPI_FIXME(&req_conv,&status_conv);
        gmax_global/=pow(h,2);
    }
    dump_domain(irank,nrank,nstep,nl,pot,n-2);
    errmax=0.0;
    for(i=1;i<nl-1;i++) {
        //printf("%14.5E%24.15E\n",(double)i*h,pot[i]);
        errmax=std::max<double>(errmax,std::abs(pot[i]-(double)(is-1+i)*h));
    }
    printf("Maximum difference relative to analytical solution=%14.5E%6d%14.5E\n",errmax,irank,gmax_global);
    delete pot;
    delete wa;
    MPI_Finalize();
    return 1;
}
