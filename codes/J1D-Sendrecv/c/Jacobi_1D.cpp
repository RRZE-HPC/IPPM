#include <stdio.h>
#include <math.h>
#include <iostream>
#include "utils.h"

int main() {
    int i, n, istep, nstep, irank, nrank, ierr, is, ie, nl;
    double h, errmax;
    double *pot, *wa;
    irank=0; //HERE
    nrank=1; //HERE
    nstep=100000;
    n=100001;
    h=1.0/(double) (n-1);
    get_share_per_proc(n-1,irank,nrank,&is,&ie);
    //is=1;
    //ie=n-2;
    nl=ie-is+3;
    printf("%8d %8d %8d\n",is,ie,nl);
    pot=new double[nl];
    wa=new double[nl];
    for(i=0;i<nl;i++) pot[i]=pow((double)(is-1+i)*h,2);
    pot[0]=0.0; //HERE
    pot[n-1]=1.0; //HERE
    wa[0]=0.0; //HERE
    wa[n-1]=1.0; //HERE
    for(istep=0;istep<=nstep;istep++) {
        for(i=1;i<nl-1;i++) wa[i]=0.5*(pot[i+1]+pot[i-1]);
        //HERE
        std::swap(pot,wa);
    }
    dump_domain(irank,nrank,nstep,nl,pot,n-2);
    errmax=0.0;
    for(i=1;i<nl-1;i++) {
        //printf("%14.5E%24.15E\n",(double)i*h,pot[i]);
        errmax=std::max(errmax,abs(pot[i]-(double)(is-1+i)*h));
      }
    printf("Maximum difference relative to analytical solution=%14.5E\n",errmax);
    delete pot;
    delete wa;
    return 1;
}
