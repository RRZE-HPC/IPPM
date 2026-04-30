#include <mpi.h>
#include <cstdio>
#include <cmath>
#include <iostream>



void get_share_per_proc(int n,int irank,int nrank,int *is,int *ie) {
    int mq, mr;
    mq=(n-1)/nrank;
    mr=(n-1)%nrank;
    is[0]=irank*mq+1;
    if(irank<mr) mq=mq+1;
    is[0]=is[0]+std::min(irank,mr);
    ie[0]=is[0]+mq-1;
}
void dump_domain(int irank,int nrank,int nstep,int nl,double *pot,int n) {
    //int 
    //int 
    char *filename=NULL;
    int ierr=asprintf(&filename,"domain-%08d.pgm",nstep);
    if(ierr==-1) {
        fprintf(stderr,"[%d] ERROR: asprintf failed.\n",irank);
#if !defined(NOMPI)
        MPI_Abort(MPI_COMM_WORLD,1);
#endif
    }
    if(irank==0) {
        FILE *f=fopen(filename,"wb");
        if(f==NULL) {
            fprintf(stderr, "[%d] ERROR: open file for writing failed.\n",irank);
#if !defined(NOMPI)
            MPI_Abort(MPI_COMM_WORLD, 1);
#endif
        }
        fprintf(f, "P5\n%d\n255\n",n);
        fclose(f);
    }
#if !defined(NOMPI)
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    for(int jrank=0;jrank<nrank;jrank++) {
        if(jrank==irank) {
            FILE *f=fopen(filename,"ab");
            if(f==NULL) {
                fprintf(stderr, "[%d] ERROR: open file for writing failed.\n",irank);
#if !defined(NOMPI)
                MPI_Abort(MPI_COMM_WORLD, 1);
#endif
            }
            for (int i=1;i<nl-1;i++) {
                    unsigned char gray = (unsigned char)(255.0*(pot[i]));
                    fwrite(&gray,sizeof(gray),1,f);
            }
            fclose(f);
        }
#if !defined(NOMPI)
        MPI_Barrier(MPI_COMM_WORLD);
#endif
    }
    free(filename);
}
