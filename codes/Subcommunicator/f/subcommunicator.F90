program subcommunicator
    use mpi_f08
    implicit none
    integer:: irank_w, nrank_w, ierr, icolor, ikey, irank_s, nrank_s, ngroup
    integer:: irank_w_of_irank_s_equal_zero
    type(MPI_COMM):: subcomm
    ngroup=3
    call MPI_INIT(ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,nrank_w,ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD,irank_w,ierr)

    icolor=10+mod(irank_w,ngroup)
    !key can be used to control of rank assignment within the subcommunicator,
    !setting the keys of a given color the same way, as done in this example, then
    !all the processes in that given color will have the relative rank order as they
    !did in their parent group.
    ikey=1
    call MPI_COMM_SPLIT(MPI_COMM_WORLD,FIXME,FIXME,FIXME,ierr)
    call MPI_COMM_RANK(subcomm,irank_s,ierr)
    call MPI_COMM_SIZE(subcomm,FIXME,ierr)

    !Let's call irank_s=0 of a subcomm the master process.
    !We want all processes of a subcomm to know
    !the irank_w of its master process.
    !The following two lines do it:
    if(irank_s==0) irank_w_of_irank_s_equal_zero=irank_w
    call MPI_BCAST(irank_w_of_irank_s_equal_zero,1,MPI_INTEGER,0,FIXME,ierr)
    write (*,'(i5,5x,a,2i4,9x,a,2i4)') &
        irank_w_of_irank_s_equal_zero, &
        'rank in WORLD: ',irank_w,nrank_w, &
        'rank in SUBCOMM: ',irank_s,nrank_s

    call MPI_COMM_FREE(FIXME,ierr)
    call MPI_FINALIZE(ierr)
end program subcommunicator
