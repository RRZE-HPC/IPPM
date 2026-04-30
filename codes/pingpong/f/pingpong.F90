program pingpong
    use mpi_f08
    implicit none
    integer:: irank, nrank, ierr
    real(kind=8):: d=0.d0;
    type(MPI_STATUS):: status
    call MPI_INIT(ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD,irank,ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,nrank,ierr)
    if(irank==0) d=100.d0
    if(irank==1) d=200.d0
    write(*,'(a,2i5,f8.1)') 'BEFORE: nrank,irank,d = ',nrank,irank,d
    if(irank==0) then
        call MPI_SEND(d,1,MPI_DOUBLE_PRECISION,1,11,MPI_COMM_WORLD,ierr)
        call MPI_RECV(d,1,MPI_DOUBLE_PRECISION,1,22,MPI_COMM_WORLD,status,ierr)
    elseif(irank==1) then
        call MPI_RECV(d,1,MPI_DOUBLE_PRECISION,0,11,MPI_COMM_WORLD,status,ierr)
        call MPI_SEND(d,1,MPI_DOUBLE_PRECISION,0,22,MPI_COMM_WORLD,ierr)
    endif
    write(*,'(a,2i5,f8.1)') 'AFTER:  nrank,irank,d = ',nrank,irank,d
    call MPI_FINALIZE(ierr)
end program pingpong
