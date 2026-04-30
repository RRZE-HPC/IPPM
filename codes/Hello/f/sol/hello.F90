program hello
    use mpi_f08
    implicit none
    integer:: irank, nrank, ierr
    !include "mpif.h"

    call MPI_INIT(ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,nrank,ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD,irank,ierr)

    write(*,'(2(a,i))') "Hello World! I am ",irank," of ",nrank

    call MPI_FINALIZE(ierr)
end program hello
