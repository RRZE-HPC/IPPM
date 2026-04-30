program hello
    use mpi_f08
    implicit none
    integer:: irank, nrank, ierr
    !include "mpif.h"

    call MPI_FIXME(ierr)
    call MPI_COMM_FIXME(MPI_COMM_WORLD,nrank,ierr)
    call MPI_COMM_RANK(FIXME,FIXME,ierr)

    write(*,'(2(a,i))') "Hello World! I am ",irank," of ",nrank

    call MPI_FIXME(ierr)
end program hello
