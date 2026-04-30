program subarray
    use mpi_f08
    implicit none
    integer:: nrank, irank
    integer:: array(0:4,0:3)
    integer:: dim_array(2), dim_subarray(2), ind_subarray(2)
    type(MPI_Datatype):: type_subarray
 
    call MPI_INIT()
    call MPI_COMM_RANK(MPI_COMM_WORLD,irank)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,nrank)
    if(nrank/=2) then
        write(*,'(a)') 'ERROR: this program must be run with 2 processes!'
        call MPI_ABORT(MPI_COMM_WORLD,-1)
    endif

    !Creating a new datatype on both processes
    dim_array=[5,4]
    dim_subarray=[3,2]
    ![1,1] is index of second row and and second column, regardless of how
    !actually indices in definition of array
    ind_subarray=[1,1]
    call MPI_TYPE_CREATE_SUBARRAY(2,dim_array,dim_subarray,ind_subarray, &
        MPI_ORDER_FORTRAN,MPI_INTEGER,type_subarray)
    call MPI_TYPE_COMMIT(type_subarray)
 
    if(irank==0) then
        array(0,0)= 1 ; array(0,1)= 6 ; array(0,2)=11 ; array(0,3)=16
        array(1,0)= 2 ; array(1,1)= 7 ; array(1,2)=12 ; array(1,3)=17
        array(2,0)= 3 ; array(2,1)= 8 ; array(2,2)=13 ; array(2,3)=18
        array(3,0)= 4 ; array(3,1)= 9 ; array(3,2)=14 ; array(3,3)=19
        array(4,0)= 5 ; array(4,1)=10 ; array(4,2)=15 ; array(4,3)=20
        write(*,'(a)') 'The array on rank=0:'
        write(*,'(4i3)') array(0,0),array(0,1),array(0,2),array(0,3)
        write(*,'(4i3)') array(1,0),array(1,1),array(1,2),array(1,3)
        write(*,'(4i3)') array(2,0),array(2,1),array(2,2),array(2,3)
        write(*,'(4i3)') array(3,0),array(3,1),array(3,2),array(3,3)
        write(*,'(4i3)') array(4,0),array(4,1),array(4,2),array(4,3)
        call MPI_SEND(array,1,type_subarray,1,22,MPI_COMM_WORLD)
    else
        array=0
        call MPI_RECV(array,1,type_subarray,0,22,MPI_COMM_WORLD,MPI_STATUS_IGNORE)
        write(*,'(a)') 'The array on rank=1:'
        write(*,'(4i3)') array(0,0),array(0,1),array(0,2),array(0,3)
        write(*,'(4i3)') array(1,0),array(1,1),array(1,2),array(1,3)
        write(*,'(4i3)') array(2,0),array(2,1),array(2,2),array(2,3)
        write(*,'(4i3)') array(3,0),array(3,1),array(3,2),array(3,3)
        write(*,'(4i3)') array(4,0),array(4,1),array(4,2),array(4,3)
    endif
    call MPI_TYPE_FREE(type_subarray)
    call MPI_FINALIZE()
end program subarray
