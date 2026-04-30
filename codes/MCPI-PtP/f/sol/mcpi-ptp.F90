program mcpi
    use mpi
    use, intrinsic:: iso_c_binding
    implicit none
    interface
        integer(c_int) function rand_r(seed) bind(c)
            import:: c_int
            integer(c_int), intent(in):: seed
        end function rand_r
        real(kind=8) function getTimeStamp() bind(c,name="getTimeStamp")
        end function getTimeStamp
    end interface
    integer(c_int):: seed
    integer(c_int), parameter:: RAND_MAX=2147483647
    real(kind=8):: x , y, time_s, time_e, pi, M_PI, count, val , one_over_rand_max
    integer:: i, nn, n_local, nrank, irank, ierr, status(MPI_STATUS_SIZE)
    !M_PI=3.14159265358979d0
    M_PI=4.d0*atan(1.d0)
    nn=1*10**9

    call MPI_Init(ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,nrank,ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD,irank,ierr)

    n_local=nn/nrank
    !write(*,*) 'nn,n_local= ',nn,n_local
    !if nrank is not a divisor (factor) of nn, then some ranks should do one more point
    if(irank<mod(nn,nrank)) n_local=n_local+1

    count=0.d0
    time_s=getTimeStamp()
    seed=2+irank
    one_over_rand_max=1.d0/real(RAND_MAX,kind=8)
    do i=1,n_local
        x=rand_r(seed)*one_over_rand_max
        y=rand_r(seed)*one_over_rand_max
        if((x*x+y*y)<1.d0) count=count+1.d0
    enddo
    !write(*,*) 'irank,count= ',irank,count
    time_e=getTimeStamp()

    if(irank==0) then
        !This is rank=0 and it receives `count' from all other ranks
        do i=1,nrank-1
            call MPI_RECV(val,1,MPI_DOUBLE_PRECISION,i,0,MPI_COMM_WORLD,status,ierr)
            count=count+val
        enddo
        pi=4.d0*count/real(nn,kind=8)
        write(*,'(a,f10.3,5x,a,es14.5)') 'Time: ',time_e-time_s,' accuracy: ',abs(M_PI-pi)/M_PI
    else 
        !Every rank except rank=0 should send its `count' to rank=0
        call MPI_SEND(count,1,MPI_DOUBLE_PRECISION,0,0,MPI_COMM_WORLD,ierr)
    endif

    call MPI_Finalize(ierr)
end program mcpi
