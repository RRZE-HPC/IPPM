program J1D_PtP_Blocking
    use mpi_f08
    use mod_utils, only: swap_real8_pointers, get_share_per_proc, dump_domain
    implicit none
    integer:: i, n, istep, nstep, irank, nrank, ierr, is, ie
    real(kind=8):: h, errmax
    real(kind=8), allocatable, target:: pot(:), wa(:)
    real(kind=8), pointer:: p1(:)=>null(), p2(:)=>null()
    integer:: irank_l, irank_r
    type(MPI_STATUS):: status
    call MPI_INIT(ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD,irank,ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,nrank,ierr)
    nstep=100000
    n=100000
    h=1.d0/real(n,kind=8)
    call get_share_per_proc(n,irank,nrank,is,ie)
    !rank of the process working on the left chunk
    irank_l=irank-1
    !rank of the process working on the right chunk
    irank_r=irank+1
    if(irank_l<0) irank_l=MPI_PROC_NULL
    if(irank_r>=nrank) irank_r=MPI_PROC_NULL
    ! pot(is-1)    pot(is)      pot(is+1) ....  pot(ie-1)    pot(ie)      pot(ie+1)
    !  left-BC  1st-interior      2nd-    ....  (ie-is)th-  (ie-is+1)th-  right-BC
    ! So each chunk contains ie-is+1 relaxing elements
    ! in the case of far left and far right chunks, pot(is-1) and pot(ie+1) are
    ! left-BC and right-BC, respectively, and in other cases they are points
    ! from a neighboring process
    allocate(pot(is-1:ie+1),wa(is-1:ie+1))
    do i=is-1,ie+1
        pot(i)=(i*h)**2
    enddo
    if(irank==0) pot(0)=0.d0
    if(irank==nrank-1) pot(n)=1.d0
    if(irank==0) wa(0)=0.d0
    if(irank==nrank-1) wa(n)=1.d0
    p1=>pot
    p2=>wa
    do istep=0,nstep
        call sub(is,ie,p1,p2)
        if(nrank>1) then
            call MPI_SEND(p2(is  ),1,MPI_DOUBLE_PRECISION,irank_l,11,MPI_COMM_WORLD,ierr)
            call MPI_RECV(p2(ie+1),1,MPI_DOUBLE_PRECISION,irank_r,11,MPI_COMM_WORLD,status,ierr)
            call MPI_SEND(p2(ie)  ,1,MPI_DOUBLE_PRECISION,irank_r,22,MPI_COMM_WORLD,ierr)
            call MPI_RECV(p2(is-1),1,MPI_DOUBLE_PRECISION,irank_l,22,MPI_COMM_WORLD,status,ierr)
        endif
        call swap_real8_pointers(p1,p2)
    enddo
    if(associated(p2,wa)) pot=wa
    call dump_domain(irank,nrank,nstep,ie-is+1,pot(is),n-1)
    errmax=0.d0
    do i=is,ie
        !write(31+irank,'(es14.5,es24.15)') i*h,pot(i)
        errmax=max(errmax,abs(pot(i)-i*h))
    enddo
    write(*,'(a,es14.5,i6)') 'Maximum difference relative to analytical solution= ',errmax,irank
    nullify(p1,p2)
    deallocate(pot,wa)
    call MPI_FINALIZE(ierr)
contains
pure subroutine sub(is,ie,p1,p2)
    implicit none
    integer, intent(in):: is, ie
    real(kind=8), intent(inout):: p1(is-1:ie+1), p2(is-1:ie+1)
    integer:: i
    do i=is,ie
        p2(i)=0.5d0*(p1(i+1)+p1(i-1))
    enddo
end subroutine sub
end program J1D_PtP_Blocking
