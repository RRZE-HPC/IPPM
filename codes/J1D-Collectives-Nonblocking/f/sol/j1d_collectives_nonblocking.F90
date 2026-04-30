program J1D_Collective_Nonblocking
    use mpi_f08
    use mod_utils, only: swap_real8_pointers, get_share_per_proc, dump_domain
    implicit none
    integer:: i, n, istep, nstep, irank, nrank, ierr, is, ie, istep_print
    real(kind=8):: h, errmax, gmax, gmax_global
    real(kind=8), allocatable, target:: pot(:), wa(:)
    real(kind=8), pointer:: p1(:)=>null(), p2(:)=>null()
    integer:: irank_l, irank_r
    logical:: flag
    type(MPI_STATUS):: statuses(4), status_conv
    TYPE(MPI_REQUEST):: req(4), req_conv
    call MPI_INIT(ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD,irank,ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,nrank,ierr)
    nstep=10**7
    n=1000
    h=1.d0/real(n,kind=8)
    call get_share_per_proc(n,irank,nrank,is,ie)
    !rank of the process working on the left chunk
    irank_l=irank-1
    !rank of the process working on the right chunk
    irank_r=irank+1
    if(irank_l<0) irank_l=MPI_PROC_NULL
    if(irank_r>=nrank) irank_r=MPI_PROC_NULL
    allocate(pot(is-1:ie+1),wa(is-1:ie+1))
    do i=is-1,ie+1
        pot(i)=(i*h)**2
    enddo
    if(irank==0) pot(0)=0.d0
    if(irank==nrank-1) pot(n)=1.d0
    if(irank==0) wa(0)=0.d0
    if(irank==nrank-1) wa(n)=1.d0
    istep_print=0
    flag=.true.
    p1=>pot
    p2=>wa
    do istep=0,nstep
        !first computation of start and end elements of the chunk on this rank
        p2(is)=0.5d0*(p1(is+1)+p1(is-1))
        p2(ie)=0.5d0*(p1(ie+1)+p1(ie-1))
        if(nrank>1) then
            call MPI_ISEND(p2(is  ),1,MPI_DOUBLE_PRECISION,irank_l,11,MPI_COMM_WORLD,req(1),ierr)
            call MPI_IRECV(p2(ie+1),1,MPI_DOUBLE_PRECISION,irank_r,11,MPI_COMM_WORLD,req(2),ierr)
            call MPI_ISEND(p2(ie  ),1,MPI_DOUBLE_PRECISION,irank_r,22,MPI_COMM_WORLD,req(3),ierr)
            call MPI_IRECV(p2(is-1),1,MPI_DOUBLE_PRECISION,irank_l,22,MPI_COMM_WORLD,req(4),ierr)
        endif
        call sub_without_endpoints(is,ie,p1,p2)
        if(nrank>1) then
            call MPI_WAITALL(4,req,statuses,ierr)
        endif
        if(mod(istep,10)==0) then
            if(flag) then
                gmax=0.d0
                do i=is,ie
                    gmax=max(gmax,abs(pot(i+1)-2.d0*pot(i)+pot(i-1)))
                enddo
                call MPI_IALLREDUCE(gmax,gmax_global,1,MPI_DOUBLE_PRECISION, &
                    MPI_MAX,MPI_COMM_WORLD,req_conv,ierr)
                flag=.false.
            else
                call MPI_WAIT(req_conv,status_conv,ierr)
                gmax_global=gmax_global/h**2
                flag=.true.
            endif
            if(irank==0 .and. flag .and. (istep-istep_print>10**5 .or. gmax_global<1.d-5)) then
                write(*,'(a,es14.5,i10,5x,l)') 'Convergence= ',gmax_global,istep,flag
                istep_print=istep
            endif
            if(flag .and. gmax_global<1.d-5) exit
        endif
        call swap_real8_pointers(p1,p2)
    enddo
    if(.not. flag) then
        call MPI_WAIT(req_conv,status_conv,ierr)
        gmax_global=gmax_global/h**2
    endif
    if(associated(p2,wa)) pot=wa
    call dump_domain(irank,nrank,nstep,ie-is+1,pot(is),n-1)
    errmax=0.d0
    do i=is,ie
        !write(31+irank,'(es14.5,es24.15)') i*h,pot(i)
        errmax=max(errmax,abs(pot(i)-i*h))
    enddo
    write(*,'(a,es14.5,i6,es14.5)') 'Maximum difference relative to analytical solution= ',errmax,irank,gmax_global
    nullify(p1,p2)
    deallocate(pot,wa)
    call MPI_FINALIZE(ierr)
contains
pure subroutine sub_without_endpoints(is,ie,p1,p2)
    implicit none
    integer, intent(in):: is, ie
    real(kind=8), intent(inout):: p1(is-1:ie+1), p2(is-1:ie+1)
    integer:: i
    do i=is+1,ie-1
        p2(i)=0.5d0*(p1(i+1)+p1(i-1))
    enddo
end subroutine sub_without_endpoints
end program J1D_Collective_Nonblocking
