module mod_utils
    implicit none
    private
    public:: swap_real8_pointers, get_share_per_proc, dump_domain
contains
pure subroutine swap_real8_pointers(p1,p2)
    implicit none
    real(kind=8), pointer, intent(inout) :: p1(:), p2(:)
    real(kind=8), pointer :: p(:)
    if(.not. associated(p1,p2)) then
        p=>p1
        p1=>p2
        p2=>p
    endif
end subroutine swap_real8_pointers
subroutine get_share_per_proc(n,irank,nrank,is,ie)
    implicit none
    integer, intent(in):: n, irank, nrank
    integer, intent(out):: is, ie
    integer:: mq, mr
    mq=(n-1)/nrank
    mr=mod(n-1,nrank)
    is=irank*mq+1
    if(irank<mr) mq=mq+1
    is=is+min(irank,mr)
    ie=is+mq-1
end subroutine get_share_per_proc
subroutine dump_domain(irank,nrank,nstep,m,pot,mt)
    use mpi_f08
    use, intrinsic :: iso_c_binding, only: c_char
    implicit none
    integer, intent(in) :: irank, nrank, nstep, m, mt
    real(kind=8), intent(in):: pot(m)
    integer:: ios, iunit, ierr
    integer:: i, jrank
    character(len=1024):: filename
    character(kind=c_char), allocatable:: gray(:)
    write(filename,'(a,i8.8,a)') 'domain-',nstep,'.pgm'
    if(irank==0) then
        open(newunit=iunit,file=trim(filename),status="replace",action='write', &
            access='stream',form='formatted',iostat=ios)
        if(ios/=0) then
            write(*,'(3a,i)') 'ERROR: could not open file ',trim(filename),', error= ',ios
#if !defined(NOMPI)
            call MPI_ABORT(MPI_COMM_WORLD,1,ierr)
#endif
            call exit(1)
        endif
        write(iunit,'(a,/,i0,/,a)') "P5",mt,"255"
        close(iunit)
    endif
#if !defined(NOMPI)
    if(nrank>1) call MPI_BARRIER(MPI_COMM_WORLD,ierr)
#endif
    do jrank=0,nrank-1
        if(jrank==irank) then
            open(newunit=iunit,file=trim(filename),status="old",position='append', &
                access='stream',action='write',iostat=ios)
            if(ios/=0) then
                write(*,'(3a,i)') 'ERROR: could not open file ',trim(filename),', error= ',ios
#if !defined(NOMPI)
                call MPI_ABORT(MPI_COMM_WORLD,1,ierr)
#endif
                call exit(1)
            endif
            allocate(gray(m))
            do i=1,m
                gray(i)=char(int(255.d0*pot(i)))
            enddo
            write(iunit) gray
            deallocate(gray)
            close(iunit)
        endif
#if !defined(NOMPI)
        if(nrank>1) call MPI_BARRIER(MPI_COMM_WORLD,ierr)
#endif
    enddo
end subroutine dump_domain
end module mod_utils
