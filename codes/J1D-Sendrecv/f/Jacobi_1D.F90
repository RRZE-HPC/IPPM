program Jacobi
    use mpi_f08
    use mod_utils, only: swap_real8_pointers, get_share_per_proc, dump_domain
    implicit none
    integer:: i, n, istep, nstep, irank, nrank, ierr, is, ie
    real(kind=8):: h, errmax
    real(kind=8), allocatable, target:: pot(:), wa(:)
    real(kind=8), pointer:: p1(:)=>null(), p2(:)=>null()
    !real(kind=8):: d=0.d0;
    irank=0 !HERE
    nrank=1 !HERE
    nstep=100000
    n=100000
    h=1.d0/real(n,kind=8)
    call get_share_per_proc(n,irank,nrank,is,ie)
    allocate(pot(is-1:ie+1),wa(is-1:ie+1))
    do i=is-1,ie+1
        pot(i)=(i*h)**2
    enddo
    pot(0)=0.d0 !HERE
    pot(n)=1.d0 !HERE
    wa(0)=0.d0 !HERE
    wa(n)=1.d0 !HERE
    p1=>pot
    p2=>wa
    do istep=0,nstep
        call sub(is,ie,p1,p2)
        !HERE
        call swap_real8_pointers(p1,p2)
    enddo
    if(associated(p2,wa)) pot=wa
    call dump_domain(irank,nrank,nstep,ie-is+1,pot(is),n-1)
    errmax=0.d0
    do i=is,ie
        !write(21,'(es14.5,es24.15)') i*h,pot(i)
        errmax=max(errmax,abs(pot(i)-i*h))
    enddo
    write(*,'(a,es14.5)') 'Maximum difference relative to analytical solution= ',errmax
    nullify(p1,p2)
    deallocate(pot,wa)
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
end program Jacobi
