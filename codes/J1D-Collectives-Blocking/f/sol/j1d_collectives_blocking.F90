program Jacobi
    use mpi_f08
    use mod_utils, only: swap_real8_pointers, get_share_per_proc, dump_domain
    implicit none
    integer:: i, n, istep, nstep, irank, nrank, ierr, is, ie
    real(kind=8):: h, errmax, tt, gnrm
    real(kind=8), allocatable, target:: pot(:), wa(:)
    real(kind=8), allocatable:: pot_global(:), errmax_all(:), gnrm_all(:)
    real(kind=8), pointer:: p1(:)=>null(), p2(:)=>null()
    integer:: irank_l, irank_r, n_loc_ext, n_tot_left, jrank, n_loc, ios
    integer, allocatable:: n_loc_all(:)
    logical:: restart
    character(len=256):: filename
    type(MPI_STATUS):: status
    call MPI_INIT(ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD,irank,ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,nrank,ierr)
    call get_input_param(restart,filename)
    if(restart) then
        if(irank==0) then
            open(unit=21,file=trim(filename),status='old',iostat=ios)
            read(21,*) n
        endif
        call MPI_BCAST(n,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    else
        n=10000
    endif
    nstep=10000000
    h=1.d0/real(n,kind=8)
    call get_share_per_proc(n,irank,nrank,is,ie)
    !write(*,*) irank,is,ie,n
    irank_l=irank-1
    irank_r=irank+1
    if(irank_l<0) irank_l=MPI_PROC_NULL
    if(irank_r>=nrank) irank_r=MPI_PROC_NULL
    allocate(pot(is-1:ie+1),wa(is-1:ie+1))
    n_loc=ie-is+1
    allocate(n_loc_all(nrank),source=0)
    call MPI_ALLGATHER(n_loc,1,MPI_INTEGER, &
        n_loc_all,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
    !n_loc_ext is the extended value of n_loc that it is the largest
    !value of n_loc across all ranks
    n_loc_ext=maxval(n_loc_all)
    !pot_global stores the potentials on all grid points
    allocate(pot_global(n_loc_ext*nrank))
    if(restart) then
        if(irank==0) then
            do jrank=0,nrank-1
                do i=1,n_loc_all(jrank+1)
                    read(21,*) tt,pot_global(jrank*n_loc_ext+i)
                enddo
            enddo
            close(21)
        endif
        !rank=0 has read the potential on all grid points from the
        !restart file and now it should scatter across all ranks
        call MPI_SCATTER(pot_global,n_loc_ext,MPI_DOUBLE_PRECISION, &
            pot(is),n_loc_ext,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        if(nrank>1) then
            call MPI_SEND(pot(is  ),1,MPI_DOUBLE_PRECISION,irank_l,11,MPI_COMM_WORLD,ierr)
            call MPI_RECV(pot(ie+1),1,MPI_DOUBLE_PRECISION,irank_r,11,MPI_COMM_WORLD,status,ierr)
            call MPI_SEND(pot(ie)  ,1,MPI_DOUBLE_PRECISION,irank_r,22,MPI_COMM_WORLD,ierr)
            call MPI_RECV(pot(is-1),1,MPI_DOUBLE_PRECISION,irank_l,22,MPI_COMM_WORLD,status,ierr)
        endif
    else
        do i=is-1,ie+1
            pot(i)=(i*h)**2
        enddo
    endif
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
    !rank=0 will write the entire potential array on all grid points
    !into a restart file and thus it should gather chunks from all ranks
    call MPI_GATHER(pot(is),n_loc_ext,MPI_DOUBLE_PRECISION, &
        pot_global,n_loc_ext,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    if(irank==0) then
        if(.not. restart) filename='starting_potential.dat'
        open(unit=21,file=trim(filename),status='replace',iostat=ios)
        write(21,'(i8)') n
        n_tot_left=0
        do jrank=0,nrank-1
            do i=1,n_loc_all(jrank+1)
                write(21,'(2es26.17)') (n_tot_left+i)*h,pot_global(jrank*n_loc_ext+i)
            enddo
        n_tot_left=n_tot_left+n_loc_all(jrank+1)
        enddo
        close(21)
    endif
    deallocate(pot_global)
    deallocate(n_loc_all)
    
    call dump_domain(irank,nrank,nstep,ie-is+1,pot(is),n-1)
    errmax=0.d0
    gnrm=0.d0
    do i=is,ie
        !write(31+irank,'(es14.5,es24.15)') i*h,pot(i)
        errmax=max(errmax,abs(pot(i)-i*h))
        gnrm=gnrm+(pot(i+1)-2.d0*pot(i)+pot(i-1))**2
    enddo
    allocate(errmax_all(nrank))
    allocate(gnrm_all(nrank))
    call MPI_GATHER(errmax,1,MPI_DOUBLE_PRECISION, &
        errmax_all,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call MPI_GATHER(gnrm,1,MPI_DOUBLE_PRECISION, &
        gnrm_all,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    gnrm=sqrt(sum(gnrm_all))
    if(irank==0) then
        write(*,'(a,2es14.5)') 'Convergence= ',maxval(errmax_all),gnrm
    endif
    deallocate(errmax_all,gnrm_all)
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
subroutine get_input_param(restart,filename)
    implicit none
    logical, intent(out):: restart
    character(len=*), intent(out):: filename
    integer:: narg, ierr
    character(len=256):: argv
    narg=command_argument_count()
    if(narg==0) then
        restart=.false.
    elseif(narg==1) then
        call get_command_argument(1,argv)
        if(irank==0) then
            if(trim(argv)=='-r') then
                write(*,'(a)') 'ERROR: missing restart filename in the argument list!'
            else
                write(*,'(a)') 'ERROR: unknown 1st argument!'
            endif
        endif
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        call MPI_ABORT(MPI_COMM_WORLD,-1,ierr)
    elseif(narg==2) then
        call get_command_argument(1,argv)
        if(trim(argv)=='-r') then
            call get_command_argument(2,filename)
            restart=.true.
        else
            if(irank==0) write(*,'(a)') 'ERROR: unknown 1st argument!'
            call MPI_ABORT(MPI_COMM_WORLD,-1,ierr)
        endif
    else
        if(irank==0) write(*,'(a)') 'ERROR: too many arguments!'
        call MPI_ABORT(MPI_COMM_WORLD,-1,ierr)
    endif
end subroutine get_input_param
end program Jacobi
