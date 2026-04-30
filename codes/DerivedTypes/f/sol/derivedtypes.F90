program derivedtypes
    use iso_c_binding
    use mpi_f08
    implicit none
    type typ_atom
        sequence
        !character(len=1):: useless1
        real(kind=8):: xat
        !character(len=1):: useless2
        character(len=3):: sat
        character(len=5):: useless3
    end type typ_atom
    integer:: rank, ssize, ierr, lengths(2), size_notresized, size_resized
    type(typ_atom):: atom(2)
    type(MPI_Datatype):: types_list(2), typ_atom_mpi_notresized, typ_atom_mpi
    type(MPI_Status):: status
    integer(kind=MPI_ADDRESS_KIND):: displacements(2), tmp_lb, extent_notresized, extent_true_notresized
    integer(kind=MPI_ADDRESS_KIND):: extent_resized, extent_true_resized, address_1, address_2
    integer(kind=MPI_ADDRESS_KIND):: extent_typ_atom
    logical:: do_it_unsafe=.false.

    call MPI_INIT(ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,ssize,ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)
    if(ssize/=2) then
        write(*,'(a)') "This exercise is supposed to be run with 2 MPI processes!"
        call MPI_ABORT(MPI_COMM_WORLD,-1,ierr)
    endif
    lengths(1)=1
    lengths(2)=3
    call MPI_Get_address(atom(1)%xat,address_1)
    call MPI_Get_address(atom(1)%sat,address_2)
    displacements(1)=int(0,kind=MPI_ADDRESS_KIND)
    displacements(2)=MPI_AINT_DIFF(address_2,address_1)
    types_list(1)=MPI_DOUBLE_PRECISION
    types_list(2)=MPI_CHARACTER
    call MPI_TYPE_CREATE_STRUCT(2,lengths,displacements,types_list,typ_atom_mpi_notresized);
    if(do_it_unsafe) then
        typ_atom_mpi=typ_atom_mpi_notresized
    else
        call MPI_GET_ADDRESS(atom(1),address_1)
        call MPI_GET_ADDRESS(atom(2),address_2)
        extent_typ_atom=MPI_AINT_DIFF(address_2,address_1)
        call MPI_TYPE_CREATE_RESIZED(typ_atom_mpi_notresized, &
            int(0,kind=MPI_ADDRESS_KIND),extent_typ_atom,typ_atom_mpi)
    endif
    call MPI_TYPE_COMMIT(typ_atom_mpi)
    if(rank==0) then
        write(*,'(a)') "------------------------------------------------"
        call MPI_TYPE_SIZE(typ_atom_mpi_notresized,size_notresized)
        call MPI_TYPE_GET_EXTENT(typ_atom_mpi_notresized,tmp_lb,extent_notresized)
        call MPI_TYPE_GET_TRUE_EXTENT(typ_atom_mpi_notresized,tmp_lb,extent_true_notresized)
        call MPI_TYPE_SIZE(typ_atom_mpi,size_resized)
        call MPI_TYPE_GET_EXTENT(typ_atom_mpi,tmp_lb,extent_resized)
        call MPI_TYPE_GET_TRUE_EXTENT(typ_atom_mpi,tmp_lb,extent_true_resized)
        write(*,'(a,i)') "notresized: size_notresized= ",size_notresized
        write(*,'(a,i)') "notresized: extent_notresized= ",extent_notresized
        write(*,'(a,i)') "notresized: extent_true_notresized= ",extent_true_notresized
        write(*,'(a)') "------------------------------------------------"
        write(*,'(a,i)') "resized: size_resized= ",size_resized
        write(*,'(a,i)') "resized: extent_resized= ",extent_resized
        write(*,'(a,i)') "resized: extent_true_resized= ",extent_true_resized
        write(*,'(a)') "------------------------------------------------"
        call MPI_GET_ADDRESS(atom(1),address_1)
        call MPI_GET_ADDRESS(atom(2),address_2)
        write(*,'(a,i)') "size_from_diff= ",MPI_AINT_DIFF(address_2,address_1)
        write(*,'(a)') "------------------------------------------------"
    endif
    if(rank==0) then
        atom(1)%xat=10.d0
        atom(1)%sat="SiA"
        atom(2)%xat=20.d0
        atom(2)%sat="TiB"
        call MPI_Send(atom(1)%xat,2,typ_atom_mpi,1,11,MPI_COMM_WORLD)
    else
        call MPI_Recv(atom(1)%xat,2,typ_atom_mpi,0,11,MPI_COMM_WORLD,status)
        write(*,'(a,f)') "atom(1): xat= ",atom(1)%xat
        write(*,'(a,a)') "atom(1): sat= ",trim(atom(1)%sat)
        write(*,'(a,f)') "atom(2): xat= ",atom(2)%xat
        write(*,'(a,a)') "atom(2): sat= ",trim(atom(2)%sat)
    endif

    call MPI_FINALIZE(ierr)
end program derivedtypes
