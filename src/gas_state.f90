! Copyright (C) 2007 Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.
!
! Gas state.

module pmc_gas_state

  type gas_state_t
     real*8, pointer :: conc(:)          ! length n_spec, concentration (ppb)
  end type gas_state_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine gas_state_alloc(gas_state, n_spec)

    ! Allocate storage for gas species.

    type(gas_state_t), intent(out) :: gas_state ! gas state to be allocated
    integer, intent(in) :: n_spec       ! number of species

    allocate(gas_state%conc(n_spec))

  end subroutine gas_state_alloc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine gas_state_free(gas_state)

    ! Free all storage.

    type(gas_state_t), intent(inout) :: gas_state ! gas state to be freed

    deallocate(gas_state%conc)

  end subroutine gas_state_free

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine gas_state_copy(from_state, to_state)

    ! Copy to an already allocated to_state.

    type(gas_state_t), intent(in) :: from_state ! existing gas state
    type(gas_state_t), intent(out) :: to_state ! must be allocated already

    integer :: n_spec

    n_spec = size(from_state%conc)
    deallocate(to_state%conc)
    allocate(to_state%conc(n_spec))
    to_state%conc = from_state%conc

  end subroutine gas_state_copy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine gas_state_scale(gas_state, alpha)

    ! Scale a gas state.

    type(gas_state_t), intent(inout) :: gas_state ! existing gas state
    real*8, intent(in) :: alpha         ! scale factor

    gas_state%conc = gas_state%conc * alpha

  end subroutine gas_state_scale

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine gas_state_add(gas_state, gas_state_delta)

    ! Adds the given gas_state_delta.

    type(gas_state_t), intent(inout) :: gas_state ! existing gas state
    type(gas_state_t), intent(in) :: gas_state_delta ! incremental state

    gas_state%conc = gas_state%conc + gas_state_delta%conc

  end subroutine gas_state_add

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine gas_state_sub(gas_state, gas_state_delta)

    ! Subtracts the given gas_state_delta.

    type(gas_state_t), intent(inout) :: gas_state ! existing gas state
    type(gas_state_t), intent(in) :: gas_state_delta ! incremental state

    gas_state%conc = gas_state%conc - gas_state_delta%conc

  end subroutine gas_state_sub

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine inout_write_gas_state(file, gas_state)
    
    ! Write full state.
    
    use pmc_inout
    
    type(inout_file_t), intent(inout) :: file ! file to write to
    type(gas_state_t), intent(in) :: gas_state ! gas_state to write

    call inout_write_real_array(file, "conc(ppb)", gas_state%conc)
    
  end subroutine inout_write_gas_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine inout_read_gas_state(file, gas_state)
    
    ! Read full state.
    
    use pmc_inout
    
    type(inout_file_t), intent(inout) :: file ! file to read from
    type(gas_state_t), intent(out) :: gas_state ! gas_state to read

    call inout_read_real_array(file, "conc(ppb)", gas_state%conc)
    
  end subroutine inout_read_gas_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine spec_read_gas_state(file, gas_data, name, gas_state)

    ! Read gas state from the file named on the line read from file.

    use pmc_inout
    use pmc_gas_data

    type(inout_file_t), intent(inout) :: file ! inout file
    type(gas_data_t), intent(in) :: gas_data ! gas data
    character(len=*), intent(in) :: name ! name of data line for filename
    type(gas_state_t), intent(out) :: gas_state ! gas data

    character(len=MAX_CHAR_LEN) :: read_name
    type(inout_file_t) :: read_file
    integer :: n_species, species, i
    character(len=MAX_CHAR_LEN), pointer :: species_name(:)
    real*8, pointer :: species_data(:,:)

    ! read the filename then read the data from that file
    call inout_read_string(file, name, read_name)
    call inout_open_read(read_name, read_file)
    call inout_read_real_named_array(read_file, 0, species_name, species_data)
    call inout_close(read_file)

    ! check the data size
    n_species = size(species_data, 1)
    if (.not. ((size(species_data, 2) == 1) .or. (n_species == 0))) then
       write(0,*) 'ERROR: each line in ', trim(read_name), &
            ' should contain exactly one data value'
       call exit(1)
    end if

    ! copy over the data
    call gas_state_alloc(gas_state, gas_data%n_spec)
    gas_state%conc = 0d0
    do i = 1,n_species
       species = gas_data_spec_by_name(gas_data, species_name(i))
       if (species == 0) then
          write(0,*) 'ERROR: unknown species ', trim(species_name(i)), &
               ' in file ', trim(read_name)
          call exit(1)
       end if
       gas_state%conc(species) = species_data(i,1)
    end do
    deallocate(species_name)
    deallocate(species_data)

  end subroutine spec_read_gas_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine gas_state_average(gas_state_vec, gas_state_avg)
    
    ! Computes the average of an array of gas_state.

    use pmc_util

    type(gas_state_t), intent(in) :: gas_state_vec(:) ! array of gas_state
    type(gas_state_t), intent(out) :: gas_state_avg   ! average of gas_state_vec

    integer :: n_spec, i_spec, i, n

    n_spec = size(gas_state_vec(1)%conc)
    call gas_state_alloc(gas_state_avg, n_spec)
    n = size(gas_state_vec)
    do i_spec = 1,n_spec
       call average_real((/(gas_state_vec(i)%conc(i_spec),i=1,n)/), &
            gas_state_avg%conc(i_spec))
    end do
    
  end subroutine gas_state_average
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function pmc_mpi_pack_gas_state_size(val)

    ! Determines the number of bytes required to pack the given value.

    use pmc_mpi

    type(gas_state_t), intent(in) :: val ! value to pack

    pmc_mpi_pack_gas_state_size = &
         + pmc_mpi_pack_real_array_size(val%conc)

  end function pmc_mpi_pack_gas_state_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine pmc_mpi_pack_gas_state(buffer, position, val)

    ! Packs the given value into the buffer, advancing position.

#ifdef PMC_USE_MPI
    use mpi
    use pmc_mpi
    use pmc_util
#endif

    character, intent(inout) :: buffer(:) ! memory buffer
    integer, intent(inout) :: position  ! current buffer position
    type(gas_state_t), intent(in) :: val ! value to pack

#ifdef PMC_USE_MPI
    integer :: prev_position

    prev_position = position
    call pmc_mpi_pack_real_array(buffer, position, val%conc)
    call assert(position - prev_position == pmc_mpi_pack_gas_state_size(val))
#endif

  end subroutine pmc_mpi_pack_gas_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine pmc_mpi_unpack_gas_state(buffer, position, val)

    ! Unpacks the given value from the buffer, advancing position.

#ifdef PMC_USE_MPI
    use mpi
    use pmc_mpi
    use pmc_util
#endif

    character, intent(inout) :: buffer(:) ! memory buffer
    integer, intent(inout) :: position  ! current buffer position
    type(gas_state_t), intent(out) :: val ! value to pack

#ifdef PMC_USE_MPI
    integer :: prev_position

    prev_position = position
    call pmc_mpi_unpack_real_array(buffer, position, val%conc)
    call assert(position - prev_position == pmc_mpi_pack_gas_state_size(val))
#endif

  end subroutine pmc_mpi_unpack_gas_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine pmc_mpi_reduce_avg_gas_state(val, val_avg)

    ! Computes the average of val across all processes, storing the
    ! result in val_avg on the root process.

    use pmc_mpi

    type(gas_state_t), intent(in) :: val ! value to average
    type(gas_state_t), intent(out) :: val_avg ! result

    call pmc_mpi_reduce_avg_real_array(val%conc, val_avg%conc)

  end subroutine pmc_mpi_reduce_avg_gas_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_gas_state
