! Copyright (C) 2007, 2008 Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_gas_state module.

!> The gas_state_t structure and associated subroutines.
module pmc_gas_state

  use pmc_util
  use pmc_inout
  use pmc_gas_data
  use pmc_mpi
#ifdef PMC_USE_MPI
  use mpi
#endif

  !> Current state of the gas concentrations in the system.
  !!
  !! The gas species are defined by the gas_data_t structure, so that
  !! \c gas_state%%conc(i) is the current concentration of the gas
  !! with name \c gas_data%%name(i), etc.
  type gas_state_t
     !> Length n_spec, concentration (ppb).
     real*8, pointer :: conc(:)
  end type gas_state_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Allocate storage for gas species.
  subroutine gas_state_alloc(gas_state, n_spec)

    !> Gas state to be allocated.
    type(gas_state_t), intent(out) :: gas_state
    !> Number of species.
    integer, intent(in) :: n_spec

    allocate(gas_state%conc(n_spec))
    call gas_state_zero(gas_state)

  end subroutine gas_state_alloc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Free all storage.
  subroutine gas_state_free(gas_state)

    !> Gas state to be freed.
    type(gas_state_t), intent(inout) :: gas_state

    deallocate(gas_state%conc)

  end subroutine gas_state_free

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Zeros the state.
  subroutine gas_state_zero(gas_state)

    !> Gas state.
    type(gas_state_t), intent(inout) :: gas_state

    gas_state%conc = 0d0

  end subroutine gas_state_zero

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Copy to an already allocated to_state.
  subroutine gas_state_copy(from_state, to_state)

    !> Existing gas state.
    type(gas_state_t), intent(in) :: from_state
    !> Must be allocated already.
    type(gas_state_t), intent(out) :: to_state

    integer :: n_spec

    n_spec = size(from_state%conc)
    deallocate(to_state%conc)
    allocate(to_state%conc(n_spec))
    to_state%conc = from_state%conc

  end subroutine gas_state_copy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Scale a gas state.
  subroutine gas_state_scale(gas_state, alpha)

    !> Existing gas state.
    type(gas_state_t), intent(inout) :: gas_state
    !> Scale factor.
    real*8, intent(in) :: alpha

    gas_state%conc = gas_state%conc * alpha

  end subroutine gas_state_scale

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Adds the given gas_state_delta.
  subroutine gas_state_add(gas_state, gas_state_delta)

    !> Existing gas state.
    type(gas_state_t), intent(inout) :: gas_state
    !> Incremental state.
    type(gas_state_t), intent(in) :: gas_state_delta

    gas_state%conc = gas_state%conc + gas_state_delta%conc

  end subroutine gas_state_add

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Subtracts the given gas_state_delta.
  subroutine gas_state_sub(gas_state, gas_state_delta)

    !> Existing gas state.
    type(gas_state_t), intent(inout) :: gas_state
    !> Incremental state.
    type(gas_state_t), intent(in) :: gas_state_delta

    gas_state%conc = gas_state%conc - gas_state_delta%conc

  end subroutine gas_state_sub

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Sets y = a * x + y.
  subroutine gas_state_axpy(alpha, gas_state_x, gas_state_y)

    !> Coefficient.
    real*8, intent(in) :: alpha
    !> Incremental state.
    type(gas_state_t), intent(in) :: gas_state_x
    !> Existing gas state.
    type(gas_state_t), intent(inout) :: gas_state_y

    gas_state_y%conc = alpha * gas_state_x%conc + gas_state_y%conc

  end subroutine gas_state_axpy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Determine the current gas_state and rate by interpolating at the
  !> current time with the lists of gas_states and rates.
  subroutine gas_state_interp_1d(gas_state_list, time_list, &
         rate_list, time, gas_state, rate)

    !> Gas states.
    type(gas_state_t), intent(in) :: gas_state_list(:)
    !> Times (s).
    real*8, intent(in) :: time_list(size(gas_state_list))
    !> Rates (s^{-1}).
    real*8, intent(in) :: rate_list(size(gas_state_list))
    !> Current time (s).
    real*8, intent(in) :: time
    !> Current gas state.
    type(gas_state_t), intent(inout) :: gas_state
    !> Current rate (s^{-1}).
    real*8, intent(out) :: rate

    integer :: n, p
    real*8 :: y, alpha

    n = size(gas_state_list)
    p = find_1d(n, time_list, time)
    if (p == 0) then
       ! before the start, just use the first state and rate
       call gas_state_copy(gas_state_list(1), gas_state)
       rate = rate_list(1)
    elseif (p == n) then
       ! after the end, just use the last state and rate
       call gas_state_copy(gas_state_list(n), gas_state)
       rate = rate_list(n)
    else
       ! in the middle, use the previous state
       call gas_state_copy(gas_state_list(p), gas_state)
       rate = rate_list(p)
    end if

  end subroutine gas_state_interp_1d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Write full state.
  subroutine inout_write_gas_state(file, gas_state)
    
    !> File to write to.
    type(inout_file_t), intent(inout) :: file
    !> Gas_state to write.
    type(gas_state_t), intent(in) :: gas_state

    call inout_write_comment(file, "begin gas_state")
    call inout_write_real_array(file, "conc(ppb)", gas_state%conc)
    call inout_write_comment(file, "end gas_state")
    
  end subroutine inout_write_gas_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read full state.
  subroutine inout_read_gas_state(file, gas_state)
    
    !> File to read from.
    type(inout_file_t), intent(inout) :: file
    !> Gas_state to read.
    type(gas_state_t), intent(out) :: gas_state

    call inout_check_comment(file, "begin gas_state")
    call inout_read_real_array(file, "conc(ppb)", gas_state%conc)
    call inout_check_comment(file, "end gas_state")
    
  end subroutine inout_read_gas_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read gas state from the file named on the line read from file.
  subroutine spec_read_gas_state(file, gas_data, name, gas_state)

    !> Inout file.
    type(inout_file_t), intent(inout) :: file
    !> Gas data.
    type(gas_data_t), intent(in) :: gas_data
    !> Name of data line for filename.
    character(len=*), intent(in) :: name
    !> Gas data.
    type(gas_state_t), intent(out) :: gas_state

    character(len=MAX_VAR_LEN) :: read_name
    type(inout_file_t) :: read_file
    integer :: n_species, species, i
    character(len=MAX_VAR_LEN), pointer :: species_name(:)
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

  !> Read an array of gas states with associated times and rates from
  !> the file named on the line read from the given file.
  subroutine spec_read_gas_states_times_rates(file, gas_data, name, &
       times, rates, gas_states)

    !> Inout file.
    type(inout_file_t), intent(inout) :: file
    !> Gas data.
    type(gas_data_t), intent(in) :: gas_data
    !> Name of data line for filename.
    character(len=*), intent(in) :: name
    !> Times (s).
    real*8, pointer :: times(:)
    !> Rates (s^{-1}).
    real*8, pointer :: rates(:)
    !> Gas states.
    type(gas_state_t), pointer :: gas_states(:)

    character(len=MAX_VAR_LEN) :: read_name
    type(inout_file_t) :: read_file
    integer :: n_lines, species, i, n_time, i_time
    character(len=MAX_VAR_LEN), pointer :: species_name(:)
    real*8, pointer :: species_data(:,:)

    ! read the filename then read the data from that file
    call inout_read_string(file, name, read_name)
    call inout_open_read(read_name, read_file)
    call inout_read_real_named_array(read_file, 0, species_name, species_data)
    call inout_close(read_file)

    ! check the data size
    n_lines = size(species_data, 1)
    if (n_lines < 2) then
       write(0,*) 'ERROR: insufficient data lines in ', trim(read_name)
       call exit(1)
    end if
    if (trim(species_name(1)) /= 'time') then
       write(0,*) 'ERROR: row 1 in ', trim(read_name), &
            ' must start with: time'
       call exit(1)
    end if
    if (trim(species_name(2)) /= 'rate') then
       write(0,*) 'ERROR: row 2 in ', trim(read_name), &
            ' must start with: rate'
       call exit(1)
    end if
    n_time = size(species_data, 2)
    if (n_time < 1) then
       write(0,*) 'ERROR: each line in ', trim(read_name), &
            ' must contain at least one data value'
       call exit(1)
    end if

    ! copy over the data
    allocate(gas_states(n_time))
    allocate(times(n_time))
    allocate(rates(n_time))
    do i_time = 1,n_time
       call gas_state_alloc(gas_states(i_time), gas_data%n_spec)
       times(i_time) = species_data(1,i_time)
       rates(i_time) = species_data(2,i_time)
    end do
    do i = 3,n_lines
       species = gas_data_spec_by_name(gas_data, species_name(i))
       if (species == 0) then
          write(0,*) 'ERROR: unknown species ', trim(species_name(i)), &
               ' in file ', trim(read_name)
          call exit(1)
       end if
       do i_time = 1,n_time
          gas_states(i_time)%conc(species) = species_data(i,i_time)
       end do
    end do
    deallocate(species_name)
    deallocate(species_data)

  end subroutine spec_read_gas_states_times_rates

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Computes the average of an array of gas_state.
  subroutine gas_state_average(gas_state_vec, gas_state_avg)

    !> Array of gas_state.
    type(gas_state_t), intent(in) :: gas_state_vec(:)
    !> Average of gas_state_vec.
    type(gas_state_t), intent(out) :: gas_state_avg

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

  !> Average val over all processes.
  subroutine gas_state_mix(val)

    !> Value to average.
    type(gas_state_t), intent(inout) :: val

#ifdef PMC_USE_MPI
    type(gas_state_t) :: val_avg

    call gas_state_alloc(val_avg, size(val%conc))
    call pmc_mpi_allreduce_average_real_array(val%conc, val_avg%conc)
    val%conc = val_avg%conc
    call gas_state_free(val_avg)
#endif

  end subroutine gas_state_mix

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Determines the number of bytes required to pack the given value.
  integer function pmc_mpi_pack_size_gas_state(val)

    !> Value to pack.
    type(gas_state_t), intent(in) :: val

    pmc_mpi_pack_size_gas_state = &
         + pmc_mpi_pack_size_real_array(val%conc)

  end function pmc_mpi_pack_size_gas_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Packs the given value into the buffer, advancing position.
  subroutine pmc_mpi_pack_gas_state(buffer, position, val)

    !> Memory buffer.
    character, intent(inout) :: buffer(:)
    !> Current buffer position.
    integer, intent(inout) :: position
    !> Value to pack.
    type(gas_state_t), intent(in) :: val

#ifdef PMC_USE_MPI
    integer :: prev_position

    prev_position = position
    call pmc_mpi_pack_real_array(buffer, position, val%conc)
    call assert(655827004, &
         position - prev_position == pmc_mpi_pack_size_gas_state(val))
#endif

  end subroutine pmc_mpi_pack_gas_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Unpacks the given value from the buffer, advancing position.
  subroutine pmc_mpi_unpack_gas_state(buffer, position, val)

    !> Memory buffer.
    character, intent(inout) :: buffer(:)
    !> Current buffer position.
    integer, intent(inout) :: position
    !> Value to pack.
    type(gas_state_t), intent(out) :: val

#ifdef PMC_USE_MPI
    integer :: prev_position

    prev_position = position
    call pmc_mpi_unpack_real_array(buffer, position, val%conc)
    call assert(520815247, &
         position - prev_position == pmc_mpi_pack_size_gas_state(val))
#endif

  end subroutine pmc_mpi_unpack_gas_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Computes the average of val across all processes, storing the
  !> result in val_avg on the root process.
  subroutine pmc_mpi_reduce_avg_gas_state(val, val_avg)

    !> Value to average.
    type(gas_state_t), intent(in) :: val
    !> Result.
    type(gas_state_t), intent(out) :: val_avg

    call pmc_mpi_reduce_avg_real_array(val%conc, val_avg%conc)

  end subroutine pmc_mpi_reduce_avg_gas_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_gas_state
