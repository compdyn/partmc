! Copyright (C) 2007-2011 Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_gas_state module.

!> The gas_state_t structure and associated subroutines.
module pmc_gas_state

  use pmc_util
  use pmc_spec_file
  use pmc_gas_data
  use pmc_mpi
  use pmc_netcdf
#ifdef PMC_USE_MPI
  use mpi
#endif

  !> Current state of the gas mixing ratios in the system.
  !!
  !! The gas species are defined by the gas_data_t structure, so that
  !! \c gas_state%%mix_rat(i) is the current mixing ratio of the gas
  !! with name \c gas_data%%name(i), etc.
  type gas_state_t
     !> Length n_spec, mixing ratio (ppb).
     real(kind=dp), pointer :: mix_rat(:)
  end type gas_state_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Allocate storage for gas species.
  subroutine gas_state_allocate(gas_state)

    !> Gas state to be allocated.
    type(gas_state_t), intent(out) :: gas_state

    allocate(gas_state%mix_rat(0))

  end subroutine gas_state_allocate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Allocate storage for gas species of the given size.
  subroutine gas_state_allocate_size(gas_state, n_spec)

    !> Gas state to be allocated.
    type(gas_state_t), intent(out) :: gas_state
    !> Number of species.
    integer, intent(in) :: n_spec

    allocate(gas_state%mix_rat(n_spec))
    call gas_state_zero(gas_state)

  end subroutine gas_state_allocate_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Free all storage.
  subroutine gas_state_deallocate(gas_state)

    !> Gas state to be freed.
    type(gas_state_t), intent(inout) :: gas_state

    deallocate(gas_state%mix_rat)

  end subroutine gas_state_deallocate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Zeros the state.
  subroutine gas_state_zero(gas_state)

    !> Gas state.
    type(gas_state_t), intent(inout) :: gas_state

    gas_state%mix_rat = 0d0

  end subroutine gas_state_zero

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Copy to an already allocated to_state.
  subroutine gas_state_copy(from_state, to_state)

    !> Existing gas state.
    type(gas_state_t), intent(in) :: from_state
    !> Must be allocated already.
    type(gas_state_t), intent(inout) :: to_state

    integer :: n_spec

    n_spec = size(from_state%mix_rat)
    deallocate(to_state%mix_rat)
    allocate(to_state%mix_rat(n_spec))
    to_state%mix_rat = from_state%mix_rat

  end subroutine gas_state_copy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Scale a gas state.
  subroutine gas_state_scale(gas_state, alpha)

    !> Existing gas state.
    type(gas_state_t), intent(inout) :: gas_state
    !> Scale factor.
    real(kind=dp), intent(in) :: alpha

    gas_state%mix_rat = gas_state%mix_rat * alpha

  end subroutine gas_state_scale

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Adds the given gas_state_delta.
  subroutine gas_state_add(gas_state, gas_state_delta)

    !> Existing gas state.
    type(gas_state_t), intent(inout) :: gas_state
    !> Incremental state.
    type(gas_state_t), intent(in) :: gas_state_delta

    gas_state%mix_rat = gas_state%mix_rat + gas_state_delta%mix_rat

  end subroutine gas_state_add

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Subtracts the given gas_state_delta.
  subroutine gas_state_sub(gas_state, gas_state_delta)

    !> Existing gas state.
    type(gas_state_t), intent(inout) :: gas_state
    !> Incremental state.
    type(gas_state_t), intent(in) :: gas_state_delta

    gas_state%mix_rat = gas_state%mix_rat - gas_state_delta%mix_rat

  end subroutine gas_state_sub

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Sets y = a * x + y.
  subroutine gas_state_axpy(alpha, gas_state_x, gas_state_y)

    !> Coefficient.
    real(kind=dp), intent(in) :: alpha
    !> Incremental state.
    type(gas_state_t), intent(in) :: gas_state_x
    !> Existing gas state.
    type(gas_state_t), intent(inout) :: gas_state_y

    gas_state_y%mix_rat = alpha * gas_state_x%mix_rat + gas_state_y%mix_rat

  end subroutine gas_state_axpy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Determine the current gas_state and rate by interpolating at the
  !> current time with the lists of gas_states and rates.
  subroutine gas_state_interp_1d(gas_state_list, time_list, &
         rate_list, time, gas_state, rate)

    !> Gas states.
    type(gas_state_t), intent(in) :: gas_state_list(:)
    !> Times (s).
    real(kind=dp), intent(in) :: time_list(size(gas_state_list))
    !> Rates (s^{-1}).
    real(kind=dp), intent(in) :: rate_list(size(gas_state_list))
    !> Current time (s).
    real(kind=dp), intent(in) :: time
    !> Current gas state.
    type(gas_state_t), intent(inout) :: gas_state
    !> Current rate (s^{-1}).
    real(kind=dp), intent(out) :: rate

    integer :: n, p
    real(kind=dp) :: y, alpha

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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read gas state from the file named on the line read from file.
  subroutine spec_file_read_gas_state(file, gas_data, gas_state)

    !> File to read gas state from.
    type(spec_file_t), intent(inout) :: file
    !> Gas data.
    type(gas_data_t), intent(in) :: gas_data
    !> Gas data to read.
    type(gas_state_t), intent(inout) :: gas_state

    integer :: n_species, species, i
    character(len=SPEC_LINE_MAX_VAR_LEN), pointer :: species_name(:)
    real(kind=dp), pointer :: species_data(:,:)

    !> \page input_format_gas_state Input File Format: Gas State
    !!
    !! A gas state input file must consist of one line per gas
    !! species, with each line having the species name followed by the
    !! species mixing ratio in ppb (parts per billion). The valid
    !! species names are those specfied by the \ref
    !! input_format_gas_data file, but not all species have to be
    !! listed. Any missing species will have mixing ratios of
    !! zero. For example, a gas state file could contain:
    !! <pre>
    !! # gas  mixing ratio (ppb)
    !! H2SO4  0
    !! HNO3   1
    !! HCl    0.7
    !! NH3    0.5
    !! </pre>
    !!
    !! See also:
    !!   - \ref spec_file_format --- the input file text format
    !!   - \ref output_format_gas_state --- the corresponding output format
    !!   - \ref input_format_gas_data --- the gas species list and
    !!     material data

    allocate(species_name(0))
    allocate(species_data(0,0))
    call spec_file_read_real_named_array(file, 0, species_name, &
         species_data)

    ! check the data size
    n_species = size(species_data, 1)
    if (.not. ((size(species_data, 2) == 1) .or. (n_species == 0))) then
       call die_msg(686719840, 'each line in ' // trim(file%name) &
            // ' must contain exactly one data value')
    end if

    ! copy over the data
    call gas_state_deallocate(gas_state)
    call gas_state_allocate_size(gas_state, gas_data%n_spec)
    gas_state%mix_rat = 0d0
    do i = 1,n_species
       species = gas_data_spec_by_name(gas_data, species_name(i))
       if (species == 0) then
          call die_msg(129794076, 'unknown species ' // &
               trim(species_name(i)) // ' in file ' // trim(file%name))
       end if
       gas_state%mix_rat(species) = species_data(i,1)
    end do
    deallocate(species_name)
    deallocate(species_data)

  end subroutine spec_file_read_gas_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read an array of gas states with associated times and rates from
  !> the file named on the line read from the given file.
  subroutine spec_file_read_gas_states_times_rates(file, gas_data, &
       times, rates, gas_states)

    !> Spec file.
    type(spec_file_t), intent(inout) :: file
    !> Gas data.
    type(gas_data_t), intent(in) :: gas_data
    !> Times (s).
    real(kind=dp), pointer :: times(:)
    !> Rates (s^{-1}).
    real(kind=dp), pointer :: rates(:)
    !> Gas states.
    type(gas_state_t), pointer :: gas_states(:)

    integer :: n_lines, species, i, n_time, i_time
    character(len=SPEC_LINE_MAX_VAR_LEN), pointer :: species_name(:)
    real(kind=dp), pointer :: species_data(:,:)

    !> \page input_format_gas_profile Input File Format: Gas Profile
    !!
    !! A gas profile input file must consist of three or more
    !! lines, consisting of:
    !! - the first line must begin with \c time and should be followed
    !!   by \f$N\f$ space-separated real scalars, giving the times (in
    !!   s after the start of the simulation) of the gas set
    !!   points --- the times must be in increasing order
    !! - the second line must begin with \c rate and should be
    !!   followed by \f$N\f$ space-separated real scalars, giving the
    !!   scalings (dimensionless) at the corresponding times
    !! - the third and subsequent lines specify gas species, one
    !!   species per line, with each line beginning with the species
    !!   name and followed by \f$N\f$ space-separated scalars giving
    !!   the gas state of that species at the corresponding times
    !!
    !! The units of the species lines depends on the type of gas
    !! profile:
    !! - emissions gas profiles have units of mol/(m^2 s) --- the
    !!   emission rate is divided by the current mixing layer
    !!   height to give a per-volume emission rate
    !! - background gas profiles have units of ppb (parts per billion)
    !!
    !! The species names must be those specified by the \ref
    !! input_format_gas_data. Any species not listed are taken to be
    !! zero.
    !!
    !! The scaling parameter is used to multiply all the gas state
    !! values at the corresponding time, giving a simple way of
    !! scaling the overall gas state.
    !!
    !! Between the specified times the gas profile is interpolated
    !! step-wise and kept constant at its last value. That is, if the
    !! times are \f$t_i\f$, the rates are \f$r_i\f$, and the gas
    !! states are \f$g_i\f$ (all with \f$i = 1,\ldots,n\f$), then
    !! between times \f$t_i\f$ and \f$t_{i+1}\f$ the gas state is
    !! constant at \f$r_i g_i\f$. Before time \f$t_1\f$ the gas state
    !! is \f$r_1 g_1\f$, while after time \f$t_n\f$ it is \f$r_n
    !! g_n\f$.
    !!
    !! Example: an emissions gas profile could be:
    !! <pre>
    !! time   0       600     1800    # time (in s) after simulation start
    !! rate   1       0.5     1       # scaling factor
    !! H2SO4  0       0       0       # emission rate in mol/(m^2 s)
    !! SO2    4e-9    5.6e-9  5e-9    # emission rate in mol/(m^2 s)
    !! </pre>
    !! Here there are no emissions of \f$\rm H_2SO_4\f$, while \f$\rm
    !! SO_2\f$ starts out being emitted at \f$4\times
    !! 10^{-9}\rm\ mol\ m^{-2}\ s^{-1}\f$ at the start of the simulation,
    !! before falling to a rate of \f$2.8\times 
    !! 10^{-9}\rm\ mol\ m^{-2}\ s^{-1}\f$ at 10&nbsp;min (note the scaling
    !! of 0.5), and then rising again to \f$5\times
    !! 10^{-9}\rm\ mol\ m^{-2}\ s^{-1}\f$ after 30&nbsp;min. Between
    !! 0&nbsp;min and 10&nbsp;min the emissions are the same as at
    !! 0&nbsp;min, while between 10&nbsp;min and 30&nbsp;min they are
    !! the same as at 10&nbsp;min. After 30&nbsp;min they are held
    !! constant at their final value.
    !!
    !! See also:
    !!   - \ref spec_file_format --- the input file text format
    !!   - \ref input_format_gas_data --- the gas species list and
    !!     material data

    ! read the data from the file
    allocate(species_name(0))
    allocate(species_data(0,0))
    call spec_file_read_real_named_array(file, 0, species_name, &
         species_data)

    ! check the data size
    n_lines = size(species_data, 1)
    if (n_lines < 2) then
       call die_msg(291542946, 'insufficient data lines in file ' &
            // trim(file%name))
    end if
    if (trim(species_name(1)) /= 'time') then
       call die_msg(398532628, 'row 1 in file ' &
            // trim(file%name) // ' must start with: time')
    end if
    if (trim(species_name(2)) /= 'rate') then
       call die_msg(398532628, 'row 2 in file ' &
            // trim(file%name) // ' must start with: rate')
    end if
    n_time = size(species_data, 2)
    if (n_time < 1) then
       call die_msg(398532628, 'each line in file ' &
            // trim(file%name) // ' must contain at least one data value')
    end if

    ! copy over the data
    do i_time = 1,size(gas_states)
       call gas_state_deallocate(gas_states(i_time))
    end do
    deallocate(gas_states)
    deallocate(times)
    deallocate(rates)
    allocate(gas_states(n_time))
    allocate(times(n_time))
    allocate(rates(n_time))
    do i_time = 1,n_time
       call gas_state_allocate_size(gas_states(i_time), gas_data%n_spec)
       times(i_time) = species_data(1,i_time)
       rates(i_time) = species_data(2,i_time)
    end do
    do i = 3,n_lines
       species = gas_data_spec_by_name(gas_data, species_name(i))
       if (species == 0) then
          call die_msg(806500079, 'unknown species ' &
               // trim(species_name(i)) // ' in file ' &
               // trim(file%name))
       end if
       do i_time = 1,n_time
          gas_states(i_time)%mix_rat(species) = species_data(i,i_time)
       end do
    end do
    deallocate(species_name)
    deallocate(species_data)

  end subroutine spec_file_read_gas_states_times_rates

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Average val over all processes.
  subroutine gas_state_mix(val)

    !> Value to average.
    type(gas_state_t), intent(inout) :: val

#ifdef PMC_USE_MPI
    type(gas_state_t) :: val_avg

    call gas_state_allocate_size(val_avg, size(val%mix_rat))
    call pmc_mpi_allreduce_average_real_array(val%mix_rat, val_avg%mix_rat)
    val%mix_rat = val_avg%mix_rat
    call gas_state_deallocate(val_avg)
#endif

  end subroutine gas_state_mix

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Average val over all processes, with the result only on the root
  !> processor.
  subroutine gas_state_reduce_avg(val)

    !> Value to average.
    type(gas_state_t), intent(inout) :: val

#ifdef PMC_USE_MPI
    type(gas_state_t) :: val_avg

    call gas_state_allocate_size(val_avg, size(val%mix_rat))
    call pmc_mpi_reduce_avg_real_array(val%mix_rat, val_avg%mix_rat)
    if (pmc_mpi_rank() == 0) then
       val%mix_rat = val_avg%mix_rat
    end if
    call gas_state_deallocate(val_avg)
#endif

  end subroutine gas_state_reduce_avg

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Determines the number of bytes required to pack the given value.
  integer function pmc_mpi_pack_size_gas_state(val)

    !> Value to pack.
    type(gas_state_t), intent(in) :: val

    pmc_mpi_pack_size_gas_state = &
         + pmc_mpi_pack_size_real_array(val%mix_rat)

  end function pmc_mpi_pack_size_gas_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
    call pmc_mpi_pack_real_array(buffer, position, val%mix_rat)
    call assert(655827004, &
         position - prev_position <= pmc_mpi_pack_size_gas_state(val))
#endif

  end subroutine pmc_mpi_pack_gas_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Unpacks the given value from the buffer, advancing position.
  subroutine pmc_mpi_unpack_gas_state(buffer, position, val)

    !> Memory buffer.
    character, intent(inout) :: buffer(:)
    !> Current buffer position.
    integer, intent(inout) :: position
    !> Value to pack.
    type(gas_state_t), intent(inout) :: val

#ifdef PMC_USE_MPI
    integer :: prev_position

    prev_position = position
    call pmc_mpi_unpack_real_array(buffer, position, val%mix_rat)
    call assert(520815247, &
         position - prev_position <= pmc_mpi_pack_size_gas_state(val))
#endif

  end subroutine pmc_mpi_unpack_gas_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Computes the average of val across all processes, storing the
  !> result in val_avg on the root process.
  subroutine pmc_mpi_reduce_avg_gas_state(val, val_avg)

    !> Value to average.
    type(gas_state_t), intent(in) :: val
    !> Result.
    type(gas_state_t), intent(inout) :: val_avg

    call pmc_mpi_reduce_avg_real_array(val%mix_rat, val_avg%mix_rat)

  end subroutine pmc_mpi_reduce_avg_gas_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Write full state.
  subroutine gas_state_output_netcdf(gas_state, ncid, gas_data)
    
    !> Gas state to write.
    type(gas_state_t), intent(in) :: gas_state
    !> NetCDF file ID, in data mode.
    integer, intent(in) :: ncid
    !> Gas data.
    type(gas_data_t), intent(in) :: gas_data

    integer :: dimid_gas_species

    !> \page output_format_gas_state Output File Format: Gas State
    !!
    !! The gas state uses the \c gas_species NetCDF dimension as specified
    !! in the \ref output_format_gas_data section.
    !!
    !! The gas state NetCDF variables are:
    !!   - \b gas_mixing_ratio (unit ppb, dim \c gas_species): current mixing
    !!     ratios of each gas species
    !!
    !! See also:
    !!   - \ref output_format_gas_data --- the gas species list and material
    !!     data output format
    !!   - \ref input_format_gas_state --- the corresponding input format

    call gas_data_netcdf_dim_gas_species(gas_data, ncid, &
         dimid_gas_species)

    call pmc_nc_write_real_1d(ncid, gas_state%mix_rat, &
         "gas_mixing_ratio", (/ dimid_gas_species /), unit="ppb", &
         long_name="mixing ratios of gas species")

  end subroutine gas_state_output_netcdf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read full state.
  subroutine gas_state_input_netcdf(gas_state, ncid, gas_data)
    
    !> Gas state to read.
    type(gas_state_t), intent(inout) :: gas_state
    !> NetCDF file ID, in data mode.
    integer, intent(in) :: ncid
    !> Gas data.
    type(gas_data_t), intent(in) :: gas_data

    call gas_state_deallocate(gas_state)
    call gas_state_allocate_size(gas_state, gas_data%n_spec)
    call pmc_nc_read_real_1d(ncid, gas_state%mix_rat, &
         "gas_mixing_ratio")

  end subroutine gas_state_input_netcdf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_gas_state
