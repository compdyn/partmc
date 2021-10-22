! Copyright (C) 2007-2012 Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_gas_state module.

!> The gas_state_t structure and associated subroutines.
module pmc_gas_state

  use pmc_util
  use pmc_spec_file
  use pmc_gas_data
  use pmc_env_state
  use pmc_mpi
  use pmc_netcdf
#ifdef PMC_USE_CAMP
  use camp_camp_state
#endif
#ifdef PMC_USE_MPI
  use mpi
#endif

  !> Current state of the gas mixing ratios in the system.
  !!
  !! The gas species are defined by the gas_data_t structure, so that
  !! \c gas_state%%mix_rat(i) is the current mixing ratio of the gas
  !! with name \c gas_data%%name(i), etc.
  !!
  !! By convention, if gas_state_is_allocated() return \c .false.,
  !! then the gas_state is treated as zero for all operations on
  !! it. This will be the case for new \c gas_state_t structures.
  type gas_state_t
     !> Length n_spec, mixing ratio (ppb).
     real(kind=dp), allocatable :: mix_rat(:)
  end type gas_state_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Determine whether the \c gas_state is correctly allocated.
  logical function gas_state_is_allocated(gas_state)

    !> Gas state to check.
    type(gas_state_t), intent(in) :: gas_state

    gas_state_is_allocated = allocated(gas_state%mix_rat)

  end function gas_state_is_allocated

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Sets the sizes of the gas state.
  subroutine gas_state_set_size(gas_state, n_spec)

    !> Gas state to be allocated.
    type(gas_state_t), intent(inout) :: gas_state
    !> Number of species.
    integer, intent(in) :: n_spec

    if (allocated(gas_state%mix_rat)) deallocate(gas_state%mix_rat)
    allocate(gas_state%mix_rat(n_spec))
    call gas_state_zero(gas_state)

  end subroutine gas_state_set_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Zeros the state.
  subroutine gas_state_zero(gas_state)

    !> Gas state.
    type(gas_state_t), intent(inout) :: gas_state

    if (gas_state_is_allocated(gas_state)) then
       gas_state%mix_rat = 0d0
    end if

  end subroutine gas_state_zero

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Scale a gas state.
  subroutine gas_state_scale(gas_state, alpha)

    !> Existing gas state.
    type(gas_state_t), intent(inout) :: gas_state
    !> Scale factor.
    real(kind=dp), intent(in) :: alpha

    if (gas_state_is_allocated(gas_state)) then
       gas_state%mix_rat = gas_state%mix_rat * alpha
    end if

  end subroutine gas_state_scale

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Adds the given gas_state_delta.
  subroutine gas_state_add(gas_state, gas_state_delta)

    !> Existing gas state.
    type(gas_state_t), intent(inout) :: gas_state
    !> Incremental state.
    type(gas_state_t), intent(in) :: gas_state_delta

    if (gas_state_is_allocated(gas_state_delta)) then
       if (gas_state_is_allocated(gas_state)) then
          gas_state%mix_rat = gas_state%mix_rat + gas_state_delta%mix_rat
       else
          gas_state%mix_rat = gas_state_delta%mix_rat
       end if
    end if

  end subroutine gas_state_add

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Adds the given \c gas_state_delta scaled by \c alpha.
  !!
  !! Does gas_state += alpha * gas_state_delta.
  subroutine gas_state_add_scaled(gas_state, gas_state_delta, alpha)

    !> Existing gas state.
    type(gas_state_t), intent(inout) :: gas_state
    !> Incremental state.
    type(gas_state_t), intent(in) :: gas_state_delta
    !> Scale factor.
    real(kind=dp), intent(in) :: alpha

    if (gas_state_is_allocated(gas_state_delta)) then
       if (gas_state_is_allocated(gas_state)) then
          gas_state%mix_rat = gas_state%mix_rat &
               + alpha * gas_state_delta%mix_rat
       else
          gas_state%mix_rat = alpha * gas_state_delta%mix_rat
       end if
    end if

  end subroutine gas_state_add_scaled

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Subtracts the given gas_state_delta.
  subroutine gas_state_sub(gas_state, gas_state_delta)

    !> Existing gas state.
    type(gas_state_t), intent(inout) :: gas_state
    !> Incremental state.
    type(gas_state_t), intent(in) :: gas_state_delta

    if (gas_state_is_allocated(gas_state_delta)) then
       if (gas_state_is_allocated(gas_state)) then
          gas_state%mix_rat = gas_state%mix_rat - gas_state_delta%mix_rat
       else
          gas_state%mix_rat = - gas_state_delta%mix_rat
       end if
    end if

  end subroutine gas_state_sub

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Set any negative values to zero.
  subroutine gas_state_ensure_nonnegative(gas_state)

    !> Gas state.
    type(gas_state_t) :: gas_state

    if (gas_state_is_allocated(gas_state)) then
       gas_state%mix_rat = max(gas_state%mix_rat, 0d0)
    end if

  end subroutine gas_state_ensure_nonnegative

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Convert (mol m^{-3}) to (ppb).
  subroutine gas_state_mole_dens_to_ppb(gas_state, env_state)

    !> Gas state.
    type(gas_state_t), intent(inout) :: gas_state
    !> Environment state.
    type(env_state_t), intent(in) :: env_state

    call gas_state_scale(gas_state, 1d9 / env_state_air_molar_den(env_state))

  end subroutine gas_state_mole_dens_to_ppb

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef PMC_USE_CAMP
  !> Set CAMP gas-phase species concentrations
  subroutine gas_state_set_camp_conc(gas_state, camp_state, gas_data)

    !> Gas state
    class(gas_state_t), intent(in) :: gas_state
    !> CAMP state
    type(camp_state_t), intent(inout) :: camp_state
    !> Gas data
    type(gas_data_t), intent(in) :: gas_data

    real(kind=dp), parameter :: t_steam = 373.15 ! steam temperature (K)
    real(kind=dp) :: a, water_vp

    camp_state%state_var(1:size(gas_state%mix_rat)) = gas_state%mix_rat(:) &
         / 1000.0d0

    ! Convert relative humidity (1) to [H2O] (ppm)
    ! From MOSAIC code - reference to Seinfeld & Pandis page 181
    ! TODO Figure out how to have consistent RH<->ppm conversions
    ! (There is only one environmental state for PartMC runs
    call assert(590005048, associated(camp_state%env_states(1)%val))
    a = 1.0 - t_steam / camp_state%env_states(1)%val%temp
    a = (((-0.1299 * a - 0.6445) * a - 1.976) * a + 13.3185) * a
    water_vp = 101325.0 * exp(a)  ! (Pa)

    camp_state%state_var(gas_data%i_camp_water) = &
         camp_state%env_states(1)%val%rel_humid * water_vp * 1.0e6 &
         / camp_state%env_states(1)%val%pressure ! (ppm)

  end subroutine gas_state_set_camp_conc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get CAMP gas-phase species concentrations
  subroutine gas_state_get_camp_conc(gas_state, camp_state)

    !> Gas state
    class(gas_state_t), intent(inout) :: gas_state
    !> CAMP state
    type(camp_state_t), intent(in) :: camp_state

    gas_state%mix_rat(:) = camp_state%state_var(1:size(gas_state%mix_rat)) &
         * 1000.0d0

  end subroutine gas_state_get_camp_conc
#endif

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
       gas_state = gas_state_list(1)
       rate = rate_list(1)
    elseif (p == n) then
       ! after the end, just use the last state and rate
       gas_state = gas_state_list(n)
       rate = rate_list(n)
    else
       ! in the middle, use the previous state
       gas_state = gas_state_list(p)
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
    character(len=SPEC_LINE_MAX_VAR_LEN), allocatable :: species_name(:)
    real(kind=dp), allocatable :: species_data(:,:)

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

    call spec_file_read_real_named_array(file, 0, species_name, &
         species_data)

    ! check the data size
    n_species = size(species_data, 1)
    if (.not. ((size(species_data, 2) == 1) .or. (n_species == 0))) then
       call die_msg(686719840, 'each line in ' // trim(file%name) &
            // ' must contain exactly one data value')
    end if

    ! copy over the data
    call gas_state_set_size(gas_state, gas_data_n_spec(gas_data))
    do i = 1,n_species
       species = gas_data_spec_by_name(gas_data, species_name(i))
       if (species == 0) then
          call die_msg(129794076, 'unknown species ' // &
               trim(species_name(i)) // ' in file ' // trim(file%name))
       end if
       gas_state%mix_rat(species) = species_data(i,1)
    end do

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
    real(kind=dp), allocatable :: times(:)
    !> Rates (s^{-1}).
    real(kind=dp), allocatable :: rates(:)
    !> Gas states.
    type(gas_state_t), allocatable :: gas_states(:)

    integer :: n_lines, species, i, n_time, i_time
    character(len=SPEC_LINE_MAX_VAR_LEN), allocatable :: species_name(:)
    real(kind=dp), allocatable :: species_data(:,:)

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
    !!   values at the corresponding times
    !! - the third and subsequent lines specify gas species, one
    !!   species per line, with each line beginning with the species
    !!   name and followed by \f$N\f$ space-separated scalars giving
    !!   the gas state of that species at the corresponding times
    !!
    !! The units and meanings of the rate and species lines depends on
    !! the type of gas profile:
    !! - emissions gas profiles have dimensionless rates that are used
    !!   to scale the species rates and species giving emission rates
    !!   with units of mol/(m^2 s) --- the emission rate is divided by
    !!   the current mixing layer height to give a per-volume emission
    !!   rate
    !! - background gas profiles have rates with units s^{-1} that are
    !!   dilution rates and species with units of ppb (parts per
    !!   billion) that are the background mixing ratios
    !!
    !! The species names must be those specified by the \ref
    !! input_format_gas_data. Any species not listed are taken to be
    !! zero.
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
    call spec_file_read_real_named_array(file, 0, species_name, &
         species_data)

    ! check the data size
    n_lines = size(species_data, 1)
    if (n_lines < 2) then
       call die_msg(291542946, 'insufficient data lines in file ' &
            // trim(file%name))
    end if
    if (trim(species_name(1)) /= 'time') then
       call die_msg(525127793, 'row 1 in file ' &
            // trim(file%name) // ' must start with: time')
    end if
    if (trim(species_name(2)) /= 'rate') then
       call die_msg(506981322, 'row 2 in file ' &
            // trim(file%name) // ' must start with: rate')
    end if
    n_time = size(species_data, 2)
    if (n_time < 1) then
       call die_msg(398532628, 'each line in file ' &
            // trim(file%name) // ' must contain at least one data value')
    end if

    ! copy over the data
    times = species_data(1,:)
    rates = species_data(2,:)
    if (allocated(gas_states)) deallocate(gas_states)
    allocate(gas_states(n_time))
    do i_time = 1,n_time
       call gas_state_set_size(gas_states(i_time), gas_data_n_spec(gas_data))
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

  end subroutine spec_file_read_gas_states_times_rates

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Average val over all processes.
  subroutine gas_state_mix(val)

    !> Value to average.
    type(gas_state_t), intent(inout) :: val

#ifdef PMC_USE_MPI
    type(gas_state_t) :: val_avg

    call gas_state_set_size(val_avg, size(val%mix_rat))
    call pmc_mpi_allreduce_average_real_array(val%mix_rat, val_avg%mix_rat)
    val%mix_rat = val_avg%mix_rat
#endif

  end subroutine gas_state_mix

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Average val over all processes, with the result only on the root
  !> process.
  subroutine gas_state_reduce_avg(val)

    !> Value to average.
    type(gas_state_t), intent(inout) :: val

#ifdef PMC_USE_MPI
    type(gas_state_t) :: val_avg

    call gas_state_set_size(val_avg, size(val%mix_rat))
    call pmc_mpi_reduce_avg_real_array(val%mix_rat, val_avg%mix_rat)
    if (pmc_mpi_rank() == 0) then
       val%mix_rat = val_avg%mix_rat
    end if
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

    call pmc_nc_read_real_1d(ncid, gas_state%mix_rat, "gas_mixing_ratio")

  end subroutine gas_state_input_netcdf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_gas_state
