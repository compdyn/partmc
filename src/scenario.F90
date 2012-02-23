! Copyright (C) 2005-2012 Nicole Riemer and Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_scenario module.

!> The scenario_t structure and associated subroutines.
module pmc_scenario

  use pmc_gas_state
  use pmc_aero_dist
  use pmc_util
  use pmc_env_state
  use pmc_spec_file
  use pmc_aero_data
  use pmc_gas_data
  use pmc_mpi
#ifdef PMC_USE_MPI
  use mpi
#endif
  
  !> Environment data.
  !!
  !! This is everything needed to support the current environment
  !! state. The environment data is not time-dependent, whereas the
  !! environment state in env_state_t is time-dependent.
  !!
  !! The temperature, emissions and background states are profiles
  !! prescribed as functions of time by giving a number of times and
  !! the corresponding data. Simple data such as temperature is
  !! linearly interpoloated between times, with constant interpolation
  !! outside of the range of times. Gases and aerosols are
  !! interpolated with gas_state_interp_1d() and
  !! aero_dist_interp_1d(), respectively.
  type scenario_t
     !> Temperature set-point times (s).
     real(kind=dp), pointer :: temp_time(:)
     !> Temperatures at set-points (K).
     real(kind=dp), pointer :: temp(:)

     !> Height set-point times (s).
     real(kind=dp), pointer :: height_time(:)
     !> Heights at set-points (m).
     real(kind=dp), pointer :: height(:)

     !> Gas emission set-point times (s).
     real(kind=dp), pointer :: gas_emission_time(:)
     !> Gas emisssion rates at set-points (s^{-1}).
     real(kind=dp), pointer :: gas_emission_rate(:)
     !> Gas emissions at set-points.
     type(gas_state_t), pointer :: gas_emission(:)

     !> Gas-background dilution set-point times (s).
     real(kind=dp), pointer :: gas_dilution_time(:)
     !> Gas-background dilution rates at set-points (s^{-1}).
     real(kind=dp), pointer :: gas_dilution_rate(:)
     !> Background gas mixing ratios at set-points.
     type(gas_state_t), pointer :: gas_background(:)

     !> Aerosol emission set-points times (s).
     real(kind=dp), pointer :: aero_emission_time(:)
     !> Aerosol emission rates at set-points (s^{-1}).
     real(kind=dp), pointer :: aero_emission_rate(:)
     !> Aerosol emissions at set-points.
     type(aero_dist_t), pointer :: aero_emission(:)

     !> Aerosol-background dilution set-point times (s).
     real(kind=dp), pointer :: aero_dilution_time(:)
     !> Aerosol-background dilution rates at set-points (s^{-1}).
     real(kind=dp), pointer :: aero_dilution_rate(:)
     !> Aerosol background at set-points.
     type(aero_dist_t), pointer :: aero_background(:)
  end type scenario_t
  
contains
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Allocate an scenario.
  subroutine scenario_allocate(scenario)

    !> Environment data.
    type(scenario_t), intent(out) :: scenario

    allocate(scenario%temp_time(0))
    allocate(scenario%temp(0))

    allocate(scenario%height_time(0))
    allocate(scenario%height(0))

    allocate(scenario%gas_emission_time(0))
    allocate(scenario%gas_emission_rate(0))
    allocate(scenario%gas_emission(0))

    allocate(scenario%gas_dilution_time(0))
    allocate(scenario%gas_dilution_rate(0))
    allocate(scenario%gas_background(0))

    allocate(scenario%aero_emission_time(0))
    allocate(scenario%aero_emission_rate(0))
    allocate(scenario%aero_emission(0))

    allocate(scenario%aero_dilution_time(0))
    allocate(scenario%aero_dilution_rate(0))
    allocate(scenario%aero_background(0))

  end subroutine scenario_allocate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Free all storage.
  subroutine scenario_deallocate(scenario)

    !> Environment data.
    type(scenario_t), intent(inout) :: scenario

    integer :: i

    deallocate(scenario%temp_time)
    deallocate(scenario%temp)

    deallocate(scenario%height_time)
    deallocate(scenario%height)

    do i = 1,size(scenario%gas_emission)
       call gas_state_deallocate(scenario%gas_emission(i))
    end do
    deallocate(scenario%gas_emission_time)
    deallocate(scenario%gas_emission_rate)
    deallocate(scenario%gas_emission)

    do i = 1,size(scenario%gas_background)
       call gas_state_deallocate(scenario%gas_background(i))
    end do
    deallocate(scenario%gas_dilution_time)
    deallocate(scenario%gas_dilution_rate)
    deallocate(scenario%gas_background)

    do i = 1,size(scenario%aero_emission)
       call aero_dist_deallocate(scenario%aero_emission(i))
    end do
    deallocate(scenario%aero_emission_time)
    deallocate(scenario%aero_emission_rate)
    deallocate(scenario%aero_emission)

    do i = 1,size(scenario%aero_background)
       call aero_dist_deallocate(scenario%aero_background(i))
    end do
    deallocate(scenario%aero_dilution_time)
    deallocate(scenario%aero_dilution_rate)
    deallocate(scenario%aero_background)

  end subroutine scenario_deallocate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Copy structure.
  subroutine scenario_copy(scenario_from, scenario_to)

    !> Source environment data.
    type(scenario_t), intent(in) :: scenario_from
    !> Destination environment data.
    type(scenario_t), intent(inout) :: scenario_to

    integer :: i

    call scenario_deallocate(scenario_to)

    allocate(scenario_to%temp_time( &
         size(scenario_from%temp_time)))
    scenario_to%temp_time = scenario_from%temp_time
    allocate(scenario_to%temp( &
         size(scenario_from%temp)))
    scenario_to%temp = scenario_from%temp

    allocate(scenario_to%height_time( &
         size(scenario_from%height_time)))
    scenario_to%height_time = scenario_from%height_time
    allocate(scenario_to%height( &
         size(scenario_from%height)))
    scenario_to%height = scenario_from%height

    allocate(scenario_to%gas_emission_time( &
         size(scenario_from%gas_emission_time)))
    scenario_to%gas_emission_time = scenario_from%gas_emission_time
    allocate(scenario_to%gas_emission_rate( &
         size(scenario_from%gas_emission_rate)))
    scenario_to%gas_emission_rate = scenario_from%gas_emission_rate
    allocate(scenario_to%gas_emission( &
         size(scenario_from%gas_emission)))
    do i = 1,size(scenario_from%gas_emission)
       call gas_state_allocate(scenario_to%gas_emission(i))
       call gas_state_copy(scenario_from%gas_emission(i), &
            scenario_to%gas_emission(i))
    end do

    allocate(scenario_to%gas_dilution_time( &
         size(scenario_from%gas_dilution_time)))
    scenario_to%gas_dilution_time = scenario_from%gas_dilution_time
    allocate(scenario_to%gas_dilution_rate( &
         size(scenario_from%gas_dilution_rate)))
    scenario_to%gas_dilution_rate = scenario_from%gas_dilution_rate
    allocate(scenario_to%gas_background( &
         size(scenario_from%gas_background)))
    do i = 1,size(scenario_from%gas_background)
       call gas_state_allocate(scenario_to%gas_background(i))
       call gas_state_copy(scenario_from%gas_background(i), &
            scenario_to%gas_background(i))
    end do

    allocate(scenario_to%aero_emission_time( &
         size(scenario_from%aero_emission_time)))
    scenario_to%aero_emission_time = scenario_from%aero_emission_time
    allocate(scenario_to%aero_emission_rate( &
         size(scenario_from%aero_emission_rate)))
    scenario_to%aero_emission_rate = scenario_from%aero_emission_rate
    allocate(scenario_to%aero_emission( &
         size(scenario_from%aero_emission)))
    do i = 1,size(scenario_from%aero_emission)
       call aero_dist_allocate(scenario_to%aero_emission(i))
       call aero_dist_copy(scenario_from%aero_emission(i), &
            scenario_to%aero_emission(i))
    end do

    allocate(scenario_to%aero_dilution_time( &
         size(scenario_from%aero_dilution_time)))
    scenario_to%aero_dilution_time = scenario_from%aero_dilution_time
    allocate(scenario_to%aero_dilution_rate( &
         size(scenario_from%aero_dilution_rate)))
    scenario_to%aero_dilution_rate = scenario_from%aero_dilution_rate
    allocate(scenario_to%aero_background( &
         size(scenario_from%aero_background)))
    do i = 1,size(scenario_from%aero_background)
       call aero_dist_allocate(scenario_to%aero_background(i))
       call aero_dist_copy(scenario_from%aero_background(i), &
            scenario_to%aero_background(i))
    end do

  end subroutine scenario_copy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Initialize the time-dependent contents of the
  !> environment. Thereafter scenario_update_state() should be used.
  subroutine scenario_init_state(scenario, env_state, time)

    !> Environment data.
    type(scenario_t), intent(in) :: scenario
    !> Environment state to update.
    type(env_state_t), intent(inout) :: env_state
    !> Current time (s).
    real(kind=dp), intent(in) :: time

    env_state%temp = interp_1d(scenario%temp_time, scenario%temp, time)
    env_state%height = interp_1d(scenario%height_time, scenario%height, time)
    env_state%elapsed_time = time

    ! init gas and aerosol emissions and background
    call gas_state_interp_1d(scenario%gas_emission, &
         scenario%gas_emission_time, scenario%gas_emission_rate, &
         time, env_state%gas_emissions, env_state%gas_emission_rate)
    call gas_state_interp_1d(scenario%gas_background, &
         scenario%gas_dilution_time, scenario%gas_dilution_rate, &
         time, env_state%gas_background, env_state%gas_dilution_rate)
    call aero_dist_interp_1d(scenario%aero_emission, &
         scenario%aero_emission_time, scenario%aero_emission_rate, &
         time, env_state%aero_emissions, env_state%aero_emission_rate)
    call aero_dist_interp_1d(scenario%aero_background, &
         scenario%aero_dilution_time, scenario%aero_dilution_rate, &
         time, env_state%aero_background, env_state%aero_dilution_rate)
    
  end subroutine scenario_init_state
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Update time-dependent contents of the environment.
  !> scenario_init_state() should have been called at the start.
  subroutine scenario_update_state(scenario, env_state, time, &
       update_rel_humid)

    !> Environment data.
    type(scenario_t), intent(in) :: scenario
    !> Environment state to update.
    type(env_state_t), intent(inout) :: env_state
    !> Current time (s).
    real(kind=dp), intent(in) :: time
    !> Whether to update the relative humidity.
    logical, intent(in) :: update_rel_humid
    
    !> Ambient water vapor pressure (Pa).
    real(kind=dp) :: pmv

    ! Update temperature and adjust relative humidity to maintain
    ! water mixing ratio. This uses the fact that we keep the total
    ! pressure constant, so ambient water vapor pressure is also
    ! constant (whatever temperature and hence volume does).
    pmv = env_state_sat_vapor_pressure(env_state) * env_state%rel_humid
    env_state%temp = interp_1d(scenario%temp_time, scenario%temp, time)
    if (update_rel_humid) then
       env_state%rel_humid = pmv / env_state_sat_vapor_pressure(env_state)
    end if

    env_state%height = interp_1d(scenario%height_time, scenario%height, time)
    env_state%elapsed_time = time

    ! update gas and aerosol emissions and background
    call gas_state_interp_1d(scenario%gas_emission, &
         scenario%gas_emission_time, scenario%gas_emission_rate, &
         time, env_state%gas_emissions, env_state%gas_emission_rate)
    call gas_state_interp_1d(scenario%gas_background, &
         scenario%gas_dilution_time, scenario%gas_dilution_rate, &
         time, env_state%gas_background, env_state%gas_dilution_rate)
    call aero_dist_interp_1d(scenario%aero_emission, &
         scenario%aero_emission_time, scenario%aero_emission_rate, &
         time, env_state%aero_emissions, env_state%aero_emission_rate)
    call aero_dist_interp_1d(scenario%aero_background, &
         scenario%aero_dilution_time, scenario%aero_dilution_rate, &
         time, env_state%aero_background, env_state%aero_dilution_rate)

  end subroutine scenario_update_state
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Whether any of the contained aerosol modes are of the given type.
  elemental logical function scenario_contains_aero_mode_type(scenario, &
       aero_mode_type)

    !> Environment data.
    type(scenario_t), intent(in) :: scenario
    !> Aerosol mode type to test for.
    integer, intent(in) :: aero_mode_type

    scenario_contains_aero_mode_type &
         = any(aero_dist_contains_aero_mode_type(scenario%aero_emission, &
         aero_mode_type)) &
         .or. any(aero_dist_contains_aero_mode_type(scenario%aero_background, &
         aero_mode_type))

  end function scenario_contains_aero_mode_type

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read environment data from an spec file.
  subroutine spec_file_read_scenario(file, gas_data, aero_data, scenario)

    !> Spec file.
    type(spec_file_t), intent(inout) :: file
    !> Gas data values.
    type(gas_data_t), intent(in) :: gas_data
    !> Aerosol data.
    type(aero_data_t), intent(inout) :: aero_data
    !> Environment data.
    type(scenario_t), intent(inout) :: scenario

    character(len=PMC_MAX_FILENAME_LEN) :: sub_filename
    type(spec_file_t) :: sub_file

    ! note that we have to hard-code the list for doxygen below

    !> \page input_format_scenario Input File Format: Environment Data
    !!
    !! The environment parameters are divided into those specified at
    !! the start of the simulation and then either held constant or
    !! computed for the rest of the simulation, and those parameters
    !! given as prescribed profiles for the entire simulation
    !! duration. The variables below are for the second type --- for
    !! the computed values see \ref input_format_env_state.
    !!
    !! The environment data parameters are:
    !! <ul>
    !! <li> \b temp_profile (string): the name of the file from which to
    !!      read the temperature profile --- the file format should be
    !!      \subpage input_format_temp_profile
    !! <li> \b height_profile (string): the name of the file from which
    !!      to read the mixing layer height profile --- the file format
    !!      should be \subpage input_format_height_profile
    !! <li> \b gas_emissions (string): the name of the file from which to
    !!      read the gas emissions profile --- the file format should be
    !!      \subpage input_format_gas_profile
    !! <li> \b gas_background (string): the name of the file from which
    !!      to read the gas background profile --- the file format should
    !!      be \subpage input_format_gas_profile
    !! <li> \b aero_emissions (string): the name of the file from which
    !!      to read the aerosol emissions profile --- the file format
    !!      should be \subpage input_format_aero_dist_profile
    !! <li> \b aero_background (string): the name of the file from which
    !!      to read the aerosol background profile --- the file format
    !!      should be \subpage input_format_aero_dist_profile
    !! </ul>
    !!
    !! See also:
    !!   - \ref spec_file_format --- the input file text format

    !> \page input_format_temp_profile Input File Format: Temperature Profile
    !!
    !! A temperature profile input file must consist of two lines:
    !! - the first line must begin with \c time and should be followed
    !!   by \f$N\f$ space-separated real scalars, giving the times (in
    !!   s after the start of the simulation) of the temperature set
    !!   points --- the times must be in increasing order
    !! - the second line must begin with \c temp and should be followed
    !!   by \f$N\f$ space-separated real scalars, giving the
    !!   temperatures (in K) at the corresponding times
    !!
    !! The temperature profile is linearly interpolated between the
    !! specified times, while before the first time it takes the first
    !! temperature value and after the last time it takes the last
    !! temperature value.
    !!
    !! Example:
    !! <pre>
    !! time  0    600  1800  # time (in s) after simulation start
    !! temp  270  290  280   # temperature (in K)
    !! </pre>
    !! Here the temperature starts at 270&nbsp;K at the start of the
    !! simulation, rises to 290&nbsp;K after 10&nbsp;min, and then
    !! falls again to 280&nbsp;K at 30&nbsp;min. Between these times
    !! the temperature is linearly interpolated, while after
    !! 30&nbsp;min it is held constant at 280&nbsp;K.
    !!
    !! See also:
    !!   - \ref spec_file_format --- the input file text format
    !!   - \ref input_format_scenario --- the environment data
    !!     containing the temperature profile

    !> \page input_format_height_profile Input File Format: Mixing Layer Height Profile
    !!
    !! A mixing layer height profile input file must consist of two
    !! lines:
    !! - the first line must begin with \c time and should be followed
    !!   by \f$N\f$ space-separated real scalars, giving the times (in
    !!   s after the start of the simulation) of the height set
    !!   points --- the times must be in increasing order
    !! - the second line must begin with \c height and should be
    !!   followed by \f$N\f$ space-separated real scalars, giving the
    !!   mixing layer heights (in m) at the corresponding times
    !!
    !! The mixing layer height profile is linearly interpolated
    !! between the specified times, while before the first time it
    !! takes the first height value and after the last time it takes
    !! the last height value.
    !!
    !! Example:
    !! <pre>
    !! time    0    600   1800  # time (in s) after simulation start
    !! height  500  1000  800   # mixing layer height (in m)
    !! </pre>
    !! Here the mixing layer height starts at 500&nbsp;m at the start
    !! of the simulation, rises to 1000&nbsp;m after 10&nbsp;min, and
    !! then falls again to 800&nbsp;m at 30&nbsp;min. Between these
    !! times the mixing layer height is linearly interpolated, while
    !! after 30&nbsp;min it is held constant at 800&nbsp;m.
    !!
    !! See also:
    !!   - \ref spec_file_format --- the input file text format
    !!   - \ref input_format_scenario --- the environment data
    !!     containing the mixing layer height profile

    ! temperature profile
    call spec_file_read_string(file, "temp_profile", sub_filename)
    call spec_file_open(sub_filename, sub_file)
    call spec_file_read_timed_real_array(sub_file, "temp", &
         scenario%temp_time, scenario%temp)
    call spec_file_close(sub_file)

    ! height profile
    call spec_file_read_string(file, "height_profile", sub_filename)
    call spec_file_open(sub_filename, sub_file)
    call spec_file_read_timed_real_array(sub_file, "height", &
         scenario%height_time, scenario%height)
    call spec_file_close(sub_file)

    ! gas emissions profile
    call spec_file_read_string(file, "gas_emissions", sub_filename)
    call spec_file_open(sub_filename, sub_file)
    call spec_file_read_gas_states_times_rates(sub_file, gas_data, &
         scenario%gas_emission_time, scenario%gas_emission_rate, &
         scenario%gas_emission)
    call spec_file_close(sub_file)

    ! gas background profile
    call spec_file_read_string(file, "gas_background", sub_filename)
    call spec_file_open(sub_filename, sub_file)
    call spec_file_read_gas_states_times_rates(sub_file, gas_data, &
         scenario%gas_dilution_time, scenario%gas_dilution_rate, &
         scenario%gas_background)
    call spec_file_close(sub_file)

    ! aerosol emissions profile
    call spec_file_read_string(file, "aero_emissions", sub_filename)
    call spec_file_open(sub_filename, sub_file)
    call spec_file_read_aero_dists_times_rates(sub_file, aero_data, &
         scenario%aero_emission_time, scenario%aero_emission_rate, &
         scenario%aero_emission)
    call spec_file_close(sub_file)

    ! aerosol background profile
    call spec_file_read_string(file, "aero_background", sub_filename)
    call spec_file_open(sub_filename, sub_file)
    call spec_file_read_aero_dists_times_rates(sub_file, aero_data, &
         scenario%aero_dilution_time, scenario%aero_dilution_rate, &
         scenario%aero_background)
    call spec_file_close(sub_file)

  end subroutine spec_file_read_scenario

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Determines the number of bytes required to pack the given value.
  integer function pmc_mpi_pack_size_scenario(val)

    !> Value to pack.
    type(scenario_t), intent(in) :: val

    integer :: total_size, i, n

    total_size = &
         pmc_mpi_pack_size_real_array(val%temp_time) &
         + pmc_mpi_pack_size_real_array(val%temp) &
         + pmc_mpi_pack_size_real_array(val%height_time) &
         + pmc_mpi_pack_size_real_array(val%height) &
         + pmc_mpi_pack_size_real_array(val%gas_emission_time) &
         + pmc_mpi_pack_size_real_array(val%gas_emission_rate) &
         + pmc_mpi_pack_size_real_array(val%gas_dilution_time) &
         + pmc_mpi_pack_size_real_array(val%gas_dilution_rate) &
         + pmc_mpi_pack_size_real_array(val%aero_emission_time) &
         + pmc_mpi_pack_size_real_array(val%aero_emission_rate) &
         + pmc_mpi_pack_size_real_array(val%aero_dilution_time) &
         + pmc_mpi_pack_size_real_array(val%aero_dilution_rate)
    do i = 1,size(val%gas_emission)
       total_size = total_size &
            + pmc_mpi_pack_size_gas_state(val%gas_emission(i))
    end do
    do i = 1,size(val%gas_background)
       total_size = total_size &
            + pmc_mpi_pack_size_gas_state(val%gas_background(i))
    end do
    do i = 1,size(val%aero_emission)
       total_size = total_size &
            + pmc_mpi_pack_size_aero_dist(val%aero_emission(i))
    end do
    do i = 1,size(val%aero_background)
       total_size = total_size &
            + pmc_mpi_pack_size_aero_dist(val%aero_background(i))
    end do

    pmc_mpi_pack_size_scenario = total_size

  end function pmc_mpi_pack_size_scenario

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Packs the given value into the buffer, advancing position.
  subroutine pmc_mpi_pack_scenario(buffer, position, val)

    !> Memory buffer.
    character, intent(inout) :: buffer(:)
    !> Current buffer position.
    integer, intent(inout) :: position
    !> Value to pack.
    type(scenario_t), intent(in) :: val

#ifdef PMC_USE_MPI
    integer :: prev_position, i

    prev_position = position
    call pmc_mpi_pack_real_array(buffer, position, val%temp_time)
    call pmc_mpi_pack_real_array(buffer, position, val%temp)
    call pmc_mpi_pack_real_array(buffer, position, val%height_time)
    call pmc_mpi_pack_real_array(buffer, position, val%height)
    call pmc_mpi_pack_real_array(buffer, position, val%gas_emission_time)
    call pmc_mpi_pack_real_array(buffer, position, val%gas_emission_rate)
    call pmc_mpi_pack_real_array(buffer, position, val%gas_dilution_time)
    call pmc_mpi_pack_real_array(buffer, position, val%gas_dilution_rate)
    call pmc_mpi_pack_real_array(buffer, position, val%aero_emission_time)
    call pmc_mpi_pack_real_array(buffer, position, val%aero_emission_rate)
    call pmc_mpi_pack_real_array(buffer, position, val%aero_dilution_time)
    call pmc_mpi_pack_real_array(buffer, position, val%aero_dilution_rate)
    do i = 1,size(val%gas_emission)
       call pmc_mpi_pack_gas_state(buffer, position, val%gas_emission(i))
    end do
    do i = 1,size(val%gas_background)
       call pmc_mpi_pack_gas_state(buffer, position, val%gas_background(i))
    end do
    do i = 1,size(val%aero_emission)
       call pmc_mpi_pack_aero_dist(buffer, position, val%aero_emission(i))
    end do
    do i = 1,size(val%aero_background)
       call pmc_mpi_pack_aero_dist(buffer, position, val%aero_background(i))
    end do
    call assert(639466930, &
         position - prev_position <= pmc_mpi_pack_size_scenario(val))
#endif

  end subroutine pmc_mpi_pack_scenario

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Unpacks the given value from the buffer, advancing position.
  subroutine pmc_mpi_unpack_scenario(buffer, position, val)

    !> Memory buffer.
    character, intent(inout) :: buffer(:)
    !> Current buffer position.
    integer, intent(inout) :: position
    !> Value to pack.
    type(scenario_t), intent(inout) :: val

#ifdef PMC_USE_MPI
    integer :: prev_position, i

    call scenario_deallocate(val)
    call scenario_allocate(val)
    prev_position = position
    call pmc_mpi_unpack_real_array(buffer, position, val%temp_time)
    call pmc_mpi_unpack_real_array(buffer, position, val%temp)
    call pmc_mpi_unpack_real_array(buffer, position, val%height_time)
    call pmc_mpi_unpack_real_array(buffer, position, val%height)
    call pmc_mpi_unpack_real_array(buffer, position, val%gas_emission_time)
    call pmc_mpi_unpack_real_array(buffer, position, val%gas_emission_rate)
    call pmc_mpi_unpack_real_array(buffer, position, val%gas_dilution_time)
    call pmc_mpi_unpack_real_array(buffer, position, val%gas_dilution_rate)
    call pmc_mpi_unpack_real_array(buffer, position, val%aero_emission_time)
    call pmc_mpi_unpack_real_array(buffer, position, val%aero_emission_rate)
    call pmc_mpi_unpack_real_array(buffer, position, val%aero_dilution_time)
    call pmc_mpi_unpack_real_array(buffer, position, val%aero_dilution_rate)
    deallocate(val%gas_emission)
    deallocate(val%gas_background)
    deallocate(val%aero_emission)
    deallocate(val%aero_background)
    allocate(val%gas_emission(size(val%gas_emission_time)))
    allocate(val%gas_background(size(val%gas_dilution_time)))
    allocate(val%aero_emission(size(val%aero_emission_time)))
    allocate(val%aero_background(size(val%aero_dilution_time)))
    do i = 1,size(val%gas_emission)
       call gas_state_allocate(val%gas_emission(i))
       call pmc_mpi_unpack_gas_state(buffer, position, val%gas_emission(i))
    end do
    do i = 1,size(val%gas_background)
       call gas_state_allocate(val%gas_background(i))
       call pmc_mpi_unpack_gas_state(buffer, position, val%gas_background(i))
    end do
    do i = 1,size(val%aero_emission)
       call aero_dist_allocate(val%aero_emission(i))
       call pmc_mpi_unpack_aero_dist(buffer, position, val%aero_emission(i))
    end do
    do i = 1,size(val%aero_background)
       call aero_dist_allocate(val%aero_background(i))
       call pmc_mpi_unpack_aero_dist(buffer, position, val%aero_background(i))
    end do
    call assert(611542570, &
         position - prev_position <= pmc_mpi_pack_size_scenario(val))
#endif

  end subroutine pmc_mpi_unpack_scenario

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
end module pmc_scenario
