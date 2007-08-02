! Copyright (C) 2005-2007 Nicole Riemer and Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.
!
! Environment parameters.
!
! The temperature profile is proscribed as a function of time by
! giving a number of times and the corresponding temperatures. Linear
! interpolation is used between the times, with constant interpolation
! outside of the range of times.

module pmc_env

  use pmc_gas_state
  use pmc_aero_dist
  
  type env_t
     real*8 :: temp                     ! temperature (K)
     real*8 :: rel_humid                ! relative humidity (1)
     real*8 :: pressure                 ! ambient pressure (Pa)
     real*8 :: air_den                  ! air density (kg m^{-3})
     real*8 :: longitude                ! longitude (degrees)
     real*8 :: latitude                 ! latitude (degrees)
     real*8 :: altitude                 ! altitude (m)
     real*8 :: start_time               ! start time (s since 00:00 UTC)
     integer :: start_day               ! start day of year (UTC)
     integer :: n_temp                  ! number of temperature set-points
     real*8, pointer :: temp_time(:)    ! times at temp set-points (s)
     real*8, pointer :: temp_set(:)     ! temps at temp set-points (K)
     type(gas_state_t) :: gas_emissions ! gas emissions
     real*8 :: gas_emission_rate        ! gas emisssion rate (s^{-1})
     type(gas_state_t) :: gas_background ! background gas concentrations
     real*8 :: gas_dilution_rate        ! gas-background dilution rate (s^{-1})
     type(aero_dist_t) :: aero_emissions ! aerosol emissions
     real*8 :: aero_emission_rate       ! aerosol emisssion rate (s^{-1})
     type(aero_dist_t) :: aero_background ! aerosol background
     real*8 :: aero_dilution_rate       ! aero-background dilution rate (s^{-1})
  end type env_t
  
contains
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine env_alloc(env)

    ! Allocate an empty environment.

    type(env_t), intent(out) :: env   ! environment

    env%temp = 0d0
    env%rel_humid = 0d0
    env%pressure = 0d0
    env%air_den = 0d0
    env%longitude = 0d0
    env%latitude = 0d0
    env%altitude = 0d0
    env%start_time = 0d0
    env%start_day = 0

    call env_temp_alloc(env, 0)
    call gas_state_alloc(env%gas_emissions, 0)
    call gas_state_alloc(env%gas_background, 0)
    env%gas_emission_rate = 0d0
    env%gas_dilution_rate = 0d0
    call aero_dist_alloc(env%aero_emissions, 0, 0, 0)
    call aero_dist_alloc(env%aero_background, 0, 0, 0)
    env%aero_emission_rate = 0d0
    env%aero_dilution_rate = 0d0

  end subroutine env_alloc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine env_free(env)

    ! Free all storage.

    type(env_t), intent(out) :: env   ! environment

    call env_temp_free(env)
    call gas_state_free(env%gas_emissions)
    call gas_state_free(env%gas_background)
    call aero_dist_free(env%aero_emissions)
    call aero_dist_free(env%aero_background)

  end subroutine env_free

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine env_temp_alloc(env, n_temp)

    ! Allocate storage for a given number of temperature set points.

    type(env_t), intent(inout) :: env   ! environment
    integer, intent(in) :: n_temp      ! number of temperature set-points

    env%n_temp = n_temp
    allocate(env%temp_time(n_temp))
    allocate(env%temp_set(n_temp))

  end subroutine env_temp_alloc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine env_temp_free(env)

    ! Free all storage.

    type(env_t), intent(inout) :: env   ! environment

    deallocate(env%temp_time)
    deallocate(env%temp_set)
    
  end subroutine env_temp_free

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine env_change_water_volume(env, aero_data, dv)
    
    ! Adds the given water volume to the water vapor and updates all
    ! environment quantities.
    
    use pmc_constants
    use pmc_aero_data
    
    type(env_t), intent(inout) :: env   ! environment state to update
    type(aero_data_t), intent(in)   :: aero_data ! aero_data constants
    real*8, intent(in) :: dv            ! conc of water added (m^3/m^3)
    
    real*8 pmv     ! ambient water vapor pressure (Pa)
    real*8 mv      ! ambient water vapor density (kg m^{-3})
                   ! pmv and mv are related by the factor molec_weight/(R*T)
    real*8 dmv     ! change of water density (kg m^{-3})
    
    dmv = dv * aero_data%density(aero_data%i_water)
    pmv = env_sat_vapor_pressure(env) * env%rel_humid
    mv = aero_data%molec_weight(aero_data%i_water)/(const%R*env%temp) * pmv
    mv = mv - dmv    
    env%rel_humid = const%R * env%temp / aero_data%molec_weight(aero_data%i_water) * mv &
         / env_sat_vapor_pressure(env)
    
  end subroutine env_change_water_volume
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine env_init(env, time)
    
    ! Initialize the time-dependent contents of the
    ! environment. Thereafter env_update() should be used.

    use pmc_util

    type(env_t), intent(inout) :: env   ! environment state to update
    real*8, intent(in) :: time          ! current time (s)

    env%temp = interp_1d(env%n_temp, env%temp_time, env%temp_set, time)
    
  end subroutine env_init
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine env_update(env, time)
    
    ! Update time-dependent contents of the environment. env_init()
    ! should have been called at the start.

    use pmc_util

    type(env_t), intent(inout) :: env   ! environment state to update
    real*8, intent(in) :: time          ! current time (s)
    
    real*8 pmv      ! ambient water vapor pressure (Pa)

    ! update temperature and relative humidity
    pmv = env_sat_vapor_pressure(env) * env%rel_humid
    env%temp = interp_1d(env%n_temp, env%temp_time, env%temp_set, time)
    env%rel_humid = pmv / env_sat_vapor_pressure(env)
    
  end subroutine env_update
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  real*8 function env_sat_vapor_pressure(env) ! Pa

    ! Computes the current saturation vapor pressure.
    
    use pmc_constants
    
    type(env_t), intent(in) :: env      ! environment state
    
    env_sat_vapor_pressure = const%p00 * 10d0**(7.45d0 * (env%temp - const%T0) &
         / (env%temp - 38d0)) ! Pa
    
  end function env_sat_vapor_pressure
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine env_update_gas_state(env, delta_t, gas_data, gas_state)

    ! Do emissions and background dilution from the environment.

    use pmc_gas_data
    use pmc_gas_state

    type(env_t), intent(in) :: env      ! current environment
    real*8, intent(in) :: delta_t       ! time increment to update over
    type(gas_data_t), intent(in) :: gas_data ! gas data values
    type(gas_state_t), intent(inout) :: gas_state ! gas state to update

    type(gas_state_t) :: emission, dilution

    call gas_state_alloc(emission, gas_data%n_spec)
    call gas_state_alloc(dilution, gas_data%n_spec)

    ! emission = delta_t * gas_emission_rate * gas_emissions
    call gas_state_copy(env%gas_emissions, emission)
    call gas_state_scale(emission, delta_t * env%gas_emission_rate)

    ! dilution = delta_t * gas_dilution_rate * (gas_background - gas_state)
    call gas_state_copy(env%gas_background, dilution)
    call gas_state_sub(dilution, gas_state)
    call gas_state_scale(dilution, delta_t * env%gas_dilution_rate)

    call gas_state_add(gas_state, emission)
    call gas_state_add(gas_state, dilution)

    call gas_state_free(emission)
    call gas_state_free(dilution)

  end subroutine env_update_gas_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine env_update_aero_state(env, delta_t, bin_grid, &
       aero_data, aero_state, aero_binned)

    ! Do emissions and background dilution from the environment for a
    ! particle aerosol distribution.

    use pmc_bin_grid
    use pmc_aero_data
    use pmc_aero_state
    use pmc_aero_binned

    type(env_t), intent(in) :: env      ! current environment
    real*8, intent(in) :: delta_t       ! time increment to update over
    type(bin_grid_t), intent(in) :: bin_grid ! bin grid
    type(aero_data_t), intent(in) :: aero_data ! aero data values
    type(aero_state_t), intent(inout) :: aero_state ! aero state to update
    type(aero_binned_t), intent(inout) :: aero_binned ! aero binned to update

    integer :: i
    real*8 :: sample_vol, sample_prop
    type(aero_state_t) :: aero_state_delta
    type(aero_binned_t) :: aero_binned_delta

    call aero_state_alloc(bin_grid%n_bin, aero_data%n_spec, aero_state_delta)
    call aero_binned_alloc(aero_binned_delta, bin_grid%n_bin, aero_data%n_spec)
    aero_state_delta%comp_vol = aero_state%comp_vol

    ! loss to background
    sample_prop = delta_t * env%aero_dilution_rate
    call aero_state_zero(aero_state_delta)
    call aero_state_sample(aero_state, aero_state_delta, sample_prop)
    call aero_state_to_binned(bin_grid, aero_data, aero_state_delta, &
         aero_binned_delta)
    call aero_binned_sub(aero_binned, aero_binned_delta)

    ! addition from background
    sample_vol = delta_t * env%aero_dilution_rate * aero_state%comp_vol
    call aero_state_zero(aero_state_delta)
    call aero_dist_sample(bin_grid, aero_data, env%aero_background, &
         sample_vol, aero_state_delta)
    call aero_state_to_binned(bin_grid, aero_data, aero_state_delta, &
         aero_binned_delta)
    call aero_state_add(aero_state, aero_state_delta)
    call aero_binned_add(aero_binned, aero_binned_delta)
    
    ! emissions
    sample_vol = delta_t * env%aero_emission_rate * aero_state%comp_vol
    call aero_state_zero(aero_state_delta)
    call aero_dist_sample(bin_grid, aero_data, env%aero_emissions, &
         sample_vol, aero_state_delta)
    call aero_state_to_binned(bin_grid, aero_data, aero_state_delta, &
         aero_binned_delta)
    call aero_state_add(aero_state, aero_state_delta)
    call aero_binned_add(aero_binned, aero_binned_delta)

    call aero_state_free(aero_state_delta)
    call aero_binned_free(aero_binned_delta)

  end subroutine env_update_aero_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine env_update_aero_binned(env, delta_t, bin_grid, &
       aero_data, aero_binned)

    ! Do emissions and background dilution from the environment for a
    ! binned aerosol distribution.

    use pmc_bin_grid
    use pmc_aero_data
    use pmc_aero_binned

    type(env_t), intent(in) :: env      ! current environment
    real*8, intent(in) :: delta_t       ! time increment to update over
    type(bin_grid_t), intent(in) :: bin_grid ! bin grid
    type(aero_data_t), intent(in) :: aero_data ! aero data values
    type(aero_binned_t), intent(inout) :: aero_binned ! aero binned to update

    type(aero_binned_t) :: emission, dilution

    call aero_binned_alloc(emission, bin_grid%n_bin, aero_data%n_spec)
    call aero_binned_alloc(dilution, bin_grid%n_bin, aero_data%n_spec)

    ! emission = delta_t * gas_emission_rate * gas_emissions
    call aero_binned_add_aero_dist(emission, bin_grid, env%aero_emissions)
    call aero_binned_scale(emission, delta_t * env%aero_emission_rate)

    ! dilution = delta_t * gas_dilution_rate * (gas_background - aero_binned)
    call aero_binned_add_aero_dist(dilution, bin_grid, env%aero_background)
    call aero_binned_sub(dilution, aero_binned)
    call aero_binned_scale(dilution, delta_t * env%aero_dilution_rate)

    call aero_binned_add(aero_binned, emission)
    call aero_binned_add(aero_binned, dilution)

    call aero_binned_free(emission)
    call aero_binned_free(dilution)

  end subroutine env_update_aero_binned

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine inout_write_env(file, env)
    
    ! Write full state.
    
    use pmc_inout
    
    type(inout_file_t), intent(inout) :: file ! file to write to
    type(env_t), intent(in) :: env      ! environment to write
    
    call inout_write_real(file, "temp(K)", env%temp)
    call inout_write_real(file, "rel_humidity(1)", env%rel_humid)
    call inout_write_real(file, "pressure(Pa)", env%pressure)
    call inout_write_real(file, "air_density(kg/m^3)", env%air_den)
    call inout_write_real(file, "longitude(deg)", env%longitude)
    call inout_write_real(file, "latitude(deg)", env%latitude)
    call inout_write_real(file, "altitude(m)", env%altitude)
    call inout_write_real(file, "start_time(s)", env%start_time)
    call inout_write_integer(file, "start_day(days)", env%start_day)
    call inout_write_integer(file, "num_temp", env%n_temp)
    call inout_write_real_array(file, "temp_time(s)", env%temp_time)
    call inout_write_real_array(file, "temp(K)", env%temp_set)
    call inout_write_gas_state(file, env%gas_emissions)
    call inout_write_real(file, "gas_emit_rate(1/s)", env%gas_emission_rate)
    call inout_write_gas_state(file, env%gas_background)
    call inout_write_real(file, "gas_dilute_rate(1/s)", env%gas_dilution_rate)
    call inout_write_aero_dist(file, env%aero_emissions)
    call inout_write_real(file, "aero_emit_rate(1/s)", env%aero_emission_rate)
    call inout_write_aero_dist(file, env%aero_background)
    call inout_write_real(file, "aero_dilute_rat(1/s)", env%aero_dilution_rate)

  end subroutine inout_write_env

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine inout_read_env(file, env)
    
    ! Read full state.
    
    use pmc_inout
    
    type(inout_file_t), intent(inout) :: file ! file to read from
    type(env_t), intent(out) :: env      ! environment to read
    
    call inout_read_real(file, "temp(K)", env%temp)
    call inout_read_real(file, "rel_humidity(1)", env%rel_humid)
    call inout_read_real(file, "pressure(Pa)", env%pressure)
    call inout_read_real(file, "air_density(kg/m^3)", env%air_den)
    call inout_read_real(file, "longitude(deg)", env%longitude)
    call inout_read_real(file, "latitude(deg)", env%latitude)
    call inout_read_real(file, "altitude(m)", env%altitude)
    call inout_read_real(file, "start_time(s)", env%start_time)
    call inout_read_integer(file, "start_day(days)", env%start_day)
    call inout_read_integer(file, "num_temp", env%n_temp)
    call inout_read_real_array(file, "temp_time(s)", env%temp_time)
    call inout_read_real_array(file, "temp(K)", env%temp_set)
    call inout_read_gas_state(file, env%gas_emissions)
    call inout_read_real(file, "gas_emit_rate(1/s)", env%gas_emission_rate)
    call inout_read_gas_state(file, env%gas_background)
    call inout_read_real(file, "gas_dilute_rate(1/s)", env%gas_dilution_rate)
    call inout_read_aero_dist(file, env%aero_emissions)
    call inout_read_real(file, "aero_emit_rate(1/s)", env%aero_emission_rate)
    call inout_read_aero_dist(file, env%aero_background)
    call inout_read_real(file, "aero_dilute_rat(1/s)", env%aero_dilution_rate)

  end subroutine inout_read_env

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine spec_read_env(file, bin_grid, gas_data, aero_data, env)

    ! Read environment specification from a inout file.

    use pmc_bin_grid
    use pmc_inout
    use pmc_aero_data
    use pmc_gas_data

    type(inout_file_t), intent(inout) :: file ! inout file
    type(bin_grid_t), intent(in) :: bin_grid ! bin grid
    type(gas_data_t), intent(in) :: gas_data ! gas data values
    type(aero_data_t), intent(in) :: aero_data ! aerosol data
    type(env_t), intent(out) :: env     ! environment data

    integer :: n_temp
    character(len=MAX_CHAR_LEN) :: read_name
    type(inout_file_t) :: read_file
    character(len=MAX_CHAR_LEN), pointer :: times_name(:), temp_name(:)
    real*8, pointer :: times_data(:,:), temp_data(:,:)

    ! read the tempurature data from the specified file
    call inout_read_string(file, 'temp_profile', read_name)
    call inout_open_read(read_name, read_file)
    call inout_read_real_named_array(read_file, 1, times_name, times_data)
    call inout_check_name(read_file, "time", times_name(1))
    ! FIXME: add a min_lines arg to inout_read_real_named_array to ensure that
    ! really got one line here
    call inout_read_real_named_array(read_file, 1, temp_name, temp_data)
    call inout_check_name(read_file, "temp", temp_name(1))
    call inout_close(read_file)

    ! check the data size
    n_temp = size(temp_data, 2)
    if (n_temp < 1) then
       write(0,*) 'ERROR: file ', trim(read_name), &
            ' must contain at least one line of data'
       call exit(1)
    end if
    if (size(times_data, 2) /= size(temp_data, 2)) then
       write(0,*) 'ERROR: file ', trim(read_name), &
            ' should contain exactly two lines with equal numbers of values'
       call exit(1)
    end if

    call env_temp_alloc(env, n_temp)
    env%temp_time = times_data(1,:)
    env%temp_set = temp_data(1,:)
    deallocate(times_name)
    deallocate(times_data)
    deallocate(temp_name)
    deallocate(temp_data)
    call inout_read_real(file, 'rel_humidity', env%rel_humid)
    call inout_read_real(file, 'pressure', env%pressure)
    call inout_read_real(file, 'air_density', env%air_den)
    call inout_read_real(file, 'latitude', env%latitude)
    call inout_read_real(file, 'longitude', env%longitude)
    call inout_read_real(file, 'altitude', env%altitude)
    call inout_read_real(file, 'start_time', env%start_time)
    call inout_read_integer(file, 'start_day', env%start_day)

    call spec_read_gas_state(file, gas_data, 'gas_emissions', env%gas_emissions)
    call inout_read_real(file, 'gas_emission_rate', env%gas_emission_rate)
    call spec_read_gas_state(file, gas_data, 'gas_background', &
         env%gas_background)
    call inout_read_real(file, 'gas_dilution_rate', env%gas_dilution_rate)
    call spec_read_aero_dist_filename(file, aero_data, bin_grid, &
         'aerosol_emissions', env%aero_emissions)
    call inout_read_real(file, 'aerosol_emission_rate', env%aero_emission_rate)
    call spec_read_aero_dist_filename(file, aero_data, bin_grid, &
         'aerosol_background', env%aero_background)
    call inout_read_real(file, 'aerosol_dilution_rate', env%aero_dilution_rate)

  end subroutine spec_read_env

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine env_average(env_vec, env_avg)
    
    ! Computes the average of an array of env.

    use pmc_util
    use pmc_gas_state
    use pmc_aero_dist
    
    type(env_t), intent(in) :: env_vec(:) ! array of env
    type(env_t), intent(out) :: env_avg   ! average of env_vec

    integer :: i_temp, n_temp, i, n

    call average_real(env_vec%temp, env_avg%temp)
    call average_real(env_vec%rel_humid, env_avg%rel_humid)
    call average_real(env_vec%pressure, env_avg%pressure)
    call average_real(env_vec%air_den, env_avg%air_den)
    call average_real(env_vec%longitude, env_avg%longitude)
    call average_real(env_vec%latitude, env_avg%latitude)
    call average_real(env_vec%altitude, env_avg%altitude)
    call average_real(env_vec%start_time, env_avg%start_time)
    call average_integer(env_vec%start_day, env_avg%start_day)
    call average_integer(env_vec%n_temp, env_avg%n_temp)
    n_temp = env_avg%n_temp
    call env_temp_alloc(env_avg, n_temp)
    n = size(env_vec)
    do i_temp = 1,n_temp
       call average_real((/(env_vec(i)%temp_time(i_temp),i=1,n)/), &
            env_avg%temp_time(i_temp))
       call average_real((/(env_vec(i)%temp_set(i_temp),i=1,n)/), &
            env_avg%temp_set(i_temp))
    end do
    call gas_state_average(env_vec%gas_emissions, env_avg%gas_emissions)
    call average_real(env_vec%gas_emission_rate, env_avg%gas_emission_rate)
    call gas_state_average(env_vec%gas_background, env_avg%gas_background)
    call average_real(env_vec%gas_dilution_rate, env_avg%gas_dilution_rate)
    call aero_dist_average(env_vec%aero_emissions, env_avg%aero_emissions)
    call average_real(env_vec%aero_emission_rate, env_avg%aero_emission_rate)
    call aero_dist_average(env_vec%aero_background, env_avg%aero_background)
    call average_real(env_vec%aero_dilution_rate, env_avg%aero_dilution_rate)
    
  end subroutine env_average
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function pmc_mpi_pack_env_size(val)

    ! Determines the number of bytes required to pack the given value.

    use pmc_mpi

    type(env_t), intent(in) :: val ! value to pack

    pmc_mpi_pack_env_size = &
         pmc_mpi_pack_real_size(val%temp) &
         + pmc_mpi_pack_real_size(val%rel_humid) &
         + pmc_mpi_pack_real_size(val%pressure) &
         + pmc_mpi_pack_real_size(val%air_den) &
         + pmc_mpi_pack_real_size(val%longitude) &
         + pmc_mpi_pack_real_size(val%latitude) &
         + pmc_mpi_pack_real_size(val%altitude) &
         + pmc_mpi_pack_real_size(val%start_time) &
         + pmc_mpi_pack_integer_size(val%start_day) &
         + pmc_mpi_pack_integer_size(val%n_temp) &
         + pmc_mpi_pack_real_array_size(val%temp_time) &
         + pmc_mpi_pack_real_array_size(val%temp_set) &
         + pmc_mpi_pack_gas_state_size(val%gas_emissions) &
         + pmc_mpi_pack_real_size(val%gas_emission_rate) &
         + pmc_mpi_pack_gas_state_size(val%gas_background) &
         + pmc_mpi_pack_real_size(val%gas_dilution_rate) &
         + pmc_mpi_pack_aero_dist_size(val%aero_emissions) &
         + pmc_mpi_pack_real_size(val%aero_emission_rate) &
         + pmc_mpi_pack_aero_dist_size(val%aero_background) &
         + pmc_mpi_pack_real_size(val%aero_dilution_rate)

  end function pmc_mpi_pack_env_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine pmc_mpi_pack_env(buffer, position, val)

    ! Packs the given value into the buffer, advancing position.

#ifdef PMC_USE_MPI
    use mpi
    use pmc_mpi
    use pmc_util
#endif

    character, intent(inout) :: buffer(:) ! memory buffer
    integer, intent(inout) :: position  ! current buffer position
    type(env_t), intent(in) :: val ! value to pack

#ifdef PMC_USE_MPI
    integer :: prev_position

    prev_position = position
    call pmc_mpi_pack_real(buffer, position, val%temp)
    call pmc_mpi_pack_real(buffer, position, val%rel_humid)
    call pmc_mpi_pack_real(buffer, position, val%pressure)
    call pmc_mpi_pack_real(buffer, position, val%air_den)
    call pmc_mpi_pack_real(buffer, position, val%longitude)
    call pmc_mpi_pack_real(buffer, position, val%latitude)
    call pmc_mpi_pack_real(buffer, position, val%altitude)
    call pmc_mpi_pack_real(buffer, position, val%start_time)
    call pmc_mpi_pack_integer(buffer, position, val%start_day)
    call pmc_mpi_pack_integer(buffer, position, val%n_temp)
    call pmc_mpi_pack_real_array(buffer, position, val%temp_time)
    call pmc_mpi_pack_real_array(buffer, position, val%temp_set)
    call pmc_mpi_pack_gas_state(buffer, position, val%gas_emissions)
    call pmc_mpi_pack_real(buffer, position, val%gas_emission_rate)
    call pmc_mpi_pack_gas_state(buffer, position, val%gas_background)
    call pmc_mpi_pack_real(buffer, position, val%gas_dilution_rate)
    call pmc_mpi_pack_aero_dist(buffer, position, val%aero_emissions)
    call pmc_mpi_pack_real(buffer, position, val%aero_emission_rate)
    call pmc_mpi_pack_aero_dist(buffer, position, val%aero_background)
    call pmc_mpi_pack_real(buffer, position, val%aero_dilution_rate)
    call assert(position - prev_position == pmc_mpi_pack_env_size(val))
#endif

  end subroutine pmc_mpi_pack_env

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine pmc_mpi_unpack_env(buffer, position, val)

    ! Unpacks the given value from the buffer, advancing position.

#ifdef PMC_USE_MPI
    use mpi
    use pmc_mpi
    use pmc_util
#endif

    character, intent(inout) :: buffer(:) ! memory buffer
    integer, intent(inout) :: position  ! current buffer position
    type(env_t), intent(out) :: val ! value to pack

#ifdef PMC_USE_MPI
    integer :: prev_position

    prev_position = position
    call pmc_mpi_unpack_real(buffer, position, val%temp)
    call pmc_mpi_unpack_real(buffer, position, val%rel_humid)
    call pmc_mpi_unpack_real(buffer, position, val%pressure)
    call pmc_mpi_unpack_real(buffer, position, val%air_den)
    call pmc_mpi_unpack_real(buffer, position, val%longitude)
    call pmc_mpi_unpack_real(buffer, position, val%latitude)
    call pmc_mpi_unpack_real(buffer, position, val%altitude)
    call pmc_mpi_unpack_real(buffer, position, val%start_time)
    call pmc_mpi_unpack_integer(buffer, position, val%start_day)
    call pmc_mpi_unpack_integer(buffer, position, val%n_temp)
    call pmc_mpi_unpack_real_array(buffer, position, val%temp_time)
    call pmc_mpi_unpack_real_array(buffer, position, val%temp_set)
    call pmc_mpi_unpack_gas_state(buffer, position, val%gas_emissions)
    call pmc_mpi_unpack_real(buffer, position, val%gas_emission_rate)
    call pmc_mpi_unpack_gas_state(buffer, position, val%gas_background)
    call pmc_mpi_unpack_real(buffer, position, val%gas_dilution_rate)
    call pmc_mpi_unpack_aero_dist(buffer, position, val%aero_emissions)
    call pmc_mpi_unpack_real(buffer, position, val%aero_emission_rate)
    call pmc_mpi_unpack_aero_dist(buffer, position, val%aero_background)
    call pmc_mpi_unpack_real(buffer, position, val%aero_dilution_rate)
    call assert(position - prev_position == pmc_mpi_pack_env_size(val))
#endif

  end subroutine pmc_mpi_unpack_env

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine pmc_mpi_reduce_average_env(val, val_avg)

    ! Computes the average of val across all processes, storing the
    ! result in val_avg on the root process.

    use pmc_mpi

    type(env_t), intent(in) :: val ! value to average
    type(env_t), intent(out) :: val_avg ! result

    call pmc_mpi_reduce_average_real(val%temp, val_avg%temp)
    call pmc_mpi_reduce_average_real(val%rel_humid, val_avg%rel_humid)
    call pmc_mpi_reduce_average_real(val%pressure, val_avg%pressure)
    call pmc_mpi_reduce_average_real(val%air_den, val_avg%air_den)

  end subroutine pmc_mpi_reduce_average_env

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
end module pmc_env
