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

module mod_environ

  use mod_gas_state
  use mod_aero_dist
  
  type environ
     real*8 :: T                        ! temperature (K)
     real*8 :: RH                       ! relative humidity (1)
     real*8 :: p                        ! ambient pressure (Pa)
     real*8 :: rho_a                    ! air density (kg m^{-3})
     real*8 :: longitude                ! longitude (degrees)
     real*8 :: latitude                 ! latitude (degrees)
     real*8 :: altitude                 ! altitude (m)
     real*8 :: start_time               ! start time (s since 00:00 UTC)
     integer :: start_day               ! start day of year (UTC)
     integer :: n_temps                 ! number of temperature set-points
     real*8, pointer :: temp_times(:)   ! times at temp set-points (s)
     real*8, pointer :: temps(:)        ! temps at temp set-points (K)
     type(gas_state_t) :: gas_emissions ! gas emissions
     real*8 :: gas_emission_rate        ! gas emisssion rate (s^{-1})
     type(gas_state_t) :: gas_background ! background gas concentrations
     real*8 :: gas_dilution_rate        ! gas-background dilution rate (s^{-1})
     type(aero_dist_t) :: aero_emissions ! aerosol emissions
     real*8 :: aero_emission_rate       ! aerosol emisssion rate (s^{-1})
     type(aero_dist_t) :: aero_background ! aerosol background
     real*8 :: aero_dilution_rate       ! aero-background dilution rate (s^{-1})
  end type environ
  
contains
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine env_alloc(env)

    ! Allocate an empty environment.

    type(environ), intent(out) :: env ! environment

    env%T = 0d0
    env%RH = 0d0
    env%p = 0d0
    env%rho_a = 0d0
    env%longitude = 0d0
    env%latitude = 0d0
    env%altitude = 0d0
    env%start_time = 0d0
    env%start_day = 0

    call environ_temps_alloc(env, 0)
    call gas_state_alloc(0, env%gas_emissions)
    call gas_state_alloc(0, env%gas_background)
    env%gas_emission_rate = 0d0
    env%gas_dilution_rate = 0d0
    call aero_dist_alloc(0, 0, 0, env%aero_emissions)
    call aero_dist_alloc(0, 0, 0, env%aero_background)
    env%aero_emission_rate = 0d0
    env%aero_dilution_rate = 0d0

  end subroutine env_alloc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine env_free(env)

    ! Free all storage.

    type(environ), intent(out) :: env ! environment

    call environ_temps_free(env)
    call gas_state_free(env%gas_emissions)
    call gas_state_free(env%gas_background)
    call aero_dist_free(env%aero_emissions)
    call aero_dist_free(env%aero_background)

  end subroutine env_free

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine environ_temps_alloc(env, n_temps)

    ! Allocate storage for a given number of temperature set points.

    type(environ), intent(inout) :: env ! environment
    integer, intent(in) :: n_temps      ! number of temperature set-points

    env%n_temps = n_temps
    allocate(env%temp_times(n_temps))
    allocate(env%temps(n_temps))

  end subroutine environ_temps_alloc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine environ_temps_free(env)

    ! Free all storage.

    type(environ), intent(inout) :: env ! environment

    deallocate(env%temp_times)
    deallocate(env%temps)
    
  end subroutine environ_temps_free

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine change_water_volume(env, aero_data, dv)
    
    ! Adds the given water volume to the water vapor and updates all
    ! environment quantities.
    
    use mod_constants
    use mod_aero_data
    
    type(environ), intent(inout) :: env ! environment state to update
    type(aero_data_t), intent(in)   :: aero_data ! aero_data constants
    real*8, intent(in) :: dv            ! conc of water added (m^3/m^3)
    
    real*8 pmv     ! ambient water vapor pressure (Pa)
    real*8 mv      ! ambient water vapor density (kg m^{-3})
                   ! pmv and mv are related by the factor M_w/(R*T)
    real*8 dmv     ! change of water density (kg m^{-3})
    
    dmv = dv * aero_data%rho(aero_data%i_water)
    pmv = sat_vapor_pressure(env) * env%RH
    mv = aero_data%M_w(aero_data%i_water)/(const%R*env%T) * pmv
    mv = mv - dmv    
    env%RH = const%R * env%T / aero_data%M_w(aero_data%i_water) * mv &
         / sat_vapor_pressure(env)
    
  end subroutine change_water_volume
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine init_environ(env, time)
    
    ! Initialize the time-dependent contents of the
    ! environment. Thereafter update_environ() should be used.

    use mod_util

    type(environ), intent(inout) :: env ! environment state to update
    real*8, intent(in) :: time          ! current time (s)

    env%T = interp_1d(env%n_temps, env%temp_times, env%temps, time)
    
  end subroutine init_environ
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine update_environ(env, time)
    
    ! Update time-dependent contents of the environment. init_environ()
    ! should have been called at the start.

    use mod_util

    type(environ), intent(inout) :: env ! environment state to update
    real*8, intent(in) :: time          ! current time (s)
    
    real*8 pmv      ! ambient water vapor pressure (Pa)

    ! update temperature and relative humidity
    pmv = sat_vapor_pressure(env) * env%RH
    env%T = interp_1d(env%n_temps, env%temp_times, env%temps, time)
    env%RH = pmv / sat_vapor_pressure(env)
    
  end subroutine update_environ
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  real*8 function sat_vapor_pressure(env) ! Pa

    ! Computes the current saturation vapor pressure.
    
    use mod_constants
    
    type(environ), intent(in) :: env    ! environment state
    
    sat_vapor_pressure = const%p00 * 10d0**(7.45d0 * (env%T - const%T0) &
         / (env%T - 38d0)) ! Pa
    
  end function sat_vapor_pressure
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine environ_update_gas_state(env, delta_t, gas_data, gas_state)

    ! Do emissions and background dilution from the environment.

    use mod_gas_data
    use mod_gas_state

    type(environ), intent(in) :: env    ! current environment
    real*8, intent(in) :: delta_t       ! time increment to update over
    type(gas_data_t), intent(in) :: gas_data ! gas data values
    type(gas_state_t), intent(inout) :: gas_state ! gas state to update

    type(gas_state_t) :: emission, dilution

    call gas_state_alloc(gas_data%n_spec, emission)
    call gas_state_alloc(gas_data%n_spec, dilution)

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

  end subroutine environ_update_gas_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine environ_update_aero_state(env, delta_t, bin_grid, &
       aero_data, aero_state, aero_binned)

    ! Do emissions and background dilution from the environment for a
    ! particle aerosol distribution.

    use mod_bin_grid
    use mod_aero_data
    use mod_aero_state
    use mod_aero_binned

    type(environ), intent(in) :: env    ! current environment
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
    call aero_binned_alloc(bin_grid%n_bin, aero_data%n_spec, aero_binned_delta)
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

  end subroutine environ_update_aero_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine environ_update_aero_binned(env, delta_t, bin_grid, &
       aero_data, aero_binned)

    ! Do emissions and background dilution from the environment for a
    ! binned aerosol distribution.

    use mod_bin_grid
    use mod_aero_data
    use mod_aero_binned

    type(environ), intent(in) :: env    ! current environment
    real*8, intent(in) :: delta_t       ! time increment to update over
    type(bin_grid_t), intent(in) :: bin_grid ! bin grid
    type(aero_data_t), intent(in) :: aero_data ! aero data values
    type(aero_binned_t), intent(inout) :: aero_binned ! aero binned to update

    type(aero_binned_t) :: emission, dilution

    call aero_binned_alloc(bin_grid%n_bin, aero_data%n_spec, emission)
    call aero_binned_alloc(bin_grid%n_bin, aero_data%n_spec, dilution)

    ! emission = delta_t * gas_emission_rate * gas_emissions
    call aero_dist_add_to_binned(bin_grid, env%aero_emissions, emission)
    call aero_binned_scale(emission, delta_t * env%aero_emission_rate)

    ! dilution = delta_t * gas_dilution_rate * (gas_background - aero_binned)
    call aero_dist_add_to_binned(bin_grid, env%aero_background, dilution)
    call aero_binned_sub(dilution, aero_binned)
    call aero_binned_scale(dilution, delta_t * env%aero_dilution_rate)

    call aero_binned_add(aero_binned, emission)
    call aero_binned_add(aero_binned, dilution)

    call aero_binned_free(emission)
    call aero_binned_free(dilution)

  end subroutine environ_update_aero_binned

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine inout_write_env(file, env)
    
    ! Write full state.
    
    use mod_inout
    
    type(inout_file_t), intent(inout) :: file ! file to write to
    type(environ), intent(in) :: env    ! environment to write
    
    call inout_write_real(file, "temp(K)", env%T)
    call inout_write_real(file, "rel_humidity(1)", env%RH)
    call inout_write_real(file, "pressure(Pa)", env%p)
    call inout_write_real(file, "air_density(kg/m^3)", env%rho_a)
    call inout_write_real(file, "longitude(deg)", env%longitude)
    call inout_write_real(file, "latitude(deg)", env%latitude)
    call inout_write_real(file, "altitude(m)", env%altitude)
    call inout_write_real(file, "start_time(s)", env%start_time)
    call inout_write_integer(file, "start_day(days)", env%start_day)
    call inout_write_integer(file, "num_temps", env%n_temps)
    call inout_write_real_array(file, "temp_times(s)", env%temp_times)
    call inout_write_real_array(file, "temps(K)", env%temps)
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
    
    use mod_inout
    
    type(inout_file_t), intent(inout) :: file ! file to read from
    type(environ), intent(out) :: env    ! environment to read
    
    call inout_read_real(file, "temp(K)", env%T)
    call inout_read_real(file, "rel_humidity(1)", env%RH)
    call inout_read_real(file, "pressure(Pa)", env%p)
    call inout_read_real(file, "air_density(kg/m^3)", env%rho_a)
    call inout_read_real(file, "longitude(deg)", env%longitude)
    call inout_read_real(file, "latitude(deg)", env%latitude)
    call inout_read_real(file, "altitude(m)", env%altitude)
    call inout_read_real(file, "start_time(s)", env%start_time)
    call inout_read_integer(file, "start_day(days)", env%start_day)
    call inout_read_integer(file, "num_temps", env%n_temps)
    call inout_read_real_array(file, "temp_times(s)", env%temp_times)
    call inout_read_real_array(file, "temps(K)", env%temps)
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

  subroutine spec_read_environ(file, bin_grid, gas_data, aero_data, env)

    ! Read environment specification from a inout file.

    use mod_bin_grid
    use mod_inout
    use mod_aero_data
    use mod_gas_data

    type(inout_file_t), intent(inout) :: file ! inout file
    type(bin_grid_t), intent(in) :: bin_grid ! bin grid
    type(gas_data_t), intent(in) :: gas_data ! gas data values
    type(aero_data_t), intent(in) :: aero_data ! aerosol data
    type(environ), intent(out) :: env   ! environment data

    integer :: n_temps
    character(len=MAX_CHAR_LEN) :: read_name
    type(inout_file_t) :: read_file
    character(len=MAX_CHAR_LEN), pointer :: times_name(:), temps_name(:)
    real*8, pointer :: times_data(:,:), temps_data(:,:)

    ! read the tempurature data from the specified file
    call inout_read_string(file, 'temp_profile', read_name)
    call inout_open_read(read_name, read_file)
    call inout_read_real_named_array(read_file, 1, times_name, times_data)
    call inout_check_name(read_file, "time", times_name(1))
    ! FIXME: add a min_lines arg to inout_read_real_named_array to ensure that
    ! really got one line here
    call inout_read_real_named_array(read_file, 1, temps_name, temps_data)
    call inout_check_name(read_file, "temp", temps_name(1))
    call inout_close(read_file)

    ! check the data size
    n_temps = size(temps_data, 2)
    if (n_temps < 1) then
       write(0,*) 'ERROR: file ', trim(read_name), &
            ' must contain at least one line of data'
       call exit(1)
    end if
    if (size(times_data, 2) /= size(temps_data, 2)) then
       write(0,*) 'ERROR: file ', trim(read_name), &
            ' should contain exactly two lines with equal numbers of values'
       call exit(1)
    end if

    call environ_temps_alloc(env, n_temps)
    env%temp_times = times_data(1,:)
    env%temps = temps_data(1,:)
    call inout_read_real(file, 'RH', env%RH)
    call inout_read_real(file, 'pressure', env%p)
    call inout_read_real(file, 'rho_a', env%rho_a)
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

  end subroutine spec_read_environ

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine average_env(env_vec, env_avg)
    
    ! Computes the average of an array of env.

    use mod_util
    use mod_gas_state
    use mod_aero_dist
    
    type(environ), intent(in) :: env_vec(:) ! array of env
    type(environ), intent(out) :: env_avg   ! average of env_vec

    integer :: i_temp, n_temps, i, n

    call average_real(env_vec%T, env_avg%T)
    call average_real(env_vec%RH, env_avg%RH)
    call average_real(env_vec%p, env_avg%p)
    call average_real(env_vec%rho_a, env_avg%rho_a)
    call average_real(env_vec%longitude, env_avg%longitude)
    call average_real(env_vec%latitude, env_avg%latitude)
    call average_real(env_vec%altitude, env_avg%altitude)
    call average_real(env_vec%start_time, env_avg%start_time)
    call average_integer(env_vec%start_day, env_avg%start_day)
    call average_integer(env_vec%n_temps, env_avg%n_temps)
    n_temps = env_avg%n_temps
    call environ_temps_alloc(env_avg, n_temps)
    n = size(env_vec)
    do i_temp = 1,n_temps
       call average_real((/(env_vec(i)%temp_times(i_temp),i=1,n)/), &
            env_avg%temp_times(i_temp))
       call average_real((/(env_vec(i)%temps(i_temp),i=1,n)/), &
            env_avg%temps(i_temp))
    end do
    call average_gas_state(env_vec%gas_emissions, env_avg%gas_emissions)
    call average_real(env_vec%gas_emission_rate, env_avg%gas_emission_rate)
    call average_gas_state(env_vec%gas_background, env_avg%gas_background)
    call average_real(env_vec%gas_dilution_rate, env_avg%gas_dilution_rate)
    call average_aero_dist(env_vec%aero_emissions, env_avg%aero_emissions)
    call average_real(env_vec%aero_emission_rate, env_avg%aero_emission_rate)
    call average_aero_dist(env_vec%aero_background, env_avg%aero_background)
    call average_real(env_vec%aero_dilution_rate, env_avg%aero_dilution_rate)
    
  end subroutine average_env
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
end module mod_environ
