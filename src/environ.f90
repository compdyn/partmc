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

  use mod_aero_state
  use mod_gas_data
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

  subroutine alloc_env(env)

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

    call allocate_environ_temps(env, 0)
    call allocate_gas_state(0, env%gas_emissions)
    call allocate_gas_state(0, env%gas_background)
    env%gas_emission_rate = 0d0
    env%gas_dilution_rate = 0d0
    call alloc_aero_dist(0, 0, 0, env%aero_emissions)
    call alloc_aero_dist(0, 0, 0, env%aero_background)
    env%aero_emission_rate = 0d0
    env%aero_dilution_rate = 0d0

  end subroutine alloc_env

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine allocate_environ_temps(env, n_temps)

    ! Allocate storage for a given number of temperature set points.

    type(environ), intent(inout) :: env ! environment
    integer, intent(in) :: n_temps      ! number of temperature set-points

    env%n_temps = n_temps
    allocate(env%temp_times(n_temps))
    allocate(env%temps(n_temps))

  end subroutine allocate_environ_temps

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

  subroutine gas_update_from_environ(env, delta_t, gas_data, gas_state)

    ! Do emissions and background dilution from the environment.

    use mod_gas_data
    use mod_gas_state

    type(environ), intent(in) :: env    ! current environment
    real*8, intent(in) :: delta_t       ! time increment to update over
    type(gas_data_t), intent(in) :: gas_data ! gas data values
    type(gas_state_t), intent(inout) :: gas_state ! gas state to update

    integer :: i
    real*8 :: emission, dilution

    do i = 1,gas_data%n_spec
       emission = delta_t * env%gas_emission_rate &
            * env%gas_emissions%conc(i)
       dilution = delta_t * env%gas_dilution_rate &
            * (env%gas_background%conc(i) - gas_state%conc(i))
       gas_state%conc(i) = gas_state%conc(i) + emission + dilution
    end do

  end subroutine gas_update_from_environ

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

  subroutine spec_read_environ(file, env)

    ! Read environment specification from a inout file.

    use mod_inout

    type(inout_file_t), intent(inout) :: file ! inout file
    type(environ), intent(out) :: env   ! environment data

    integer :: n_temps
    character(len=MAX_CHAR_LEN) :: read_name
    type(inout_file_t) :: read_file
    character(len=MAX_CHAR_LEN), pointer :: times_name(:), temps_name(:)
    real*8, pointer :: times_data(:,:), temps_data(:,:)
    integer :: times_data_shape(2), temps_data_shape(2)

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
    times_data_shape = shape(times_data)
    temps_data_shape = shape(temps_data)
    n_temps = temps_data_shape(2)
    if (n_temps < 1) then
       write(0,*) 'ERROR: file ', trim(read_name), &
            ' must contain at least one line of data'
       call exit(1)
    end if
    if (times_data_shape(2) /= temps_data_shape(2)) then
       write(0,*) 'ERROR: file ', trim(read_name), &
            ' should contain exactly two lines with equal numbers of values'
       call exit(1)
    end if

    call allocate_environ_temps(env, n_temps)
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

  end subroutine spec_read_environ

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
end module mod_environ
