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

  use mod_array
  use mod_gas
  use mod_init_dist
  
  type environ
     real*8 :: T                        ! temperature (K)
     real*8 :: RH                       ! relative humidity (1)
     real*8 :: V_comp                   ! computational volume (m^3)
     real*8 :: p                        ! ambient pressure (Pa)
     real*8 :: rho_a                    ! air density (kg m^{-3})
     real*8 :: longitude                ! longitude (degrees)
     real*8 :: latitude                 ! latitude (degrees)
     real*8 :: altitude                 ! altitude (m)
     real*8 :: start_time               ! start time (s since 00:00 UTC)
     integer :: start_day               ! start day of year (UTC)
     integer :: n_temps                 ! number of temperature set-points
     real*8, dimension(:), pointer :: temp_times ! times at temp set-points
     real*8, dimension(:), pointer :: temps      ! temps at temp set-points
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

  subroutine allocate_environ_temps(env, n_temps)

    ! Allocate storage for a given number of temperature set points.

    type(environ), intent(inout) :: env ! environment
    integer, intent(in) :: n_temps      ! number of temperature set-points

    env%n_temps = n_temps
    allocate(env%temp_times(n_temps))
    allocate(env%temps(n_temps))

  end subroutine allocate_environ_temps

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine change_water_volume(env, mat, dv)
    
    ! Adds the given water volume to the water vapor and updates all
    ! environment quantities.
    
    use mod_constants
    use mod_material
    
    type(environ), intent(inout) :: env ! environment state to update
    type(material), intent(in)   :: mat ! material constants
    real*8, intent(in) :: dv            ! volume of water added (m^3)
    
    real*8 pmv     ! ambient water vapor pressure (Pa)
    real*8 mv      ! ambient water vapor density (kg m^{-3})
                   ! pmv and mv are related by the factor M_w/(R*T)
    real*8 dmv     ! change of water density (kg m^{-3})
    
    dmv = dv * mat%rho(mat%i_water) / env%V_comp
    pmv = sat_vapor_pressure(env) * env%RH
    mv = mat%M_w(mat%i_water)/(const%R*env%T) * pmv
    mv = mv - dmv    
    env%RH = const%R * env%T / mat%M_w(mat%i_water) * mv &
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

    use mod_gas

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
  
end module mod_environ
