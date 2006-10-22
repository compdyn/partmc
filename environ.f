! -*- mode: f90; -*-
! Copyright (C) 2005,2006 Nicole Riemer and Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.
! -*- mode: f90;-*-

module mod_environ

  type environ
     real*8 :: T    ! temperature (K)
     real*8 :: RH   ! relative humidity (1)
     real*8 :: V_comp ! computational volume (m^3)
     real*8 :: p    ! ambient pressure (Pa)
     real*8 :: dTdt ! change in temperature due to updraft/subsidence (K s^{-1})
  end type environ

contains

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
  
  subroutine change_temp(env, dt)
    
    type(environ), intent(inout) :: env ! environment state to update
    real*8, intent(in) :: dt            ! time step (s)
    
    real*8 pmv      ! ambient water vapor pressure (Pa)
    
    pmv = sat_vapor_pressure(env) * env%RH
    env%T = env%T + env%dTdt * dt
    env%RH = pmv / sat_vapor_pressure(env)
    
  end subroutine change_temp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real*8 function sat_vapor_pressure(env) ! Pa

    use mod_constants

    type(environ), intent(in) :: env ! environment state

    sat_vapor_pressure = const%p00 * 10d0**(7.45d0 * (env%T - const%T0) &
         / (env%T - 38d0)) ! Pa

  end function sat_vapor_pressure

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module mod_environ
