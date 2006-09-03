
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

  subroutine add_water_volume(env, dv)

  use mod_constants

    ! Adds the given water volume to the water vapor and updates all
    ! environment quantities.

    type(environ), intent(inout) :: env ! environment state to update
    real*8, intent(in) :: dv            ! volume of water added (m^3)

    real*8 pmv                          ! ambient water vapor pressure (Pa)
    real*8 mv                           ! ambient water vapor density (kg m^{-3})
                                        ! pmv and mv are related by the factor M_w/(R*T)
    real*8 dmv                          ! change of water density (kg m^{-3})

    dmv = dv * rho(i_water) / V_comp    

    pmv = sat_vapor_pressure * RH

    mv = M_w(i_water)/(R*T) * pmv

    mv = mv - dmv    

    RH = R * T / M_w(i_water) * mv / sat_vapor_pressure

  end subroutine add_water_volume

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine change_temp(env, dt)

  type(environ), intent(inout) :: env ! environment state to update
  real*8 intent(in) :: dt
  real*8 pmv                          ! ambient water vapor pressure (Pa)

  pmv = sat_vapor_pressure * RH

  T = T + dTdt * dt

  RH = pmv / sat_vapor_pressure

  end subroutine change_temp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real*8 function sat_vapor_pressure(env) ! Pa

    use mod_constants

    type(environ), intent(in) :: env ! environment state

    sat_vapor_pressure = p00 * 10d0**(7.45d0 * (env%T - T0) &
         / (env%T - 38d0)) ! Pa

  end function sat_vapor_pressure

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module mod_environ
