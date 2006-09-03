
module mod_environ

  type environ
     real*8 :: T    ! temperature (K)
     real*8 :: RH   ! relative humidity (1)
     real*8 :: V_comp ! computational volume (m^3)
  end type environ

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine add_water_volume(env, dv)

    ! Adds the given water volume to the water vapor and updates all
    ! enironment quantities.

    type(environ), intent(inout) :: env ! environment state to update
    real*8, intent(in) :: dv            ! volume of water added (m^3)

    ! FIXME: change RH according to dv
    
  end subroutine add_water_volume

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real*8 function sat_vapor_pressure(env) ! Pa

    use mod_constants

    type(environ), intent(in) :: env ! environment state

    sat_vapor_pressure = p00 * 10d0**(7.45d0 * (env%T - T0) &
         / (env%T - 38d0)) ! Pa

  end function sat_vapor_pressure

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module mod_environ
