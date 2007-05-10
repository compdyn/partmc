! Copyright (C) 2005-2007 Nicole Riemer and Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.
!
! Brownian coagulation kernel.
! See Seinfeld, Atmospheric chemistry and physics of air pollution,
! page 394 (equation 10.18)
! This expression is based on the assumption that the continuum regime applies.
    
module mod_kernel_brown
contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine kernel_brown(v1, v2, env, k)

    ! Compute the Brownian coagulation kernel.

    use mod_environ
    use mod_constants

    real*8, intent(in) :: v1            ! volume of first particle (m^3)
    real*8, intent(in) :: v2            ! volume of second particle (m^3)
    type(environ), intent(in) :: env    ! environment state
    real*8, intent(out) :: k            ! kernel k(a,b) (m^3/s)

    real*8 c_1, a_third, b_third

    c_1 = 2d0 * const%k_b * env%T / (3.d0 * const%mu)
    a_third = v1**(1.d0/3.d0)
    b_third = v2**(1.d0/3.d0)
    
    k = c_1 * (a_third + b_third) *(1d0/a_third + 1d0/b_third) 
    
  end subroutine kernel_brown
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module mod_kernel_brown

