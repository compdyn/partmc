! Copyright (C) 2011-2012 Jian Tian
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details. 

module pmc_newton_mod_Rme2Rmec

  use pmc_constants
  
  public :: f_Rmec, df_Rmec
  
  contains
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> JT Convert fractal geometric radius (m) to volume (m^3).
  real(kind=dp) function r2vol(r)

    !> Radius (m).
    real(kind=dp), intent(in) :: r
    
    r2vol = 4d0*const%pi*const%R0**3d0*(r/const%R0)**const%df/3d0/const%f
    
  end function r2vol

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> JT Slip correction function from continuum to free molecular regime
  real(kind=dp) function C_Reff(r)
  
    !> radius (m)
    real(kind=dp), intent(in) :: r
 
    C_Reff = 1d0 + const%A_slip*const%l/r + const%Q*const%l/r*exp(-const%b*r/const%l)

  end function C_Reff
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function f_Rmec(x,r) result (y)
    real(kind=dp), intent(in) :: x
    real(kind=dp), intent(in) :: r
    real(kind=dp) :: y

    real(kind=dp) :: C_Rme, phi, h_KR
    C_Rme = C_Reff(r)
    h_KR = -0.06483d0*const%df**2 + 0.6353d0*const%df - 0.4898d0
    phi = const%R0**(2d0-const%df*const%v_gamma)/(const%f**const%v_gamma*  &
            h_KR**(const%df*const%v_gamma))

    y = C_Rme*x - const%A_slip*const%l*r/phi*x**(1d0-const%df*const%v_gamma) &
            - const%Q*const%l*r/phi*x**(1d0-const%df*const%v_gamma) *        &
            exp(-const%b*phi/const%l*x**(const%df*const%v_gamma-1d0))-r 

  end function f_Rmec

  function df_Rmec(x,r) result (y)
    real(kind=dp), intent(in) :: x
    real(kind=dp), intent(in) :: r
    real(kind=dp) :: y

    real(kind=dp) :: C_Rme, phi, h_KR
    C_Rme = C_Reff(r)
    h_KR = -0.06483d0*const%df**2 + 0.6353d0*const%df - 0.4898d0
    phi = const%R0**(2d0-const%df*const%v_gamma)/(const%f**const%v_gamma*  &
            h_KR**(const%df*const%v_gamma))

    y = C_Rme - const%l*r/phi*(1d0-const%df*const%v_gamma)*x**(-const%df*  &
            const%v_gamma)*(const%A_slip+const%Q*exp(-const%b*phi/const%l  &
            *x**(const%df*const%v_gamma-1d0))) - const%Q*const%b*r*        &
            (1d0-const%df*const%v_gamma)/x*exp(-const%b*phi/const%l        &
            *x**(const%df*const%v_gamma-1d0))
 
  end function df_Rmec

!  subroutine newton_Rme2Rmec(f_Rmec,df_Rmec,x,r)
!    
!    interface
!      elemental function f_Rmec(x,r) result (y)
!        integer, parameter :: dp = kind(0.d0) 
!        real(kind=dp), intent(in) :: x 
!        real(kind=dp), intent(in) :: r
!        real(kind=dp) :: y
!      end function f_Rmec
!    end interface
!
!    interface
!      elemental function df_Rmec(x,r) result (y)
!        integer, parameter :: dp = kind(0.d0)
!        real(kind=dp), intent(in) :: x
!        real(kind=dp), intent(in) :: r
!        real(kind=dp) :: y
!      end function df_Rmec
!    end interface
!
!    real(kind=dp), intent(inout) :: x
!    real(kind=dp), intent(in) :: r
!    integer :: iter
!    iter = 1
!
!    do
!      x = x - f_Rmec(x,r)/df_Rmec(x,r)
!      if (iter>MAX_ITERATIONS) then
!         exit
!      end if 
!      iter=iter+1
!    end do 
!  
!  end subroutine newton_Rme2Rmec

end module pmc_newton_mod_Rme2Rmec
