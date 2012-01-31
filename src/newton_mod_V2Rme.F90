!> Jian Tian 
!> Convert volume to mobility equivalent radius

module pmc_newton_mod_V2Rme

  use pmc_constants
  
  public :: f_Rme, df_Rme
  
  contains

 !> JT Convert volume (m^3) to geometric radius (m) for fractal particles
  real(kind=dp) function v2rad(v)
  
    !> Volume (m^3)
    real(kind=dp), intent(in) :: v

    v2rad = (3d0*v*const%f/4d0/const%pi/const%R0**3d0)**(1d0/const%df)*const%R0    
    
  end function v2rad

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> JT Convert volume (m^3) to continuum regime mobility equivalent radius
  real(kind=dp) function vol2R_me_c(v)

    !> Volume (m^3)
    real(kind=dp), intent(in) :: v
    !> Geometric radius (m)
    real(kind=dp) :: Rgeo
    
    Rgeo = v2rad(v)
    vol2R_me_c = (-0.06483d0*const%df**2d0+0.6353d0*const%df-0.4898d0)*Rgeo

  end function vol2R_me_c

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
!> JT Convert volume to accessible particle surface
  real(kind=dp) function vol2S_acc(v)

    !> Volume (m^3)
    real(kind=dp), intent(in) :: v

    real(kind=dp) :: N, ds
    
    !> Surface fractal dimension
    if (const%df .le. 2d0) then
       ds = 3d0
    elseif ((const%df .gt. 2d0) .and. (const%df .le. 3d0)) then
       ds = 6d0 / const%df
    end if

    N = vol2N(v)
    vol2S_acc = 4d0*const%pi*const%R0**2*N**(ds/3d0)*     &
          ((ds-2d0)*(const%z/N)**(1d0-const%v_gamma)-ds+3d0)

  end function vol2S_acc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> JT Convert volume to effective radius
  real(kind=dp) function vol2R_eff(v)

    !> Volume (m^3)
    real(kind=dp), intent(in) :: v
    
    real(kind=dp) :: R_me_c, S_acc
    R_me_c = vol2R_me_c(v)
    S_acc = vol2S_acc(v)
    vol2R_eff = S_acc / 4d0 / const%pi /R_me_c

  end function vol2R_eff

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> JT Slip correction function from continuum to free molecular regime
  real(kind=dp) function C_r(r)
  
    !> radius (m)
    real(kind=dp), intent(in) :: r
 
    C_r = 1d0 + const%A_slip*const%l/r + const%Q*const%l/r*exp(-const%b*r/const%l)

  end function C_r
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
   function f_Rme(x,v) result (y)
      real(kind=dp), intent(in) :: x
      real(kind=dp), intent(in) :: v
      real(kind=dp) :: y

      real(kind=dp) :: R_me_c, R_eff, C_Reff
      R_me_c = vol2R_me_c(v)
      R_eff = vol2R_eff(v)
      C_Reff = C_r(R_eff)

      y=C_Reff*x**2-R_me_c*x-R_me_c*const%Q*const%l*exp(-const%b*x/    &
           const%l)-R_me_c*const%A_slip*const%l
    end function f_Rme

    function df_Rme(x,v) result (y)
      real(kind=dp), intent(in) :: x
      real(kind=dp), intent(in) :: v
      real(kind=dp) :: y

      real(kind=dp) :: R_me_c, R_eff, C_Reff
      R_me_c = vol2R_me_c(v)
      R_eff = vol2R_eff(v)
      C_Reff = C_r(R_eff) 

      y = 2d0*C_Reff*x - R_me_c + R_me_c*const%Q*const%b*exp(-const%b*x/const%l)
    end function df_Rme
    
  !subroutine newton_V2Rme(f_Rme,df_Rme,x,v)
    
  !  interface
  !    function f_Rme(x,v) result (y)
  !      integer, parameter :: dp = kind(0.d0) 
  !      real(kind=dp), intent(in) :: x 
  !      real(kind=dp), intent(in) :: v
  !      real(kind=dp) :: y
  !    end function f_Rme
  !  end interface

  !  interface
  !    function df_Rme(x,v) result (y)
  !      integer, parameter :: dp = kind(0.d0)
  !      real(kind=dp), intent(in) :: x
  !      real(kind=dp), intent(in) :: v
  !      real(kind=dp) :: y
  !    end function df_Rme
  !  end interface

  !  real(kind=dp), intent(inout) :: x
  !  real(kind=dp), intent(in) :: v
  !  integer :: iter
  !  iter = 1

  !  do
  !    x = x - f_Rme(x,v)/df_Rme(x,v)
  !    if (iter>MAX_ITERATIONS) then
  !       exit
  !    end if 
  !    iter=iter+1
  !  end do 
  
  !end subroutine newton_V2Rme
  
end module pmc_newton_mod_V2Rme
