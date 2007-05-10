! -*- mode: f90; -*-
! Copyright (C) 2005-2007 Nicole Riemer and Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.
!
! Golovin coagulation kernel.

module mod_kernel_golovin
contains
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine kernel_golovin(v1, v2, env, k)

    ! Golovin coagulation kernel.

    use mod_environ
    
    real*8, intent(in) :: v1            ! volume of first particle
    real*8, intent(in) :: v2            ! volume of second particle
    type(environ), intent(in) :: env    ! environment state
    real*8, intent(out) :: k            ! coagulation kernel
    
    real*8, parameter :: beta_1 = 1000d0
    
    k = beta_1 * (v1 + v2)
    
  end subroutine kernel_golovin
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine soln_golovin_exp(n_bin, bin_v, bin_g_den, bin_n_den, &
       time, num_conc, mean_vol, rho_p, env)

    ! Exact solution with the Golovin coagulation kernel and
    ! exponential initial condition.

    use mod_environ
    use mod_util
    use mod_constants
    
    integer, intent(in) :: n_bin        ! number of bins
    real*8, intent(in) :: bin_v(n_bin)  ! volume of particles in bins
    real*8, intent(out) :: bin_g_den(n_bin) ! volume density in bins
    real*8, intent(out) :: bin_n_den(n_bin) ! number density in bins
    
    real*8, intent(in) :: time          ! current time
    real*8, intent(in) :: num_conc      ! particle number concentration (#/m^3)
    real*8, intent(in) :: mean_vol      ! mean init volume (m^3)
    real*8, intent(in) :: rho_p         ! particle density (kg/m^3)
    type(environ), intent(in) :: env    ! environment state
    
    real*8 beta_1, tau, T, rat_v, nn, b, x
    integer k
    
    call kernel_golovin(1d0, 0d0, env, beta_1)
    
    if (time .eq. 0d0) then
       do k = 1,n_bin
          bin_n_den(k) = const%pi/2d0 * (2d0*vol2rad(bin_v(k)))**3 &
               * num_conc / mean_vol * exp(-(bin_v(k)/mean_vol))
       end do
    else
       tau = num_conc * mean_vol * beta_1 * time
       T = 1d0 - exp(-tau)
       do k = 1,n_bin
          rat_v = bin_v(k) / mean_vol
          x = 2d0 * rat_v * sqrt(T)
          if (x .lt. 500d0) then
             call bessi1(x, b)
          else
             b = 0d0
          end if
          nn = num_conc / bin_v(k) * (1d0 - T) / sqrt(T) &
               * exp(-((1d0 + T) * rat_v)) * b
          bin_n_den(k) = const%pi/2d0 * (2d0*vol2rad(bin_v(k)))**3 * nn
       end do
    end if

    do k = 1,n_bin
       bin_g_den(k) = const%pi/6d0 * (2d0*vol2rad(bin_v(k)))**3 * bin_n_den(k)
    end do
    
  end subroutine soln_golovin_exp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine bessi1(x, r)

    ! $I_1(x)$
    ! Modified Bessel function of the first kind

    ! Bessel function.
    ! This looks like it was taken from Numerical Recipes.
    ! FIXME: Where did this come from? What license does it have?
    
    real*8, intent(in) :: x             ! function argument
    real*8, intent(out) :: r            ! function value
    
    real*8 ax
    real*8 p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,q8,q9,y
    data p1,p2,p3,p4,p5,p6,p7/0.5d0,0.87890594d0,0.51498869d0, &
         0.15084934d0,0.2658733d-1,0.301532d-2,0.32411d-3/
    data q1,q2,q3,q4,q5,q6,q7,q8,q9/0.39894228d0,-0.3988024d-1, &
         -0.362018d-2,0.163801d-2,-0.1031555d-1,0.2282967d-1, &
         -0.2895312d-1,0.1787654d-1,-0.420059d-2/
    
    if(abs(x) .lt. 3.75d0) then
       y = (x / 3.75d0)**2
       r = x*(p1+y*(p2+y*(p3+y*(p4+y*(p5+y*(p6+y*p7))))))
    else
       ax = abs(x)
       y = 3.75d0 / ax
       r = (exp(ax)/sqrt(ax))*(q1+y*(q2+y*(q3+y*(q4+ &
            y*(q5+y*(q6+y*(q7+y*(q8+y*q9))))))))
       if (x .lt. 0d0) r = -r
    end if
    
  end subroutine bessi1
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
end module mod_kernel_golovin
