! -*- mode: f90; -*-
! Copyright (C) 2005,2006 Nicole Riemer and Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.
!
! Golovin coagulation kernel.

module mod_kernel_golovin
contains
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine kernel_golovin(v1, v2, env, k)

    use mod_environ
    
    real*8, intent(in) :: v1  !  volume of first particle
    real*8, intent(in) :: v2  !  volume of second particle
    real*8, intent(out) :: k  !  coagulation kernel
    
    real*8 beta_1
    parameter (beta_1 = 1000d0)
    type(environ), intent(in) :: env  ! environment state
    
    k = beta_1 * (v1 + v2)
    
  end subroutine kernel_golovin
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine soln_golovin_exp(n_bin, bin_v, bin_r, &
       bin_g, bin_n, dlnr, &
       time, N_0, V_0, rho_p, V_comp, env)

    use mod_environ
    
    integer, intent(in) :: n_bin        !  number of bins
    real*8, intent(in) :: bin_v(n_bin)  !  volume of particles in bins
    real*8, intent(in) :: bin_r(n_bin)  !  radius of particles in bins
    real*8, intent(out) :: bin_g(n_bin)  !  volume in bins
    integer, intent(out) :: bin_n(n_bin) !  number in bins
    real*8, intent(in) :: dlnr          !  bin scale factor
    
    real*8, intent(in) :: time          !  cubin_rent time
    real*8, intent(in) :: N_0           !  particle number concentration (#/m^3)
    real*8, intent(in) :: V_0           ! 
    real*8, intent(in) :: rho_p         !  particle density (kg/m^3)
    real*8, intent(in) :: V_comp        !  computational volume

    type(environ), intent(in) :: env  ! environment state
    
    real*8 beta_1, tau, T, rat_v, nn, b, x
    integer k
    
    real*8 pi
    parameter (pi = 3.14159265358979323846d0)
    
    call kernel_golovin(1d0, 0d0, env, beta_1)
    
    if (time .eq. 0d0) then
       do k = 1,n_bin
          bin_n(k) = int(pi/2d0 * (2d0*bin_r(k))**3 * N_0/V_0 &
               * exp(-(bin_v(k)/V_0)))
       enddo
    else
       tau = N_0 * V_0 * beta_1 * time
       T = 1d0 - exp(-tau)
       do k = 1,n_bin
          rat_v = bin_v(k) / V_0
          x = 2d0 * rat_v * sqrt(T)
          if (x .lt. 500d0) then
             call bessi1(x, b)
          else
             b = 0d0
          endif
          nn = N_0/bin_v(k) * (1d0 - T) / sqrt(T) &
               * exp(-((1d0 + T) * rat_v)) * b
          bin_n(k) = int(pi/2d0 * (2d0*bin_r(k))**3 * nn)
       enddo
    endif
    
    do k = 1,n_bin
       bin_g(k) = pi/6d0 * rho_p * (2d0*bin_r(k))**3 * dble(bin_n(k))
    enddo
    
    do k = 1,n_bin
       bin_g(k) = bin_g(k) * dlnr * V_comp
       bin_n(k) = int(dble(bin_n(k)) * dlnr * V_comp)
    enddo
    
  end subroutine soln_golovin_exp
  
  ! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  
  subroutine bessi1(x, r)
    ! bessel function
    
    real*8, intent(in) :: x   !  function argument
    real*8, intent(out) :: r   !  function value
    
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
    endif
    
  end subroutine bessi1
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
end module mod_kernel_golovin
