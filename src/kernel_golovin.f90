! Copyright (C) 2005-2007 Nicole Riemer and Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.
!
! Golovin coagulation kernel.

module pmc_kernel_golovin

  use pmc_bin_grid
  use pmc_env
  use pmc_util
  use pmc_constants
  use pmc_aero_binned
  use pmc_aero_data
  use pmc_aero_dist
  
contains
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine kernel_golovin(v1, v2, env, k)

    ! Golovin coagulation kernel.

    real*8, intent(in) :: v1            ! volume of first particle
    real*8, intent(in) :: v2            ! volume of second particle
    type(env_t), intent(in) :: env      ! environment state
    real*8, intent(out) :: k            ! coagulation kernel
    
    real*8, parameter :: beta_1 = 1000d0
    
    k = beta_1 * (v1 + v2)
    
  end subroutine kernel_golovin
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine soln_golovin_exp(bin_grid, aero_data, time, num_den, &
       mean_radius, rho_p, aero_dist_init, env, aero_binned)

    ! Exact solution with the Golovin coagulation kernel and
    ! exponential initial condition.

    type(bin_grid_t), intent(in) :: bin_grid ! bin grid
    type(aero_data_t), intent(in) :: aero_data ! aerosol data
    real*8, intent(in) :: time          ! current time
    real*8, intent(in) :: num_den       ! particle number concentration (#/m^3)
    real*8, intent(in) :: mean_radius   ! mean init radius (m)
    real*8, intent(in) :: rho_p         ! particle density (kg/m^3)
    type(aero_dist_t), intent(in) :: aero_dist_init ! initial distribution
    type(env_t), intent(in) :: env      ! environment state
    type(aero_binned_t), intent(out) :: aero_binned ! output state
    
    real*8 :: beta_1, tau, T, rat_v, nn, b, x, mean_vol
    integer :: k
    
    call kernel_golovin(1d0, 0d0, env, beta_1)

    mean_vol = rad2vol(mean_radius)
    if (time .eq. 0d0) then
       do k = 1,bin_grid%n_bin
          aero_binned%num_den(k) = const%pi/2d0 &
               * (2d0*vol2rad(bin_grid%v(k)))**3 * num_den / mean_vol &
               * exp(-(bin_grid%v(k)/mean_vol))
       end do
    else
       tau = num_den * mean_vol * beta_1 * time
       T = 1d0 - exp(-tau)
       do k = 1,bin_grid%n_bin
          rat_v = bin_grid%v(k) / mean_vol
          x = 2d0 * rat_v * sqrt(T)
          if (x .lt. 500d0) then
             call bessi1(x, b)
          else
             b = 0d0
          end if
          nn = num_den / bin_grid%v(k) * (1d0 - T) / sqrt(T) &
               * exp(-((1d0 + T) * rat_v)) * b
          aero_binned%num_den(k) = const%pi/2d0 &
               * (2d0*vol2rad(bin_grid%v(k)))**3 * nn
       end do
    end if

    aero_binned%vol_den = 0d0
    do k = 1,bin_grid%n_bin
       aero_binned%vol_den(k,1) = const%pi/6d0 &
            * (2d0*vol2rad(bin_grid%v(k)))**3 * aero_binned%num_den(k)
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
  
end module pmc_kernel_golovin
