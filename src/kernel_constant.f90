! Copyright (C) 2005-2007 Nicole Riemer and Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.
!
! Constant coagulation kernel.

module pmc_kernel_constant

  use pmc_env_state
  use pmc_bin_grid
  use pmc_util
  use pmc_constants
  use pmc_aero_binned
  use pmc_aero_data
  use pmc_aero_dist
  
contains
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine kernel_constant(v1, v2, env_state, k)

    ! Constant coagulation kernel.

    real*8, intent(in) :: v1            ! volume of first particle
    real*8, intent(in) :: v2            ! volume of second particle
    type(env_state_t), intent(in) :: env_state      ! environment state
    real*8, intent(out) :: k            ! coagulation kernel
    
    real*8, parameter :: beta_0 = 0.25d0 / (60d0 * 2d8)
    
    k = beta_0
    
  end subroutine kernel_constant
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine soln_constant_exp_cond(bin_grid, aero_data, time, num_den, &
       mean_radius, rho_p, aero_dist_init, env_state, aero_binned)

    ! Exact solution with a constant coagulation kernel and an
    ! exponential initial condition.

    type(bin_grid_t), intent(in) :: bin_grid ! bin grid
    type(aero_data_t), intent(in) :: aero_data ! aerosol data
    real*8, intent(in) :: time          ! current time
    real*8, intent(in) :: num_den       ! particle number concentration (#/m^3)
    real*8, intent(in) :: mean_radius   ! mean init radius (m)
    real*8, intent(in) :: rho_p         ! particle density (kg/m^3)
    type(aero_dist_t), intent(in) :: aero_dist_init ! initial distribution
    type(env_state_t), intent(in) :: env_state      ! environment state
    type(aero_binned_t), intent(out) :: aero_binned ! output state
    
    real*8 :: beta_0, tau, T, rat_v, nn, b, x, sigma, mean_vol
    integer :: k
    
    real*8, parameter :: lambda = 1d0   ! FIXME: what is this?
    
    call kernel_constant(1d0, 1d0, env_state, beta_0)

    mean_vol = rad2vol(mean_radius)
    if (time .eq. 0d0) then
       do k = 1,bin_grid%n_bin
          aero_binned%num_den(k) = const%pi/2d0 &
               * (2d0*vol2rad(bin_grid%v(k)))**3 * num_den / mean_vol &
               * exp(-(bin_grid%v(k)/mean_vol))
       end do
    else
       tau = num_den * beta_0 * time
       do k = 1,bin_grid%n_bin
          rat_v = bin_grid%v(k) / mean_vol
          x = 2d0 * rat_v / (tau + 2d0)
          nn = 4d0 * num_den / (mean_vol * ( tau + 2d0 ) ** 2d0) &
               * exp(-2d0*rat_v/(tau+2d0)*exp(-lambda*tau)-lambda*tau)
          aero_binned%num_den(k) = const%pi/2d0 &
               * (2d0*vol2rad(bin_grid%v(k)))**3d0 * nn
       end do
    end if
    
    aero_binned%vol_den = 0d0
    do k = 1,bin_grid%n_bin
       aero_binned%vol_den(k,1) = const%pi/6d0 &
            * (2d0*vol2rad(bin_grid%v(k)))**3d0 * aero_binned%num_den(k)
    end do
    
  end subroutine soln_constant_exp_cond
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
end module pmc_kernel_constant
