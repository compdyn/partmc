! Copyright (C) 2005-2009 Nicole Riemer and Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_kernel_constant module.

!> Constant coagulation kernel.
module pmc_kernel_constant

  use pmc_env_state
  use pmc_bin_grid
  use pmc_util
  use pmc_constants
  use pmc_aero_binned
  use pmc_aero_data
  use pmc_aero_dist
  use pmc_aero_data
  use pmc_aero_particle
  
contains
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Constant coagulation kernel.
  subroutine kernel_constant(aero_particle_1, aero_particle_2, &
       aero_data, env_state, k)

    !> First particle.
    type(aero_particle_t), intent(in) :: aero_particle_1
    !> Second particle.
    type(aero_particle_t), intent(in) :: aero_particle_2
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Environment state.
    type(env_state_t), intent(in) :: env_state
    !> Coagulation kernel.
    real*8, intent(out) :: k
    
    call kernel_constant_max(aero_particle_volume(aero_particle_1), &
         aero_particle_volume(aero_particle_2), aero_data, env_state, k)
    
  end subroutine kernel_constant
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Maximum value of the constant coagulation kernel.
  subroutine kernel_constant_max(v1, v2, aero_data, env_state, k_max)

    !> Volume of first particle.
    real*8, intent(in) :: v1
    !> Volume of second particle.
    real*8, intent(in) :: v2
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Environment state.
    type(env_state_t), intent(in) :: env_state
    !> Coagulation kernel maximum value.
    real*8, intent(out) :: k_max

    real*8, parameter :: beta_0 = 0.25d0 / (60d0 * 2d8)

    k_max = beta_0
    
  end subroutine kernel_constant_max
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Exact solution with a constant coagulation kernel and an
  !> exponential initial condition.
  subroutine soln_constant_exp_cond(bin_grid, aero_data, time, num_conc, &
       mean_radius, rho_p, aero_dist_init, env_state, aero_binned)

    !> Bin grid.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Current time.
    real*8, intent(in) :: time
    !> Particle number concentration (#/m^3).
    real*8, intent(in) :: num_conc
    !> Mean init radius (m).
    real*8, intent(in) :: mean_radius
    !> Particle density (kg/m^3).
    real*8, intent(in) :: rho_p
    !> Initial distribution.
    type(aero_dist_t), intent(in) :: aero_dist_init
    !> Environment state.
    type(env_state_t), intent(in) :: env_state
    !> Output state.
    type(aero_binned_t), intent(out) :: aero_binned
    
    real*8 :: beta_0, tau, T, rat_v, nn, b, x, sigma, mean_vol
    integer :: k
    
    real*8, parameter :: lambda = 1d0 ! FIXME: what is this?

    call kernel_constant_max(1d0, 1d0, aero_data, env_state, beta_0)
    
    mean_vol = rad2vol(mean_radius)
    if (time .eq. 0d0) then
       do k = 1,bin_grid%n_bin
          aero_binned%num_conc(k) = const%pi/2d0 &
               * (2d0*vol2rad(bin_grid%v(k)))**3 * num_conc / mean_vol &
               * exp(-(bin_grid%v(k)/mean_vol))
       end do
    else
       tau = num_conc * beta_0 * time
       do k = 1,bin_grid%n_bin
          rat_v = bin_grid%v(k) / mean_vol
          x = 2d0 * rat_v / (tau + 2d0)
          nn = 4d0 * num_conc / (mean_vol * ( tau + 2d0 ) ** 2d0) &
               * exp(-2d0*rat_v/(tau+2d0)*exp(-lambda*tau)-lambda*tau)
          aero_binned%num_conc(k) = const%pi/2d0 &
               * (2d0*vol2rad(bin_grid%v(k)))**3d0 * nn
       end do
    end if
    
    aero_binned%vol_conc = 0d0
    do k = 1,bin_grid%n_bin
       aero_binned%vol_conc(k,1) = const%pi/6d0 &
            * (2d0*vol2rad(bin_grid%v(k)))**3d0 * aero_binned%num_conc(k)
    end do
    
  end subroutine soln_constant_exp_cond
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
end module pmc_kernel_constant
