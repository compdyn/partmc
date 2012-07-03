! Copyright (C) 2005-2012 Nicole Riemer and Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_coag_kernel_constant module.

!> Constant coagulation kernel.
module pmc_coag_kernel_constant

  use pmc_env_state
  use pmc_bin_grid
  use pmc_util
  use pmc_constants
  use pmc_aero_binned
  use pmc_aero_data
  use pmc_aero_dist
  use pmc_aero_data
  use pmc_aero_particle

  !> Coefficient for constant kernel.
  real(kind=dp), parameter :: beta_0 = 0.25d0 / (60d0 * 2d8)

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
    real(kind=dp), intent(out) :: k

    k = beta_0

  end subroutine kernel_constant

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Minimum and maximum values of the constant coagulation kernel.
  subroutine kernel_constant_minmax(v1, v2, aero_data, env_state, k_min, k_max)

    !> Volume of first particle.
    real(kind=dp), intent(in) :: v1
    !> Volume of second particle.
    real(kind=dp), intent(in) :: v2
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Environment state.
    type(env_state_t), intent(in) :: env_state
    !> Coagulation kernel minimum value.
    real(kind=dp), intent(out) :: k_min
    !> Coagulation kernel maximum value.
    real(kind=dp), intent(out) :: k_max

    k_min = beta_0
    k_max = beta_0

  end subroutine kernel_constant_minmax

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Exact solution with a constant coagulation kernel and an
  !> exponential initial condition.
  !!
  !! Given input paramaters \f$R\f$ and \f$N_0\f$ we let the mean
  !! volume be \f$v_\mu = \frac{4\pi}{3} R^3\f$ and define the
  !! rescaled time \f$\tau = N_0 \beta_0 t\f$, where \f$\beta_0\f$ is
  !! the fixed constant kernel value. We also set the parameter
  !! \f$\lambda = 1\f$. Then the solution is
  !! \f[
  !!     n(D,t) \ {\rm d}\ln D
  !!     = \frac{\pi}{2} D^3 \frac{N_0}{v_\mu} \frac{4}{(\tau + 2)^2}
  !!     \exp\left(-\frac{v}{v_\mu} \frac{2}{\tau + 2}
  !!     \exp(-\lambda \tau) - \lambda \tau\right) {\rm d}\ln D
  !! \f]
  !! This thus has initial condition
  !! \f[
  !!     n(D,t) \ {\rm d}\ln D
  !!     = \frac{\pi}{2} D^3 \frac{N_0}{v_\mu}
  !!     \exp\left(-\frac{v}{v_\mu}\right) {\rm d}\ln D
  !! \f]
  subroutine soln_constant_exp(bin_grid, aero_data, time, num_conc, &
       radius_at_mean_vol, env_state, aero_binned)

    !> Bin grid.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Current time.
    real(kind=dp), intent(in) :: time
    !> Particle number concentration (#/m^3).
    real(kind=dp), intent(in) :: num_conc
    !> Mean init radius (m).
    real(kind=dp), intent(in) :: radius_at_mean_vol
    !> Environment state.
    type(env_state_t), intent(in) :: env_state
    !> Output state.
    type(aero_binned_t), intent(inout) :: aero_binned

    real(kind=dp) :: tau, T, rat_v, nn, b, sigma, mean_vol
    integer :: k

    real(kind=dp), parameter :: lambda = 1d0

    call aero_binned_set_sizes(aero_binned, bin_grid_size(bin_grid), &
         aero_data_n_spec(aero_data))

    mean_vol = aero_data_rad2vol(aero_data, radius_at_mean_vol)
    if (time .eq. 0d0) then
       do k = 1,bin_grid_size(bin_grid)
          aero_binned%num_conc(k) = const%pi / 2d0 &
               * (2d0 * bin_grid%centers(k))**3 * num_conc / mean_vol &
               * exp(-(aero_data_rad2vol(aero_data, bin_grid%centers(k)) &
               / mean_vol))
       end do
    else
       tau = num_conc * beta_0 * time
       do k = 1,bin_grid_size(bin_grid)
          rat_v = aero_data_rad2vol(aero_data, bin_grid%centers(k)) / mean_vol
          nn = 4d0 * num_conc / (mean_vol * ( tau + 2d0 ) ** 2d0) &
               * exp(-2d0*rat_v/(tau+2d0)*exp(-lambda*tau)-lambda*tau)
          aero_binned%num_conc(k) = const%pi / 2d0 &
               * (2d0 * bin_grid%centers(k))**3d0 * nn
       end do
    end if

    aero_binned%vol_conc = 0d0
    do k = 1,bin_grid_size(bin_grid)
       aero_binned%vol_conc(k,1) = aero_data_rad2vol(aero_data, &
            bin_grid%centers(k)) * aero_binned%num_conc(k)
    end do

  end subroutine soln_constant_exp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_coag_kernel_constant
