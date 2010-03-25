! Copyright (C) 2007-2010 Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_kernel_zero module.

!> Constant kernel equal to zero.
!!
!! This is only of interest for the exact solution to the
!! no-coagulation, no-condensation case that can be used to test
!! emissions and background dilution.
module pmc_kernel_zero

  use pmc_bin_grid
  use pmc_env_state
  use pmc_util
  use pmc_aero_binned
  use pmc_aero_dist
  use pmc_aero_data
  use pmc_aero_data
  use pmc_aero_particle
  
contains
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Zero coagulation kernel.
  subroutine kernel_zero(aero_particle_1, aero_particle_2, &
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
    
    k = 0d0
    
  end subroutine kernel_zero
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Zero coagulation kernel.
  subroutine kernel_zero_max(v1, v2, aero_data, env_state, k_max)

    !> Volume of first particle.
    real(kind=dp), intent(in) :: v1
    !> Volume of second particle.
    real(kind=dp), intent(in) :: v2
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Environment state.
    type(env_state_t), intent(in) :: env_state
    !> Coagulation kernel maximum value.
    real(kind=dp), intent(out) :: k_max
    
    k_max = 0d0
    
  end subroutine kernel_zero_max
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Exact solution with the zero coagulation kernel. Only useful for
  !> testing emissions and background dilution.
  !!
  !! With only constant-rate emissions and dilution the number
  !! distribution \f$n(D,t)\f$ at diameter \f$D\f$ and time \f$t\f$
  !! satisfies:
  !! \f[
  !!     \frac{d n(D,t)}{dt} = k_{\rm emit} n_{\rm emit}(D)
  !!                           + k_{\rm dilute} (n_{\rm back}(D) - n(D,t))
  !! \f]
  !! together with the initial condition \f$ n(D,0) = n_0(D) \f$. Here
  !! \f$n_{\rm emit}(D)\f$ and \f$n_{\rm back}(D)\f$ are emission and
  !! background size distributions, with corresponding rates \f$k_{\rm
  !! emit}\f$ and \f$k_{\rm dilute}\f$. All values are taken at time
  !! \f$t = 0\f$ and held constant, so there is no support for
  !! time-varying emissions or background dilution.
  !!
  !! This is a family of ODEs parameterized by \f$D\f$ with
  !! solution:
  !! \f[
  !!     n(D,t) = n_{\infty}(D)
  !!              + (n_0(D) - n_{\infty}(D)) \exp(-k_{\rm dilute} t)
  !! \f]
  !! where the steady state limit is:
  !! \f[
  !!     n_{\infty}(D) = n(D,\infty)
  !!                   = n_{\rm back}(D)
  !!                     + \frac{k_{\rm emit}}{k_{\rm dilute}} n_{\rm emit}(D)
  !! \f]
  subroutine soln_zero(bin_grid, aero_data, time, num_conc, mean_vol, &
       rho_p, aero_dist_init, env_state, aero_binned)

    !> Bin grid.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Current time (s).
    real(kind=dp), intent(in) :: time
    !> Particle number concentration (#/m^3).
    real(kind=dp), intent(in) :: num_conc
    !> Mean init volume (m^3).
    real(kind=dp), intent(in) :: mean_vol
    !> Particle density (kg/m^3).
    real(kind=dp), intent(in) :: rho_p
    !> Initial distribution.
    type(aero_dist_t), intent(in) :: aero_dist_init
    !> Environment state.
    type(env_state_t), intent(in) :: env_state
    !> Output state.
    type(aero_binned_t), intent(out) :: aero_binned

    type(aero_binned_t) :: aero_binned_limit

    if (env_state%aero_dilution_rate == 0d0) then
       call aero_binned_zero(aero_binned)
       call aero_binned_add_aero_dist(aero_binned, bin_grid, aero_data, &
            env_state%aero_emissions)
       call aero_binned_scale(aero_binned, &
            env_state%aero_emission_rate * time / env_state%height)
    else
       ! calculate the limit steady state distribution
       call aero_binned_allocate_size(aero_binned_limit, bin_grid%n_bin, &
            aero_data%n_spec)
       call aero_binned_add_aero_dist(aero_binned_limit, bin_grid, &
            aero_data, env_state%aero_emissions)
       call aero_binned_scale(aero_binned_limit, &
            env_state%aero_emission_rate / env_state%height &
            / env_state%aero_dilution_rate)
       call aero_binned_add_aero_dist(aero_binned_limit, bin_grid, &
            aero_data, env_state%aero_background)

       ! calculate the current state
       call aero_binned_zero(aero_binned)
       call aero_binned_add_aero_dist(aero_binned, bin_grid, aero_data, &
            aero_dist_init)
       call aero_binned_sub(aero_binned, aero_binned_limit)
       call aero_binned_scale(aero_binned, &
            exp(-env_state%aero_dilution_rate * time))
       call aero_binned_add(aero_binned, aero_binned_limit)

       call aero_binned_deallocate(aero_binned_limit)
    end if

  end subroutine soln_zero

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
end module pmc_kernel_zero
