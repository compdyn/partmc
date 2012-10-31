! Copyright (C) 2007-2012 Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_coag_kernel_zero module.

!> Constant kernel equal to zero.
!!
!! This is only of interest for the exact solution to the
!! no-coagulation, no-condensation case that can be used to test
!! emissions and background dilution.
module pmc_coag_kernel_zero

  use pmc_bin_grid
  use pmc_scenario
  use pmc_env_state
  use pmc_util
  use pmc_aero_binned
  use pmc_aero_dist
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

  !> Minimum and maximum of the zero coagulation kernel.
  subroutine kernel_zero_minmax(v1, v2, aero_data, env_state, k_min, k_max)

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
    
    k_min = 0d0
    k_max = 0d0
    
  end subroutine kernel_zero_minmax
  
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
  subroutine soln_zero(bin_grid, aero_data, time, aero_dist_init, &
       scenario, env_state, loss_function_type, aero_binned)

    !> Bin grid.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Current time (s).
    real(kind=dp), intent(in) :: time
    !> Initial distribution.
    type(aero_dist_t), intent(in) :: aero_dist_init
    !> Scenario.
    type(scenario_t), intent(in) :: scenario
    !> Environment state.
    type(env_state_t), intent(in) :: env_state
    !> Particle loss function type.
    integer, intent(in) :: loss_function_type
    !> Output state.
    type(aero_binned_t), intent(inout) :: aero_binned
    
    
    !FIXME: add functionality for loss_function_type
    

    real(kind=dp) :: emission_rate_scale, dilution_rate, p
    type(aero_dist_t) :: emissions, background
    type(aero_binned_t) :: background_binned, aero_binned_limit

    call aero_dist_allocate(emissions)
    call aero_dist_allocate(background)
    call aero_binned_allocate_size(background_binned, bin_grid%n_bin, &
         aero_data%n_spec)

    call aero_dist_interp_1d(scenario%aero_emission, &
         scenario%aero_emission_time, scenario%aero_emission_rate_scale, &
         env_state%elapsed_time, emissions, emission_rate_scale)
    call aero_dist_interp_1d(scenario%aero_background, &
         scenario%aero_dilution_time, scenario%aero_dilution_rate, 0d0, &
         background, dilution_rate)
    call aero_binned_add_aero_dist(background_binned, bin_grid, aero_data, &
         background)

    if (dilution_rate == 0d0) then
       call aero_binned_zero(aero_binned)
       call aero_binned_add_aero_dist(aero_binned, bin_grid, aero_data, &
            emissions)
       call aero_binned_scale(aero_binned, &
            emission_rate_scale * time / env_state%height)
    else
       ! calculate the limit steady state distribution
       call aero_binned_allocate_size(aero_binned_limit, bin_grid%n_bin, &
            aero_data%n_spec)
       call aero_binned_add_aero_dist(aero_binned_limit, bin_grid, &
            aero_data, emissions)
       call aero_binned_scale(aero_binned_limit, &
            emission_rate_scale / env_state%height / dilution_rate)
       call aero_binned_add_aero_dist(aero_binned_limit, bin_grid, &
            aero_data, background)

       ! calculate the current state
       call aero_binned_zero(aero_binned)
       call aero_binned_add_aero_dist(aero_binned, bin_grid, aero_data, &
            aero_dist_init)
       call aero_binned_sub(aero_binned, aero_binned_limit)
       call aero_binned_scale(aero_binned, &
            exp(-dilution_rate * time))
       call aero_binned_add(aero_binned, aero_binned_limit)

       call aero_binned_deallocate(aero_binned_limit)
    end if

    call aero_dist_deallocate(emissions)
    call aero_dist_deallocate(background)
    call aero_binned_deallocate(background_binned)

  end subroutine soln_zero

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
end module pmc_coag_kernel_zero
