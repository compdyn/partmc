! Copyright (C) 2007 Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.
!
! Zero-valued kernel. This is only of interest for the exact solution
! to the no-coagulation, no-condensation case that can be used to test
! emissions and background dilution.

module pmc_kernel_zero

  use pmc_bin_grid
  use pmc_env
  use pmc_util
  use pmc_constants
  use pmc_aero_binned
  use pmc_aero_dist
  use pmc_aero_data
  
contains
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine kernel_zero(v1, v2, env, k)

    ! Zero coagulation kernel.

    real*8, intent(in) :: v1            ! volume of first particle
    real*8, intent(in) :: v2            ! volume of second particle
    type(env_t), intent(in) :: env      ! environment state
    real*8, intent(out) :: k            ! coagulation kernel
    
    k = 0d0
    
  end subroutine kernel_zero
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine soln_zero(bin_grid, aero_data, time, num_conc, mean_vol, &
       rho_p, aero_dist_init, env, aero_binned)

    ! Exact solution with the zero coagulation kernel. Only useful for
    ! testing emissions and background dilution.

    type(bin_grid_t), intent(in) :: bin_grid ! bin grid
    type(aero_data_t), intent(in) :: aero_data ! aerosol data
    real*8, intent(in) :: time          ! current time (s)
    real*8, intent(in) :: num_conc      ! particle number concentration (#/m^3)
    real*8, intent(in) :: mean_vol      ! mean init volume (m^3)
    real*8, intent(in) :: rho_p         ! particle density (kg/m^3)
    type(aero_dist_t), intent(in) :: aero_dist_init ! initial distribution
    type(env_t), intent(in) :: env      ! environment state
    type(aero_binned_t), intent(out) :: aero_binned ! output state

    type(aero_binned_t) :: aero_binned_limit

    ! With only emissions and dilution the number distribution n(r,t)
    ! satisfies:
    !
    !  d n(r,t)
    ! ---------- = k_emit * n_emit(r) + k_dilute * (n_back(r) - n(r,t))
    !     dt
    !
    ! n(r,0) = n_init(r)
    !
    ! This is a family of ODEs parameterized by r with solution:
    !
    ! n(r,t) = (n_init(r) - n_lim(r)) * exp(-k_dilute * t) + n_lim(r)
    !
    ! where the steady state limit is:
    !
    !                                     k_emit
    ! n(r,inf) = n_lim(r) = n_back(r) + ---------- n_emit(r)
    !                                    k_dilute

    if (env%aero_dilution_rate == 0d0) then
       call aero_binned_zero(aero_binned)
       call aero_binned_add_aero_dist(aero_binned, bin_grid, env%aero_emissions)
       call aero_binned_scale(aero_binned, env%aero_emission_rate * time &
            / env%height)
    else
       ! calculate the limit steady state distribution
       call aero_binned_alloc(aero_binned_limit, bin_grid%n_bin, &
            aero_data%n_spec)
       call aero_binned_add_aero_dist(aero_binned_limit, bin_grid, &
            env%aero_emissions)
       call aero_binned_scale(aero_binned_limit, &
            env%aero_emission_rate / env%height / env%aero_dilution_rate)
       call aero_binned_add_aero_dist(aero_binned_limit, bin_grid, &
            env%aero_background)

       ! calculate the current state
       call aero_binned_zero(aero_binned)
       call aero_binned_add_aero_dist(aero_binned, bin_grid, aero_dist_init)
       call aero_binned_sub(aero_binned, aero_binned_limit)
       call aero_binned_scale(aero_binned, exp(-env%aero_dilution_rate * time))
       call aero_binned_add(aero_binned, aero_binned_limit)

       call aero_binned_free(aero_binned_limit)
    end if

  end subroutine soln_zero

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
end module pmc_kernel_zero
