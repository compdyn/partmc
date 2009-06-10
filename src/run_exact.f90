! Copyright (C) 2005-2009 Nicole Riemer and Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_run_exact module.

!> Exact solution simulation.
module pmc_run_exact

  use pmc_aero_dist
  use pmc_bin_grid
  use pmc_aero_state
  use pmc_env_data
  use pmc_env_state
  use pmc_aero_data
  use pmc_output_sectional
  use pmc_aero_binned
  use pmc_gas_data
  use pmc_gas_state

  !> Options controlling the execution of run_exact().
  type run_exact_opt_t
     ! FIXME: following few items depend on kernel/soln choice
     !> Particle number concentration (#/m^3).
     real*8 :: num_den
     !> Mean init radius (m).
     real*8 :: mean_radius
     !> Aerosol initial distribution.
     type(aero_dist_t) :: aero_dist_init
     !> Particle density (kg/m^3).
     real*8 :: rho_p
     !> Total simulation time.
     real*8 :: t_max
     !> Interval to output info (s).
     real*8 :: t_output
     !> Output prefix.
     character(len=300) :: prefix
  end type run_exact_opt_t

contains
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Run an exact simulation.
  !!
  !! FIXME: num_den and mean_radius are really parameters for the
  !! initial value of the particle distribution. They should be
  !! replaced by a n_param, params() pair.
  subroutine run_exact(bin_grid, env_data, env_state, aero_data, &
       exact_opt, soln)

    !> Bin grid.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Environment data.
    type(env_data_t), intent(in) :: env_data
    !> Environment state.
    type(env_state_t), intent(inout) :: env_state
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Options.
    type(run_exact_opt_t), intent(in) :: exact_opt
    
    integer :: i_time, n_time, ncid
    type(aero_binned_t) :: aero_binned
    real*8 :: time
    type(gas_data_t) :: gas_data
    type(gas_state_t) :: gas_state
    
    interface
       subroutine soln(bin_grid, aero_data, time, num_den, mean_radius, &
            rho_p, aero_dist_init, env_state, aero_binned)

         use pmc_bin_grid
         use pmc_env_state
         use pmc_aero_binned
         use pmc_aero_data

         !> Bin grid.
         type(bin_grid_t), intent(in) :: bin_grid
         !> Aerosol data.
         type(aero_data_t), intent(in) :: aero_data
         !> Current time.
         real*8, intent(in) :: time
         !> Particle number conc (#/m^3).
         real*8, intent(in) :: num_den
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
       end subroutine soln
    end interface

    call aero_binned_alloc(aero_binned, bin_grid%n_bin, aero_data%n_spec)
    call gas_data_alloc(gas_data, 0)
    call gas_state_alloc(gas_state, 0)

    n_time = nint(exact_opt%t_max / exact_opt%t_output)
    do i_time = 0,n_time
       time = dble(i_time) / dble(n_time) * exact_opt%t_max
       call env_data_update_state(env_data, env_state, time)
       call soln(bin_grid, aero_data, time, exact_opt%num_den, &
            exact_opt%mean_radius, exact_opt%rho_p, &
            exact_opt%aero_dist_init, env_state, aero_binned)
       call output_sectional(exact_opt%prefix, bin_grid, aero_data, &
            aero_binned, gas_data, gas_state, env_state, i_time + 1, &
            time, exact_opt%t_output)
    end do

    call gas_data_free(gas_data)
    call gas_state_free(gas_state)
    call aero_binned_free(aero_binned)
    
  end subroutine run_exact
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
end module pmc_run_exact
