! Copyright (C) 2005-2016 Nicole Riemer and Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_run_exact module.

!> Exact solution simulation.
module pmc_run_exact

  use pmc_aero_dist
  use pmc_bin_grid
  use pmc_aero_state
  use pmc_scenario
  use pmc_env_state
  use pmc_aero_data
  use pmc_output
  use pmc_aero_binned
  use pmc_gas_data
  use pmc_gas_state
  use pmc_exact_soln

  !> Options controlling the execution of run_exact().
  type run_exact_opt_t
     !> Total simulation time.
     real(kind=dp) :: t_max
     !> Interval to output info (s).
     real(kind=dp) :: t_output
     !> Output prefix.
     character(len=300) :: prefix
     !> Whether to do coagulation.
     logical :: do_coagulation
     !> Type of coagulation kernel.
     integer :: coag_kernel_type
     !> UUID of the simulation.
     character(len=PMC_UUID_LEN) :: uuid
  end type run_exact_opt_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Run an exact simulation.
  subroutine run_exact(bin_grid, scenario, env_state, aero_data, &
       aero_dist_init, gas_data, run_exact_opt)

    !> Bin grid.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Environment data.
    type(scenario_t), intent(in) :: scenario
    !> Environment state.
    type(env_state_t), intent(inout) :: env_state
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Initial aerosol distribution.
    type(aero_dist_t), intent(in) :: aero_dist_init
    !> Gas data.
    type(gas_data_t), intent(in) :: gas_data
    !> Options.
    type(run_exact_opt_t), intent(in) :: run_exact_opt

    integer :: i_time, n_time, ncid
    type(aero_binned_t) :: aero_binned
    real(kind=dp) :: time
    type(gas_state_t) :: gas_state

    call check_time_multiple("t_max", run_exact_opt%t_max, &
         "t_output", run_exact_opt%t_output)

    call gas_state_set_size(gas_state, gas_data_n_spec(gas_data))

    n_time = nint(run_exact_opt%t_max / run_exact_opt%t_output)
    do i_time = 0,n_time
       time = real(i_time, kind=dp) / real(n_time, kind=dp) &
            * run_exact_opt%t_max
       call scenario_update_env_state(scenario, env_state, time)
       call exact_soln(bin_grid, aero_data, run_exact_opt%do_coagulation, &
            run_exact_opt%coag_kernel_type, aero_dist_init, scenario, &
            env_state, time, aero_binned)
       call output_sectional(run_exact_opt%prefix, bin_grid, aero_data, &
            aero_binned, gas_data, gas_state, env_state, i_time + 1, &
            time, run_exact_opt%t_output, run_exact_opt%uuid)
    end do

  end subroutine run_exact

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_run_exact
