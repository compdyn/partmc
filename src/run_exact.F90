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

  !> Read the specification for a run_exact simulation from a spec file.
  subroutine spec_file_read_run_exact(file, run_exact_opt, aero_data, &
       bin_grid, gas_data, env_state, aero_dist_init, scenario)

    !> Spec file.
    type(spec_file_t), intent(inout) :: file
    !> Options controlling the operation of run_exact().
    type(run_exact_opt_t), intent(inout) :: run_exact_opt
    !> Aerosol data.
    type(aero_data_t), intent(out) :: aero_data
    !> Bin grid.
    type(bin_grid_t), intent(out) :: bin_grid
    !> Initial aerosol state.
    type(aero_dist_t), intent(out) :: aero_dist_init
    !> Scenario data.
    type(scenario_t), intent(out) :: scenario
    !> Environmental state.
    type(env_state_t), intent(out) :: env_state
    !> Gas data.
    type(gas_data_t), intent(out) :: gas_data

    character(len=PMC_MAX_FILENAME_LEN) :: sub_filename
    type(spec_file_t) :: sub_file

    call spec_file_read_string(file, 'output_prefix', run_exact_opt%prefix)

    call spec_file_read_real(file, 't_max', run_exact_opt%t_max)
    call spec_file_read_real(file, 't_output', run_exact_opt%t_output)

    call spec_file_read_radius_bin_grid(file, bin_grid)

    call spec_file_read_string(file, 'gas_data', sub_filename)
    call spec_file_open(sub_filename, sub_file)
    call spec_file_read_gas_data(sub_file, gas_data)
    call spec_file_close(sub_file)

    call spec_file_read_string(file, 'aerosol_data', sub_filename)
    call spec_file_open(sub_filename, sub_file)
    call spec_file_read_aero_data(sub_file, aero_data)
    call spec_file_close(sub_file)

    call spec_file_read_fractal(file, aero_data%fractal)

    call spec_file_read_string(file, 'aerosol_init', sub_filename)
    call spec_file_open(sub_filename, sub_file)
    call spec_file_read_aero_dist(sub_file, aero_data, .false., aero_dist_init)
    call spec_file_close(sub_file)

    call spec_file_read_scenario(file, gas_data, aero_data, .false., scenario)
    call spec_file_read_env_state(file, env_state)

    call spec_file_read_logical(file, 'do_coagulation', &
         run_exact_opt%do_coagulation)
    if (run_exact_opt%do_coagulation) then
       call spec_file_read_coag_kernel_type(file, &
            run_exact_opt%coag_kernel_type)
       if (run_exact_opt%coag_kernel_type == COAG_KERNEL_TYPE_ADDITIVE) then
          call spec_file_read_real(file, 'additive_kernel_coeff', &
               env_state%additive_kernel_coefficient)
       end if
    else
       run_exact_opt%coag_kernel_type = COAG_KERNEL_TYPE_INVALID
    end if

    call spec_file_close(file)

    ! finished reading .spec data, now do the run

    call pmc_srand(0, 0)

    call uuid4_str(run_exact_opt%uuid)

    call scenario_init_env_state(scenario, env_state, 0d0)

  end subroutine spec_file_read_run_exact

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_run_exact
