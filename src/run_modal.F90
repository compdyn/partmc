!> \file
!> The pmc_run_modal module.

!> 1D modal simulation.

module pmc_run_modal

    use pmc_util
    use pmc_aero_dist
    use pmc_scenario
    use pmc_env_state
    use pmc_aero_data
    use pmc_output
    use pmc_bin_grid
    use pmc_aero_binned
    use pmc_gas_data
    use pmc_gas_state
    use pmc_constants

    !> Options controlling the operation of run_modal()
    type run_modal_opt_t
      !> Final time (s).
      real(kind=dp) :: t_max
      !> Timestep (s).
      real(kind=dp) :: del_t
      !> Output interval (0 disables) (s).
      real(kind=dp) :: t_output
      !> Progress interval (0 disables) (s).
      real(kind=dp) :: t_progress
      !> Output prefix.
      character(len=300) :: prefix
      !> UUID of the simulation.
      character(len=PMC_UUID_LEN) :: uuid
    end type run_modal_opt_t

contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !> Run a modal simulation.
    subroutine run_modal(aero_data, aero_dist, scenario, env_state, &
        gas_data, bin_grid, run_modal_opt)

      !> Aerosol data.
      type(aero_data_t), intent(in) :: aero_data
      !> Aerosol distribution.
      type(aero_dist_t), intent(inout) :: aero_dist
      !> Environment data.
      type(scenario_t), intent(inout) :: scenario
      !> Environment state.
      type(env_state_t), intent(inout) :: env_state
      !> Gas data.
      type(gas_data_t), intent(in) :: gas_data
      !> Bin grid.
      type(bin_grid_t), intent(in) :: bin_grid
      !> Modal options.
      type(run_modal_opt_t), intent(in) :: run_modal_opt

      type(env_state_t) :: old_env_state
      type(gas_state_t) :: gas_state
      type(aero_binned_t) :: aero_binned

      real(kind=dp) time, last_output_time, last_progress_time

      integer i, i_time, n_time, i_summary, i_mode, mode
      logical do_output, do_progress

      call check_time_multiple("t_max", run_modal_opt%t_max, &
                               "del_t", run_modal_opt%del_t)
      call check_time_multiple("t_output", run_modal_opt%t_output, &
                               "del_t", run_modal_opt%del_t)
      call check_time_multiple("t_progress", run_modal_opt%t_progress, &
                               "del_t", run_modal_opt%del_t)

      if (aero_data_n_spec(aero_data) /= 1) then
        call die_msg(927384615, &
            'run_modal() can only use one aerosol species')
      end if

      ! output data structure
      call gas_state_set_size(gas_state, gas_data_n_spec(gas_data))   

      ! Initialize time.
      time = 0d0
      last_output_time = 0d0
      last_progress_time = 0d0
      i_summary = 1

      ! initial output
      call check_event(time, run_modal_opt%del_t, run_modal_opt%t_output, &
                       last_output_time, do_output)
      if (do_output) then
        call aero_binned_add_aero_dist(aero_binned, bin_grid, aero_data, &
                                       aero_dist)
        call output_modal(run_modal_opt%prefix, aero_binned, aero_dist, aero_data, &
                          env_state, gas_data, gas_state, bin_grid, scenario, i_summary, time, &
                          run_modal_opt%del_t, run_modal_opt%uuid)
        call aero_binned_zero(aero_binned)
      end if

      ! Main time-stepping loop
      n_time = nint(run_modal_opt%t_max / run_modal_opt%del_t)
      do i_time = 1,n_time

        time = run_modal_opt%t_max * real(i_time, kind=dp) &
             / real(n_time, kind=dp)

        old_env_state = env_state

        call scenario_update_env_state(scenario, env_state, time)
        call scenario_update_gas_state(scenario, run_modal_opt%del_t, &
             env_state, old_env_state, gas_data, gas_state)
        call scenario_update_aero_modes(aero_dist, run_modal_opt%del_t, &
             env_state, aero_data%density(1), scenario)

        call check_event(time, run_modal_opt%del_t, run_modal_opt%t_output, &
             last_output_time, do_output)
        if (do_output) then
          i_summary = i_summary + 1
          call aero_binned_add_aero_dist(aero_binned, bin_grid, aero_data, &
               aero_dist)
          call output_modal(run_modal_opt%prefix, aero_binned, aero_dist, aero_data, &
               env_state, gas_data, gas_state, bin_grid, scenario,i_summary, time, &
               run_modal_opt%del_t, run_modal_opt%uuid)
          call aero_binned_zero(aero_binned)
        end if

        call check_event(time, run_modal_opt%del_t, run_modal_opt%t_progress, &
             last_progress_time, do_progress)
        if (do_progress) then
          write(*, '(a6,a8)') 'step', 'time'
          write(*, '(i6,f8.1)') i_time, time
        end if
      end do

    end subroutine run_modal

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_run_modal
