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
     !> Whether to run CAMP.
     logical :: do_camp_chem
     !> Whether to run TChem.
     logical :: do_tchem
     !> Whether to do coagulation.
     logical :: do_coagulation
     !> Whether to do condensation.
     logical :: do_condensation
     !> Whether to run MOSAIC.
     logical :: do_mosaic
     !> Whether to compute optical properties.
     logical :: do_optical
     !> Whether to do nucleation.
     logical :: do_nucleation
     !> Whether to do immersion freezing.
     logical :: do_immersion_freezing
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
    type(scenario_t), intent(in) :: scenario
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

    integer i_time, n_time, i_summary
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
       call output_modal(run_modal_opt%prefix, aero_binned, aero_dist, &
            aero_data, env_state, gas_data, gas_state, bin_grid, scenario, &
            i_summary, time, run_modal_opt%del_t, run_modal_opt%uuid)
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
          call output_modal(run_modal_opt%prefix, aero_binned, aero_dist, &
               aero_data, env_state, gas_data, gas_state, bin_grid, &
               scenario, i_summary, time, run_modal_opt%del_t, &
               run_modal_opt%uuid)
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

  !> Read the specification for a run_modal simulation from a spec file.
  subroutine spec_file_read_run_modal(file, run_modal_opt, aero_data, &
       bin_grid, gas_data, env_state, aero_dist_init, scenario)

    !> Spec file.
    type(spec_file_t), intent(inout) :: file
    !> Options controlling the operation of run_modal().
    type(run_modal_opt_t), intent(inout) :: run_modal_opt
    !> Aerosol data.
    type(aero_data_t), intent(out) :: aero_data
    !> Bin grid.
    type(bin_grid_t), intent(out) :: bin_grid
    !> Gas data.
    type(gas_data_t), intent(out) :: gas_data
    !> Environment state.
    type(env_state_t), intent(out) :: env_state
    !> Initial aerosol distribution.
    type(aero_dist_t), intent(out) :: aero_dist_init
    !> Scenario.
    type(scenario_t), intent(out) :: scenario

    character(len=PMC_MAX_FILENAME_LEN) :: sub_filename
    type(spec_file_t) :: sub_file

    call spec_file_read_string(file, 'output_prefix', run_modal_opt%prefix)

    call spec_file_read_real(file, 't_max', run_modal_opt%t_max)
    call spec_file_read_real(file, 'del_t', run_modal_opt%del_t)
    call spec_file_read_real(file, 't_output', run_modal_opt%t_output)
    call spec_file_read_real(file, 't_progress', run_modal_opt%t_progress)

    call spec_file_read_radius_bin_grid(file, bin_grid)

    call spec_file_read_logical(file, 'do_camp_chem', run_modal_opt%do_camp_chem)
    if (run_modal_opt%do_camp_chem) then
       call spec_file_die_msg(263948175, file, &
            "modal run does not support CAMP chemistry")
    end if

    call spec_file_read_logical(file, 'do_tchem', run_modal_opt%do_tchem)
    if (run_modal_opt%do_tchem) then
       call spec_file_die_msg(195837264, file, &
            "modal run does not support TChem chemistry")
    end if

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

    call spec_file_read_logical(file, 'do_coagulation', run_modal_opt%do_coagulation)
    if (run_modal_opt%do_coagulation) then
       call spec_file_die_msg(473829156, file, &
            "modal run does not support coagulation")
    end if

    call spec_file_read_logical(file, 'do_condensation', run_modal_opt%do_condensation)
    if (run_modal_opt%do_condensation) then
       call spec_file_die_msg(612938475, file, &
            "modal run does not support condensation")
    end if

    call spec_file_read_logical(file, 'do_mosaic', run_modal_opt%do_mosaic)
    if (run_modal_opt%do_mosaic) then
       call spec_file_die_msg(584729163, file, &
            "modal run does not support MOSAIC chemistry")
    end if

    call spec_file_read_logical(file, 'do_optical', run_modal_opt%do_optical)
    if (run_modal_opt%do_optical) then
       call spec_file_die_msg(527436819, file, &
            "modal run does not support optical properties calculation")
    end if

    call spec_file_read_logical(file, 'do_nucleation', run_modal_opt%do_nucleation)
    if (run_modal_opt%do_nucleation) then
       call spec_file_die_msg(391847265, file, &
            "modal run does not support nucleation")
    end if

    call spec_file_read_logical(file, 'do_immersion_freezing', &
         run_modal_opt%do_immersion_freezing)
    if (run_modal_opt%do_immersion_freezing) then
       call spec_file_die_msg(748291635, file, &
            "modal run does not support immersion freezing")
    end if

    call spec_file_close(file)

    ! Modal runs do not support aerosol emissions, background dilution,
    ! or loss functions other than drydep/none.
    if (size(scenario%aero_emission_rate_scale) > 0) then
       call assert_msg(583920417, &
            all(scenario%aero_emission_rate_scale == 0.0d0), &
            "modal run does not support aerosol emissions: " &
            // "all aero_emission rates must be zero")
    end if
    if (size(scenario%aero_dilution_rate) > 0) then
       call assert_msg(742916380, &
            all(scenario%aero_dilution_rate == 0.0d0), &
            "modal run does not support aerosol background dilution: " &
            // "all aero_dilution rates must be zero")
    end if
    call assert_msg(631804952, &
         scenario%loss_function_type == SCENARIO_LOSS_FUNCTION_NONE &
         .or. scenario%loss_function_type == SCENARIO_LOSS_FUNCTION_DRYDEP, &
         "modal run only supports loss_function none or drydep")

    ! All data from spec file read. Do the modal run.

    call pmc_srand(0,0)

    call uuid4_str(run_modal_opt%uuid)

    call scenario_init_env_state(scenario, env_state, 0d0)

  end subroutine spec_file_read_run_modal

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_run_modal
