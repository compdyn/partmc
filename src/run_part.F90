! Copyright (C) 2005-2017 Nicole Riemer and Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_run_part module.

!> Monte Carlo simulation.
module pmc_run_part

  use pmc_constants
  use pmc_util
  use pmc_aero_state
  use pmc_scenario
  use pmc_env_state
  use pmc_aero_data
  use pmc_gas_data
  use pmc_gas_state
  use pmc_output
  use pmc_mosaic
  use pmc_coagulation
  use pmc_coagulation_dist
  use pmc_coag_kernel
  use pmc_nucleate
  use pmc_mpi
  use pmc_camp_interface
  use pmc_photolysis
#ifdef PMC_USE_CAMP
  use camp_camp_core
  use camp_camp_state
  use camp_env_state, only: camp_env_state_t => env_state_t
#endif
#ifdef PMC_USE_SUNDIALS
  use pmc_condense
#endif
#ifdef PMC_USE_MPI
  use mpi
#endif

  !> Type code for undefined or invalid parallel coagulation method.
  integer, parameter :: PARALLEL_COAG_TYPE_INVALID = 0
  !> Type code for local parallel coagulation.
  integer, parameter :: PARALLEL_COAG_TYPE_LOCAL   = 1
  !> Type code for distributed parallel coagulation.
  integer, parameter :: PARALLEL_COAG_TYPE_DIST    = 2

  !> Options controlling the execution of run_part().
  type run_part_opt_t
     !> Final time (s).
     real(kind=dp) :: t_max
     !> Output interval (0 disables) (s).
     real(kind=dp) :: t_output
     !> Progress interval (0 disables) (s).
     real(kind=dp) :: t_progress
     !> Timestep for coagulation.
     real(kind=dp) :: del_t
     !> Prefix for output files.
     character(len=300) :: output_prefix
     !> Type of coagulation kernel.
     integer :: coag_kernel_type
     !> Type of nucleation.
     integer :: nucleate_type
     !> Source of nucleation.
     integer :: nucleate_source
     !> Whether to do coagulation.
     logical :: do_coagulation
     !> Whether to do nucleation.
     logical :: do_nucleation
     !> Allow doubling if needed.
     logical :: allow_doubling
     !> Allow halving if needed.
     logical :: allow_halving
     !> Whether to do condensation.
     logical :: do_condensation
     !> Whether to do MOSAIC.
     logical :: do_mosaic
     !> Whether to compute optical properties.
     logical :: do_optical
     !> Whether to have explicitly selected weighting.
     logical :: do_select_weighting
     !> Type of particle weighting scheme.
     integer :: weighting_type
     !> Weighting exponent for power weighting scheme.
     real(kind=dp) :: weighting_exponent
     !> Repeat number of run.
     integer :: i_repeat
     !> Total number of repeats.
     integer :: n_repeat
     !> Cpu_time() of start.
     real(kind=dp) :: t_wall_start
     !> Whether to record particle removal information.
     logical :: record_removals
     !> Whether to run in parallel.
     logical :: do_parallel
     !> Parallel output type.
     integer :: output_type
     !> Mixing timescale between processes (s).
     real(kind=dp) :: mix_timescale
     !> Whether to average gases each timestep.
     logical :: gas_average
     !> Whether to average environment each timestep.
     logical :: env_average
     !> Parallel coagulation method type.
     integer :: parallel_coag_type
     !> Whether to run CAMP
     logical :: do_camp_chem
     !> UUID for this simulation.
     character(len=PMC_UUID_LEN) :: uuid
  end type run_part_opt_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Do a particle-resolved Monte Carlo simulation.
#ifdef PMC_USE_CAMP
  subroutine run_part(scenario, env_state, aero_data, aero_state, gas_data, &
       gas_state, run_part_opt, camp_core, photolysis)
#else
  subroutine run_part(scenario, env_state, aero_data, aero_state, gas_data, &
       gas_state, run_part_opt)
#endif
    !> Environment state.
    type(scenario_t), intent(in) :: scenario
    !> Environment state.
    type(env_state_t), intent(inout) :: env_state
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Aerosol state.
    type(aero_state_t), intent(inout) :: aero_state
    !> Gas data.
    type(gas_data_t), intent(in) :: gas_data
    !> Gas state.
    type(gas_state_t), intent(inout) :: gas_state
    !> Monte Carlo options.
    type(run_part_opt_t), intent(in) :: run_part_opt
#ifdef PMC_USE_CAMP
    !> CAMP chemistry core
    type(camp_core_t), pointer, intent(inout), optional :: camp_core
    !> Photolysis calculator
    type(photolysis_t), pointer, intent(inout), optional :: photolysis

    type(camp_state_t), pointer :: camp_state
    type(camp_state_t), pointer :: camp_pre_aero_state, camp_post_aero_state
#endif

    real(kind=dp) :: time, pre_time, pre_del_t, prop_done
    real(kind=dp) :: last_output_time, last_progress_time
    integer :: rank, n_proc, pre_index, ncid
    integer :: pre_i_repeat
    integer :: n_samp, n_coag, n_emit, n_dil_in, n_dil_out, n_nuc
    integer :: progress_n_samp, progress_n_coag
    integer :: progress_n_emit, progress_n_dil_in, progress_n_dil_out
    integer :: progress_n_nuc, n_part_before
    integer :: global_n_part, global_n_samp, global_n_coag
    integer :: global_n_emit, global_n_dil_in, global_n_dil_out
    integer :: global_n_nuc
    logical :: do_output, do_state, do_state_netcdf, do_progress, did_coag
    real(kind=dp) :: t_start, t_wall_now, t_wall_elapsed, t_wall_remain
    type(env_state_t) :: old_env_state
    integer :: n_time, i_time, i_time_start, pre_i_time
    integer :: i_state, i_state_netcdf, i_output
#ifdef PMC_USE_CAMP
    type(camp_env_state_t) :: camp_env_state
#endif

    rank = pmc_mpi_rank()
    n_proc = pmc_mpi_size()

    i_time = 0
    i_output = 1
    i_state = 1
    i_state_netcdf = 1
    time = 0d0
    progress_n_samp = 0
    progress_n_coag = 0
    progress_n_emit = 0
    progress_n_dil_in = 0
    progress_n_dil_out = 0
    progress_n_nuc = 0

    call check_time_multiple("t_max", run_part_opt%t_max, &
         "del_t", run_part_opt%del_t)
    call check_time_multiple("t_output", run_part_opt%t_output, &
         "del_t", run_part_opt%del_t)
    call check_time_multiple("t_progress", run_part_opt%t_progress, &
         "del_t", run_part_opt%del_t)
#ifdef PMC_USE_CAMP
    if (run_part_opt%do_camp_chem) then
       camp_env_state%temp = env_state%temp
       camp_env_state%pressure = env_state%pressure
       camp_state =>  camp_core%new_state_one_cell(camp_env_state)
       camp_pre_aero_state => camp_core%new_state_one_cell(camp_env_state)
       camp_post_aero_state => camp_core%new_state_one_cell(camp_env_state)
    end if
#endif

    if (run_part_opt%do_mosaic) then
       call mosaic_init(env_state, aero_data, run_part_opt%del_t, &
            run_part_opt%do_optical)
       if (run_part_opt%do_optical) then
          call mosaic_aero_optical_init(env_state, aero_data, &
            aero_state, gas_data, gas_state)
       end if
    end if

    if (run_part_opt%t_output > 0d0) then
       call output_state(run_part_opt%output_prefix, &
            run_part_opt%output_type, aero_data, aero_state, gas_data, &
            gas_state, env_state, i_state, time, run_part_opt%del_t, &
            run_part_opt%i_repeat, run_part_opt%record_removals, &
            run_part_opt%do_optical, run_part_opt%uuid)
       call aero_info_array_zero(aero_state%aero_info_array)
    end if

    call aero_state_rebalance(aero_state, aero_data, &
         run_part_opt%allow_doubling, &
         run_part_opt%allow_halving, initial_state_warning=.true.)

    t_start = env_state%elapsed_time
    last_output_time = time
    last_progress_time = time
    n_time = nint(run_part_opt%t_max / run_part_opt%del_t)
    i_time_start = nint(time / run_part_opt%del_t) + 1

    global_n_part = aero_state_total_particles_all_procs(aero_state)
    if (rank == 0) then
       ! progress only printed from root process
       if (run_part_opt%i_repeat == 1) then
          t_wall_elapsed = 0d0
          t_wall_remain = 0d0
       else
          call cpu_time(t_wall_now)
          prop_done = real(run_part_opt%i_repeat - 1, kind=dp) &
               / real(run_part_opt%n_repeat, kind=dp)
          t_wall_elapsed = t_wall_now - run_part_opt%t_wall_start
          t_wall_remain = (1d0 - prop_done) / prop_done &
               * t_wall_elapsed
       end if
       call print_part_progress(run_part_opt%i_repeat, time, &
            global_n_part, 0, 0, 0, 0, 0, t_wall_elapsed, t_wall_remain)
    end if

    do i_time = i_time_start,n_time

       time = real(i_time, kind=dp) * run_part_opt%del_t

       old_env_state = env_state
       call scenario_update_env_state(scenario, env_state, time + t_start)

       if (run_part_opt%do_nucleation) then
          n_part_before = aero_state_total_particles(aero_state)
          call nucleate(run_part_opt%nucleate_type, &
               run_part_opt%nucleate_source, env_state, gas_data, aero_data, &
               aero_state, gas_state, run_part_opt%del_t, &
               run_part_opt%allow_doubling, run_part_opt%allow_halving)
          n_nuc = aero_state_total_particles(aero_state) &
               - n_part_before
          progress_n_nuc = progress_n_nuc + n_nuc
       end if

       if (run_part_opt%do_coagulation) then
          if (run_part_opt%parallel_coag_type &
               == PARALLEL_COAG_TYPE_LOCAL) then
             call mc_coag(run_part_opt%coag_kernel_type, env_state, &
                  aero_data, aero_state, run_part_opt%del_t, n_samp, n_coag)
          elseif (run_part_opt%parallel_coag_type &
               == PARALLEL_COAG_TYPE_DIST) then
             call mc_coag_dist(run_part_opt%coag_kernel_type, env_state, &
                  aero_data, aero_state, run_part_opt%del_t, n_samp, n_coag)
          else
             call die_msg(323011762, "unknown parallel coagulation type: " &
                  // trim(integer_to_string(run_part_opt%parallel_coag_type)))
          end if
          progress_n_samp = progress_n_samp + n_samp
          progress_n_coag = progress_n_coag + n_coag
       end if

#ifdef PMC_USE_SUNDIALS
       if (run_part_opt%do_condensation) then
          call condense_particles(aero_state, aero_data, old_env_state, &
               env_state, run_part_opt%del_t)
       end if
#endif

       call scenario_update_gas_state(scenario, run_part_opt%del_t, &
            env_state, old_env_state, gas_data, gas_state)
       call scenario_update_aero_state(scenario, run_part_opt%del_t, &
            env_state, old_env_state, aero_data, aero_state, n_emit, &
            n_dil_in, n_dil_out, run_part_opt%allow_doubling, &
            run_part_opt%allow_halving)
       progress_n_emit = progress_n_emit + n_emit
       progress_n_dil_in = progress_n_dil_in + n_dil_in
       progress_n_dil_out = progress_n_dil_out + n_dil_out

       if (run_part_opt%do_camp_chem) then
#ifdef PMC_USE_CAMP
          call pmc_camp_interface_solve(camp_core, camp_state, &
               camp_pre_aero_state, camp_post_aero_state, env_state, &
               aero_data, aero_state, gas_data, gas_state, photolysis, &
               run_part_opt%del_t)
#endif
       end if

       if (run_part_opt%do_mosaic) then
          call mosaic_timestep(env_state, aero_data, aero_state, gas_data, &
               gas_state, run_part_opt%do_optical)
       end if

       if (run_part_opt%mix_timescale > 0d0) then
          call aero_state_mix(aero_state, run_part_opt%del_t, &
               run_part_opt%mix_timescale, aero_data)
       end if
       if (run_part_opt%gas_average) then
          call gas_state_mix(gas_state)
       end if
       if (run_part_opt%gas_average) then
          call env_state_mix(env_state)
       end if

       call aero_state_rebalance(aero_state, aero_data, &
            run_part_opt%allow_doubling, &
            run_part_opt%allow_halving, initial_state_warning=.false.)

       ! DEBUG: enable to check consistency
       ! call aero_state_check(aero_state, aero_data)
       ! DEBUG: end

       if (run_part_opt%t_output > 0d0) then
          call check_event(time, run_part_opt%del_t, run_part_opt%t_output, &
               last_output_time, do_output)
          if (do_output) then
             i_output = i_output + 1
             call output_state(run_part_opt%output_prefix, &
                  run_part_opt%output_type, aero_data, aero_state, gas_data, &
                  gas_state, env_state, i_output, time, run_part_opt%del_t, &
                  run_part_opt%i_repeat, run_part_opt%record_removals, &
                  run_part_opt%do_optical, run_part_opt%uuid)
             call aero_info_array_zero(aero_state%aero_info_array)
          end if
       end if

       if (.not. run_part_opt%record_removals) then
          ! If we are not recording removals then we can zero them as
          ! often as possible to minimize the cost of maintaining
          ! them.
          call aero_info_array_zero(aero_state%aero_info_array)
       end if

       if (run_part_opt%t_progress > 0d0) then
          call check_event(time, run_part_opt%del_t, &
               run_part_opt%t_progress, last_progress_time, do_progress)
          if (do_progress) then
             global_n_part = aero_state_total_particles_all_procs(aero_state)
             call pmc_mpi_reduce_sum_integer(progress_n_samp, global_n_samp)
             call pmc_mpi_reduce_sum_integer(progress_n_coag, global_n_coag)
             call pmc_mpi_reduce_sum_integer(progress_n_emit, global_n_emit)
             call pmc_mpi_reduce_sum_integer(progress_n_dil_in, &
                  global_n_dil_in)
             call pmc_mpi_reduce_sum_integer(progress_n_dil_out, &
                  global_n_dil_out)
             call pmc_mpi_reduce_sum_integer(progress_n_nuc, global_n_nuc)
             if (rank == 0) then
                ! progress only printed from root process
                call cpu_time(t_wall_now)
                prop_done = (real(run_part_opt%i_repeat - 1, kind=dp) &
                     + time / run_part_opt%t_max) &
                     / real(run_part_opt%n_repeat, kind=dp)
                t_wall_elapsed = t_wall_now - run_part_opt%t_wall_start
                t_wall_remain = (1d0 - prop_done) / prop_done &
                     * t_wall_elapsed
                call print_part_progress(run_part_opt%i_repeat, time, &
                     global_n_part, global_n_coag, global_n_emit, &
                     global_n_dil_in, global_n_dil_out, global_n_nuc, &
                     t_wall_elapsed, t_wall_remain)
             end if
             ! reset counters so they show information since last
             ! progress display
             progress_n_samp = 0
             progress_n_coag = 0
             progress_n_emit = 0
             progress_n_dil_in = 0
             progress_n_dil_out = 0
             progress_n_nuc = 0
          end if
       end if

    end do

    if (run_part_opt%do_mosaic) then
       call mosaic_cleanup()
    end if

  end subroutine run_part

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Print the current simulation progress to the screen.
  subroutine print_part_progress(i_repeat, t_sim_elapsed, n_part, n_coag, &
       n_emit, n_dil_in, n_dil_out, n_nuc, t_wall_elapsed, t_wall_remain)

    !> Repeat number of simulation.
    integer, intent(in) :: i_repeat
    !> Elapsed simulation time (s).
    real(kind=dp), intent(in) :: t_sim_elapsed
    !> Number of particles.
    integer, intent(in) :: n_part
    !> Number of coagulated particles since last progress printing.
    integer, intent(in) :: n_coag
    !> Number of emitted particles since last progress printing.
    integer, intent(in) :: n_emit
    !> Number of diluted-in particles since last progress printing.
    integer, intent(in) :: n_dil_in
    !> Number of diluted-out particles since last progress printing.
    integer, intent(in) :: n_dil_out
    !> Number of nucleated particles since last progress printing.
    integer, intent(in) :: n_nuc
    !> Elapsed wall time (s).
    real(kind=dp), intent(in) :: t_wall_elapsed
    !> Estimated remaining wall time (s).
    real(kind=dp), intent(in) :: t_wall_remain

    write(*,'(a6,a1,a6,a1,a7,a1,a7,a1,a7,a1,a8,a1,a9,a1,a7,a1,a6,a1,a6)') &
         "repeat", " ", "t_sim", " ", "n_part", " ", "n_coag", " ", &
         "n_emit", " ", "n_dil_in", " ", "n_dil_out", " ", "n_nuc", " ", &
         "t_wall", " ", "t_rem"
    write(*,'(a6,a1,a6,a1,a7,a1,a7,a1,a7,a1,a8,a1,a9,a1,a7,a1,a6,a1,a6)') &
         trim(integer_to_string_max_len(i_repeat, 6)), " ", &
         trim(time_to_string_max_len(t_sim_elapsed, 6)), " ", &
         trim(integer_to_string_max_len(n_part, 7)), " ", &
         trim(integer_to_string_max_len(n_coag, 7)), " ", &
         trim(integer_to_string_max_len(n_emit, 7)), " ", &
         trim(integer_to_string_max_len(n_dil_in, 7)), " ", &
         trim(integer_to_string_max_len(n_dil_out, 7)), " ", &
         trim(integer_to_string_max_len(n_nuc, 7)), " ", &
         trim(time_to_string_max_len(t_wall_elapsed, 6)), " ", &
         trim(time_to_string_max_len(t_wall_remain, 6))

  end subroutine print_part_progress

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Determines the number of bytes required to pack the given value.
  integer function pmc_mpi_pack_size_run_part_opt(val)

    !> Value to pack.
    type(run_part_opt_t), intent(in) :: val

    pmc_mpi_pack_size_run_part_opt = &
         pmc_mpi_pack_size_real(val%t_max) &
         + pmc_mpi_pack_size_real(val%t_output) &
         + pmc_mpi_pack_size_real(val%t_progress) &
         + pmc_mpi_pack_size_real(val%del_t) &
         + pmc_mpi_pack_size_string(val%output_prefix) &
         + pmc_mpi_pack_size_integer(val%coag_kernel_type) &
         + pmc_mpi_pack_size_integer(val%nucleate_type) &
         + pmc_mpi_pack_size_integer(val%nucleate_source) &
         + pmc_mpi_pack_size_logical(val%do_coagulation) &
         + pmc_mpi_pack_size_logical(val%do_nucleation) &
         + pmc_mpi_pack_size_logical(val%allow_doubling) &
         + pmc_mpi_pack_size_logical(val%allow_halving) &
         + pmc_mpi_pack_size_logical(val%do_condensation) &
         + pmc_mpi_pack_size_logical(val%do_mosaic) &
         + pmc_mpi_pack_size_logical(val%do_optical) &
         + pmc_mpi_pack_size_logical(val%do_select_weighting) &
         + pmc_mpi_pack_size_integer(val%weighting_type) &
         + pmc_mpi_pack_size_real(val%weighting_exponent) &
         + pmc_mpi_pack_size_integer(val%i_repeat) &
         + pmc_mpi_pack_size_integer(val%n_repeat) &
         + pmc_mpi_pack_size_real(val%t_wall_start) &
         + pmc_mpi_pack_size_logical(val%record_removals) &
         + pmc_mpi_pack_size_logical(val%do_parallel) &
         + pmc_mpi_pack_size_integer(val%output_type) &
         + pmc_mpi_pack_size_real(val%mix_timescale) &
         + pmc_mpi_pack_size_logical(val%gas_average) &
         + pmc_mpi_pack_size_logical(val%env_average) &
         + pmc_mpi_pack_size_integer(val%parallel_coag_type) &
         + pmc_mpi_pack_size_string(val%uuid)

  end function pmc_mpi_pack_size_run_part_opt

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Packs the given value into the buffer, advancing position.
  subroutine pmc_mpi_pack_run_part_opt(buffer, position, val)

    !> Memory buffer.
    character, intent(inout) :: buffer(:)
    !> Current buffer position.
    integer, intent(inout) :: position
    !> Value to pack.
    type(run_part_opt_t), intent(in) :: val

#ifdef PMC_USE_MPI
    integer :: prev_position

    prev_position = position
    call pmc_mpi_pack_real(buffer, position, val%t_max)
    call pmc_mpi_pack_real(buffer, position, val%t_output)
    call pmc_mpi_pack_real(buffer, position, val%t_progress)
    call pmc_mpi_pack_real(buffer, position, val%del_t)
    call pmc_mpi_pack_string(buffer, position, val%output_prefix)
    call pmc_mpi_pack_integer(buffer, position, val%coag_kernel_type)
    call pmc_mpi_pack_integer(buffer, position, val%nucleate_type)
    call pmc_mpi_pack_integer(buffer, position, val%nucleate_source)
    call pmc_mpi_pack_logical(buffer, position, val%do_coagulation)
    call pmc_mpi_pack_logical(buffer, position, val%do_nucleation)
    call pmc_mpi_pack_logical(buffer, position, val%allow_doubling)
    call pmc_mpi_pack_logical(buffer, position, val%allow_halving)
    call pmc_mpi_pack_logical(buffer, position, val%do_condensation)
    call pmc_mpi_pack_logical(buffer, position, val%do_mosaic)
    call pmc_mpi_pack_logical(buffer, position, val%do_optical)
    call pmc_mpi_pack_logical(buffer, position, val%do_select_weighting)
    call pmc_mpi_pack_integer(buffer, position, val%weighting_type)
    call pmc_mpi_pack_real(buffer, position, val%weighting_exponent)
    call pmc_mpi_pack_integer(buffer, position, val%i_repeat)
    call pmc_mpi_pack_integer(buffer, position, val%n_repeat)
    call pmc_mpi_pack_real(buffer, position, val%t_wall_start)
    call pmc_mpi_pack_logical(buffer, position, val%record_removals)
    call pmc_mpi_pack_logical(buffer, position, val%do_parallel)
    call pmc_mpi_pack_integer(buffer, position, val%output_type)
    call pmc_mpi_pack_real(buffer, position, val%mix_timescale)
    call pmc_mpi_pack_logical(buffer, position, val%gas_average)
    call pmc_mpi_pack_logical(buffer, position, val%env_average)
    call pmc_mpi_pack_integer(buffer, position, val%parallel_coag_type)
    call pmc_mpi_pack_string(buffer, position, val%uuid)
    call assert(946070052, &
         position - prev_position <= pmc_mpi_pack_size_run_part_opt(val))
#endif

  end subroutine pmc_mpi_pack_run_part_opt

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Unpacks the given value from the buffer, advancing position.
  subroutine pmc_mpi_unpack_run_part_opt(buffer, position, val)

    !> Memory buffer.
    character, intent(inout) :: buffer(:)
    !> Current buffer position.
    integer, intent(inout) :: position
    !> Value to pack.
    type(run_part_opt_t), intent(inout) :: val

#ifdef PMC_USE_MPI
    integer :: prev_position

    prev_position = position
    call pmc_mpi_unpack_real(buffer, position, val%t_max)
    call pmc_mpi_unpack_real(buffer, position, val%t_output)
    call pmc_mpi_unpack_real(buffer, position, val%t_progress)
    call pmc_mpi_unpack_real(buffer, position, val%del_t)
    call pmc_mpi_unpack_string(buffer, position, val%output_prefix)
    call pmc_mpi_unpack_integer(buffer, position, val%coag_kernel_type)
    call pmc_mpi_unpack_integer(buffer, position, val%nucleate_type)
    call pmc_mpi_unpack_integer(buffer, position, val%nucleate_source)
    call pmc_mpi_unpack_logical(buffer, position, val%do_coagulation)
    call pmc_mpi_unpack_logical(buffer, position, val%do_nucleation)
    call pmc_mpi_unpack_logical(buffer, position, val%allow_doubling)
    call pmc_mpi_unpack_logical(buffer, position, val%allow_halving)
    call pmc_mpi_unpack_logical(buffer, position, val%do_condensation)
    call pmc_mpi_unpack_logical(buffer, position, val%do_mosaic)
    call pmc_mpi_unpack_logical(buffer, position, val%do_optical)
    call pmc_mpi_unpack_logical(buffer, position, val%do_select_weighting)
    call pmc_mpi_unpack_integer(buffer, position, val%weighting_type)
    call pmc_mpi_unpack_real(buffer, position, val%weighting_exponent)
    call pmc_mpi_unpack_integer(buffer, position, val%i_repeat)
    call pmc_mpi_unpack_integer(buffer, position, val%n_repeat)
    call pmc_mpi_unpack_real(buffer, position, val%t_wall_start)
    call pmc_mpi_unpack_logical(buffer, position, val%record_removals)
    call pmc_mpi_unpack_logical(buffer, position, val%do_parallel)
    call pmc_mpi_unpack_integer(buffer, position, val%output_type)
    call pmc_mpi_unpack_real(buffer, position, val%mix_timescale)
    call pmc_mpi_unpack_logical(buffer, position, val%gas_average)
    call pmc_mpi_unpack_logical(buffer, position, val%env_average)
    call pmc_mpi_unpack_integer(buffer, position, val%parallel_coag_type)
    call pmc_mpi_unpack_string(buffer, position, val%uuid)
    call assert(480118362, &
         position - prev_position <= pmc_mpi_pack_size_run_part_opt(val))
#endif

  end subroutine pmc_mpi_unpack_run_part_opt

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read the specification for a parallel coagulation type from a spec file.
  subroutine spec_file_read_parallel_coag_type(file, parallel_coag_type)

    !> Spec file.
    type(spec_file_t), intent(inout) :: file
    !> Kernel type.
    integer, intent(out) :: parallel_coag_type

    character(len=SPEC_LINE_MAX_VAR_LEN) :: parallel_coag_type_name

    !> \page input_format_parallel_coag Input File Format: Parallel Coagulation Type
    !!
    !! The output type is specified by the parameter:
    !!   - \b parallel_coag (string): type of parallel coagulation ---
    !!     must be one of: \c local for only within-process
    !!     coagulation or \c dist to have all processes perform
    !!     coagulation globally, requesting particles from other
    !!     processes as needed
    !!
    !! See also:
    !!   - \ref spec_file_format --- the input file text format

    call spec_file_read_string(file, 'parallel_coag', &
         parallel_coag_type_name)
    if (trim(parallel_coag_type_name) == 'local') then
       parallel_coag_type = PARALLEL_COAG_TYPE_LOCAL
    elseif (trim(parallel_coag_type_name) == 'dist') then
       parallel_coag_type = PARALLEL_COAG_TYPE_DIST
    else
       call spec_file_die_msg(494684716, file, &
            "Unknown parallel coagulation type: " &
            // trim(parallel_coag_type_name))
    end if

  end subroutine spec_file_read_parallel_coag_type

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_run_part
