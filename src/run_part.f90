! Copyright (C) 2005-2010 Nicole Riemer and Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_run_part module.

!> Monte Carlo simulation.
module pmc_run_part

  use pmc_util
  use pmc_aero_state
  use pmc_bin_grid
  use pmc_env_data
  use pmc_env_state
  use pmc_aero_data
  use pmc_gas_data
  use pmc_gas_state
  use pmc_output
  use pmc_mosaic
  use pmc_coagulation
  use pmc_coagulation_mpi_controlled
  use pmc_coagulation_mpi_equal
  use pmc_kernel
  use pmc_nucleate
  use pmc_mpi
#ifdef PMC_USE_SUNDIALS
  use pmc_condense
#endif
#ifdef PMC_USE_MPI
  use mpi
#endif

  integer, parameter :: RUN_PART_OPT_CHAR_LEN = 100

  !> Options controlling the execution of run_part().
  type run_part_opt_t
     !> Maximum number of particles.
    integer :: n_part_max
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
    !> Type of nucleation.
    integer :: nucleate_type
    !> Whether to do coagulation.
    logical :: do_coagulation
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
    !> Loop number of run.
    integer :: i_loop
    !> Total number of loops.
    integer :: n_loop
    !> Cpu_time() of start.
    real(kind=dp) :: t_wall_start
    !> Whether to record particle removal information.
    logical :: record_removals
    !> Whether to run in parallel.
    logical :: do_parallel
    !> Parallel output type (central/dist/single).
    character(len=RUN_PART_OPT_CHAR_LEN) :: output_type
    !> Mixing timescale between processes (s).
    real(kind=dp) :: mix_timescale
    !> Whether to average gases each timestep.
    logical :: gas_average
    !> Whether to average environment each timestep.
    logical :: env_average
    !> Parallel coagulation method (local/collect/central/dist).
    character(len=RUN_PART_OPT_CHAR_LEN) :: coag_method
 end type run_part_opt_t
  
contains
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Do a particle-resolved Monte Carlo simulation.
  subroutine run_part(kernel, kernel_max, bin_grid, env_data, env_state, &
       aero_data, aero_weight, aero_state, gas_data, gas_state, part_opt)
    
    !> Bin grid.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Environment state.
    type(env_data_t), intent(in) :: env_data
    !> Environment state.
    type(env_state_t), intent(inout) :: env_state
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Aerosol weight.
    type(aero_weight_t), intent(in) :: aero_weight
    !> Aerosol state.
    type(aero_state_t), intent(inout) :: aero_state
    !> Gas data.
    type(gas_data_t), intent(in) :: gas_data
    !> Gas state.
    type(gas_state_t), intent(inout) :: gas_state
    !> Monte Carlo options.
    type(run_part_opt_t), intent(in) :: part_opt

#ifndef DOXYGEN_SKIP_DOC
    interface
       subroutine kernel(aero_particle_1, aero_particle_2, aero_data, &
            env_state, k)
         use pmc_aero_particle
         use pmc_aero_data
         use pmc_env_state
         type(aero_particle_t), intent(in) :: aero_particle_1
         type(aero_particle_t), intent(in) :: aero_particle_2
         type(aero_data_t), intent(in) :: aero_data
         type(env_state_t), intent(in) :: env_state  
         real(kind=dp), intent(out) :: k
       end subroutine kernel
       subroutine kernel_max(v1, v2, aero_data, env_state, k_max)
         use pmc_aero_data
         use pmc_env_state
         real(kind=dp), intent(in) :: v1
         real(kind=dp), intent(in) :: v2
         type(aero_data_t), intent(in) :: aero_data
         type(env_state_t), intent(in) :: env_state  
         real(kind=dp), intent(out) :: k_max
       end subroutine kernel_max
    end interface
#endif
    
    real(kind=dp) :: time, pre_time, pre_del_t
    real(kind=dp) :: last_output_time, last_progress_time
    real(kind=dp) :: k_max(bin_grid%n_bin, bin_grid%n_bin)
    integer :: tot_n_samp, tot_n_coag, rank, n_proc, pre_index, ncid, pre_i_loop
    integer :: progress_n_samp, progress_n_coag
    integer :: global_n_part, global_n_samp, global_n_coag
    logical :: do_output, do_state, do_state_netcdf, do_progress, did_coag
    logical :: update_rel_humid
    real(kind=dp) :: t_start, t_wall_now, t_wall_elapsed, t_wall_remain, prop_done
    type(env_state_t) :: old_env_state
    integer :: n_time, i_time, i_time_start, pre_i_time
    integer :: i_state, i_state_netcdf, i_output
    character*100 :: filename
  
    rank = pmc_mpi_rank()
    n_proc = pmc_mpi_size()

    i_time = 0
    i_output = 1
    i_state = 1
    i_state_netcdf = 1
    time = 0d0
    tot_n_samp = 0
    tot_n_coag = 0
    progress_n_samp = 0
    progress_n_coag = 0

    call env_state_allocate(old_env_state)

    call est_k_max_binned(bin_grid, kernel_max, aero_data, aero_weight, &
         env_state, k_max)

    if (part_opt%do_mosaic) then
       call mosaic_init(bin_grid, env_state, part_opt%del_t, &
            part_opt%do_optical)
    end if

    if (part_opt%t_output > 0d0) then
       call output_state(part_opt%output_prefix, part_opt%output_type, &
            bin_grid, aero_data, aero_weight, aero_state, gas_data, &
            gas_state, env_state, i_state, time, part_opt%del_t, &
            part_opt%i_loop, part_opt%record_removals, part_opt%do_optical)
       call aero_info_array_zero(aero_state%aero_info_array)
    end if
    
    ! Do an initial double/halve test. This shouldn't happen (except
    ! for restart with inconsistent number) so issue a warning.
    if (part_opt%allow_doubling) then
       do while ((aero_state_total_particles(aero_state) &
            < part_opt%n_part_max / 2) &
            .and. (aero_state_total_particles(aero_state) > 0))
          call warn_msg(716882783, "doubling particles in initial condition")
          call aero_state_double(aero_state)
       end do
    end if
    if (part_opt%allow_halving) then
       do while (aero_state_total_particles(aero_state) &
            > part_opt%n_part_max * 2)
          call warn_msg(661936373, "halving particles in initial condition")
          call aero_state_halve(aero_state, bin_grid)
       end do
    end if
    
    t_start = env_state%elapsed_time
    last_output_time = time
    last_progress_time = time
    n_time = nint(part_opt%t_max / part_opt%del_t)
    i_time_start = nint(time / part_opt%del_t) + 1
    do i_time = i_time_start,n_time

       time = real(i_time, kind=dp) * part_opt%del_t

       call env_state_copy(env_state, old_env_state)
       update_rel_humid = .not. part_opt%do_condensation
       call env_data_update_state(env_data, env_state, time + t_start, &
            update_rel_humid)
       call nucleate(part_opt%nucleate_type, bin_grid, env_state, gas_data, &
            aero_data, aero_weight, aero_state, gas_state, part_opt%del_t)

       if (part_opt%do_coagulation) then
          if (part_opt%coag_method == "local") then
             call mc_coag(kernel, bin_grid, env_state, aero_data, &
                  aero_weight, aero_state, part_opt, k_max, tot_n_samp, &
                  tot_n_coag)
          elseif (part_opt%coag_method == "collect") then
             call mc_coag_mpi_centralized(kernel, bin_grid, env_state, aero_data, &
                  aero_weight, aero_state, part_opt, k_max, tot_n_samp, tot_n_coag)
          elseif (part_opt%coag_method == "central") then
             call mc_coag_mpi_controlled(kernel, bin_grid, env_state, aero_data, &
                  aero_weight, aero_state, part_opt%del_t, k_max, tot_n_samp, &
                  tot_n_coag)
          elseif (part_opt%coag_method == "dist") then
             call mc_coag_mpi_equal(kernel, bin_grid, env_state, aero_data, &
                  aero_weight, aero_state, part_opt%del_t, k_max, tot_n_samp, &
                  tot_n_coag)
          else
             call die_msg(323011762, "unknown coag_method: " &
                  // trim(part_opt%coag_method))
          end if
          progress_n_samp = progress_n_samp + tot_n_samp
          progress_n_coag = progress_n_coag + tot_n_coag
       end if

       call env_state_update_gas_state(env_state, part_opt%del_t, &
            old_env_state, gas_data, gas_state)
       call env_state_update_aero_state(env_state, part_opt%del_t, &
            old_env_state, bin_grid, aero_data, aero_weight, aero_state)

#ifdef PMC_USE_SUNDIALS
       if (part_opt%do_condensation) then
          call condense_particles(bin_grid, env_state, &
               env_data, aero_data, aero_weight, aero_state, part_opt%del_t)
       end if
#endif

       if (part_opt%do_mosaic) then
          call mosaic_timestep(bin_grid, env_state, aero_data, &
               aero_weight, aero_state, gas_data, gas_state, &
               part_opt%do_optical)
       end if

       if (part_opt%mix_timescale > 0d0) then
          call aero_state_mix(aero_state, part_opt%del_t, &
               part_opt%mix_timescale, aero_data, bin_grid)
       end if
       if (part_opt%gas_average) then
          call gas_state_mix(gas_state)
       end if
       if (part_opt%gas_average) then
          call env_state_mix(env_state)
       end if

       ! if we have less than half the maximum number of particles then
       ! double until we fill up the array
       if (part_opt%allow_doubling) then
          do while ((aero_state_total_particles(aero_state) &
               < part_opt%n_part_max / 2) &
               .and. (aero_state_total_particles(aero_state) > 0))
             call aero_state_double(aero_state)
          end do
       end if
       ! same for halving if we have too many particles
       if (part_opt%allow_halving) then
          do while (aero_state_total_particles(aero_state) &
               > part_opt%n_part_max * 2)
             call aero_state_halve(aero_state, bin_grid)
          end do
       end if
    
       ! DEBUG: enable to check array handling
       ! call aero_state_check(bin_grid, aero_data, aero_state)
       ! DEBUG: end
       
       if (part_opt%t_output > 0d0) then
          call check_event(time, part_opt%del_t, part_opt%t_output, &
               last_output_time, do_output)
          if (do_output) then
             i_output = i_output + 1
             call output_state(part_opt%output_prefix, &
                  part_opt%output_type, bin_grid, aero_data, aero_weight, &
                  aero_state, gas_data, gas_state, env_state, i_output, &
                  time, part_opt%del_t, part_opt%i_loop, &
                  part_opt%record_removals, part_opt%do_optical)
             call aero_info_array_zero(aero_state%aero_info_array)
          end if
       end if

       if (.not. part_opt%record_removals) then
          ! If we are not recording removals then we can zero them as
          ! often as possible to minimize the cost of maintaining
          ! them.
          call aero_info_array_zero(aero_state%aero_info_array)
       end if

       if (part_opt%t_progress > 0d0) then
          call check_event(time, part_opt%del_t, part_opt%t_progress, &
               last_progress_time, do_progress)
          if (do_progress) then
             call pmc_mpi_reduce_sum_integer(&
                  aero_state_total_particles(aero_state), global_n_part)
             call pmc_mpi_reduce_sum_integer(progress_n_samp, global_n_samp)
             call pmc_mpi_reduce_sum_integer(progress_n_coag, global_n_coag)
             if (rank == 0) then
                ! progress only printed from root process
                call cpu_time(t_wall_now)
                prop_done = (real(part_opt%i_loop - 1, kind=dp) &
                     + time / part_opt%t_max) / real(part_opt%n_loop, kind=dp)
                t_wall_elapsed = t_wall_now - part_opt%t_wall_start
                t_wall_remain = (1d0 - prop_done) / prop_done &
                     * t_wall_elapsed
                write(*,'(a6,a9,a11,a12,a12,a13,a12)') 'loop', 'time(s)', &
                     'n_particle', 'n_samples', 'n_coagulate', &
                     't_elapsed(s)', 't_remain(s)'
                write(*,'(i6,f9.1,i11,i12,i12,f13.0,f12.0)') &
                     part_opt%i_loop, time, global_n_part, global_n_samp, &
                     global_n_coag, t_wall_elapsed, t_wall_remain
             end if
             ! reset counters so they show information since last
             ! progress display
             progress_n_samp = 0
             progress_n_coag = 0
          end if
       end if
       
    enddo

    if (part_opt%do_mosaic) then
       call mosaic_cleanup()
    end if

    call env_state_deallocate(old_env_state)

  end subroutine run_part
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Do coagulation for time del_t.
  subroutine mc_coag(kernel, bin_grid, env_state, aero_data, &
       aero_weight, aero_state, part_opt, k_max, tot_n_samp, tot_n_coag)

    !> Bin grid.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Environment state.
    type(env_state_t), intent(in) :: env_state
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Aerosol weight.
    type(aero_weight_t), intent(in) :: aero_weight
    !> Aerosol state.
    type(aero_state_t), intent(inout) :: aero_state
    !> Monte Carlo options.
    type(run_part_opt_t), intent(in) :: part_opt
    !> Maximum kernel.
    real(kind=dp), intent(in) :: k_max(bin_grid%n_bin,bin_grid%n_bin)
    !> Total number of samples tested.
    integer, intent(out) :: tot_n_samp
    !> Number of coagulation events.
    integer, intent(out) :: tot_n_coag

#ifndef DOXYGEN_SKIP_DOC
    interface
       subroutine kernel(aero_particle_1, aero_particle_2, aero_data, &
            env_state, k)
         use pmc_aero_particle
         use pmc_aero_data
         use pmc_env_state
         type(aero_particle_t), intent(in) :: aero_particle_1
         type(aero_particle_t), intent(in) :: aero_particle_2
         type(aero_data_t), intent(in) :: aero_data
         type(env_state_t), intent(in) :: env_state  
         real(kind=dp), intent(out) :: k
       end subroutine kernel
    end interface
#endif
    
    logical did_coag
    integer i, j, n_samp, i_samp
    real(kind=dp) n_samp_real

    tot_n_samp = 0
    tot_n_coag = 0
    do i = 1,bin_grid%n_bin
       do j = 1,bin_grid%n_bin
          call compute_n_samp(aero_state%bin(i)%n_part, &
               aero_state%bin(j)%n_part, i == j, k_max(i,j), &
               aero_state%comp_vol, part_opt%del_t, n_samp_real)
          ! probabalistically determine n_samp to cope with < 1 case
          n_samp = prob_round(n_samp_real)
          tot_n_samp = tot_n_samp + n_samp
          do i_samp = 1,n_samp
             ! check we still have enough particles to coagulate
             if ((aero_state%bin(i)%n_part < 1) &
                  .or. (aero_state%bin(j)%n_part < 1) &
                  .or. ((i == j) .and. (aero_state%bin(i)%n_part < 2))) then
                exit
             end if
             call maybe_coag_pair(bin_grid, env_state, &
                  aero_data, aero_weight, aero_state, i, j, part_opt%del_t, &
                  k_max(i,j), kernel, did_coag)
             if (did_coag) tot_n_coag = tot_n_coag + 1
          enddo
       enddo
    enddo

  end subroutine mc_coag

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Do coagulation for time del_t in parallel by centralizing on node 0.
  subroutine mc_coag_mpi_centralized(kernel, bin_grid, env_state, &
       aero_data, aero_weight, aero_state, part_opt, k_max, tot_n_samp, &
       tot_n_coag)

    !> Bin grid.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Environment state.
    type(env_state_t), intent(in) :: env_state
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Aerosol weight.
    type(aero_weight_t), intent(in) :: aero_weight
    !> Aerosol state.
    type(aero_state_t), intent(inout) :: aero_state
    !> Monte Carlo options.
    type(run_part_opt_t), intent(in) :: part_opt
    !> Maximum kernel.
    real(kind=dp), intent(in) :: k_max(bin_grid%n_bin,bin_grid%n_bin)
    !> Total number of samples tested.
    integer, intent(out) :: tot_n_samp
    !> Number of coagulation events.
    integer, intent(out) :: tot_n_coag

#ifndef DOXYGEN_SKIP_DOC
    interface
       subroutine kernel(aero_particle_1, aero_particle_2, aero_data, &
            env_state, k)
         use pmc_aero_particle
         use pmc_aero_data
         use pmc_env_state
         type(aero_particle_t), intent(in) :: aero_particle_1
         type(aero_particle_t), intent(in) :: aero_particle_2
         type(aero_data_t), intent(in) :: aero_data
         type(env_state_t), intent(in) :: env_state  
         real(kind=dp), intent(out) :: k
       end subroutine kernel
    end interface
#endif

    type(aero_state_t) :: aero_state_total
    
#ifdef PMC_USE_MPI

    call aero_state_allocate(aero_state_total)
    call aero_state_mpi_gather(aero_state, aero_state_total)
    if (pmc_mpi_rank() == 0) then
       call mc_coag(kernel, bin_grid, env_state, aero_data, &
            aero_weight, aero_state_total, part_opt, k_max, tot_n_samp, &
            tot_n_coag)
    end if
    call aero_state_mpi_scatter(aero_state_total, aero_state, aero_data)
    call aero_state_deallocate(aero_state_total)
    
#endif

  end subroutine mc_coag_mpi_centralized
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Determines the number of bytes required to pack the given value.
  integer function pmc_mpi_pack_size_part_opt(val)

    !> Value to pack.
    type(run_part_opt_t), intent(in) :: val

    pmc_mpi_pack_size_part_opt = &
         pmc_mpi_pack_size_integer(val%n_part_max) &
         + pmc_mpi_pack_size_real(val%t_max) &
         + pmc_mpi_pack_size_real(val%t_output) &
         + pmc_mpi_pack_size_real(val%t_progress) &
         + pmc_mpi_pack_size_real(val%del_t) &
         + pmc_mpi_pack_size_string(val%output_prefix) &
         + pmc_mpi_pack_size_logical(val%do_coagulation) &
         + pmc_mpi_pack_size_logical(val%allow_doubling) &
         + pmc_mpi_pack_size_logical(val%do_condensation) &
         + pmc_mpi_pack_size_logical(val%do_mosaic) &
         + pmc_mpi_pack_size_logical(val%do_optical) &
         + pmc_mpi_pack_size_integer(val%i_loop) &
         + pmc_mpi_pack_size_integer(val%n_loop) &
         + pmc_mpi_pack_size_real(val%t_wall_start) &
         + pmc_mpi_pack_size_logical(val%record_removals) &
         + pmc_mpi_pack_size_logical(val%do_parallel) &
         + pmc_mpi_pack_size_string(val%output_type) &
         + pmc_mpi_pack_size_real(val%mix_timescale) &
         + pmc_mpi_pack_size_logical(val%gas_average) &
         + pmc_mpi_pack_size_logical(val%env_average) &
         + pmc_mpi_pack_size_string(val%coag_method)

  end function pmc_mpi_pack_size_part_opt

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Packs the given value into the buffer, advancing position.
  subroutine pmc_mpi_pack_part_opt(buffer, position, val)

    !> Memory buffer.
    character, intent(inout) :: buffer(:)
    !> Current buffer position.
    integer, intent(inout) :: position
    !> Value to pack.
    type(run_part_opt_t), intent(in) :: val

#ifdef PMC_USE_MPI
    integer :: prev_position

    prev_position = position
    call pmc_mpi_pack_integer(buffer, position, val%n_part_max)
    call pmc_mpi_pack_real(buffer, position, val%t_max)
    call pmc_mpi_pack_real(buffer, position, val%t_output)
    call pmc_mpi_pack_real(buffer, position, val%t_progress)
    call pmc_mpi_pack_real(buffer, position, val%del_t)
    call pmc_mpi_pack_string(buffer, position, val%output_prefix)
    call pmc_mpi_pack_logical(buffer, position, val%do_coagulation)
    call pmc_mpi_pack_logical(buffer, position, val%allow_doubling)
    call pmc_mpi_pack_logical(buffer, position, val%do_condensation)
    call pmc_mpi_pack_logical(buffer, position, val%do_mosaic)
    call pmc_mpi_pack_logical(buffer, position, val%do_optical)
    call pmc_mpi_pack_integer(buffer, position, val%i_loop)
    call pmc_mpi_pack_integer(buffer, position, val%n_loop)
    call pmc_mpi_pack_real(buffer, position, val%t_wall_start)
    call pmc_mpi_pack_logical(buffer, position, val%record_removals)
    call pmc_mpi_pack_logical(buffer, position, val%do_parallel)
    call pmc_mpi_pack_string(buffer, position, val%output_type)
    call pmc_mpi_pack_real(buffer, position, val%mix_timescale)
    call pmc_mpi_pack_logical(buffer, position, val%gas_average)
    call pmc_mpi_pack_logical(buffer, position, val%env_average)
    call pmc_mpi_pack_string(buffer, position, val%coag_method)
    call assert(946070052, &
         position - prev_position == pmc_mpi_pack_size_part_opt(val))
#endif

  end subroutine pmc_mpi_pack_part_opt

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Unpacks the given value from the buffer, advancing position.
  subroutine pmc_mpi_unpack_part_opt(buffer, position, val)

    !> Memory buffer.
    character, intent(inout) :: buffer(:)
    !> Current buffer position.
    integer, intent(inout) :: position
    !> Value to pack.
    type(run_part_opt_t), intent(out) :: val

#ifdef PMC_USE_MPI
    integer :: prev_position

    prev_position = position
    call pmc_mpi_unpack_integer(buffer, position, val%n_part_max)
    call pmc_mpi_unpack_real(buffer, position, val%t_max)
    call pmc_mpi_unpack_real(buffer, position, val%t_output)
    call pmc_mpi_unpack_real(buffer, position, val%t_progress)
    call pmc_mpi_unpack_real(buffer, position, val%del_t)
    call pmc_mpi_unpack_string(buffer, position, val%output_prefix)
    call pmc_mpi_unpack_logical(buffer, position, val%do_coagulation)
    call pmc_mpi_unpack_logical(buffer, position, val%allow_doubling)
    call pmc_mpi_unpack_logical(buffer, position, val%do_condensation)
    call pmc_mpi_unpack_logical(buffer, position, val%do_mosaic)
    call pmc_mpi_unpack_logical(buffer, position, val%do_optical)
    call pmc_mpi_unpack_integer(buffer, position, val%i_loop)
    call pmc_mpi_unpack_integer(buffer, position, val%n_loop)
    call pmc_mpi_unpack_real(buffer, position, val%t_wall_start)
    call pmc_mpi_unpack_logical(buffer, position, val%record_removals)
    call pmc_mpi_unpack_logical(buffer, position, val%do_parallel)
    call pmc_mpi_unpack_string(buffer, position, val%output_type)
    call pmc_mpi_unpack_real(buffer, position, val%mix_timescale)
    call pmc_mpi_unpack_logical(buffer, position, val%gas_average)
    call pmc_mpi_unpack_logical(buffer, position, val%env_average)
    call pmc_mpi_unpack_string(buffer, position, val%coag_method)
    call assert(480118362, &
         position - prev_position == pmc_mpi_pack_size_part_opt(val))
#endif

  end subroutine pmc_mpi_unpack_part_opt

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
end module pmc_run_part
