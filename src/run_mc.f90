! Copyright (C) 2005-2007 Nicole Riemer and Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.
!
! Monte Carlo with fixed timestep and particles stored per-bin.

module pmc_run_mc

  use pmc_inout
  use pmc_process_spec

  type run_mc_opt_t
    integer :: n_part_max               ! maximum number of particles
    real*8 :: t_max                     ! final time (s)
    real*8 :: t_output                  ! output interval (0 disables) (s)
    real*8 :: t_state                   ! state output interval (0 disables) (s)
    real*8 :: t_progress                ! progress interval (0 disables) (s)
    real*8 :: del_t                     ! timestep for coagulation
    character(len=300) :: state_prefix  ! prefix for state files
    logical :: do_coagulation           ! whether to do coagulation
    logical :: allow_double             ! allow doubling if needed
    logical :: do_condensation          ! whether to do condensation
    logical :: do_mosaic                ! whether to do MOSAIC
    logical :: do_restart               ! whether to restart from state
    character(len=300) :: restart_name  ! name of state to restart from
    integer :: i_loop                   ! loop number of run
    integer :: n_loop                   ! total number of loops
    real*8 :: t_wall_start              ! cpu_time() of start
    real*8 :: mix_rate                  ! mix rate for parallel states (0 to 1)
 end type run_mc_opt_t
  
contains
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine run_mc(kernel, bin_grid, aero_binned, env_data, env, &
       aero_data, aero_state, gas_data, gas_state, mc_opt, summary_file, &
       process_spec_list)

    ! Do a particle-resolved Monte Carlo simulation.
    
    use pmc_util
    use pmc_aero_state
    use pmc_bin_grid 
    use pmc_aero_binned
    use pmc_condensation
    use pmc_env_data
    use pmc_env
    use pmc_aero_data
    use pmc_gas_data
    use pmc_gas_state
    use pmc_output_state
    use pmc_mosaic
    use pmc_coagulation
    use pmc_kernel
    use pmc_output_summary
    use pmc_mpi

    type(bin_grid_t), intent(in) :: bin_grid ! bin grid
    type(aero_binned_t), intent(out) :: aero_binned ! binned distributions
    type(env_data_t), intent(in) :: env_data ! environment state
    type(env_t), intent(inout) :: env   ! environment state
    type(aero_data_t), intent(in) :: aero_data ! aerosol data
    type(aero_state_t), intent(inout) :: aero_state ! aerosol state
    type(gas_data_t), intent(in) :: gas_data ! gas data
    type(gas_state_t), intent(inout) :: gas_state ! gas state
    type(run_mc_opt_t), intent(in) :: mc_opt ! Monte Carlo options
    type(inout_file_t), intent(inout) :: summary_file ! summary output file
    type(process_spec_t), intent(in) :: process_spec_list(:) ! processing spec

    ! FIXME: can we shift this to a module? pmc_kernel presumably
    interface
       subroutine kernel(v1, v2, env, k)
         use pmc_env
         real*8, intent(in) :: v1
         real*8, intent(in) :: v2
         type(env_t), intent(in) :: env   
         real*8, intent(out) :: k
       end subroutine kernel
    end interface
    
    real*8 time, pre_time
    real*8 last_output_time, last_state_time, last_progress_time
    real*8 k_max(bin_grid%n_bin, bin_grid%n_bin)
    integer n_coag, tot_n_samp, tot_n_coag, rank, pre_index
    logical do_output, do_state, do_progress, did_coag
    real*8 t_start, t_wall_now, t_wall_est, prop_done, old_height
    integer n_time, i_time, i_time_start, pre_i_time, i_state, i_summary
    character*100 filename
    type(bin_grid_t) :: restart_bin_grid
    type(aero_data_t) :: restart_aero_data
    type(gas_data_t) :: restart_gas_data

    rank = pmc_mpi_rank() ! MPI process rank (0 is root)

    i_time = 0
    i_summary = 0
    i_state = 0
    time = 0d0
    call env_data_init_state(env_data, env, time)
    n_coag = 0
    tot_n_samp = 0
    tot_n_coag = 0
    
    if (mc_opt%do_restart) then
#ifdef PMC_USE_MPI
       if (rank == 0) then
          write(0,*) 'ERROR: restarting not currently supported with MPI'
          call pmc_mpi_abort(1)
       end if
#endif
       call inout_read_state(mc_opt%restart_name, restart_bin_grid, &
            restart_aero_data, aero_state, restart_gas_data, gas_state, &
            env, time, pre_index)
       ! FIXME: should we check whether bin_grid == restart_bin_grid, etc?
       i_time = nint(time / mc_opt%del_t)
       if (mc_opt%allow_double) then
          do while (total_particles(aero_state) .lt. mc_opt%n_part_max / 2)
             call aero_state_double(aero_state)
          end do
       end if
       ! write data into output file so that it will look correct
       pre_time = 0d0
       last_output_time = pre_time
       call aero_state_to_binned(bin_grid, aero_data, aero_state, aero_binned)
       do pre_i_time = 0,(i_time - 1)
          call env_data_update_state(env_data, env, pre_time)
          if (mc_opt%t_output > 0d0) then
             call check_event(pre_time, mc_opt%del_t, mc_opt%t_output, &
                  last_output_time, do_output)
             if (do_output) call output_summary(summary_file, &
                  pre_time, bin_grid, aero_data, aero_binned, &
                  gas_data, gas_state, env, mc_opt%i_loop)
          end if
          pre_time = pre_time + mc_opt%del_t
       end do
    end if

    call aero_state_to_binned(bin_grid, aero_data, aero_state, aero_binned)
    
    call est_k_max_binned(bin_grid, kernel, env, k_max)

    if (mc_opt%do_mosaic) then
       call mosaic_init(bin_grid, env, mc_opt%del_t)
    end if

    if (mc_opt%t_output > 0d0) then
       call output_summary(summary_file, &
            time, bin_grid, aero_data, aero_binned, &
            gas_data, gas_state, env, mc_opt%i_loop)
       call output_processed(mc_opt%state_prefix, process_spec_list, &
            bin_grid, aero_data, aero_state, gas_data, gas_state, &
            env, i_summary, time, mc_opt%i_loop)
    end if

    if (mc_opt%t_state > 0d0) then
       call inout_write_state(mc_opt%state_prefix, bin_grid, &
            aero_data, aero_state, gas_data, gas_state, env, i_state, &
            time, mc_opt%i_loop)
    end if

    t_start = time
    last_progress_time = time
    last_state_time = time
    last_output_time = time
    n_time = nint(mc_opt%t_max / mc_opt%del_t)
    i_time_start = nint(time / mc_opt%del_t) + 1
    do i_time = i_time_start,n_time

       time = dble(i_time) * mc_opt%del_t

       old_height = env%height
       call env_data_update_state(env_data, env, time)
       call env_update_gas_state(env, mc_opt%del_t, old_height, gas_data, &
            gas_state)
       call env_update_aero_state(env, mc_opt%del_t, old_height, bin_grid, &
            aero_data, aero_state, aero_binned)

       if (mc_opt%do_coagulation) then
          call mc_coag(kernel, bin_grid, aero_binned, env, aero_data, &
               aero_state, mc_opt, k_max, tot_n_samp, n_coag)
       end if
       tot_n_coag = tot_n_coag + n_coag

       if (mc_opt%do_condensation) then
          call condense_particles(bin_grid, aero_binned, env, aero_data, &
               aero_state, mc_opt%del_t)
       end if

       if (mc_opt%do_mosaic) then
          call mosaic_timestep(bin_grid, env, aero_data, &
               aero_state, aero_binned, gas_data, gas_state, time)
       end if

       call mc_mix(aero_data, aero_state, gas_data, gas_state, &
            aero_binned, env, bin_grid, mc_opt%mix_rate)
       
       ! if we have less than half the maximum number of particles then
       ! double until we fill up the array, and the same for halving
       if (mc_opt%allow_double) then
          do while (total_particles(aero_state) < mc_opt%n_part_max / 2)
             call aero_state_double(aero_state)
          end do
          do while (total_particles(aero_state) > mc_opt%n_part_max * 2)
             call aero_state_halve(aero_state, aero_binned, bin_grid)
          end do
       end if
    
       ! DEBUG: enable to check array handling
       ! call aero_state_check(bin_grid, aero_binned, aero_data, aero_state)
       ! DEBUG: end
       
       if (mc_opt%t_output > 0d0) then
          call check_event(time, mc_opt%del_t, mc_opt%t_output, &
               last_output_time, do_output)
          if (do_output) then
             i_summary = i_summary + 1
             call output_summary(summary_file, &
                  time, bin_grid, aero_data, aero_binned, &
                  gas_data, gas_state, env, mc_opt%i_loop)
             call output_processed(mc_opt%state_prefix, process_spec_list, &
                  bin_grid, aero_data, aero_state, gas_data, gas_state, &
                  env, i_summary, time, mc_opt%i_loop)
          end if
       end if

       if (mc_opt%t_state > 0d0) then
          call check_event(time, mc_opt%del_t, mc_opt%t_state, &
               last_state_time, do_state)
          if (do_state) then
             i_state = i_state + 1
             call inout_write_state(mc_opt%state_prefix, bin_grid, &
                  aero_data, aero_state, gas_data, gas_state, env, i_state, &
                  time, mc_opt%i_loop)
          end if
       end if

       if (mc_opt%t_progress > 0d0) then
          if (rank == 0) then
             ! progress only printed from root process
             call check_event(time, mc_opt%del_t, mc_opt%t_progress, &
                  last_progress_time, do_progress)
             if (do_progress) then
                call cpu_time(t_wall_now)
                prop_done = (dble(mc_opt%i_loop - 1) + (time - t_start) &
                     / (mc_opt%t_max - t_start)) / dble(mc_opt%n_loop)
                t_wall_est = (1d0 - prop_done) / prop_done &
                     * (t_wall_now - mc_opt%t_wall_start)
                write(6,'(a6,a8,a9,a11,a9,a11,a10)') 'loop', 'time', &
                     'n_part', 'tot_n_samp', 'n_coag', 'tot_n_coag', 't_est'
                write(6,'(i6,f8.1,i9,i11,i9,i11,f10.0)') mc_opt%i_loop, time, &
                     total_particles(aero_state), tot_n_samp, n_coag, &
                     tot_n_coag, t_wall_est
             end if
          end if
       end if
       
    enddo
    
  end subroutine run_mc
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine mc_coag(kernel, bin_grid, aero_binned, env, aero_data, &
       aero_state, mc_opt, k_max, tot_n_samp, n_coag)

    ! Do coagulation for time del_t.

    use pmc_util
    use pmc_aero_state
    use pmc_bin_grid
    use pmc_aero_binned
    use pmc_env
    use pmc_aero_data
    use pmc_coagulation

    type(bin_grid_t), intent(in) :: bin_grid ! bin grid
    type(aero_binned_t), intent(out) :: aero_binned ! binned distributions
    type(env_t), intent(inout) :: env   ! environment state
    type(aero_data_t), intent(in) :: aero_data ! aerosol data
    type(aero_state_t), intent(inout) :: aero_state ! aerosol state
    type(run_mc_opt_t), intent(in) :: mc_opt ! Monte Carlo options
    real*8, intent(in) :: k_max(bin_grid%n_bin,bin_grid%n_bin) ! maximum kernel
    integer, intent(out) :: tot_n_samp  ! total number of samples tested
    integer, intent(out) :: n_coag      ! number of coagulation events

    interface
       subroutine kernel(v1, v2, env, k)
         use pmc_env
         real*8, intent(in) :: v1
         real*8, intent(in) :: v2
         type(env_t), intent(in) :: env   
         real*8, intent(out) :: k
       end subroutine kernel
    end interface
    
    logical did_coag
    integer i, j, n_samp, i_samp, M
    real*8 n_samp_real

    tot_n_samp = 0
    n_coag = 0
    do i = 1,bin_grid%n_bin
       do j = 1,bin_grid%n_bin
          call compute_n_samp(aero_state%bins(i)%n_part, &
               aero_state%bins(j)%n_part, i == j, k_max(i,j), &
               aero_state%comp_vol, mc_opt%del_t, n_samp_real)
          ! probabalistically determine n_samp to cope with < 1 case
          n_samp = prob_round(n_samp_real)
          tot_n_samp = tot_n_samp + n_samp
          do i_samp = 1,n_samp
             M = total_particles(aero_state)
             ! check we still have enough particles to coagulate
             if ((aero_state%bins(i)%n_part < 1) &
                  .or. (aero_state%bins(j)%n_part < 1) &
                  .or. ((i == j) .and. (aero_state%bins(i)%n_part < 2))) then
                exit
             end if
             call maybe_coag_pair(bin_grid, aero_binned, env, aero_data, &
                  aero_state, i, j, mc_opt%del_t, k_max(i,j), kernel, &
                  did_coag)
             if (did_coag) n_coag = n_coag + 1
          enddo
       enddo
    enddo

  end subroutine mc_coag

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine compute_n_samp(ni, nj, same_bin, k_max, comp_vol, &
       del_t, n_samp_real)
  
    ! Compute the number of samples required for the pair of bins.

    use pmc_env

    integer, intent(in) :: ni           ! number particles in first bin 
    integer, intent(in) :: nj           ! number particles in second bin
    logical, intent(in) :: same_bin     ! whether first bin is second bin
    real*8, intent(in) :: k_max         ! maximum kernel value
    real*8, intent(in) :: comp_vol      ! computational volume (m^3)
    real*8, intent(in) :: del_t         ! timestep (s)
    real*8, intent(out) :: n_samp_real  ! number of samples per timestep
    
    real*8 r_samp
    real*8 n_possible ! use real*8 to avoid integer overflow
    ! FIXME: should use integer*8 or integer(kind = 8)
    
    if (same_bin) then
       n_possible = dble(ni) * (dble(nj) - 1d0) / 2d0
    else
       n_possible = dble(ni) * dble(nj) / 2d0
    endif
    
    r_samp = k_max / comp_vol * del_t
    n_samp_real = r_samp * n_possible
    
  end subroutine compute_n_samp
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine mc_mix(aero_data, aero_state, gas_data, gas_state, &
       aero_binned, env, bin_grid, mix_rate)

    ! Mix data between processes.

    use pmc_util
    use pmc_aero_data
    use pmc_aero_state
    use pmc_gas_data
    use pmc_gas_state
    use pmc_aero_binned
    use pmc_bin_grid
    use pmc_env

    type(aero_data_t), intent(in) :: aero_data ! aerosol data
    type(aero_state_t), intent(inout) :: aero_state ! aerosol state
    type(gas_data_t), intent(in) :: gas_data ! gas data
    type(gas_state_t), intent(inout) :: gas_state ! gas state
    type(aero_binned_t), intent(inout) :: aero_binned ! binned aerosol data
    type(env_t), intent(inout) :: env   ! environment
    type(bin_grid_t), intent(in) :: bin_grid ! bin grid
    real*8, intent(in) :: mix_rate      ! amount to mix (0 to 1)

    call assert(173605827, (mix_rate >= 0d0) .and. (mix_rate <= 1d0))
    if (mix_rate == 0d0) return

    call aero_state_mix(aero_state, mix_rate, &
         aero_binned, aero_data, bin_grid)
    call gas_state_mix(gas_state)
    call env_mix(env)
    
  end subroutine mc_mix

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function pmc_mpi_pack_size_mc_opt(val)

    ! Determines the number of bytes required to pack the given value.

    use pmc_mpi

    type(run_mc_opt_t), intent(in) :: val ! value to pack

    pmc_mpi_pack_size_mc_opt = &
         pmc_mpi_pack_size_integer(val%n_part_max) &
         + pmc_mpi_pack_size_real(val%t_max) &
         + pmc_mpi_pack_size_real(val%t_output) &
         + pmc_mpi_pack_size_real(val%t_state) &
         + pmc_mpi_pack_size_real(val%t_progress) &
         + pmc_mpi_pack_size_real(val%del_t) &
         + pmc_mpi_pack_size_string(val%state_prefix) &
         + pmc_mpi_pack_size_logical(val%do_coagulation) &
         + pmc_mpi_pack_size_logical(val%allow_double) &
         + pmc_mpi_pack_size_logical(val%do_condensation) &
         + pmc_mpi_pack_size_logical(val%do_mosaic) &
         + pmc_mpi_pack_size_logical(val%do_restart) &
         + pmc_mpi_pack_size_string(val%restart_name) &
         + pmc_mpi_pack_size_integer(val%i_loop) &
         + pmc_mpi_pack_size_integer(val%n_loop) &
         + pmc_mpi_pack_size_real(val%t_wall_start) &
         + pmc_mpi_pack_size_real(val%mix_rate)

  end function pmc_mpi_pack_size_mc_opt

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine pmc_mpi_pack_mc_opt(buffer, position, val)

    ! Packs the given value into the buffer, advancing position.

#ifdef PMC_USE_MPI
    use mpi
    use pmc_mpi
    use pmc_util
#endif

    character, intent(inout) :: buffer(:) ! memory buffer
    integer, intent(inout) :: position  ! current buffer position
    type(run_mc_opt_t), intent(in) :: val ! value to pack

#ifdef PMC_USE_MPI
    integer :: prev_position

    prev_position = position
    call pmc_mpi_pack_integer(buffer, position, val%n_part_max)
    call pmc_mpi_pack_real(buffer, position, val%t_max)
    call pmc_mpi_pack_real(buffer, position, val%t_output)
    call pmc_mpi_pack_real(buffer, position, val%t_state)
    call pmc_mpi_pack_real(buffer, position, val%t_progress)
    call pmc_mpi_pack_real(buffer, position, val%del_t)
    call pmc_mpi_pack_string(buffer, position, val%state_prefix)
    call pmc_mpi_pack_logical(buffer, position, val%do_coagulation)
    call pmc_mpi_pack_logical(buffer, position, val%allow_double)
    call pmc_mpi_pack_logical(buffer, position, val%do_condensation)
    call pmc_mpi_pack_logical(buffer, position, val%do_mosaic)
    call pmc_mpi_pack_logical(buffer, position, val%do_restart)
    call pmc_mpi_pack_string(buffer, position, val%restart_name)
    call pmc_mpi_pack_integer(buffer, position, val%i_loop)
    call pmc_mpi_pack_integer(buffer, position, val%n_loop)
    call pmc_mpi_pack_real(buffer, position, val%t_wall_start)
    call pmc_mpi_pack_real(buffer, position, val%mix_rate)
    call assert(946070052, position - prev_position == pmc_mpi_pack_size_mc_opt(val))
#endif

  end subroutine pmc_mpi_pack_mc_opt

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine pmc_mpi_unpack_mc_opt(buffer, position, val)

    ! Unpacks the given value from the buffer, advancing position.

#ifdef PMC_USE_MPI
    use mpi
    use pmc_mpi
    use pmc_util
#endif

    character, intent(inout) :: buffer(:) ! memory buffer
    integer, intent(inout) :: position  ! current buffer position
    type(run_mc_opt_t), intent(out) :: val ! value to pack

#ifdef PMC_USE_MPI
    integer :: prev_position

    prev_position = position
    call pmc_mpi_unpack_integer(buffer, position, val%n_part_max)
    call pmc_mpi_unpack_real(buffer, position, val%t_max)
    call pmc_mpi_unpack_real(buffer, position, val%t_output)
    call pmc_mpi_unpack_real(buffer, position, val%t_state)
    call pmc_mpi_unpack_real(buffer, position, val%t_progress)
    call pmc_mpi_unpack_real(buffer, position, val%del_t)
    call pmc_mpi_unpack_string(buffer, position, val%state_prefix)
    call pmc_mpi_unpack_logical(buffer, position, val%do_coagulation)
    call pmc_mpi_unpack_logical(buffer, position, val%allow_double)
    call pmc_mpi_unpack_logical(buffer, position, val%do_condensation)
    call pmc_mpi_unpack_logical(buffer, position, val%do_mosaic)
    call pmc_mpi_unpack_logical(buffer, position, val%do_restart)
    call pmc_mpi_unpack_string(buffer, position, val%restart_name)
    call pmc_mpi_unpack_integer(buffer, position, val%i_loop)
    call pmc_mpi_unpack_integer(buffer, position, val%n_loop)
    call pmc_mpi_unpack_real(buffer, position, val%t_wall_start)
    call pmc_mpi_unpack_real(buffer, position, val%mix_rate)
    call assert(480118362, position - prev_position == pmc_mpi_pack_size_mc_opt(val))
#endif

  end subroutine pmc_mpi_unpack_mc_opt

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
end module pmc_run_mc
