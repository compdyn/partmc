! Copyright (C) 2005-2007 Nicole Riemer and Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.
!
! Monte Carlo with fixed timestep and particles stored per-bin.

module mod_run_mc

  use mod_inout

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
 end type run_mc_opt_t
  
contains
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine run_mc(kernel, bin_grid, aero_binned, env, aero_data, &
    aero_state, gas_data, gas_state, mc_opt, summary_file)

    ! Do a particle-resolved Monte Carlo simulation.
    
    use mod_util
    use mod_aero_state
    use mod_bin_grid 
    use mod_aero_binned
    use mod_condensation
    use mod_environ
    use mod_aero_data
    use mod_gas_data
    use mod_gas_state
    use mod_output_state
    use mod_mosaic
    use mod_coagulation
    use mod_kernel
    use mod_output_summary

    type(bin_grid_t), intent(in) :: bin_grid ! bin grid
    type(aero_binned_t), intent(out) :: aero_binned ! binned distributions
    type(environ), intent(inout) :: env ! environment state
    type(aero_data_t), intent(in) :: aero_data ! aerosol data
    type(aero_state_t), intent(inout) :: aero_state ! aerosol state
    type(gas_data_t), intent(in) :: gas_data ! gas data
    type(gas_state_t), intent(inout) :: gas_state ! gas state
    type(run_mc_opt_t), intent(in) :: mc_opt ! Monte Carlo options
    type(inout_file_t), intent(inout) :: summary_file ! summary output file

    ! FIXME: can we shift this to a module? mod_kernel presumably
    interface
       subroutine kernel(v1, v2, env, k)
         use mod_environ
         real*8, intent(in) :: v1
         real*8, intent(in) :: v2
         type(environ), intent(in) :: env 
         real*8, intent(out) :: k
       end subroutine kernel
    end interface
    
    real*8 time, pre_time
    real*8 last_output_time, last_state_time, last_progress_time
    real*8 k_max(bin_grid%n_bin, bin_grid%n_bin)
    integer n_coag, tot_n_samp, tot_n_coag
    logical do_output, do_state, do_progress, did_coag
    real*8 t_start, t_wall_now, t_wall_est, prop_done
    integer i_time, pre_i_time
    character*100 filename
    type(bin_grid_t) :: restart_bin_grid
    type(aero_data_t) :: restart_aero_data
    type(gas_data_t) :: restart_gas_data

    i_time = 0
    time = 0d0
    call init_environ(env, time)
    n_coag = 0
    tot_n_samp = 0
    tot_n_coag = 0
    
    if (mc_opt%do_restart) then
       call inout_read_state(mc_opt%restart_name, restart_bin_grid, &
            restart_aero_data, aero_state, restart_gas_data, gas_state, &
            env, time)
       ! FIXME: should we check whether bin_grid == restart_bin_grid, etc?
       i_time = nint(time / mc_opt%del_t)
       if (mc_opt%allow_double) then
          do while (total_particles(aero_state) .lt. mc_opt%n_part_max / 2)
             call double(aero_state)
          end do
       end if
       ! write data into output file so that it will look correct
       pre_time = 0d0
       last_output_time = pre_time
       call aero_state_to_binned(bin_grid, aero_data, aero_state, aero_binned)
       do pre_i_time = 0,(i_time - 1)
          call update_environ(env, pre_time)
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

    if (mc_opt%t_output > 0d0) then
       call output_summary(summary_file, &
            time, bin_grid, aero_data, aero_binned, &
            gas_data, gas_state, env, mc_opt%i_loop)
    end if

    if (mc_opt%t_state > 0d0) then
       call inout_write_state(mc_opt%state_prefix, bin_grid, &
            aero_data, aero_state, gas_data, gas_state, env, i_time, &
            time, mc_opt%i_loop)
    end if
    
    t_start = time
    last_progress_time = time
    last_state_time = time
    last_output_time = time
    do while (time < mc_opt%t_max)
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
          call singlestep_mosaic(bin_grid, env, aero_data, &
               aero_state, gas_data, gas_state, time, mc_opt%del_t)
       end if
       
       ! DEBUG: enable to check array handling
       ! call check_aero_state(bin_grid, aero_binned, aero_data, aero_state)
       ! DEBUG: end
       
       i_time = i_time + 1
       time = time + mc_opt%del_t

       call update_environ(env, time)
       call environ_update_gas_state(env, mc_opt%del_t, gas_data, gas_state)
       call environ_update_aero_state(env, mc_opt%del_t, bin_grid, &
            aero_data, aero_state, aero_binned)

       if (mc_opt%t_output > 0d0) then
          call check_event(time, mc_opt%del_t, mc_opt%t_output, &
               last_output_time, do_output)
          if (do_output) call output_summary(summary_file, &
               time, bin_grid, aero_data, aero_binned, &
               gas_data, gas_state, env, mc_opt%i_loop)
       end if

       if (mc_opt%t_state > 0d0) then
          call check_event(time, mc_opt%del_t, mc_opt%t_state, &
               last_state_time, do_state)
          if (do_state) call inout_write_state(mc_opt%state_prefix, bin_grid, &
            aero_data, aero_state, gas_data, gas_state, env, i_time, &
            time, mc_opt%i_loop)
       end if
       
       if (mc_opt%t_progress > 0d0) then
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
       
    enddo
    
  end subroutine run_mc
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine mc_coag(kernel, bin_grid, aero_binned, env, aero_data, &
       aero_state, mc_opt, k_max, tot_n_samp, n_coag)

    ! Do coagulation for time del_t.

    use mod_util
    use mod_aero_state
    use mod_bin_grid
    use mod_aero_binned
    use mod_environ
    use mod_aero_data
    use mod_coagulation

    type(bin_grid_t), intent(in) :: bin_grid ! bin grid
    type(aero_binned_t), intent(out) :: aero_binned ! binned distributions
    type(environ), intent(inout) :: env ! environment state
    type(aero_data_t), intent(in) :: aero_data ! aerosol data
    type(aero_state_t), intent(inout) :: aero_state ! aerosol state
    type(run_mc_opt_t), intent(in) :: mc_opt ! Monte Carlo options
    real*8, intent(in) :: k_max(bin_grid%n_bin,bin_grid%n_bin) ! maximum kernel
    integer, intent(out) :: tot_n_samp  ! total number of samples tested
    integer, intent(out) :: n_coag      ! number of coagulation events

    interface
       subroutine kernel(v1, v2, env, k)
         use mod_environ
         real*8, intent(in) :: v1
         real*8, intent(in) :: v2
         type(environ), intent(in) :: env 
         real*8, intent(out) :: k
       end subroutine kernel
    end interface
    
    logical did_coag, bin_change
    integer i, j, n_samp, i_samp, M
    real*8 n_samp_real

    tot_n_samp = 0
    n_coag = 0
    do i = 1,bin_grid%n_bin
       do j = 1,bin_grid%n_bin
          call compute_n_samp(aero_state%n(i), aero_state%n(j), i == j, &
               k_max(i,j), aero_state%comp_vol, mc_opt%del_t, n_samp_real)
          ! probabalistically determine n_samp to cope with < 1 case
          n_samp = int(n_samp_real)
          if (util_rand() .lt. mod(n_samp_real, 1d0)) then
             n_samp = n_samp + 1
          endif
          tot_n_samp = tot_n_samp + n_samp
          do i_samp = 1,n_samp
             M = total_particles(aero_state)
             ! check we still have enough particles to coagulate
             if ((aero_state%n(i) < 1) .or. (aero_state%n(j) < 1) &
                  .or. ((i == j) .and. (aero_state%n(i) < 2))) then
                exit
             end if
             call maybe_coag_pair(bin_grid, aero_binned, env, aero_data, &
                  aero_state, i, j, mc_opt%del_t, k_max(i,j), kernel, &
                  did_coag, bin_change)
             if (did_coag) n_coag = n_coag + 1
          enddo
       enddo
    enddo

    ! if we have less than half the maximum number of particles
    ! then double until we fill up the array
    if (mc_opt%allow_double) then
       do while (total_particles(aero_state) .lt. mc_opt%n_part_max / 2)
          call double(aero_state)
       end do
    end if
    
  end subroutine mc_coag

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine compute_n_samp(ni, nj, same_bin, k_max, comp_vol, &
       del_t, n_samp_real)
  
    ! Compute the number of samples required for the pair of bins.

    use mod_environ

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
  
end module mod_run_mc
