! -*- mode: f90; -*-
! Copyright (C) 2005-2007 Nicole Riemer and Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.
!
! Monte Carlo with fixed timestep and particles stored per-bin.

module mod_run_mc
contains
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine run_mc(MM, M, n_spec, n_bin, MH, VH, bin_v, bin_g, &
       bin_gs, bin_n, dlnr, kernel, t_max, t_output, t_state, &
       t_progress, del_t, output_unit, state_unit, state_name, &
       do_coagulation, do_condensation, do_restart, restart_name, &
       i_loop, n_loop, t_wall_start, env, mat)
    
    use mod_util
    use mod_array
    use mod_bin 
    use mod_condensation
    use mod_environ
    use mod_material
    use mod_state
    
    integer, intent(in) :: MM                ! maximum number of particles
    integer, intent(inout) :: M              ! actual number of particles
    integer, intent(in) :: n_spec            ! number of species
    integer, intent(in) :: n_bin             ! number of bins
    integer, intent(out) :: MH(n_bin)        ! number of particles per bin
    type(bin_p), intent(inout) :: VH(n_bin)  ! particle volumes (m^3)
    
    real*8, intent(in) :: bin_v(n_bin)       ! volume of particles in bins (m^3)
    real*8, intent(out) :: bin_g(n_bin)      ! volume in bins  
    real*8, intent(out) :: bin_gs(n_bin,n_spec) ! species volume in bins
    integer, intent(out) :: bin_n(n_bin)     ! number in bins
    real*8, intent(in) :: dlnr               ! bin scale factor
    
    real*8, intent(in) :: t_max              ! final time (seconds)
    real*8, intent(in) :: t_output           ! interval to output data, or
                                             ! zero to not output (seconds)
    real*8, intent(in) :: t_state            ! interval to output state, or
                                             ! zero to not output (seconds)
    real*8, intent(in) :: t_progress         ! interval to print progress, or
                                             ! zero to not print (seconds)
    real*8, intent(in) :: del_t              ! timestep for coagulation
    integer, intent(in) :: output_unit       ! unit number to output to
    integer, intent(in) :: state_unit        ! unit number for state files
    character(len=*), intent(in) :: state_name ! name for state files
    
    logical, intent(in) :: do_coagulation    ! whether to do coagulation
    logical, intent(in) :: do_condensation   ! whether to do condensation
    logical, intent(in) :: do_restart        ! whether to restart from state
    character(len=*), intent(in) :: restart_name ! name of state to restart from
    integer, intent(in) :: i_loop            ! loop number of run
    integer, intent(in) :: n_loop            ! total number of loops
    real*8, intent(in) :: t_wall_start       ! cpu_time() of start
    
    type(environ), intent(inout) :: env      ! environment state
    type(material), intent(in) :: mat        ! material properties
    
    interface
       subroutine kernel(v1, v2, env, k)
         use mod_environ
         real*8, intent(in) :: v1
         real*8, intent(in) :: v2
         type(environ), intent(in) :: env 
         real*8, intent(out) :: k
       end subroutine kernel
    end interface
    
    real*8 time, last_output_time, last_state_time, last_progress_time
    real*8 k_max(n_bin, n_bin)
    integer n_coag, tot_n_samp, tot_n_coag
    logical do_output, do_state, do_progress, did_coag
    real*8 t_start, t_wall_now, t_wall_est, prop_done
    integer i_time
    character*100 filename

    i_time = 0
    time = 0d0
    tot_n_coag = 0
    
    if (do_restart) then
       call read_state(restart_name, n_bin, n_spec, MH, VH, bin_v, dlnr, &
              env, time)
       i_time = nint(time / del_t)
       M = sum(MH)
       do while (M .lt. MM / 2)
          call double(M, n_bin, MH, VH, n_spec, &
               bin_v, bin_g, bin_gs, bin_n, dlnr, env)
       end do
    end if
    
    call moments(n_bin, n_spec, MH, VH, bin_v, &
         bin_g, bin_gs, bin_n, dlnr)
    
    call est_k_max_binned(n_bin, bin_v, kernel, env, k_max)

    if (t_output > 0d0) then
       call output_info(output_unit, time, n_bin, n_spec, bin_v, &
            bin_g, bin_gs, bin_n, dlnr, env, mat)
    end if

    if (t_state > 0d0) then
       call write_state(state_unit, state_name, n_bin, n_spec, &
            MH, VH, bin_v, dlnr, env, i_time, time)
    end if
    
    t_start = time
    last_progress_time = time
    last_state_time = time
    last_output_time = time
    do while (time < t_max)
       if (do_coagulation) then
          call mc_coag(MM, M, n_spec, n_bin, MH, VH, &
               bin_v, bin_g, bin_gs, bin_n, dlnr, kernel, k_max, del_t, &
               env, mat, tot_n_samp, n_coag)
       end if

       tot_n_coag = tot_n_coag + n_coag

       if (do_condensation) then
          call condense_particles(n_bin, n_spec, MH, VH, del_t, &
               bin_v, bin_g, bin_gs, bin_n, dlnr, env, mat)
       end if
       
       ! DEBUG: enable to check array handling
       ! call check_array(M, n_bin, n_spec, MH, VH, bin_v, &
       !     bin_g, bin_gs, bin_n, dlnr)
       ! DEBUG: end
       
       i_time = i_time + 1
       time = time + del_t

       ! FIXME: change to linear interpolation
       call change_temp(env, del_t)
       if (time .ge. 1200d0) then
          env%dTdt = 0d0
       endif

       if (t_output > 0d0) then
          call check_event(time, del_t, t_output, last_output_time, &
               do_output)
          if (do_output) call output_info(output_unit, time, n_bin, n_spec, &
               bin_v, bin_g, bin_gs, bin_n, dlnr, env, mat)
       end if

       if (t_state > 0d0) then
          call check_event(time, del_t, t_state, last_state_time, &
               do_state)
          if (do_state) call write_state(state_unit, state_name, n_bin, &
               n_spec, MH, VH, bin_v, dlnr, env, i_time, time)
       end if
       
       if (t_progress > 0d0) then
          call check_event(time, del_t, t_progress, last_progress_time, &
               do_progress)
          if (do_progress) then
             call cpu_time(t_wall_now)
             prop_done = (dble(i_loop - 1) + (time - t_start) &
                  / (t_max - t_start)) / dble(n_loop)
             t_wall_est = (1d0 - prop_done) / prop_done &
                  * (t_wall_now - t_wall_start)
             write(6,'(a6,a8,a9,a11,a9,a11,a10)') 'loop', 'time', 'M', &
                  'tot_n_samp', 'n_coag', 'tot_n_coag', 't_est'
             write(6,'(i6,f8.1,i9,i11,i9,i11,f10.0)') i_loop, time, M, &
                  tot_n_samp, n_coag, tot_n_coag, t_wall_est
          end if
       end if
       
    enddo
    
  end subroutine run_mc
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine mc_coag(MM, M, n_spec, n_bin, MH, VH, bin_v, &
       bin_g, bin_gs, bin_n, dlnr, kernel, k_max, del_t, env, mat, &
       tot_n_samp, n_coag)

    use mod_util
    use mod_array
    use mod_bin 
    use mod_environ
    use mod_material
    
    integer, intent(in) :: MM                ! maximum number of particles
    integer, intent(inout) :: M              ! actual number of particles
    integer, intent(in) :: n_spec            ! number of species
    integer, intent(in) :: n_bin             ! number of bins
    integer, intent(out) :: MH(n_bin)        ! number of particles per bin
    type(bin_p), intent(out) :: VH(n_bin)    ! particle volumes (m^3)
    real*8, intent(in) :: bin_v(n_bin)       ! volume of particles in bins (m^3)
    real*8, intent(out) :: bin_g(n_bin)      ! volume in bins  
    real*8, intent(out) :: bin_gs(n_bin,n_spec) ! species volume in bins
    integer, intent(out) :: bin_n(n_bin)     ! number in bins
    real*8, intent(in) :: dlnr               ! bin scale factor
    real*8, intent(in) :: k_max(n_bin,n_bin) ! maximum kernel values
    real*8, intent(in) :: del_t              ! timestep for coagulation
    type(environ), intent(inout) :: env      ! environment state
    type(material), intent(in) :: mat        ! material properties
    integer, intent(out) :: tot_n_samp       ! total number of samples tested
    integer, intent(out) :: n_coag           ! number of coagulation events

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
    integer i, j, n_samp, i_samp
    real*8 n_samp_real

    tot_n_samp = 0
    n_coag = 0
    do i = 1,n_bin
       do j = 1,n_bin
          call compute_n_samp(n_bin, MH, i, j, &
               k_max, del_t, env, n_samp_real)
          ! probabalistically determine n_samp to cope with < 1 case
          n_samp = int(n_samp_real)
          if (util_rand() .lt. mod(n_samp_real, 1d0)) then
             n_samp = n_samp + 1
          endif
          tot_n_samp = tot_n_samp + n_samp
          do i_samp = 1,n_samp
             call maybe_coag_pair(M, n_bin, MH, VH, &
                  n_spec, bin_v, bin_g, bin_gs, &
                  bin_n, dlnr, i, j, del_t, k_max(i,j), kernel, &
                  env, did_coag, bin_change)
             if (did_coag) n_coag = n_coag + 1
          enddo
       enddo
    enddo
    
    ! if we have less than half the maximum number of particles
    ! then double until we fill up the array
    do while (M .lt. MM / 2)
       call double(M, n_bin, MH, VH,  n_spec, &
            bin_v, bin_g, bin_gs, bin_n, dlnr, env)
    end do
    
  end subroutine mc_coag

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine compute_n_samp(n_bin, MH, i, j, k_max, &
       del_t, env, n_samp_real)

    use mod_environ
    
    integer, intent(in) :: n_bin         ! number of bins
    integer, intent(in) :: MH(n_bin)     ! number particles per bin
    integer, intent(in) :: i             ! first bin 
    integer, intent(in) :: j             ! second bin
    real*8, intent(in) :: k_max(n_bin,n_bin) ! maximum kernel values
    real*8, intent(in) :: del_t          ! timestep (s)
    type(environ), intent(in) :: env        ! environment state
    real*8, intent(out) :: n_samp_real   ! number of samples per timestep
                                         ! for bin-i to bin-j events
    
    real*8 r_samp
    real*8 n_possible ! use real*8 to avoid integer overflow
    
    if (i .eq. j) then
       n_possible = dble(MH(i)) * (dble(MH(j)) - 1d0) / 2d0
    else
       n_possible = dble(MH(i)) * dble(MH(j)) / 2d0
    endif
    
    r_samp = k_max(i,j) * 1d0/env%V_comp * del_t
    n_samp_real = r_samp * n_possible
    
  end subroutine compute_n_samp
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
end module mod_run_mc
