! -*- mode: f90; -*-
! Copyright (C) 2005,2006 Nicole Riemer and Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.
!
! Monte Carlo with fixed timestep and a hybrid array.

module mod_mc_fix_hybrid
contains
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine mc_fix_hybrid(MM, M, n_spec, V, n_bin, TDV, &
       MH, VH, &
       bin_v, bin_r, bin_g, bin_gs, bin_n, dlnr, &
       kernel, t_max, t_print, &
       t_progress, del_t, loop, env, mat)
    
    use mod_util
    use mod_array
    use mod_array_hybrid
    use mod_bin 
    use mod_condensation
    use mod_environ
    use mod_material
    use mod_state
    
    integer, intent(in) :: MM                !  physical dimension of V
    integer, intent(inout) :: M                 !  logical dimension of V
    integer, intent(in) :: n_spec            !  number of species
    real*8, intent(inout) :: V(MM,n_spec)       !  particle volumes (m^3)
    integer, intent(in) :: n_bin             !  number of bins
    integer, intent(in) :: TDV               !  trailing dimension of VH
    integer, intent(out) :: MH(n_bin)         !  number of particles per bin
    real*8, intent(out) :: VH(n_bin,TDV,n_spec) !  particle volumes (m^3)
    
    real*8, intent(in) :: bin_v(n_bin)       !  volume of particles in bins (m^3)
    real*8, intent(in) :: bin_r(n_bin)       !  radius of particles in bins (m)
    real*8, intent(out) :: bin_g(n_bin)       !  mass in bins  
    real*8, intent(out) :: bin_gs(n_bin,n_spec) !  species mass in bins
    integer, intent(out) :: bin_n(n_bin)      !  number in bins
    real*8, intent(in) :: dlnr               !  bin scale factor
    
    real*8, intent(in) :: t_max              !  final time (seconds)
    real*8, intent(in) :: t_print            !  interval to output data (seconds)
    real*8, intent(in) :: t_progress         !  interval to print progress (seconds)
    real*8, intent(in) :: del_t              !  timestep for coagulation
    
    integer, intent(in) :: loop              !  loop number of run
    
    type(environ), intent(inout) :: env  ! environment state
    type(material), intent(in) :: mat    ! material properties
    
    interface
       subroutine kernel(v1, v2, env, k)
         use mod_environ
         real*8, intent(in) :: v1
         real*8, intent(in) :: v2
         type(environ), intent(in) :: env 
         real*8, intent(out) :: k
       end subroutine kernel
    end interface
    
    real*8 time, last_print_time, last_progress_time
    real*8 k_max(n_bin, n_bin), n_samp_real
    integer n_samp, i_samp, n_coag, i, j, tot_n_samp, tot_n_coag, k
    logical do_print, do_progress, did_coag, bin_change
    real*8 t_start, t_wall_start, t_wall_now, t_wall_est
    integer i_time
    character*100 filename
    
    i_time = 0
    time = 0d0
    tot_n_coag = 0
    call array_to_hybrid(MM, M, V, n_spec, n_bin, bin_v, TDV, MH, VH)
    
    ! RESTART
    !      filename = 'start_state_0800_2e8.d'
    !      i_time = 800
    !      call read_state(filename, n_bin, TDV, n_spec, MH, VH, env, time)
    !      M = sum(MH)
    !      do while (M .lt. MM / 2)
    !         call double_hybrid(M, n_bin, TDV, MH, VH, env%V_comp, n_spec
    !     $        ,bin_v,bin_r, bin_g, bin_gs, bin_n, dlnr)
    !      enddo
    ! RESTART
    
    call moments_hybrid(n_bin, TDV, n_spec, MH, VH, bin_v, &
         bin_r, bin_g, bin_gs, bin_n, dlnr)
    
    call est_k_max_binned(n_bin, bin_v, kernel, env, k_max)
    
    call print_info(time, env%V_comp, n_spec, n_bin, bin_v, &
         bin_r,bin_g, bin_gs, bin_n, dlnr, env, mat)
    call write_state_hybrid(n_bin, TDV, n_spec, MH, VH, env, i_time, &
         time)
    
    call cpu_time(t_wall_start)
    t_start = time
    last_progress_time = time
    last_print_time = time
    do while (time < t_max)
       tot_n_samp = 0
       n_coag = 0
       do i = 1,n_bin
          do j = 1,n_bin
             call compute_n_samp_hybrid(n_bin, MH, i, j, env%V_comp, &
                  k_max, del_t, n_samp_real)
             ! probabalistically determine n_samp to cope with < 1 case
             n_samp = int(n_samp_real)
             if (util_rand() .lt. mod(n_samp_real, 1d0)) then
                n_samp = n_samp + 1
             endif
             tot_n_samp = tot_n_samp + n_samp
             do i_samp = 1,n_samp
                call maybe_coag_pair_hybrid(M, n_bin, TDV, MH, VH, &
                     env%V_comp, n_spec, bin_v, bin_r, bin_g, bin_gs, &
                     bin_n, dlnr, i, j, del_t, k_max(i,j), kernel, &
                     env, did_coag, bin_change)
                if (did_coag) n_coag = n_coag + 1
             enddo
          enddo
       enddo
       
       tot_n_coag = tot_n_coag + n_coag
       if (M .lt. MM / 2) then
          call double_hybrid(M, n_bin, TDV, MH, VH, env%V_comp, n_spec &
               ,bin_v,bin_r, bin_g, bin_gs, bin_n, dlnr)
       endif
       
       ! NO CONDENSATION IN RESTART RUN
       !call condense_particles(n_bin, TDV, n_spec, MH, VH, del_t, &
       !     bin_v, bin_r, bin_g, bin_gs, bin_n, dlnr, env, mat)
       ! NO CONDENSATION IN RESTART RUN
       
       ! DEBUG
       !         call check_hybrid(M, n_bin, n_spec, TDV, MH, VH, bin_v, bin_r,
       !     &        bin_g, bin_gs, bin_n, dlnr)
       ! DEBUG
       
       i_time = i_time + 1
       time = time + del_t
       call change_temp(env, del_t)
       if (time .ge. 1200d0) then
          env%dTdt = 0d0
       endif
       
       call check_event(time, del_t, t_print, last_print_time, &
            do_print)
       if (do_print) call print_info(time, env%V_comp, n_spec, n_bin, &
            bin_v, bin_r, bin_g, bin_gs, bin_n, dlnr, env, mat)
       if (do_print) call write_state_hybrid(n_bin, TDV, n_spec, MH, &
            VH, env, i_time, time)
       
       call check_event(time, del_t, t_progress, last_progress_time, &
            do_progress)
       if (do_progress) then
          call cpu_time(t_wall_now)
          t_wall_est = (t_max - time) * (t_wall_now - t_wall_start) &
               / (time - t_start)
          write(6,'(a6,a8,a9,a11,a9,a11,a10)') 'loop', 'time', 'M', &
               'tot_n_samp', 'n_coag', 'tot_n_coag', 't_est'
          write(6,'(i6,f8.1,i9,i11,i9,i11,f10.0)') loop, time, M, &
               tot_n_samp, n_coag, tot_n_coag, t_wall_est
       endif
       
    enddo
    
  end subroutine mc_fix_hybrid
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine compute_n_samp_hybrid(n_bin, MH, i, j, V_comp, k_max, &
       del_t, n_samp_real)
    
    integer, intent(in) :: n_bin              !  number of bins
    integer, intent(in) :: MH(n_bin)          !  number particles per bin
    integer, intent(in) :: i                  !  first bin 
    integer, intent(in) :: j                  !  second bin
    real*8, intent(in) :: V_comp              !  computational volume
    real*8, intent(in) :: k_max(n_bin,n_bin)  !  maximum kernel values
    real*8, intent(in) :: del_t               !  timestep (s)
    real*8, intent(out) :: n_samp_real         !  number of samples per timestep
    !         for bin-i to bin-j events
    
    real*8 r_samp
    real*8 n_possible ! use real*8 to avoid integer overflow
    
    if (i .eq. j) then
       n_possible = dble(MH(i)) * (dble(MH(j)) - 1d0) / 2d0
    else
       n_possible = dble(MH(i)) * dble(MH(j)) / 2d0
    endif
    
    r_samp = k_max(i,j) * 1d0/V_comp * del_t
    n_samp_real = r_samp * n_possible
    
  end subroutine compute_n_samp_hybrid
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
end module mod_mc_fix_hybrid
