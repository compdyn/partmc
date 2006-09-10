C     Monte Carlo with fixed timestep and a hybrid array.

      module mod_mc_fix_hybrid
      contains

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine mc_fix_hybrid(MM, M, V, n_spec, n_bin, TDV, 
     $     MH, VH, V_comp,
     $     bin_v, i_water,
     $     bin_r, bin_g, bin_gs, bin_n, dlnr, 
     $     kernel, t_max, t_print,
     $     t_progress, del_t, loop, env, mat)

      use mod_array
      use mod_array_hybrid
      use mod_bin 
      use mod_condensation
      use mod_environ
      use mod_material
      use mod_state

      integer MM                ! INPUT: physical dimension of V
      integer M                 ! INPUT/OUTPUT: logical dimension of V
      integer n_spec            ! INPUT: number of species
      real*8 V(MM,n_spec)       ! INPUT/OUTPUT: particle volumes (m^3)
      integer n_bin             ! INPUT: number of bins
      integer TDV               ! INPUT: trailing dimension of VH
      integer MH(n_bin)         ! OUTPUT: number of particles per bin
      real*8 VH(n_bin,TDV,n_spec) ! OUTPUT: particle volumes (m^3)
      real*8 V_comp             ! INPUT/OUTPUT: computational volume (m^3)
      integer i_water           ! INPUT: water species number
      
      real*8 bin_v(n_bin)       ! INPUT: volume of particles in bins (m^3)
      real*8 bin_r(n_bin)       ! INPUT: radius of particles in bins (m)
      real*8 bin_g(n_bin)       ! OUTPUT: mass in bins  
      real*8 bin_gs(n_bin,n_spec) ! OUTPUT: species mass in bins             
      integer bin_n(n_bin)      ! OUTPUT: number in bins
      real*8 dlnr               ! INPUT: bin scale factor
      
      real*8 t_max              ! INPUT: final time (seconds)
      real*8 t_print            ! INPUT: interval to output data (seconds)
      real*8 t_progress         ! INPUT: interval to print progress (seconds)
      real*8 del_t              ! INPUT: timestep for coagulation
      
      integer loop              ! INPUT: loop number of run

      type(environ), intent(inout) :: env  ! environment state
      type(material), intent(in) :: mat    ! material properties

      interface
         subroutine kernel(v1, v2, k)
         real*8, intent(in) :: v1
         real*8, intent(in) :: v2
         real*8, intent(out) :: k
         end subroutine
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
      filename = 'start_state1150.d'
      i_time = 1150
      call read_state(filename, n_bin, TDV, n_spec, MH, VH, env, time)
      M = sum(MH)
! RESTART

      call moments_hybrid(n_bin, TDV, n_spec, MH, VH, bin_v,
     &     bin_r, bin_g, bin_gs, bin_n, dlnr)
      
      call est_k_max_binned(n_bin, bin_v, kernel, k_max)

      call print_info(time, V_comp, n_spec, n_bin, bin_v,
     $     bin_r,bin_g, bin_gs, bin_n, dlnr, env, mat)
      call write_state_hybrid(n_bin, TDV, n_spec, MH, VH, env, i_time,
     $     time)

      call cpu_time(t_wall_start)
      t_start = time
      last_progress_time = time
      last_print_time = time
      do while (time < t_max)
         tot_n_samp = 0
         n_coag = 0
         do i = 1,n_bin
            do j = 1,n_bin
               call compute_n_samp_hybrid(n_bin, MH, i, j, V_comp,
     $              k_max, del_t, n_samp_real)
               ! probabalistically determine n_samp to cope with < 1 case
               n_samp = int(n_samp_real)
               if (dble(rand()) .lt. mod(n_samp_real, 1d0)) then
                  n_samp = n_samp + 1
               endif
               tot_n_samp = tot_n_samp + n_samp
               do i_samp = 1,n_samp
                  call maybe_coag_pair_hybrid(M, n_bin, TDV, MH, VH,
     $                 V_comp, n_spec, bin_v, bin_r, bin_g, bin_gs,
     $                 bin_n, dlnr, i, j, del_t, k_max(i,j), kernel,
     $                 did_coag, bin_change)
                  if (did_coag) n_coag = n_coag + 1
               enddo
            enddo
         enddo

         tot_n_coag = tot_n_coag + n_coag
         if (M .lt. MM / 2) then
            call double_hybrid(M, n_bin, TDV, MH, VH, V_comp, n_spec
     $           ,bin_v,bin_r, bin_g, bin_gs, bin_n, dlnr)
         endif

         call condense_particles(n_bin, TDV, n_spec, MH, VH, del_t,
     $        bin_v, bin_r, bin_g, bin_gs, bin_n, dlnr, env, mat)

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

         call check_event(time, del_t, t_print, last_print_time,
     &        do_print)
         if (do_print) call print_info(time, V_comp, n_spec, n_bin,
     $        bin_v, bin_r, bin_g, bin_gs, bin_n, dlnr, env, mat)
         if (do_print) call write_state_hybrid(n_bin, TDV, n_spec, MH,
     $        VH, env, i_time, time)

         call check_event(time, del_t, t_progress, last_progress_time,
     $        do_progress)
         if (do_progress) then
            call cpu_time(t_wall_now)
            t_wall_est = (t_max - time) * (t_wall_now - t_wall_start)
     &           / (time - t_start)
            write(6,'(a6,a8,a9,a11,a9,a11,a10)') 'loop', 'time', 'M',
     $           'tot_n_samp', 'n_coag', 'tot_n_coag', 't_est'
            write(6,'(i6,f8.1,i9,i11,i9,i11,f10.0)') loop, time, M,
     $           tot_n_samp, n_coag, tot_n_coag, t_wall_est
         endif

      enddo

      end subroutine

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine compute_n_samp_hybrid(n_bin, MH, i, j, V_comp, k_max,
     $     del_t, n_samp_real)

      integer n_bin              ! INPUT: number of bins
      integer MH(n_bin)          ! INPUT: number particles per bin
      integer i                  ! INPUT: first bin 
      integer j                  ! INPUT: second bin
      real*8 V_comp              ! INPUT: computational volume
      real*8 k_max(n_bin,n_bin)  ! INPUT: maximum kernel values
      real*8 del_t               ! INPUT: timestep (s)
      real*8 n_samp_real         ! OUTPUT: number of samples per timestep
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

      end subroutine

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      end module
