C     Monte Carlo with fixed timestep and a hybrid array.

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine mc_fix_hybrid(MM, M, V, V_comp, n_bin, bin_v, bin_r,
     $     bin_g, bin_n, dlnr, kernel, t_max, t_print, t_progress ,
     $     del_t, loop)

      integer MM           ! INPUT: physical dimension of V
      integer M            ! INPUT/OUTPUT: logical dimension of V
      real*8 V(MM)         ! INPUT/OUTPUT: particle volumes (m^3)
      real*8 V_comp        ! INPUT/OUTPUT: computational volume (m^3)

      integer n_bin        ! INPUT: number of bins
      real*8 bin_v(n_bin)  ! INPUT: volume of particles in bins (m^3)
      real*8 bin_r(n_bin)  ! INPUT: radius of particles in bins (m)
      real*8 bin_g(n_bin)  ! OUTPUT: mass in bins               
      integer bin_n(n_bin) ! OUTPUT: number in bins
      real*8 dlnr          ! INPUT: bin scale factor

      external kernel      ! INPUT: kernel function
      real*8 t_max         ! INPUT: final time (seconds)
      real*8 t_print       ! INPUT: interval to output data (seconds)
      real*8 t_progress    ! INPUT: interval to print progress (seconds)
      real*8 del_t         ! INPUT: timestep
      integer loop         ! INPUT: loop number of run

      real*8 VH(n_bin, MM)
      integer MH(n_bin)
      real*8 time, last_print_time, last_progress_time
      real*8 k_max(n_bin, n_bin)
      integer n_samp, i_samp, n_coag, i, j
      logical do_print, do_progress, did_coag, bin_change
      real*8 t_start, t_end, t_est

      last_progress_time = 0d0
      time = 0d0
      n_coag = 0
      call moments(MM, M, V, V_comp, n_bin, bin_v, bin_r, bin_g, bin_n,
     $     dlnr)
      call check_event(time, t_print, last_print_time, do_print)
      if (do_print) call print_info(time, V_comp, n_bin, bin_v, bin_r,
     $     bin_g, bin_n, dlnr)
      call array_to_hybrid(MM, M, V, n_bin, bin_v, MH, VH)
      call est_k_max_binned(n_bin, bin_v, kernel, k_max)

      call cpu_time(t_start)
      do while (time < t_max)
         do i = 1,n_bin
            do j = 1,n_bin
               call compute_n_samp_hybrid(n_bin, MH, i, j, V_comp,
     $              k_max(i,j), del_t, n_samp_real)
               ! probabalistically determine n_samp to cope with < 1 case
               n_samp = int(n_samp_real)
               if (rand() .lt. mod(n_samp_real, 1)) n_samp = n_samp + 1
               do i_samp = 1,n_samp
                  call maybe_coag_pair(MM, M, n_bin, MH, VH, V_comp,
     $                 bin_v,bin_r, bin_g, bin_n, dlnr, del_t, k_max(i,j
     $                 ), kernel, did_coag,bin_change)
                  if (did_coag) n_coag = n_coag + 1
                  ! FIXME: just store up bin_change and do once outside loop
                  if (bin_change) call est_k_max_split(n_bin, s_bin,
     $                 bin_v,bin_n, kernel, k_max_small, k_max_big)
               enddo
            enddo
         enddo
         if (M .lt. MM / 2) then
            call double_hybrid(MM, n_bin, M, V, V_comp, bin_v, bin_r,
     $           bin_g, bin_n, dlnr)
         endif

         time = time + del_t

         call check_event(time, t_print, last_print_time, do_print)
         if (do_print) call print_info(time, V_comp, n_bin, bin_v, bin_r
     $        , bin_g, bin_n, dlnr)

         call check_event(time, t_progress, last_progress_time,
     $        do_progress)
         if (do_progress) then
            call cpu_time(t_end)
            t_est = (t_max - time) / time * (t_end - t_start)
            write(6,'(a6,a8,a9,a9,a10)') 'loop', 'time', 'M',
     $           'n_coag', 't_est'
            write(6,'(i6,f8.1,i9,i9,f10)') loop, time, M,
     $           n_coag, t_est
         endif
      enddo

      return
      end

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
         n_possible = dble(MH(i)) * (dble(MH(j)) - 1d0)
      else
         n_possible = dble(MH(i)) * dble(MH(j))
      endif

      r_samp = k_max(i,j) * 1d0/V_comp * del_t
      n_samp_real = r_samp_real

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
