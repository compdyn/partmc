C Monte Carlo with variable timestep.

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine mc_var(MM, M, V, V_comp,
     &     n_bin, bin_v, bin_r, bin_g, bin_n, dlnr,
     &     kernel, t_max, t_print, t_k_avg,
     &     k_avg_samp, loop)

      integer MM           ! INPUT: physical dimension of V
      integer M            ! INPUT/OUTPUT: logical dimension of V
      real*8 V(MM)         ! INPUT/OUTPUT: particle volumes
      real*8 V_comp        ! INPUT/OUTPUT: computational volume

      integer n_bin        ! INPUT: number of bins
      real*8 bin_v(n_bin)  ! INPUT: volume of particles in bins
      real*8 bin_r(n_bin)  ! INPUT: radius of particles in bins
      real*8 bin_g(n_bin)  ! OUTPUT: mass in bins
      integer bin_n(n_bin) ! OUTPUT: number in bins
      real*8 dlnr          ! INPUT: bin scale factor

      external kernel    ! INPUT: kernel function
      real*8 t_max       ! INPUT: final time (seconds)
      real*8 t_print     ! INPUT: interval to print info (seconds)
      real*8 t_k_avg     ! INPUT: interval to update k_avg (seconds)
      integer k_avg_samp ! INPUT: number of samples to estimate k_avg
      integer loop       ! INPUT: loop number of run

      real*8 del_t, k_max, k_avg
      real*8 time, last_print_time, last_k_avg_time
      real*8 t_progress, last_progress_time
      logical do_print, do_k_avg, do_progress, bin_change
      integer s1, s2, n_coag

      t_progress = 1d0 ! how often to print progress (simulation seconds)
      n_coag = 0

      time = 0d0
      call moments(MM, M, V, V_comp,
     &     n_bin, bin_v, bin_r, bin_g, bin_n, dlnr)
      call est_k_max(n_bin, bin_v, bin_n, kernel, k_max)
      call check_event(time, t_print, last_print_time, do_print)
      if (do_print) call print_info(time, V_comp,
     &     n_bin, bin_v, bin_r, bin_g, bin_n, dlnr)
      call check_event(time, t_k_avg, last_k_avg_time, do_k_avg)
      if (do_k_avg) call kernel_avg(MM, M, V, kernel, k_avg_samp,
     &     k_avg)

      do while (time < t_max)
         ! coagulate a pair and increment time
         call find_rand_pair_acc_rej(MM, M, V, k_max,
     &        kernel, s1, s2)
         call coagulate(MM, M, V, V_comp,
     &        n_bin, bin_v, bin_r, bin_g, bin_n, dlnr,
     &        s1, s2, bin_change)
         if (bin_change) call est_k_max(n_bin, bin_v, bin_n,
     &        kernel, k_max)
         del_t = V_comp / (k_avg * (M*(M-1)/2d0))
         time = time + del_t
         n_coag = n_coag + 1

         ! update things if required
         call check_event(time, t_print, last_print_time, do_print)
         call check_event(time, t_k_avg, last_k_avg_time, do_k_avg)
         if (do_print) call print_info(time, V_comp,
     &        n_bin, bin_v, bin_r, bin_g, bin_n, dlnr)
          if (do_k_avg) call kernel_avg(MM, M, V, kernel,
     &        k_avg_samp, k_avg)

         ! print progress
         call check_event(time, t_progress, last_progress_time,
     &        do_progress)
         if (do_progress) then
            write(6,'(a6,a6,a6,a6,a10,a9,a11)')
     &           'loop', 'time','del_t','M','k_max','n_coag'
            write(6,'(i6,f6.1,f6.3,i6,e10.3,i9)')
     &           loop, time, del_t, M, k_max, n_coag
         endif

         ! if we are running low on particles then top-up
         if (M < MM / 2) then
            call double(MM, M, V, V_comp,
     &           n_bin, bin_v, bin_r, bin_g, bin_n, dlnr)
         endif
      enddo

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
