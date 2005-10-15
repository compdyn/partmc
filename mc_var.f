C Monte Carlo with variable timestep.

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine mc_var(MM, M, M_comp, V, V_comp, kernel, n_bin, vv,
     &     rr, g, n_ln, dlnr, t_max, t_print, t_k_max, t_k_avg,
     &     k_avg_samp)

      integer MM         ! INPUT: physical dimension of V
      integer M          ! INPUT/OUTPUT: number of particles
      integer M_comp     ! INPUT/OUTPUT: logical dimension of V
      real*8 V(MM)       ! INPUT/OUTPUT: particle volumes
      real*8 V_comp      ! INPUT/OUTPUT: computational volume
      external kernel    ! INPUT: kernel function
      integer n_bin      ! INPUT: number of bins
      real*8 vv(n_bin)   ! INPUT: volume of particles in bins
      real*8 rr(n_bin)   ! INPUT: radius of particles in bins
      real*8 g(n_bin)    ! OUTPUT: mass in bins
      real*8 n_ln(n_bin) ! OUTPUT: number in bins
      real*8 dlnr        ! INPUT: scale factor
      real*8 t_max       ! INPUT: final time (seconds)
      real*8 t_print     ! INPUT: interval to print info (seconds)
      real*8 t_k_max     ! INPUT: interval to update k_max (seconds)
      real*8 t_k_avg     ! INPUT: interval to update k_avg (seconds)
      integer k_avg_samp ! INPUT: number of samples to estimate k_avg

      real*8 del_t, k_max, k_avg
      real*8 time, last_print_time, last_k_max_time, last_k_avg_time
      real*8 t_progress, last_progress_time
      logical do_print, do_k_max, do_k_avg, do_progress
      integer s1, s2, n_coag

      t_progress = 1d0 ! how often to print progress (simulation seconds)
      n_coag = 0
      last_progress_time = 0d0

      time = 0d0
      call moments(MM, V, n_bin, M_comp, V_comp, vv, dlnr, g, n_ln)
      call check_event(time, t_print, last_print_time, do_print)
      if (do_print) call print_info(n_bin, time, rr, g, n_ln)
      call check_event(time, t_k_max, last_k_max_time, do_k_max)
      if (do_k_max) call est_k_max(n_bin, rr, n_ln, dlnr, kernel, k_max)
      call check_event(time, t_k_avg, last_k_avg_time, do_k_avg)
      if (do_k_avg) call kernel_avg(MM, M_comp, V, kernel, k_avg_samp,
     &     k_avg)

      do while (time < t_max)
         ! coagulate a pair and increment time
         call find_rand_pair_acc_rej(MM, V, M_comp, k_max,
     &        kernel, s1, s2)
         call coagulate(MM, M, V, s1, s2)
         del_t = V_comp / (k_avg * (M*(M-1)/2d0))
         time = time + del_t
         n_coag = n_coag + 1

         ! update things if required
         call check_event(time, t_print, last_print_time, do_print)
         call check_event(time, t_k_max, last_k_max_time, do_k_max)
         call check_event(time, t_k_avg, last_k_avg_time, do_k_avg)
         if (do_print .or. do_k_max) call moments(MM, V, n_bin, M_comp,
     &        V_comp, vv, dlnr, g, n_ln)
         if (do_print) call print_info(n_bin, time, rr, g, n_ln)
         if (do_k_max) call est_k_max(n_bin, rr, n_ln, dlnr,
     &        kernel, k_max)
         if (do_k_avg) call kernel_avg(MM, M_comp, V, kernel,
     &        k_avg_samp, k_avg)

         ! print progress
         call check_event(time, t_progress, last_progress_time,
     &        do_progress)
         if (do_progress) then
            write(6,'(a6,a6,a6,a7,a10,a9,a11)')
     &           'time','del_t','M','M_comp','k_max','n_coag'
            write(6,'(f6.1,f6.3,i6,i7,e10.3,i9)')
     &           time, del_t, M, M_comp, k_max, n_coag
         endif

         ! if we are running low on particles then top-up
         if (M < MM / 2) then
            call double(MM, M, M_comp, V, V_comp)
            write(6,*)'double'
         endif
      enddo

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
