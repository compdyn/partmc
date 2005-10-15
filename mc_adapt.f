C Monte Carlo with adaptive timestep.

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine mc_adapt(MM, M, M_comp, V, V_comp, kernel, n_bin, vv,
     &     rr, g, n_ln, dlnr, t_max, t_print,
     &     p_max, r_samp_max, del_t_max)

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
      real*8 p_max       ! INPUT: maximum coagulation probability
      real*8 r_samp_max  ! INPUT: maximum sampling ratio per timestep
      real*8 del_t_max   ! INPUT: maximum timestep

      real*8 del_t, time, last_print_time, k_max
      integer n_samp, i_samp, n_coag
      logical do_print, did_coag
      real*8 t_start, t_end, t_loop, t_per_samp

      time = 0
      n_coag = 0
      call moments(MM, V, n_bin, M_comp, V_comp, vv, dlnr, g, n_ln)
      call check_event(time, t_print, last_print_time, do_print)
      if (do_print) call print_info(n_bin, time, rr, g, n_ln)

      do while (time < t_max)
         call cpu_time(t_start)

         call est_k_max(n_bin, rr, n_ln, dlnr, kernel, k_max)
         call compute_n_samp_del_t(M, V_comp, k_max, p_max, r_samp_max,
     &        del_t_max, n_samp, del_t)
         do i_samp = 1,n_samp
            call maybe_coag_pair(MM, V, M, M_comp, V_comp,
     &           del_t, n_samp, kernel, did_coag)
            if (did_coag) n_coag = n_coag + 1
            if (M .lt. MM / 2) then
               call double(MM, M, M_comp, V, V_comp)
               write(6,*)'double'
            endif
         enddo

         time = time + del_t

         call moments(MM, V, n_bin, M_comp, V_comp, vv, dlnr, g, n_ln)
         call check_event(time, t_print, last_print_time, do_print)
         if (do_print) call print_info(n_bin, time, rr, g, n_ln)

         call cpu_time(t_end)
         t_loop = t_end - t_start
         t_per_samp = t_loop / n_samp
         write(6,'(a6,a6,a6,a7,a10,a9,a11,a9)')
     &        'time', 'del_t', 'M', 'M_comp', 'k_max',
     &        'n_samp', 't_per_samp', 'n_coag'
         write(6,'(f6.1,f6.3,i6,i7,e10.3,i9,e11.3,i9)')
     &        time, del_t, M, M_comp, k_max, n_samp,
     &        t_per_samp, n_coag
      enddo

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine compute_n_samp_del_t(M, V_comp, k_max, p_max,
     &     r_samp_max, del_t_max, n_samp, del_t)

      integer M            ! INPUT: number of particles
      real*8 V_comp        ! INPUT: computational volume
      real*8 k_max         ! INPUT: maximum kernel value
      real*8 p_max         ! INPUT: maximum coagulation probability
      real*8 r_samp_max    ! INPUT: maximum sampling ratio per timestep
      real*8 del_t_max     ! INPUT: maximum timestep
      integer n_samp       ! OUTPUT: number of samples per timestep
      real*8 del_t         ! OUTPUT: timestep

      real*8 r_samp, c

      c = - (k_max * 1d0/V_comp / log(1 - p_max))
      del_t = r_samp_max / c
      if (del_t .gt. del_t_max) del_t = del_t_max
      r_samp = del_t * c
      n_samp = int(r_samp * M*(M-1)/2)

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
