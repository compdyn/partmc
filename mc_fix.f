C Monte Carlo with fixed timestep.

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine mc_fix(MM, M, V, V_comp,
     &     n_bin, bin_v, bin_r, bin_g, bin_n, dlnr,
     &     kernel, t_max, del_t, p_max, t_print, loop)

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

      external kernel    ! INPUT: kernel procedure
      real*8 t_max       ! INPUT: total simulation time
      real*8 del_t       ! INPUT: timestep
      real*8 p_max       ! INPUT: maximum coagulation probability
      real*8 t_print     ! INPUT: interval to print info (seconds)
      integer loop       ! INPUT: loop number of run

      integer i_top, nt, n_samp, i_samp, n_print, n_coag
      real*8 k_max, time
      logical did_coag, bin_change
      real*8 t_start, t_end, t_loop, t_per_samp

      time = 0d0
      n_coag = 0

      call moments(MM, M, V, V_comp,
     &     n_bin, bin_v, bin_r, bin_g, bin_n, dlnr)
      call print_info(time, V_comp,
     &     n_bin, bin_v, bin_r, bin_g, bin_n, dlnr)
      call est_k_max(n_bin, bin_v, bin_n, kernel, k_max)

      nt = int(dble(t_max) / del_t)
      n_print = int(dble(t_print) / del_t)
      do i_top = 1,nt
         call cpu_time(t_start)

         time = dble(i_top) / dble(nt) * dble(t_max)
         
         call compute_n_samp(M, k_max, V_comp, del_t, p_max, n_samp)
         do i_samp = 1,n_samp
            call maybe_coag_pair(MM, M, V, V_comp,
     &           n_bin, bin_v, bin_r, bin_g, bin_n, dlnr,
     &           del_t, n_samp, kernel, did_coag, bin_change)
            if (did_coag) n_coag = n_coag + 1
            if (bin_change) call est_k_max(n_bin, bin_v, bin_n,
     &           kernel, k_max)
            if (M .lt. MM / 2) then
               call double(MM, M, V, V_comp,
     &              n_bin, bin_v, bin_r, bin_g, bin_n, dlnr)
            endif
         enddo

         if (mod(i_top, n_print) .eq. 0) then
            call print_info(time, V_comp,
     &           n_bin, bin_v, bin_r, bin_g, bin_n, dlnr)
         endif

         call cpu_time(t_end)
         t_loop = t_end - t_start
         t_per_samp = t_loop / n_samp
         write(6,'(a6,a6,a6,a6,a6,a10,a9,a11,a9)')
     &        'loop', 'i_top', 'time', 'del_t', 'M',
     &        'k_max', 'n_samp', 't_per_samp', 'n_coag'
         write(6,'(i6,i6,f6.1,f6.3,i6,e10.3,i9,e11.3,i9)')
     &        loop, i_top, time, del_t, M, k_max, n_samp,
     &        t_per_samp, n_coag
      enddo

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine compute_n_samp(M, k_max, V_comp,
     &     del_T, p_max, n_samp)

      integer M            ! INPUT: number of particles
      real*8 k_max         ! INPUT: maximum kernel value
      real*8 V_comp        ! INPUT: computational volume
      real*8 del_T         ! INPUT: timestep
      real*8 p_max         ! INPUT: maximum coagulation probability
      integer n_samp       ! OUTPUT: number of samples per timestep

      real*8 r_samp

      r_samp = - (k_max * 1d0/V_comp * del_T / log(1 - p_max))
      n_samp = int(r_samp * M*(M-1)/2)

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
