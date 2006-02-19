C     Monte Carlo with fixed timestep and a superparticle array.

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine mc_fix_super(MM, M, n_bin, n_fact, TDV, fac_base, MS,
     $     VS, min_fill, V_comp, bin_v, bin_r, bin_g, bin_n, dlnr,
     $     kernel, t_max, t_print, t_progress, del_t, loop)

      integer MM                  ! INPUT: maximum number of phys. particles
      integer M                   ! INPUT/OUTPUT: number of phys. particles
      integer n_bin               ! INPUT: number of bins
      integer n_fact              ! INPUT: number of allowable factors
      integer TDV                 ! INPUT: trailing dimension of VS
      integer fac_base            ! INPUT: factor base of a superparticle
      integer MS(n_bin,n_fact)    ! INPUT/OUTPUT: number of superparticles
      real*8 VS(n_bin,n_fact,TDV) ! INPUT/OUTPUT: volume of physical particles
      integer min_fill            ! INPUT: minimum comp. part. per bin
      real*8 V_comp               ! INPUT/OUTPUT: computational volume

      real*8 bin_v(n_bin)         ! INPUT: volume of particles in bins
      real*8 bin_r(n_bin)         ! INPUT: radius of particles in bins
      real*8 bin_g(n_bin)         ! INPUT/OUTPUT: mass in bins
      integer bin_n(n_bin)        ! INPUT/OUTPUT: number in bins
      real*8 dlnr                 ! INPUT: bin scale factor

      external kernel             ! INPUT: kernel function
      real*8 t_max                ! INPUT: final time (seconds)
      real*8 t_print              ! INPUT: interval to output data (seconds)
      real*8 t_progress           ! INPUT: interval to print progress (seconds)
      real*8 del_t                ! INPUT: timestep
      integer loop                ! INPUT: loop number of run

      real*8 time, last_print_time, last_progress_time
      real*8 k_max(n_bin, n_fact, n_bin, n_fact), n_samp_real
      integer n_samp, i_samp, n_coag, b1, b2, f1, f2
      integer tot_n_samp, tot_n_coag, max_f, coag_factor
      logical do_print, do_progress, did_coag, bin_change
      real*8 t_start, t_end, t_est
! DEBUG
      integer M_comp
! DEBUG

      last_progress_time = 0d0
      time = 0d0
      tot_n_coag = 0
      call init_to_super(n_bin, n_fact, TDV, fac_base, MS, VS, bin_v,
     $     bin_n, min_fill)
      call check_event(time, t_print, last_print_time, do_print)
      if (do_print) call print_info(time, V_comp, n_bin, bin_v, bin_r,
     $     bin_g, bin_n, dlnr)
      call est_k_max_super(n_bin, n_fact, fac_base, bin_v, kernel,
     $     k_max)

! DEBUG
      call check_super(M, n_bin, n_fact, TDV, fac_base, MS, VS, bin_v,
     $     bin_r)
! DEBUG

      call cpu_time(t_start)
      do while (time < t_max)
         tot_n_samp = 0
         n_coag = 0
         do b1 = 1,n_bin
            do f1 = 1,n_fact
               do b2 = 1,n_bin
                  do f2 = 1,n_fact
                     max_f = max(f1, f2)
                     coag_factor = fac_base**(max_f - 1)
                     call compute_n_samp_super(n_bin, n_fact, MS, b1, f1
     $                    , b2, f2, V_comp, k_max, del_t, n_samp_real)
                     ! probabalistically determine n_samp to cope with < 1 case
                     n_samp = int(n_samp_real)
                     if (dble(rand()) .lt. mod(n_samp_real, 1d0)) then
                        n_samp = n_samp + 1
                     endif
                     tot_n_samp = tot_n_samp + n_samp
                     do i_samp = 1,n_samp
                        call maybe_coag_pair_super(M, n_bin, n_fact, TDV
     $                       , fac_base, MS, VS, min_fill, V_comp, bin_v
     $                       , bin_r, bin_g, bin_n, dlnr, b1, b2, f1, f2
     $                       , del_t, k_max(b1, f1, b2, f2), kernel,
     $                       did_coag, bin_change)
                        if (did_coag) n_coag = n_coag + coag_factor
                     enddo
                  enddo
               enddo
            enddo
         enddo
         tot_n_coag = tot_n_coag + n_coag
         if (M .lt. MM / 2) then
            call double_super(M, n_bin, n_fact, TDV, MS, VS, V_comp,
     $           bin_v, bin_r, bin_g, bin_n, dlnr)
         endif

! DEBUG
c         call check_super(M, n_bin, n_fact, TDV, fac_base, MS, VS, bin_v
c     $        , bin_r)
! DEBUG

         time = time + del_t

         call check_event(time, t_print, last_print_time, do_print)
         if (do_print) call print_info(time, V_comp, n_bin, bin_v, bin_r
     $        , bin_g, bin_n, dlnr)

         call check_event(time, t_progress, last_progress_time,
     $        do_progress)
         if (do_progress) then
            call cpu_time(t_end)
            t_est = (t_max - time) / time * (t_end - t_start)
            write(6,'(a6,a8,a9,a11,a9,a11,a10)') 'loop', 'time', 'M',
     $           'tot_n_samp', 'n_coag', 'tot_n_coag', 't_est'
            write(6,'(i6,f8.1,i9,i11,i9,i11,f10.0)') loop, time, M,
     $           tot_n_samp, n_coag, tot_n_coag, t_est
! DEBUG
            call sum_int_2d(n_bin, n_fact, MS, M_comp)
            write(6,'(a9)'), 'M_comp'
            write(6,'(i9)'), M_comp
! DEBUG
         endif
      enddo

      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine compute_n_samp_super(n_bin, n_fact, MS, b1, f1, b2, f2
     $     , V_comp, k_max, del_t, n_samp_real)

      integer n_bin               ! INPUT: number of bins
      integer n_fact              ! INPUT: number of allowable factors
      integer MS(n_bin,n_fact)    ! INPUT: number of superparticles
      integer b1                  ! INPUT: first particle (bin number)
      integer f1                  ! INPUT: first particle (factor step)
      integer b2                  ! INPUT: second particle (bin number)
      integer f2                  ! INPUT: second particle (factor step)
      real*8 V_comp               ! INPUT: computational volume
      real*8 k_max(n_bin,n_fact,n_bin,n_fact) ! INPUT: maximum kernel values
      real*8 del_t                ! INPUT: timestep (s)
      real*8 n_samp_real          ! OUTPUT: number of samples per timestep
                                  !         for bin-i to bin-j events
      
      real*8 r_samp
      real*8 n_possible ! use real*8 to avoid integer overflow

      if ((b1 .eq. b2) .and. (f1 .eq. f2)) then
         if (f1 .eq. 1) then
            n_possible = dble(MS(b1, f1)) * (dble(MS(b1, f1)) - 1d0) /
     $           2d0
         else
            ! FIXME: this is only an approximation (better for larger f)
            n_possible = dble(MS(b1, f1)) * dble(MS(b1, f1)) / 2d0
         endif
      else
         n_possible = dble(MS(b1, f1)) * dble(MS(b2, f2)) / 2d0
      endif

      r_samp = k_max(b1, f1, b2, f2) * 1d0/V_comp * del_t
      n_samp_real = r_samp * n_possible

      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
