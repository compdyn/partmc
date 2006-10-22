C     Monte Carlo with fixed timestep and a split array.

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine mc_fix_split(MM, M, MS, V, V_comp, n_bin, s_bin, bin_v,
     $     bin_r, bin_g, bin_n, dlnr, kernel, t_max, t_print, t_progress
     $     , del_t, loop)

      integer MM           ! INPUT: physical dimension of V
      integer M            ! INPUT/OUTPUT: logical dimension of V
      integer MS           ! INPUT/OUTPUT: split point in V
      real*8 V(MM)         ! INPUT/OUTPUT: particle volumes (m^3)
      real*8 V_comp        ! INPUT/OUTPUT: computational volume (m^3)

      integer n_bin        ! INPUT: number of bins
      integer s_bin        ! INPUT: split point in bin grid
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

      real*8 time, last_print_time, last_progress_time
      real*8 k_max_small, k_max_big
      integer n_samp_small, n_samp_big, i_samp, n_coag
      logical do_print, do_progress, did_coag, bin_change
      real*8 t_start, t_end, t_est
! DEBUG
      integer call_small, call_big
      call_small = 0
      call_big = 0
! DEBUG

      last_progress_time = 0d0
      time = 0d0
      n_coag = 0
      call moments(MM, M, V, V_comp, n_bin, bin_v, bin_r, bin_g, bin_n,
     $     dlnr)
      call check_event(time, t_print, last_print_time, do_print)
      if (do_print) call print_info(time, V_comp, n_bin, bin_v, bin_r,
     $     bin_g, bin_n, dlnr)
      call est_k_max_split(n_bin, s_bin, bin_v, bin_n, kernel,
     $     k_max_small, k_max_big)

      call cpu_time(t_start)
      do while (time < t_max)
         call compute_n_samp_split(M, MS, V_comp, k_max_small, k_max_big
     $        , del_t, n_samp_small, n_samp_big)
         do i_samp = 1,n_samp_small
            call maybe_coag_pair_small(MM, M, MS, V, V_comp, n_bin,
     $           s_bin, bin_v, bin_r, bin_g, bin_n, dlnr, del_t,
     $           n_samp_small, kernel, did_coag, bin_change)
            if (did_coag) n_coag = n_coag + 1
            ! FIXME: just store up bin_change and do once outside loop
            if (bin_change) call est_k_max_split(n_bin, s_bin, bin_v,
     $           bin_n, kernel, k_max_small, k_max_big)
! DEBUG
            call_small = call_small + 1
! DEBUG
         enddo
         do i_samp = 1,n_samp_big
            call maybe_coag_pair_big(MM, M, MS, V, V_comp, n_bin, s_bin,
     $           bin_v, bin_r, bin_g, bin_n, dlnr, del_t, n_samp_big,
     $           kernel, did_coag, bin_change)
            if (did_coag) n_coag = n_coag + 1
            if (bin_change) call est_k_max_split(n_bin, s_bin, bin_v,
     $           bin_n, kernel, k_max_small, k_max_big)
! DEBUG
            call_big = call_big + 1
! DEBUG
         enddo
         if (M .lt. MM / 2) then
            call double_split(MM, M, MS, V, V_comp, n_bin, s_bin, bin_v,
     $           bin_r, bin_g, bin_n, dlnr)
         endif


! DEBUG
c         write(*,*)'call_small, call_big = ', call_small, call_big
c            call cpu_time(t_end)
c            t_est = (t_max - time) / time * (t_end - t_start)
c            write(6,'(a6,a8,a8,a8,a10,a10,a9,a10)') 'loop', 'time', 'M',
c     $           'M - MS', 'k_max_sml', 'k_max_big', 'n_coag', 't_est'
c            write(6,'(i6,f8.1,i8,i8,e10.3,e10.3,i9,f10)') loop, time, M,
c     $           M - MS, k_max_small, k_max_big, n_coag, t_est
c            write(*,'(a12,a12,a12)') 'n_samp_small', 'n_samp_big',
c     $           'n_samp_eff'
c            write(*,'(i12,i12,i12)') n_samp_small, n_samp_big,
c     $           n_samp_big + int(dble(n_samp_small) * max(k_max_small,
c     $           k_max_big) / k_max_small)
c         call exit(2)

c         call check_split(MM, M, MS, V, n_bin, s_bin, bin_v, bin_r)
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
            write(6,'(a6,a8,a8,a8,a10,a10,a9,a10)') 'loop', 'time', 'M',
     $           'M - MS', 'k_max_sml', 'k_max_big', 'n_coag', 't_est'
            write(6,'(i6,f8.1,i8,i8,e10.3,e10.3,i9,f10)') loop, time, M,
     $           M - MS, k_max_small, k_max_big, n_coag, t_est
            write(*,'(a12,a12,a12)') 'n_samp_small', 'n_samp_big',
     $           'n_samp_eff'
            write(*,'(i12,i12,i12)') n_samp_small, n_samp_big,
     $           n_samp_big + int(dble(n_samp_small) * max(k_max_small,
     $           k_max_big) / k_max_small)
         endif
      enddo

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine compute_n_samp_split(M, MS, V_comp, k_max_small,
     $     k_max_big, del_t, n_samp_small, n_samp_big)

      integer M            ! INPUT: number of particles
      integer MS           ! INPUT/OUTPUT: split point in V
      real*8 V_comp        ! INPUT: computational volume (m^3)
      real*8 k_max_small   ! INPUT: maximum small kernel value (m^3/s)
      real*8 k_max_big     ! INPUT: maximum big kernel value (m^3/s)
      real*8 del_t         ! INPUT: timestep (s)
      integer n_samp_small ! OUTPUT: number of small samples per timestep
      integer n_samp_big   ! OUTPUT: number of big samples per timestep

      real*8 r_samp
      real*8 M_small, M_big, n_possible ! use real*8 to avoid integer overflow

      M_small = MS
      M_big = M - MS

      ! small-small collisions
      r_samp = k_max_small * 1d0/V_comp * del_t
      n_possible = M_small * (M_small - 1d0) / 2d0
      n_samp_small = int(r_samp * n_possible)

      ! big-big and big-small collisions
      r_samp = k_max_big * 1d0/V_comp * del_t
      n_possible = M_big * (M_big - 1d0) / 2d0 + M_big * M_small
      n_samp_big = int(r_samp * n_possible)

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
