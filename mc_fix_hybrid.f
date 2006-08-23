C     Monte Carlo with fixed timestep and a hybrid array.

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine mc_fix_hybrid(MM, M, V, n_spec, n_bin, TDV, 
     &     MH, VH, V_comp,
     $     bin_v, rho_p, bin_r, bin_g, bin_gs, bin_n, dlnr, 
     &     kernel, t_max, t_print,
     $     t_progress, del_t, loop)

      integer MM           ! INPUT: physical dimension of V
      integer M            ! INPUT/OUTPUT: logical dimension of V
      integer n_spec       ! INPUT: number of species
      real*8 V(MM,n_spec)  ! INPUT/OUTPUT: particle volumes (m^3)
      integer n_bin        ! INPUT: number of bins
      integer TDV          ! INPUT: trailing dimension of VH
      integer MH(n_bin)    ! OUTPUT: number of particles per bin
      real*8 VH(n_bin,TDV,n_spec) ! OUTPUT: particle volumes
      real*8 V_comp        ! INPUT/OUTPUT: computational volume (m^3)
      real*8 rho_p(n_spec) ! INPUT: density of species (kg m^{-3})

      real*8 bin_v(n_bin)  ! INPUT: volume of particles in bins (m^3)
      real*8 bin_r(n_bin)  ! INPUT: radius of particles in bins (m)
      real*8 bin_g(n_bin)  ! OUTPUT: mass in bins  
      real*8 bin_gs(n_bin,n_spec) !OUTPUT: species mass in bins             
      integer bin_n(n_bin) ! OUTPUT: number in bins
      real*8 dlnr          ! INPUT: bin scale factor

      external kernel      ! INPUT: kernel function
      real*8 t_max         ! INPUT: final time (seconds)
      real*8 t_print       ! INPUT: interval to output data (seconds)
      real*8 t_progress    ! INPUT: interval to print progress (seconds)
      real*8 del_t         ! INPUT: timestep
      integer loop         ! INPUT: loop number of run

      real*8 time, last_print_time, last_progress_time
      real*8 k_max(n_bin, n_bin), n_samp_real
      integer n_samp, i_samp, n_coag, i, j, tot_n_samp, tot_n_coag, k
      logical do_print, do_progress, did_coag, bin_change
      real*8 t_start, t_end, t_est

      last_progress_time = 0d0
      time = 0d0
      tot_n_coag = 0
      call moments(MM, M, V, V_comp, n_spec, n_bin, bin_v, bin_r, bin_g,
     $     bin_gs, bin_n, dlnr)
      call check_event(time, t_print, last_print_time, do_print)
      write(6,*)'do_print ',do_print
      if (do_print) call print_info(time, V_comp, n_spec, n_bin, bin_v,
     $     bin_r,bin_g, bin_gs, bin_n, dlnr)

      call array_to_hybrid(MM, M, V, n_spec, n_bin, bin_v, TDV, MH, VH)
      call est_k_max_binned(n_bin, bin_v, kernel, k_max)

      call cpu_time(t_start)
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

! DEBUG
c         call check_hybrid(MM, M, n_bin, MH, VH, bin_v, bin_r)
! DEBUG

         time = time + del_t

         call check_event(time, t_print, last_print_time, do_print)
         if (do_print) call print_info(time, V_comp, n_spec, n_bin,
     $        bin_v, bin_r, bin_g, bin_gs, bin_n, dlnr)

         call check_event(time, t_progress, last_progress_time,
     $        do_progress)
         if (do_progress) then
            call cpu_time(t_end)
            t_est = (t_max - time) / time * (t_end - t_start)
            write(6,'(a6,a8,a9,a11,a9,a11,a10)') 'loop', 'time', 'M',
     $           'tot_n_samp', 'n_coag', 'tot_n_coag', 't_est'
            write(6,'(i6,f8.1,i9,i11,i9,i11,f10.0)') loop, time, M,
     $           tot_n_samp, n_coag, tot_n_coag, t_est
         endif

         write(6,*)'vor condensation'
         call condensation(n_bin, TDV, n_spec, MH, VH, rho_p)

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
         n_possible = dble(MH(i)) * (dble(MH(j)) - 1d0) / 2d0
      else
         n_possible = dble(MH(i)) * dble(MH(j)) / 2d0
      endif

      r_samp = k_max(i,j) * 1d0/V_comp * del_t
      n_samp_real = r_samp * n_possible

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
