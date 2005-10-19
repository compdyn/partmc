C     Utility functions for handling V array of particle volumes.
C
C     There are two different representations of particle size
C     distributions used throughout this code: a sectional
C     representation and an explicit particle representation.
C
C     The sectional representation stores the number and mass of
C     particles in bins, which are logarithmicly spaced. The bins are
C     described by the bin_v(n_bin) and bin_r(n_bin) arrays, which store the
C     volume and radius of the centerpoint of each bin. The variable
C     dlnr ... FIXME

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine make_grid(n_bin, scal, rho_p, bin_v, bin_r, dlnr)

      integer n_bin        ! INPUT: number of bins
      integer scal         ! INPUT: scale factor
      real*8 rho_p         ! INPUT: density
      real*8 bin_v(n_bin)  ! OUTPUT: volume of particles in bins
      real*8 bin_r(n_bin)  ! OUTPUT: radius of particles in bins
      real*8 dlnr          ! OUTPUT: scale factor

      integer i
      real*8 ax, e(n_bin), r(n_bin), emin

      real*8 pi
      parameter (pi = 3.14159265358979323846d0)

      dlnr = dlog(2d0) / (3d0 * scal)
      ax = 2d0**(1d0 / scal)
      emin = 1d-15

      ! FIXME: rewrite in a sane way
      do i = 1,n_bin
         ! mass (FIXME: units?)
         e(i) = emin * 0.5d0 * (ax + 1d0) * ax**(i - 1)
         ! radius (um)
         ! FIXME: following line assumes rho_p = 1000
         r(i) = 1000d0 * dexp(dlog(3d0 * e(i) / (4d0 * pi)) / 3d0)
         ! volume (FIXME: units?)
         bin_v(i) = 1d-6 * e(i) / rho_p
         ! radius (m)
         bin_r(i) = 1d-6 * r(i)
      enddo

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine compute_volumes(n_bin, MM, n_ini, bin_r, dlnr, V, M)

      integer n_bin        ! INPUT: number of bins
      integer MM           ! INPUT: physical size of V
      real*8 n_ini(n_bin)  ! INPUT: initial number distribution
      real*8 bin_r(n_bin)  ! INPUT: diameter of particles in bins
      real*8 dlnr          ! INPUT: scale factor
      real*8 V(MM)         ! OUTPUT: particle volumes
      integer M            ! OUTPUT: logical dimension of V

      integer k, i, sum_e, sum_a, delta_n

      real*8 pi
      parameter (pi = 3.14159265358979323846d0)

      sum_e = 0
      do k = 1,n_bin
         delta_n = int(n_ini(k) * dlnr)
         sum_a = sum_e + 1
         sum_e = sum_e + delta_n
         do i = sum_a,sum_e
            V(i) = 4d0/3d0 * pi * bin_r(k)**3
         enddo
      enddo

      M = sum_e

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
     
      subroutine find_rand_pair(M, s1, s2)
      
      integer M       ! INPUT: number of particles
      integer s1, s2  ! OUTPUT: s1 and s2 are not equal, random
                      !         particles with (1 <= s1,s2 <= M)

 100  s1 = int(rand() * M) + 1
      if ((s1 .lt. 1) .or. (s1 .gt. M)) goto 100
 101  s2 = int(rand() * M) + 1
      if ((s2 .lt. 1) .or. (s2 .gt. M)) goto 101
      if (s1 .eq. s2) goto 101

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
     
      subroutine find_rand_pair_acc_rej(MM, M, V, max_k, kernel,
     &     s1, s2)
      
      integer MM      ! INPUT: physical dimension of V
      integer M       ! INPUT: logical dimension of V
      real*8 V(MM)    ! INPUT: array of particle volumes
      real*8 max_k    ! INPUT: maximum value of the kernel
      external kernel ! INPUT: kernel function
      integer s1, s2  ! OUTPUT: s1 and s2 are not equal, random
                      !         particles with V(s1/s2) != 0

      real*8 k, p

 200  continue
      call find_rand_pair(M, s1, s2) ! test particles s1, s2
      call kernel(V(s1), V(s2), k)
      p = k / max_k     ! collision probability   
      if (dble(rand()) .gt. p ) goto 200

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine coagulate(MM, M, V, V_comp,
     &        n_bin, bin_v, bin_r, bin_g, bin_n, dlnr,
     &        s1, s2, bin_change)

      integer MM           ! INPUT: physical dimension of V
      integer M            ! INPUT/OUTPUT: logical dimension of V
      real*8 V(MM)         ! INPUT/OUTPUT: particle volumes
      real*8 V_comp        ! INPUT: computational volume

      integer n_bin        ! INPUT: number of bins
      real*8 bin_v(n_bin)  ! INPUT: volume of particles in bins
      real*8 bin_r(n_bin)  ! INPUT: radius of particles in bins
      real*8 bin_g(n_bin)  ! INPUT/OUTPUT: mass in bins
      integer bin_n(n_bin) ! INPUT/OUTPUT: number in bins
      real*8 dlnr          ! INPUT: bin scale factor

      integer s1           ! INPUT: first particle to coagulate
      integer s2           ! INPUT: second particle to coagulate
      logical bin_change   ! OUTPUT: whether an empty bin filled,
                           !         or a filled bin became empty

      integer k1, k2, kn

      bin_change = .false.

      ! remove s1 and s2 from bins
      call particle_in_bin(V(s1), n_bin, bin_v, k1)
      call particle_in_bin(V(s2), n_bin, bin_v, k2)
      bin_n(k1) = bin_n(k1) - 1
      bin_n(k2) = bin_n(k2) - 1
      bin_g(k1) = bin_g(k1) - V(s1)
      bin_g(k2) = bin_g(k2) - V(s2)

      if ((bin_n(k1) .lt. 0) .or. (bin_n(k2) .lt. 0)) then
         write(*,*)'ERROR: invalid bin_n'
         call exit(2)
      endif

      V(s1) = V(s1) + V(s2) ! add particle 2 onto particle 1
      V(s2) = V(M)          ! shift the last particle into empty slot
      M = M - 1             ! shorten array

      ! add new particle to bins
      call particle_in_bin(V(s1), n_bin, bin_v, kn)
      bin_n(kn) = bin_n(kn) + 1
      bin_g(kn) = bin_g(kn) + V(s1)

      if ((bin_n(k1) .eq. 0) .or. (bin_n(k2) .eq. 0))
     &     bin_change = .true.
      if ((bin_n(kn) .eq. 1) .and. (kn .ne. k1) .and. (kn .ne. k2))
     &     bin_change = .true.

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
     
      subroutine maybe_coag_pair(MM, M, V, V_comp,
     &           n_bin, bin_v, bin_r, bin_g, bin_n, dlnr,
     &           del_t, n_samp, kernel, did_coag, bin_change)

      integer MM           ! INPUT: physical dimension of V
      integer M            ! INPUT/OUTPUT: logical dimension of V
      real*8 V(MM)         ! INPUT/OUTPUT: particle volumes
      real*8 V_comp        ! INPUT: computational volume

      integer n_bin        ! INPUT: number of bins
      real*8 bin_v(n_bin)  ! INPUT: volume of particles in bins
      real*8 bin_r(n_bin)  ! INPUT: radius of particles in bins
      real*8 bin_g(n_bin)  ! INPUT/OUTPUT: mass in bins
      integer bin_n(n_bin) ! INPUT/OUTPUT: number in bins
      real*8 dlnr          ! INPUT: bin scale factor
      
      real*8 del_t         ! INPUT: timestep
      integer n_samp       ! INPUT: number of samples per timestep
      external kernel      ! INPUT: kernel function
      logical did_coag     ! OUTPUT: whether a coagulation occured
      logical bin_change   ! OUTPUT: whether bin structure changed

      integer s1, s2
      real*8 expo, p, k

      call find_rand_pair(M, s1, s2) ! test particles s1, s2
      call kernel(V(s1), V(s2), k)
      expo = k * 1d0/V_comp * del_t * (M*(M-1)/2d0) / n_samp
      p = 1d0 - exp(-expo) ! probability of coagulation
      bin_change = .false.
      if (dble(rand()) .lt. p) then
         call coagulate(MM, M, V, V_comp,
     &        n_bin, bin_v, bin_r, bin_g, bin_n, dlnr,
     &        s1, s2, bin_change)
         did_coag = .true.
      else
         did_coag = .false.
      endif

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine kernel_avg(MM, M, V, kernel, n_samp, k_avg)
      ! FIXME: use binned data instead

      integer MM      ! INPUT: physical dimension of V
      integer M       ! INPUT: logical dimension of V
      real*8 V(MM)    ! INPUT: array of particle volumes
      external kernel ! INPUT: kernel function
      integer n_samp  ! INPUT: number of samples to use (squared)
      real*8 k_avg    ! OUTPUT: estimated average of kernel values

      integer i, s1, s2
      real*8 k, k_sum

      k_sum = 0d0
      do i = 1,(n_samp**2)
         call find_rand_pair(M, s1, s2)
         call kernel(V(s1), V(s2), k)
         k_sum = k_sum + k
      enddo
      k_avg  = k_sum / (n_samp**2)

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine double(MM, M, V, V_comp,
     &     n_bin, bin_v, bin_r, bin_g, bin_n, dlnr)

      integer MM           ! INPUT: physical dimension of V
      integer M            ! INPUT/OUTPUT: logical dimension of V
      real*8 V(MM)         ! INPUT/OUTPUT: particle volumes
      real*8 V_comp        ! INPUT/OUTPUT: computational volume

      integer n_bin        ! INPUT: number of bins
      real*8 bin_v(n_bin)  ! INPUT: volume of particles in bins
      real*8 bin_r(n_bin)  ! INPUT: radius of particles in bins
      real*8 bin_g(n_bin)  ! INPUT/OUTPUT: mass in bins
      integer bin_n(n_bin) ! INPUT/OUTPUT: number in bins
      real*8 dlnr          ! INPUT: bin scale factor

      integer i

      ! only double if we have enough space to do so
      if (M .gt. MM / 2) then
         write(*,*)'ERROR: double without enough space'
         call exit(2)
      endif
      
      ! double V and associated structures
      do i = 1,M
         V(i + M) = V(i)
      enddo
      M = 2 * M
      V_comp = 2d0 * V_comp
      
      ! double bin structures
      do i = 1,n_bin
         bin_g(i) = bin_g(i) * 2d0
         bin_n(i) = bin_n(i) * 2
      enddo

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine est_k_max(n_bin, bin_v, bin_n, kernel, k_max)
      
      integer n_bin         ! INPUT: number of bins
      real*8 bin_v(n_bin)   ! INPUT: volume of particles in bins
      integer bin_n(n_bin)  ! INPUT: number in each bin
      external kernel       ! INPUT: kernel function
      real*8 k_max          ! OUTPUT: maximum kernel value

      real*8 k
      integer i, j
      logical use_bin(n_bin)

      ! use_bin starts as non-empty bins
      do i = 1,n_bin
         use_bin(i) = (bin_n(i) .gt. 0)
      enddo

      ! add all bins downstream of non-empty bins
      do i = 2,n_bin
         if (use_bin(i)) use_bin(i-1) = .true.
      enddo

      ! add all bins upstream of non-empty bins
      do i = (n_bin-1),1,-1
         if (use_bin(i)) use_bin(i+1) = .true.
      enddo
      
      k_max = 0d0
      do i = 1,n_bin
         if (use_bin(i)) then
            do j = 1,i
               if (use_bin(j)) then
                  call kernel(bin_v(i), bin_v(j), k)
                  if (k .gt. k_max) then
                     k_max = k
                  endif
               endif
            enddo
         endif
      enddo

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine est_k_avg(n_bin, bin_v, bin_n, kernel, k_avg)
      
      integer n_bin         ! INPUT: number of bins
      real*8 bin_v(n_bin)   ! INPUT: volume of particles in bins
      integer bin_n(n_bin)  ! INPUT: number in each bin
      external kernel       ! INPUT: kernel function
      real*8 k_avg          ! OUTPUT: average kernel value

      real*8 k
      integer i, j, div
      
      k_avg = 0d0
      div = 0
      do i = 1,n_bin
         if (bin_n(i) .gt. 0) then
            do j = 1,i
               if (bin_n(j) .gt. 0) then
                  call kernel(bin_v(i), bin_v(j), k)
                  k_avg = k_avg + k *  bin_n(i) * bin_n(j)
                  div = div + bin_n(i) * bin_n(j)
               endif
            enddo
         endif
      enddo

      k_avg = k_avg / div

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine particle_in_bin(v, n_bin, bin_v, k)
      ! FIXME: for log-spaced bins we can do this without search

      real*8 v             ! INPUT: volume of particle
      integer n_bin        ! INPUT: number of bins
      real*8 bin_v(n_bin)  ! INPUT: volume of particles in bins
      integer k            ! OUTPUT: bin number containing particle

      k = 0
 300  k = k + 1
c      write(*,*)'k,bin_v(k) = ', k, bin_v(k)
      if ((k .lt. n_bin) .and.
     &     (v .gt. (bin_v(k) + bin_v(k+1)) / 2d0)) goto 300

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine moments(MM, M, V, V_comp,
     &     n_bin, bin_v, bin_r, bin_g, bin_n, dlnr)

      integer MM           ! INPUT: physical dimension of V
      integer M            ! INPUT: logical dimension of V
      real*8 V(MM)         ! INPUT: particle volumes
      real*8 V_comp        ! INPUT: computational volume

      integer n_bin        ! INPUT: number of bins
      real*8 bin_v(n_bin)  ! INPUT: volume of particles in bins
      real*8 bin_r(n_bin)  ! INPUT: radius of particles in bins
      real*8 bin_g(n_bin)  ! OUTPUT: mass in bins
      integer bin_n(n_bin) ! OUTPUT: number in bins
      real*8 dlnr          ! INPUT: bin scale factor
      
      integer i, k

      do k = 1,n_bin
         bin_g(k) = 0d0
         bin_n(k) = 0
      enddo
      do i = 1,M
         call particle_in_bin(V(i), n_bin, bin_v, k)
         bin_g(k) = bin_g(k) + V(i)
         bin_n(k) = bin_n(k) + 1
      enddo

      return
      end
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine check_event(time, interval, last_time, do_event)

      real*8 time       ! INPUT: cubin_rent time
      real*8 interval   ! INPUT: how often the event should be done
      real*8 last_time  ! INPUT/OUTPUT: when the event was last done
      logical do_event  ! OUTPUT: whether the event should be done

      real*8 interval_below

      if (time .eq. 0d0) then
         last_time = 0d0
         do_event = .true.
      else
         interval_below = aint(time / interval) * interval
         if (last_time .lt. interval_below) then
            last_time = time
            do_event = .true.
         else
            do_event = .false.
         endif
      endif

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      
      subroutine print_header(n_loop, n_bin, n_time)

      integer n_loop  ! INPUT: number of loops
      integer n_bin   ! INPUT: number of bins
      integer n_time  ! INPUT: number of times

      write(30,'(a10,i10)') 'n_loop', n_loop
      write(30,'(a10,i10)') 'n_bin', n_bin
      write(30,'(a10,i10)') 'n_time', n_time

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      
      subroutine print_info(time, V_comp,
     &     n_bin, bin_v, bin_r, bin_g, bin_n, dlnr)

      real*8 time          ! INPUT: cubin_rent simulation time
      real*8 V_comp        ! INPUT: computational volume
      
      integer n_bin        ! INPUT: number of bins
      real*8 bin_v(n_bin)  ! INPUT: volume of particles in bins
      real*8 bin_r(n_bin)  ! INPUT: radius of particles in bins
      real*8 bin_g(n_bin)  ! INPUT: mass in bins
      integer bin_n(n_bin) ! INPUT: number in bins
      real*8 dlnr          ! INPUT: bin scale factor

      integer k

      write(30,'(a10,e14.5)') 'time', time
      do k = 1,n_bin
         write(30, '(i8,3e14.5)') k, bin_r(k), bin_n(k) / V_comp / dlnr,
     &        bin_g(k) / V_comp / dlnr
      enddo

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
