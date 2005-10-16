C     Utility functions for handling V array of particle volumes.
C
C     There are two different representations of particle size
C     distributions used throughout this code: a sectional
C     representation and an explicit particle representation.
C
C     The sectional representation stores the number and mass of
C     particles in bins, which are logarithmicly spaced. The bins are
C     described by the vv(n_bin) and rr(n_bin) arrays, which store the
C     volume and radius of the centerpoint of each bin. The variable
C     dlnr ... FIXME

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine make_grid(n_bin, scal, rho_p, vv, rr, dlnr)

      integer n_bin     ! INPUT: number of bins
      integer scal      ! INPUT: scale factor
      real*8 rho_p      ! INPUT: density
      real*8 vv(n_bin)  ! OUTPUT: volume of particles in bins
      real*8 rr(n_bin)  ! OUTPUT: radius of particles in bins
      real*8 dlnr       ! OUTPUT: scale factor

      integer i
      real*8 ax, e(n_bin), r(n_bin), emin

      real*8 pi
      parameter (pi = 3.14159265358979323846d0)

      dlnr = dlog(2d0) / (3d0 * scal)
      ax = 2d0**(1d0 / scal)
      emin = 1d-15

      do i = 1,n_bin
         ! mass (?)
         e(i) = emin * 0.5d0 * (ax + 1d0) * ax**(i - 1)
         ! radius (um)
         ! FIXME: following line assumes rho_p = 1000
         r(i) = 1000d0 * dexp(dlog(3d0 * e(i) / (4d0 * pi)) / 3d0)
         ! volume (?)
         vv(i) = 1d-6 * e(i) / rho_p
         ! radius (m)
         rr(i) = 1d-6 * r(i)
      enddo

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine compute_volumes(n_bin, MM, n_ini, rr, dlnr, V, M_comp)

      integer n_bin        ! INPUT: number of bins
      integer MM           ! INPUT: physical size of V
      real*8 n_ini(n_bin)  ! INPUT: initial number distribution
      real*8 rr(n_bin)     ! INPUT: diameter of particles in bins
      real*8 dlnr          ! INPUT: scale factor
      real*8 V(MM)         ! OUTPUT: particle volumes
      integer M_comp       ! OUTPUT: logical dimension of V

      integer k, i, sum_e, sum_a, delta_n

      real*8 pi
      parameter (pi = 3.14159265358979323846d0)

      sum_e = 0
      do k = 1,n_bin
         delta_n = int(n_ini(k) * dlnr)
         sum_a = sum_e + 1
         sum_e = sum_e + delta_n
         do i = sum_a,sum_e
            V(i) = pi/6d0 * (2d0*rr(k))**3
         enddo
      enddo

      M_comp = sum_e

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
     
      subroutine find_rand_pair(MM, V, M_comp, s1, s2)
      
      integer MM      ! INPUT: dimension of V
      real*8 V(MM)    ! INPUT: array of particle volumes
      integer M_comp  ! INPUT: maximum index of non-zero entries in V
      integer s1, s2  ! OUTPUT: s1 and s2 are not equal, random
                      !         particles with V(s1/s2) != 0

 100  s1 = int(rand() * M_comp) + 1
      if ((s1 .gt. M_comp) .or. (s1 .lt. 1) .or. (V(s1) .eq. 0d0))
     &     goto 100
 101  s2 = int(rand() * M_comp) + 1
      if ((s2 .gt. M_comp) .or. (s2 .lt. 1) .or. (V(s2) .eq. 0d0))
     &     goto 101
      if (s1 .eq. s2) goto 101

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
     
      subroutine find_rand_pair_acc_rej(MM, V, M_comp, max_k, kernel,
     &     s1, s2)
      
      integer MM      ! INPUT: dimension of V
      real*8 V(MM)    ! INPUT: array of particle volumes
      integer M_comp  ! INPUT: maximum index of non-zero entries in V
      real*8 max_k    ! INPUT: maximum value of the kernel
      external kernel ! INPUT: kernel function
      integer s1, s2  ! OUTPUT: s1 and s2 are not equal, random
                      !         particles with V(s1/s2) != 0

      real*8 k, p

 200  continue
      call find_rand_pair(MM, V, M_comp, s1, s2) ! test particles s1, s2
      call kernel(V(s1), V(s2), k)
      p = k / max_k     ! collision probability   
      if (dble(rand()) .gt. p ) goto 200

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine coagulate(MM, M, V, s1, s2)

      integer MM    ! INPUT: physical dimension of V
      integer M     ! INPUT/OUTPUT: number of particles
      real*8 V(MM)  ! INPUT/OUTPUT: array of particle sizes
      integer s1    ! INPUT: first particle to coagulate
      integer s2    ! INPUT: second particle to coagulate

      V(s1) = V(s1) + V(s2)          
      V(s2) = 0d0
      M = M - 1

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
     
      subroutine maybe_coag_pair(MM, V, M, M_comp, V_comp,
     &     del_T, n_samp, kernel, did_coag)
      
      integer MM       ! INPUT: physical dimension of V
      real*8 V(MM)     ! INPUT/OUTPUT: particle volumes
      integer M        ! INPUT/OUTPUT: number of particles
      integer M_comp   ! INPUT: logical dimension of V
      real*8 V_comp    ! INPUT: computational volume
      real*8 del_T     ! INPUT: timestep
      integer n_samp   ! INPUT: number of samples per timestep
      external kernel  ! INPUT: kernel function
      logical did_coag ! OUTPUT: whether a coagulation occured

      integer s1, s2
      real*8 expo, p, k

      call find_rand_pair(MM, V, M_comp, s1, s2) ! test particles s1, s2
      call kernel(V(s1), V(s2), k)
      expo = k * 1d0/V_comp * del_T * (M*(M-1)/2d0) / n_samp
      p = 1d0 - exp(-expo) ! probability of coagulation
      if (dble(rand()) .lt. p) then
         call coagulate(MM, M, V, s1, s2)
         did_coag = .true.
      else
         did_coag = .false.
      endif

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine kernel_avg(MM, M_comp, V, kernel, n_samp, k_avg)

      integer MM      ! INPUT: physical dimension of V
      integer M_comp  ! INPUT: logical dimension of V
      real*8 V(MM)    ! INPUT: array of particle volumes
      external kernel ! INPUT: kernel function
      integer n_samp  ! INPUT: number of samples to use (squared)
      real*8 k_avg    ! OUTPUT: estimated average of kernel values

      integer i, s1, s2
      real*8 k, k_sum

      k_sum = 0d0
      do i = 1,(n_samp**2)
         call find_rand_pair(MM, V, M_comp, s1, s2)
         call kernel(V(s1), V(s2), k)
         k_sum = k_sum + k
      enddo
      k_avg  = k_sum / (n_samp**2)

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine compress(MM, M, M_comp, V)

      integer MM      ! INPUT: physical dimension of V
      integer M       ! INPUT/OUTPUT: number of particles
      integer M_comp  ! INPUT/OUTPUT: logical dimension of V
      real*8 V(MM)    ! INPUT/OUTPUT: on exit, all non-zero entries are
                      !               at the beginning, followed by all
                      !               zeros.

      integer i, i_w, i_v

      ! group non-zero entries at the start of V
      i_w = 1
      do i_v = 1,MM
         if (V(i_v) .ne. 0d0) then
            V(i_w) = V(i_v)
            i_w = i_w + 1
         endif
      enddo
      M = i_w
      M_comp = i_w

      ! set remaining entries to zero
      do i = (i_w + 1),MM
         V(i) = 0d0
      enddo

      return
      end


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine double(MM, M, M_comp, V, V_comp)

      integer MM      ! INPUT: physical size of V
      integer M       ! INPUT/OUTPUT: number of particles
      integer M_comp  ! INPUT/OUTPUT: logical size of V
      real*8 V(MM)    ! INPUT/OUTPUT: particle volume array
      real*8 V_comp   ! INPUT/OUTPUT: computational volume

      integer i

      ! only double if we have enough space to do so
      if (2 * M .le. MM) then
         call compress(MM, M, M_comp, V)
         do i = 1,M_comp
            V(i + M_comp) = V(i)
         enddo
         M_comp = 2 * M_comp
         V_comp = 2 * V_comp
         M = 2 * M
      endif

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine est_k_max(n_bin, rr, n_ln, dlnr, kernel, k_max)
      
      integer n_bin        ! INPUT: number of bins
      real*8 rr(n_bin)     ! INPUT: radii of bins
      real*8 n_ln(n_bin)   ! INPUT: number in each bin
      real*8 dlnr          ! INPUT: scale factor
      external kernel      ! INPUT: kernel function
      real*8 k_max         ! OUTPUT: maximum kernel value

      real*8 V_bin(n_bin), cck
      integer ll, k

      real*8 pi
      parameter (pi = 3.14159265358979323846d0)
      
      do k = 1,n_bin
         V_bin(k) = 4d0 / 3d0 * pi * rr(k)**3
      enddo
      
      k_max = 0d0
      do k = 1,n_bin
         if (n_ln(k) * dlnr .ge. 1d0) then
            do ll = 1,k
               if (n_ln(ll) * dlnr .ge. 1d0) then
                  call kernel(V_bin(k), V_bin(ll), cck)
                  if (cck .gt. k_max) then
                     k_max = cck
                  endif
               endif
            enddo
         endif
      enddo
      
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      
      subroutine moments(MM, V, n_bin, M_comp, V_comp,
     &     vv, dlnr, g, n_ln)
      
      integer MM           ! INPUT: dimension of V
      real*8 V(MM)         ! INPUT: particle volume array
      integer n_bin        ! INPUT: number of bins
      integer M_comp       ! INPUT: maximum index of particle in V
      real*8 V_comp        ! INPUT: computational volume
      real*8 vv(n_bin)     ! INPUT: volumes of particles in bins
      real*8 dlnr          ! INPUT: scale factor
      real*8 g(n_bin)      ! OUTPUT: total mass in each bin
      real*8 n_ln(n_bin)   ! OUTPUT: total number in each bin (log scaled)

      real*8 nv_conc, vv_cnt, vv_conc, vv_low, vv_high
      integer NN_cnt, k, i

      do k = 1,n_bin
         NN_cnt = 0
         vv_cnt = 0d0
         if (k .eq. 1) then
            vv_low = vv(k) - 1d100
            vv_high = (vv(k) + vv(k+1)) / 2d0
         elseif (k .eq. n_bin) then
            vv_low = (vv(k-1) + vv(k)) / 2d0
            vv_high = vv(k) + 1d100
         else
            vv_low = (vv(k-1) + vv(k)) / 2d0
            vv_high = (vv(k) + vv(k+1)) / 2d0
         endif
         do i = 1,M_comp
            if ((V(i) .ge. vv_low) .and. (V(i) .lt. vv_high)) then
               NN_cnt = NN_cnt + 1
               vv_cnt = vv_cnt + V(i)
            endif
         enddo
         nv_conc = NN_cnt / V_comp
         vv_conc = vv_cnt / V_comp
         n_ln(k) = nv_conc / dlnr
         g(k) = vv_conc / dlnr
      enddo

      return
      end
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine check_event(time, interval, last_time, do_event)

      real*8 time       ! INPUT: current time
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
      
      subroutine print_info(n_bin, time, rr, g, n_ln)

      integer n_bin        ! INPUT: number of bins
      real*8 time          ! INPUT: current simulation time
      real*8 rr(n_bin)     ! INPUT: radius of particles in bins
      real*8 g(n_bin)      ! OUTPUT: total mass in each bin
      real*8 n_ln(n_bin)   ! OUTPUT: total number in each bin (log scaled)

      integer k

      write(30,'(a10,e14.5)') 'time', time
      do k = 1,n_bin
         write(30, '(i8,3e14.5)') k, rr(k), n_ln(k), g(k)
      enddo
      
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
