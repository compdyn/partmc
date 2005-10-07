C particle_array.f
C
C utility functions for handling V array of particle volumes

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine make_grid(n_bin, scal, rho_p, vv, dp, rr, dlnr)

      integer n_bin     ! INPUT: number of bins
      integer scal      ! INPUT: scale factor
      real*8 rho_p      ! INPUT: density
      real*8 vv(n_bin)  ! OUTPUT: volume of particles in bins
      real*8 dp(n_bin)  ! OUTPUT: diameter of particles in bins
      real*8 rr(n_bin)  ! OUTPUT: radius of particles in bins
      real*8 dlnr       ! OUTPUT: scale factor

      integer i
      real*8 ax, e(n_bin), r(n_bin), emin

      real*8 pi
      parameter (pi = 3.14159265358979323846)

      dlnr=dlog(2.d0)/(3.*scal)
      ax=2.d0**(1.0/scal)
      emin = 1.e-15

      e(1)=emin*0.5*(ax+1.)
      r(1)=1000.*dexp(dlog(3*e(1)/(4.*pi))/3.)
      vv(1)=1.e-06*e(1)/rho_p
      dp(1) = 1.e-06*r(1)
      rr(1) = dp(1)/2.
      
      do i=2,n_bin
         e(i)=ax*e(i-1)
         r(i)=1000.*dexp(dlog(3.*e(i)/(4.*pi))/3.)
         vv(i)=1.e-06*e(i)/rho_p
         dp(i)=1.e-06*2.*r(i)
         rr(i) = dp(i)/2.
      enddo

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine compute_volumes(n_bin, MM, n_ini, dp, dlnr, V)

      integer n_bin        ! INPUT: number of bins
      integer MM           ! INPUT: physical size of V
      real*8 n_ini(n_bin)  ! INPUT: initial number distribution
      real*8 dp(n_bin)     ! INPUT: diameter of particles in bins
      real*8 dlnr          ! INPUT: scale factor
      real*8 V(MM)         ! OUTPUT: particle volumes

      integer k, i
      real*8 sum_e, delta_n, sum_a

      real*8 pi
      parameter (pi = 3.14159265358979323846)

      sum_e = 0.
      do k = 1,n_bin
         delta_n = int(n_ini(k) * dlnr)
         sum_a = sum_e + 1
         sum_e = sum_e + delta_n
         do i = sum_a,sum_e
            V(i) = pi/6. * dp(k)**3.
         enddo
      enddo

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
      s2 = int(rand() * M_comp) + 1
      if ((s1 .gt. M_comp) .or. (V(s1) .eq. 0) .or.
     &     (s2 .gt. M_comp) .or. (V(s2) .eq. 0) .or.
     &     (s1 .eq. s2)) then
         goto 100
      endif

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
      if (rand() .gt. p ) goto 200

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
      V(s2) = 0.
      M = M - 1

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

      k_sum = 0.
      do i = 1,(n_samp**2)
         call find_rand_pair(MM, V, M_comp, s1, s2)
         call kernel(V(s1), V(s2), k)
         k_sum = k_sum + k
      enddo
      k_avg  = k_sum / (n_samp * (n_samp - 1))

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine compress(MM, M_comp, V)

      integer MM      ! INPUT: physical dimension of V
      integer M_comp  ! INPUT/OUTPUT: logical dimension of V
      real*8 V(MM)    ! INPUT/OUTPUT: on exit, all non-zero entries are
                      !               at the beginning, followed by all
                      !               zeros.

      integer i, i_w, i_v

      i_w = 1
      do i_v = 1,MM
         if (V(i_v) .ne. 0.) then
            V(i_w) = V(i_v)
            i_w = i_w + 1
         endif
      enddo
      M_comp = i_w

      do i = (i_w + 1),MM
         V(i) = 0.
      enddo

      return
      end


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine double(MM, M_comp, V, V_comp)

      integer MM      ! INPUT: physical size of V
      integer M_comp  ! INPUT/OUTPUT: logical size of V
      real*8 V(MM)    ! INPUT/OUTPUT: particle volume array
      real*8 V_comp   ! INPUT/OUTPUT: computational volume

      integer i

      call compress(MM, M_comp, V)
      do i = 1,M_comp
         V(i + M_comp) = V(i)
      enddo
      M_comp = 2 * M_comp
      V_comp = 2 * V_comp

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
      parameter (pi = 3.14159265358979323846)
      
      do k = 1,n_bin
         V_bin(k) = 4./3.*pi*rr(k)**3.
      enddo
      
      k_max = 0.
      do k = 1,n_bin
         if (n_ln(k)*dlnr .ge. 1.) then
            do ll = 1,k
               if (n_ln(ll)*dlnr .ge. 1.) then
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

      real*8 nv_conc
      real*8 V_0, d_0, vv_cnt, vv_conc
      integer NN_cnt, k, i

      real*8 pi, rho_w
      parameter (pi = 3.14159265358979323846)
      parameter (rho_w = 1000.)
      
      V_0 = 1.e-12
      d_0 = (6*V_0/pi)**(1./3.)
      
      do k=1,n_bin
         NN_cnt = 0
         vv_cnt = 0.
         do i=1,M_comp
            if ((V(i).ge. vv(k-1)) .and. (V(i) .lt. vv(k))) then
               NN_cnt = NN_cnt +1
               vv_cnt = vv_cnt + V(i)
            endif
         enddo
         nv_conc = NN_cnt/V_comp
         vv_conc = vv_cnt/V_comp
         n_ln(k) = nv_conc/dlnr
         g(k) = vv_conc/dlnr
      enddo
      
      return
      end
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      
      subroutine print_info(n_bin, time, dp, g, n_ln)

      integer n_bin        ! INPUT: number of bins
      real*8 time          ! INPUT: current simulation time
      real*8 dp(n_bin)     ! INPUT: diameter of particles in bins
      real*8 g(n_bin)      ! OUTPUT: total mass in each bin
      real*8 n_ln(n_bin)   ! OUTPUT: total number in each bin (log scaled)

      integer k
      real*8 seconds

      seconds = time - anint(time / 60) * 60
      if (abs(seconds) < 1.0e-10) then
         write(30,*)'Time = ', time
         do k = 1,n_bin
            write(30, '(i4,6e14.5)'), k, dp(k)/2., n_ln(k), g(k)
         enddo
      endif
      
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
