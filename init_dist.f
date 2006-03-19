C Initial size distributions.

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine init_bidisperse(MM, n_bin, n_ini)

      integer MM           ! INPUT: physical dimension of V
      integer n_bin        ! INPUT: number of bins
      integer n_ini(n_bin) ! OUTPUT: initial number distribution

      integer k

      if (MM .lt. 126) then
         write(*,*)'ERROR: MM too small for bidisperse'
         call exit(2)
      endif

      do k = 1,n_bin
         n_ini(k) = 0
      enddo
      n_ini(97) = MM - 1
      n_ini(126) = 1

      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine dist_to_n(N, dlnr, n_bin, bin_v, bin_r, n_den, bin_n)

C     Convert a number density (in ln(r)) to actual number of particles
C     in each bin.

      integer N            ! INPUT: total number of particles (approximate)
      real*8 dlnr          ! INPUT: bin scale factor
      integer n_bin        ! INPUT: number of bins
      real*8 bin_v(n_bin)  ! INPUT: volume of particles in bins (m^3)
      real*8 bin_r(n_bin)  ! INPUT: radius of particles in bins (m)
      real*8 n_den(n_bin)  ! INPUT: initial number density (#(ln(r))d(ln(r)))
      integer bin_n(n_bin) ! OUTPUT: number distribution

      integer k

      do k = 1,n_bin
         bin_n(k) = int(N * n_den(k) * dlnr)
      enddo

      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine init_exp(V_0, n_bin, bin_v, bin_r, n_den)

C     Exponential distribution in volume (or mass)
C     n(v) = N_0 / V_0 exp(- v / V_0)

      real*8 V_0           ! INPUT: mean volume of initial distribution (m^3)
      integer n_bin        ! INPUT: number of bins
      real*8 bin_v(n_bin)  ! INPUT: volume of particles in bins (m^3)
      real*8 bin_r(n_bin)  ! INPUT: radius of particles in bins (m)
      real*8 n_den(n_bin)  ! OUTPUT: initial number density (#(ln(r))d(ln(r)))

      integer k
      real*8 n_den_vol

      do k = 1,n_bin
         n_den_vol = 1d0 / V_0 * exp(-(bin_v(k) / V_0))
         call vol_to_lnr(bin_r(k), n_den_vol, n_den(k))
      enddo

      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
