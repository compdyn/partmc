C Initial size distributions.

      module mod_init_dist
      contains

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

      return
      end subroutine

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine init_exp(MM, V_0, dlnr, n_bin, bin_v, bin_r, n_ini)

      integer MM           ! INPUT: physical dimension of V
      real*8 V_0           ! INPUT: mean volume of initial distribution (m^3)
      real*8 dlnr          ! INPUT: bin scale factor
      integer n_bin        ! INPUT: number of bins
      real*8 bin_v(n_bin)  ! INPUT: volume of particles in bins (m^3)
      real*8 bin_r(n_bin)  ! INPUT: radius of particles in bins (m)
      integer n_ini(n_bin) ! OUTPUT: initial number distribution

      integer k

      real*8 pi
      parameter (pi = 3.14159265358979323846d0)

      do k = 1,n_bin
         n_ini(k) = int(4d0 * pi * bin_r(k)**3 * dble(MM)/V_0
     &        * exp(-(bin_v(k) / V_0)) * dlnr)
      enddo

      end subroutine

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine init_log_normal(MM, d_mean, log_sigma, dlnr, n_bin,
     &     bin_v, bin_r, n_ini)

      integer MM           ! INPUT: physical dimension of V
      real*8 d_mean        ! INPUT: mean diameter of initial distribution (m)
      real*8 log_sigma     ! INPUT: log_e of the geometric standard
                           ! deviation of initial distribution (1)
      real*8 dlnr          ! INPUT: bin scale factor
      integer n_bin        ! INPUT: number of bins
      real*8 bin_v(n_bin)  ! INPUT: volume of particles in bins (m^3)
      real*8 bin_r(n_bin)  ! INPUT: radius of particles in bins (m)
      integer n_ini(n_bin) ! OUTPUT: initial number distribution

      real*8 pi
      parameter (pi = 3.14159265358979323846d0)

      integer k

      do k = 1,n_bin
          n_ini(k) = dble(MM) / (sqrt(2d0 * pi) * log_sigma) *
     &        dexp(-(dlog10(bin_r(k)) - dlog10(d))**2d0
     &        / (2d0 * log_sigma**2d0)) * dlnr / dlog(10)
      enddo

      ! The formula above was originally for a distribution in
      ! log_10(r), while we are using log_e(r). The division by dlog(10)
      ! at the end corrects for this.

      end subroutine

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      end module
