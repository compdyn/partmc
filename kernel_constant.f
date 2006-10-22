C Constant coagulation kernel.

      module mod_kernel_constant
      contains

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine kernel_constant(a, b, k)

      real*8, intent(in) :: a  ! volume of first particle
      real*8, intent(in) :: b  ! volume of second particle
      real*8, intent(out) :: k ! coagulation kernel

      real*8, parameter :: beta_0 = 0.25d0 / (60d0 * 2d8)

      k = beta_0

      end subroutine

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
 
      subroutine soln_constant_exp_cond(n_bin, bin_v, bin_r,
     &     bin_g, bin_n, dlnr,
     &     time, N_0, V_0, rho_p, V_comp)

      integer n_bin        ! INPUT: number of bins
      real*8 bin_v(n_bin)  ! INPUT: volume of particles in bins
      real*8 bin_r(n_bin)  ! INPUT: radius of particles in bins
      real*8 bin_g(n_bin)  ! OUTPUT: mass in bins
      integer bin_n(n_bin) ! OUTPUT: number in bins
      real*8 dlnr          ! INPUT: bin scale factor

      real*8 time          ! INPUT: cubin_rent time
      real*8 N_0           ! INPUT: particle number concentration (#/m^3)
      real*8 V_0           ! INPUT:
      real*8 rho_p         ! INPUT: particle density (kg/m^3)
      real*8 V_comp        ! INPUT: computational volume

      real*8 beta_0, tau, T, rat_v, nn, b, x, lambda, sigma
      integer k

      real*8 pi
      parameter (pi = 3.14159265358979323846d0)
      parameter (lambda = 1d0)

      call kernel_constant(1d0, 1d0, beta_0)

      if (time .eq. 0d0) then
         do k = 1,n_bin
            bin_n(k) = int(pi/2d0 * (2d0*bin_r(k))**3 * N_0/V_0
     &           * exp(-(bin_v(k)/V_0)))
         enddo
      else
         tau = N_0 * beta_0 * time
         do k = 1,n_bin
            rat_v = bin_v(k) / V_0
            x = 2d0 * rat_v / (tau + 2d0)
            nn = 4d0 * N_0 / (V_0 * ( tau + 2d0 ) ** 2d0)
     &           * exp(-2d0*rat_v/(tau+2d0)*exp(-lambda*tau)-lambda*tau)
            bin_n(k) = int(pi/2d0 * (2d0*bin_r(k))**3d0 * nn)
         enddo
      endif

      do k = 1,n_bin
         bin_g(k) = pi/6d0 * rho_p * (2d0*bin_r(k))**3d0
     &        * dble(bin_n(k))
      enddo

      do k = 1,n_bin
         bin_g(k) = bin_g(k) * dlnr * V_comp
         bin_n(k) = int(dble(bin_n(k)) * dlnr * V_comp)
      enddo

      end subroutine

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      end module
