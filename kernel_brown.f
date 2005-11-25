C Brownian coagulation kernel.
C See Seinfeld, Atmospheric chemistry and physics of air pollution,
C page 396.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine kernel_brown(a, b, k)

      real*8 a  ! INPUT: volume of first particle
      real*8 b  ! INPUT: volume of second particle
      real*8 k  ! OUTPUT: coagulation kernel

      real*8 c_1, c_2, k_b, T, rho_p
      parameter (k_b = 1.38d-23) ! Boltzmann constant in J K-1
      parameter ( T = 298. )     ! Temperature of gas medium in K
      parameter (rho_p = 1800. ) ! particle density in kg m-3
      
      real*8 pi
      parameter (pi = 3.14159265358979323846d0)


      c_1 = (3./(4.*pi))**(1./6.)
      c_2 = (6*k_b*T/rho_p)**0.5

      k = c_1*c_2*(1/a+1/b)**0.5*(a**(1./3.)+b**(1./3))**2

      return
      end
