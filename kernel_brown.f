C Brownian coagulation kernel.
C See Seinfeld, Atmospheric chemistry and physics of air pollution,
C page 396.
    
      module mod_kernel_brown
      contains

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine kernel_brown(a, b, k)

      real*8 a  ! INPUT: volume of first particle
      real*8 b  ! INPUT: volume of second particle
      real*8 k  ! OUTPUT: coagulation kernel

      real*8 c_1,  k_b, T, mu, a_third, b_third
      parameter (k_b = 1.38d-23)   ! Boltzmann constant in J K^{-1}
      parameter (mu = 1.78d-5 )    ! dynamic viscosity of air in kg m^{-1} s^{-1}
      parameter (T = 288.d0)       ! Temperature of gas medium in K
      
      real*8 pi
      parameter (pi = 3.14159265358979323846d0)


      c_1 = 2d0 * k_b * T / (3.d0 * mu)
      a_third = a**(1.d0/3.d0)
      b_third = b**(1.d0/3.d0)

      k = c_1 * (a_third + b_third) *(1d0/a_third + 1d0/b_third) 

      end subroutine

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      end module

