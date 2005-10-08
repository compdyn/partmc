C Golovin coagulation kernel.

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine kernel_golovin(a, b, k)

      real*8 a  ! INPUT: volume of first particle
      real*8 b  ! INPUT: volume of second particle
      real*8 k  ! OUTPUT: coagulation kernel

      real*8 beta_1, beta_0
      parameter (beta_1 = 1000.)
      parameter (beta_0 = 1.d-11)

      k = beta_1 * (a + b)

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
