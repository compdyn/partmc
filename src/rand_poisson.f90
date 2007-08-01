! Copyright (C) 2007 Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.
!
! Routine to generate random numbers from a Poisson distribution.
!
! The method used at present is rather inefficient and inaccurate
! (brute force for below 10 and normal approximation above that
! point).
!
! The best known method appears to be due to Ahrens and Dieter (ACM
! Trans. Math. Software, 8(2), 163-179, 1982) and is available (in
! various forms) from:
!     http://www.netlib.org/toms/599
!     http://www.netlib.org/random/ranlib.f.tar.gz
!     http://users.bigpond.net.au/amiller/random/random.f90
!     http://www.netlib.org/random/random.f90
! Unfortunately this code is under the non-free license:
!     http://www.acm.org/pubs/copyright_policy/softwareCRnotice.html

module pmc_rand_poisson
contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function rand_poisson(mean)

    ! Generate a Poisson-distributed random number with the given
    ! mean. See http://en.wikipedia.org/wiki/Poisson_distribution for
    ! information on the method.

    use pmc_util

    real*8, intent(in) :: mean          ! mean of the distribution

    real*8 :: L, p
    integer :: k

    if (mean <= 10d0) then
       ! exact method due to Knuth
       L = exp(-mean)
       k = 0
       p = 1d0
       do
          k = k + 1
          p = p * util_rand()
          if (p < L) exit
       end do
       rand_poisson = k - 1
    else
       ! normal approximation with a continuity correction
       k = nint(rand_normal(mean - 0.5d0, sqrt(mean)))
       rand_poisson = max(k, 0)
    end if

  end function rand_poisson

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real*8 function rand_normal(mean, stddev)

    ! Generates a normally distributed random number with the given
    ! mean and standard deviation.

    use pmc_util
    use pmc_constants

    real*8, intent(in) :: mean          ! mean of the distribution
    real*8, intent(in) :: stddev        ! standard deviation of the distribution

    real*8 :: u1, u2, r, theta, z0, z1

    ! Uses the Box-Muller transform
    ! http://en.wikipedia.org/wiki/Box-Muller_transform
    u1 = util_rand()
    u2 = util_rand()
    r = sqrt(-2d0 * log(u1))
    theta = 2d0 * const%pi * u2
    z0 = r * cos(theta)
    z1 = r * sin(theta)
    ! z0 and z1 are now independent N(0,1) random variables
    ! We through away z1, but we could use a SAVE variable to only do
    ! the computation on every second call of this function.
    rand_normal = stddev * z0 + mean

  end function rand_normal

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_rand_poisson
