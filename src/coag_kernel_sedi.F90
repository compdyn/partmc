! Copyright (C) 2005-2011 Nicole Riemer and Matthew West
! Copyright (C) Andreas Bott
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_coag_kernel_sedi module.
!!
!! Contains code based on \c coad1d.f by Andreas Bott
!!     - http://www.meteo.uni-bonn.de/mitarbeiter/ABott/
!!     - Released under the GPL to Nicole Riemer (personal communication)
!!     - A. Bott, A flux method for the numerical solution of the
!!       stochastic collection equation, J. Atmos. Sci. 55, 2284-2293,
!!       1998.

!> Gravitational sedimentation coagulation kernel.
module pmc_coag_kernel_sedi

  use pmc_env_state
  use pmc_constants
  use pmc_aero_data
  use pmc_aero_particle

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Sedimentation coagulation kernel.
  subroutine kernel_sedi(aero_particle_1, aero_particle_2, &
       aero_data, env_state, k)

    !> First particle.
    type(aero_particle_t), intent(in) :: aero_particle_1
    !> Second particle.
    type(aero_particle_t), intent(in) :: aero_particle_2
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Environment state.
    type(env_state_t), intent(in) :: env_state
    !> Kernel \c k(a,b) (m^3/s).
    real(kind=dp), intent(out) :: k

    call kernel_sedi_helper(aero_particle_volume(aero_particle_1), &
         aero_particle_volume(aero_particle_2), aero_data, env_state%temp, &
         env_state%pressure, k)

  end subroutine kernel_sedi

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Minimum and maximum values of the sedimentation coagulation.
  subroutine kernel_sedi_minmax(v1, v2, aero_data, env_state, k_min, k_max)

    !> Volume of first particle (m^3).
    real(kind=dp), intent(in) :: v1
    !> Volume of second particle (m^3).
    real(kind=dp), intent(in) :: v2
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Environment state.
    type(env_state_t), intent(in) :: env_state
    !> Minimum kernel \c k(a,b) (m^3/s).
    real(kind=dp), intent(out) :: k_min
    !> Maximum kernel \c k(a,b) (m^3/s).
    real(kind=dp), intent(out) :: k_max

    call kernel_sedi_helper(v1, v2, aero_data, env_state%temp, &
         env_state%pressure, k_min)
    k_max = k_min

  end subroutine kernel_sedi_minmax

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Helper function that does the actual sedimentation kernel computation.
  !!
  !! Helper function. Do not call directly. Instead use kernel_sedi().
  subroutine kernel_sedi_helper(v1, v2, aero_data, temp, pressure, k)

    !> Volume of first particle (m^3).
    real(kind=dp), intent(in) :: v1
    !> Volume of second particle (m^3).
    real(kind=dp), intent(in) :: v2
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Temperature (K).
    real(kind=dp), intent(in) :: temp
    !> Pressure (Pa).
    real(kind=dp), intent(in) :: pressure
    !> Kernel k(a,b) (m^3/s).
    real(kind=dp), intent(out) :: k

    real(kind=dp) r1, r2, winf1, winf2, ec

    r1 = aero_data_vol2rad(aero_data, v1) ! m
    r2 = aero_data_vol2rad(aero_data, v2) ! m
    call fall_g(aero_data_vol_to_mobility_rad(aero_data, v1, temp, pressure), &
         winf1) ! winf1 in m/s
    call fall_g(aero_data_vol_to_mobility_rad(aero_data, v2, temp, pressure), &
         winf2) ! winf2 in m/s
    call effic(r1, r2, ec) ! ec is dimensionless
    k = ec * const%pi * (r1 + r2)**2 * abs(winf1 - winf2)

  end subroutine kernel_sedi_helper

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Finds the terminal velocity of a particle based on its size.
  subroutine fall_g(r, w_inf)

    !> Particle mobility radius (m).
    real(kind=dp), intent(in) :: r
    !> Terminal velocity (m/s).
    real(kind=dp), intent(out) :: w_inf

    ! terminal velocity of falling drops
    real(kind=dp) eta, xlamb, rhow, rhoa, grav, cunh, t0, sigma
    real(kind=dp) stok, stb, phy, py, rr, x, y, xrey, bond
    integer i
    real(kind=dp) b(7),c(6)
    data b /-0.318657d1,0.992696d0,-0.153193d-2,-0.987059d-3, &
         -0.578878d-3,0.855176d-4,-0.327815d-5/
    data c /-0.500015d1,0.523778d1,-0.204914d1,0.475294d0, &
         -0.542819d-1,0.238449d-2/

    eta = 1.818d-4
    xlamb = 6.62d-6
    rhow = 1d0
    rhoa = 1.225d-3
    grav = 980.665d0
    cunh = 1.257d0 * xlamb
    t0 = 273.15d0
    sigma = 76.1d0 - 0.155d0 * (293.15d0 - t0)
    stok = 2d0 * grav * (rhow - rhoa) / (9d0 * eta)
    stb = 32d0 * rhoa * (rhow - rhoa) * grav / (3d0 * eta * eta)
    phy = sigma * sigma * sigma * rhoa * rhoa  &
         / (eta**4 * grav * (rhow - rhoa))
    py = phy**(1d0/6d0)

    ! rr: radius in cm-units
    rr = r * 1d2

    if (rr .le. 1d-3) then
       w_inf = stok * (rr * rr + cunh * rr)
    elseif (rr .gt. 1d-3 .and. rr .le. 5.35d-2) then
       x = log(stb * rr * rr * rr)
       y = 0d0
       do i = 1,7
          y = y + b(i) * (x**(i - 1))
       end do
       xrey = (1d0 + cunh/rr) * exp(y)
       w_inf = xrey * eta / (2d0 * rhoa * rr)
    elseif (rr .gt. 5.35d-2) then
       bond = grav * (rhow - rhoa) * rr**2 / sigma
       if (rr .gt. 0.35d0) then
          bond = grav * (rhow - rhoa) * 0.35d0**2 / sigma
       end if
       x = log(16d0 * bond * py / 3d0)
       y = 0d0
       do i = 1,6
          y = y + c(i) * (x**(i - 1))
       end do
       xrey = py * exp(y)
       w_inf = xrey * eta / (2d0 * rhoa * rr)
       if (rr .gt. 0.35d0) then
          w_inf = xrey * eta / (2d0 * rhoa * 0.35d0)
       end if
    end if
    w_inf = w_inf / 100d0

  end subroutine fall_g

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Coagulation efficiency.
  !!
  !! Determines the chance that two particles will actually coagulate,
  !! given that they approach close enough to do so.
  subroutine effic(r1, r2, ec)

    !> Geometric radius of first particle (m).
    real(kind=dp), intent(in) :: r1
    !> Geometric radius of second particle (m).
    real(kind=dp), intent(in) :: r2
    !> Collision efficiency (dimensionless).
    real(kind=dp), intent(out) :: ec

    real(kind=dp) :: r_small, r_big, rq, p, q, ek
    integer :: k, ir, kk, iq
    ! collision efficiencies of hall kernel
    real(kind=dp) :: rat(21),r0(15),ecoll(15,21)

    data r0 /6.0d0,8.0d0,10.0d0,15.0d0,20.0d0,25.0d0,30.0d0,40.0d0 &
         ,50.0d0,60.0d0,70.0d0,100.0d0,150.0d0,200.0d0,300.0d0/
    data rat /0.0d0,0.05d0,0.1d0,0.15d0,0.2d0,0.25d0,0.3d0,0.35d0 &
         ,0.4d0,0.45d0,0.5d0,0.55d0,0.6d0,0.65d0,0.7d0,0.75d0,0.8d0 &
         ,0.85d0,0.9d0,0.95d0,1.0d0/
    ! two-dimensional linear interpolation of the collision efficiency
    data ecoll /0.001d0,0.001d0,0.001d0,0.001d0,0.001d0,0.001d0 &
         ,0.001d0,0.001d0,0.001d0,0.001d0 ,0.001d0,0.001d0,0.001d0 &
         ,0.001d0,0.001d0,0.003d0,0.003d0,0.003d0,0.004d0,0.005d0 &
         ,0.005d0,0.005d0,0.010d0,0.100d0,0.050d0,0.200d0,0.500d0 &
         ,0.770d0,0.870d0,0.970d0 ,0.007d0,0.007d0,0.007d0,0.008d0 &
         ,0.009d0,0.010d0,0.010d0,0.070d0,0.400d0,0.430d0 ,0.580d0 &
         ,0.790d0,0.930d0,0.960d0,1.000d0,0.009d0,0.009d0,0.009d0 &
         ,0.012d0,0.015d0 ,0.010d0,0.020d0,0.280d0,0.600d0,0.640d0 &
         ,0.750d0,0.910d0,0.970d0,0.980d0,1.000d0 ,0.014d0,0.014d0 &
         ,0.014d0,0.015d0,0.016d0,0.030d0,0.060d0,0.500d0,0.700d0 &
         ,0.770d0 ,0.840d0,0.950d0,0.970d0,1.000d0,1.000d0,0.017d0 &
         ,0.017d0,0.017d0,0.020d0,0.022d0 ,0.060d0,0.100d0,0.620d0 &
         ,0.780d0,0.840d0,0.880d0,0.950d0,1.000d0,1.000d0,1.000d0 &
         ,0.030d0,0.030d0,0.024d0,0.022d0,0.032d0,0.062d0,0.200d0 &
         ,0.680d0,0.830d0,0.870d0 ,0.900d0,0.950d0,1.000d0,1.000d0 &
         ,1.000d0,0.025d0,0.025d0,0.025d0,0.036d0,0.043d0 ,0.130d0 &
         ,0.270d0,0.740d0,0.860d0,0.890d0,0.920d0,1.000d0,1.000d0 &
         ,1.000d0,1.000d0 ,0.027d0,0.027d0,0.027d0,0.040d0,0.052d0 &
         ,0.200d0,0.400d0,0.780d0,0.880d0,0.900d0 ,0.940d0,1.000d0 &
         ,1.000d0,1.000d0,1.000d0,0.030d0,0.030d0,0.030d0,0.047d0 &
         ,0.064d0 ,0.250d0,0.500d0,0.800d0,0.900d0,0.910d0,0.950d0 &
         ,1.000d0,1.000d0,1.000d0,1.000d0 ,0.040d0,0.040d0,0.033d0 &
         ,0.037d0,0.068d0,0.240d0,0.550d0,0.800d0,0.900d0,0.910d0 &
         ,0.950d0,1.000d0,1.000d0,1.000d0,1.000d0,0.035d0,0.035d0 &
         ,0.035d0,0.055d0,0.079d0 ,0.290d0,0.580d0,0.800d0,0.900d0 &
         ,0.910d0,0.950d0,1.000d0,1.000d0,1.000d0,1.000d0 ,0.037d0 &
         ,0.037d0,0.037d0,0.062d0,0.082d0,0.290d0,0.590d0,0.780d0 &
         ,0.900d0,0.910d0 ,0.950d0,1.000d0,1.000d0,1.000d0,1.000d0 &
         ,0.037d0,0.037d0,0.037d0,0.060d0,0.080d0 ,0.290d0,0.580d0 &
         ,0.770d0,0.890d0,0.910d0,0.950d0,1.000d0,1.000d0,1.000d0 &
         ,1.000d0 ,0.037d0,0.037d0,0.037d0,0.041d0,0.075d0,0.250d0 &
         ,0.540d0,0.760d0,0.880d0,0.920d0 ,0.950d0,1.000d0,1.000d0 &
         ,1.000d0,1.000d0,0.037d0,0.037d0,0.037d0,0.052d0,0.067d0 &
         ,0.250d0,0.510d0,0.770d0,0.880d0,0.930d0,0.970d0,1.000d0 &
         ,1.000d0,1.000d0,1.000d0 ,0.037d0,0.037d0,0.037d0,0.047d0 &
         ,0.057d0,0.250d0,0.490d0,0.770d0,0.890d0,0.950d0 ,1.000d0 &
         ,1.000d0,1.000d0,1.000d0,1.000d0,0.036d0,0.036d0,0.036d0 &
         ,0.042d0,0.048d0 ,0.230d0,0.470d0,0.780d0,0.920d0,1.000d0 &
         ,1.020d0,1.020d0,1.020d0,1.020d0,1.020d0 ,0.040d0,0.040d0 &
         ,0.035d0,0.033d0,0.040d0,0.112d0,0.450d0,0.790d0,1.010d0 &
         ,1.030d0 ,1.040d0,1.040d0,1.040d0,1.040d0,1.040d0,0.033d0 &
         ,0.033d0,0.033d0,0.033d0,0.033d0 ,0.119d0,0.470d0,0.950d0 &
         ,1.300d0,1.700d0,2.300d0,2.300d0,2.300d0,2.300d0,2.300d0 &
         ,0.027d0,0.027d0,0.027d0,0.027d0,0.027d0,0.125d0,0.520d0 &
         ,1.400d0,2.300d0,3.000d0 ,4.000d0,4.000d0,4.000d0,4.000d0 &
         ,4.000d0/

    r_small = min(r1 * 1d6, r2 * 1d6) ! um
    r_big = max(r1 * 1d6, r2 * 1d6) ! um
    rq = r_small / r_big

    ir = 1
    do k = 1, 15
       if (r_big .gt. r0(k)) then
          ir = k + 1
       end if
    end do

    iq = 1
    do kk = 1,21
       if (rq .gt. rat(kk)) then
          iq = kk + 1
       end if
    end do

    if (ir .lt. 16) then
       if (ir .ge. 2) then
          p = (r_big - r0(ir - 1)) / (r0(ir) - r0(ir - 1))
          q = (rq - rat(iq - 1)) / (rat(iq) - rat(iq - 1))
          ec = (1d0 - p) * (1d0 - q) * ecoll(ir - 1, iq - 1) &
               + p * (1d0 - q) * ecoll(ir, iq - 1) &
               + q * (1d0 - p) * ecoll(ir - 1, iq) &
               + p * q * ecoll(ir, iq)
       else
          q = (rq - rat(iq - 1)) / (rat(iq) - rat(iq - 1))
          ec = (1d0 - q) * ecoll(1, iq - 1) + q * ecoll(1, iq)
       end if
    else
       q = (rq - rat(iq - 1)) / (rat(iq) - rat(iq - 1))
       ek = (1d0 - q) * ecoll(15, iq - 1) + q * ecoll(15, iq)
       ec = min(ek, 1d0)
    end if

    if (ec .lt. 1d-20) stop 99

  end subroutine effic

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_coag_kernel_sedi
