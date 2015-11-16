! Copyright (C) 2005-2011 Nicole Riemer and Matthew West
! Copyright (C) 2007 Richard Easter
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_coag_kernel_brown module.

!> Brownian coagulation kernel.
module pmc_coag_kernel_brown

  use pmc_env_state
  use pmc_constants
  use pmc_util
  use pmc_aero_particle
  use pmc_aero_data

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Compute the Brownian coagulation kernel.
  !!
  !! Uses equation (15.33) of M. Z. Jacobson, Fundamentals of Atmospheric
  !! Modeling Second Edition, Cambridge University Press, 2005.
  subroutine kernel_brown(aero_particle_1, aero_particle_2, &
       aero_data, env_state, k)

    !> First particle.
    type(aero_particle_t), intent(in) :: aero_particle_1
    !> Second particle.
    type(aero_particle_t), intent(in) :: aero_particle_2
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Environment state.
    type(env_state_t), intent(in) :: env_state
    !> Kernel k(a,b) (m^3/s).
    real(kind=dp), intent(out) :: k

    real(kind=dp) :: v1, v2, d1, d2

    v1 = aero_particle_volume(aero_particle_1)
    v2 = aero_particle_volume(aero_particle_2)
    d1 = aero_particle_density(aero_particle_1, aero_data)
    d2 = aero_particle_density(aero_particle_2, aero_data)

    call kernel_brown_helper(v1, d1, v2, d2, aero_data, env_state%temp, &
         env_state%pressure, k)

  end subroutine kernel_brown

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Compute the minimum and maximum Brownian coagulation kernel.
  !!
  !! Finds the minimum and maximum kernel values between particles of
  !! volumes v1 and v2, by sampling over possible densities.
  subroutine kernel_brown_minmax(v1, v2, aero_data, env_state, k_min, k_max)

    !> Volume of first particle (m^3).
    real(kind=dp), intent(in) :: v1
    !> Volume of second particle (m^3).
    real(kind=dp), intent(in) :: v2
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Environment state.
    type(env_state_t), intent(in) :: env_state
    !> Minimum kernel value (m^3/s).
    real(kind=dp), intent(out) :: k_min
    !> Maximum kernel value (m^3/s).
    real(kind=dp), intent(out) :: k_max

    !> Number of density sample points.
    integer, parameter :: n_sample = 3

    real(kind=dp) :: d1, d2, d_min, d_max, k
    integer :: i, j
    logical :: first

    d_min = minval(aero_data%density)
    d_max = maxval(aero_data%density)

    first = .true.
    do i = 1,n_sample
       do j = 1,n_sample
          d1 = interp_linear_disc(d_min, d_max, n_sample, i)
          d2 = interp_linear_disc(d_min, d_max, n_sample, j)
          call kernel_brown_helper(v1, d1, v2, d2, aero_data, &
               env_state%temp, env_state%pressure, k)
          if (first) then
             first = .false.
             k_min = k
             k_max = k
          else
             k_min = min(k_min, k)
             k_max = max(k_max, k)
          end if
       end do
    end do

  end subroutine kernel_brown_minmax

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Helper function that does the actual Brownian kernel computation.
  !!
  !! Helper function. Do not call directly. Instead use kernel_brown().
  !!
  !! Uses equation (15.33) of M. Z. Jacobson, Fundamentals of Atmospheric
  !! Modeling Second Edition, Cambridge University Press, 2005.
  subroutine kernel_brown_helper(vol_i, den_i, vol_j, den_j, aero_data, &
       temp, pressure, bckernel)

    !> Volume of first particle (m^3).
    real(kind=dp), intent(in) :: vol_i
    !> Density of first particle (kg/m^3).
    real(kind=dp), intent(in) :: den_i
    !> Volume of second particle (m^3).
    real(kind=dp), intent(in) :: vol_j
    !> Density of second particle (kg/m^3).
    real(kind=dp), intent(in) :: den_j
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Temperature (K).
    real(kind=dp), intent(in) :: temp
    !> Pressure (Pa).
    real(kind=dp), intent(in) :: pressure
    !> Kernel k(a,b) (m^3/s).
    real(kind=dp), intent(out) :: bckernel

    real(kind=dp) :: cunning, deltasq_i, &
         deltasq_j, diffus_i, diffus_j, diffus_sum, &
         freepath, gasfreepath, gasspeed, knud, rad_i, rad_j, &
         rad_sum, rhoair, speedsq_i, speedsq_j, tmp1, tmp2, &
         viscosd, viscosk, Rme_i, Rme_j

    ! rhoair  = air density (kg/m^3)
    ! viscosd = air dynamic viscosity (kg/m/s)
    ! viscosk = air kinematic viscosity (m^2/s)
    ! gasspeed    = air molecule mean thermal velocity (m/s)
    ! gasfreepath = air molecule mean free path (m)

    rhoair = (pressure * const%air_molec_weight) / &
         (const%univ_gas_const * temp)

    viscosd = 1.8325d-05 * (416.16d0 / (temp + 120d0)) * &
         (temp / 296.16d0)**1.5d0
    viscosk = viscosd / rhoair
    gasspeed = sqrt((8.0d0 * const%boltzmann * temp * const%avagadro) / &
         (const%pi * const%air_molec_weight))
    gasfreepath = 2d0 * viscosk / gasspeed

    ! coagulation kernel from eqn 15.33 of Jacobson (2005) text
    !
    ! diffus_i/j  = particle brownian diffusion coefficient  (m^2/s)
    ! speedsq_i/j = square of particle mean thermal velocity (m/s)
    ! freepath    = particle mean free path (m)
    ! cunning     = cunningham slip-flow correction factor
    ! deltasq_i/j = square of "delta_i" in eqn 15.34
    !
    ! bckernel   = brownian coagulation kernel (m3/s)

    rad_i     = aero_data_vol2rad(aero_data, vol_i)
    Rme_i     = aero_data_vol_to_mobility_rad(aero_data, vol_i, &
         temp, pressure)

    knud      = gasfreepath/Rme_i
    cunning   = 1d0 + knud*(1.249d0 + 0.42d0*exp(-0.87d0/knud))
    diffus_i  = (const%boltzmann * temp * cunning) / &
         (6.0d0 * const%pi * Rme_i * viscosd)
    speedsq_i = 8d0 * const%boltzmann * temp / (const%pi * den_i * vol_i)
    freepath  = 8d0*diffus_i/(const%pi*sqrt(speedsq_i))
    tmp1      = (2d0*Rme_i + freepath)**3
    tmp2      = (4d0*Rme_i*Rme_i + freepath*freepath)**1.5d0
    deltasq_i = ( (tmp1-tmp2)/(6d0*Rme_i*freepath) - 2d0*Rme_i )**2

    rad_j     = aero_data_vol2rad(aero_data, vol_j)
    Rme_j     = aero_data_vol_to_mobility_rad(aero_data, vol_j, &
         temp, pressure)

    knud      = gasfreepath/Rme_j
    cunning   = 1d0 + knud*(1.249d0 + 0.42d0*exp(-0.87d0/knud))
    diffus_j  = (const%boltzmann * temp * cunning) / &
         (6.0d0 * const%pi * Rme_j * viscosd)
    speedsq_j = 8d0 * const%boltzmann * temp / (const%pi * den_j * vol_j)
    freepath  = 8d0*diffus_j/(const%pi*sqrt(speedsq_j))
    tmp1      = (2d0*Rme_j + freepath)**3
    tmp2      = (4d0*Rme_j*Rme_j + freepath*freepath)**1.5d0
    deltasq_j = ( (tmp1-tmp2)/(6d0*Rme_j*freepath) - 2d0*Rme_j )**2

    rad_sum    = rad_i + rad_j
    diffus_sum = diffus_i + diffus_j
    tmp1       = rad_sum/(rad_sum + sqrt(deltasq_i + deltasq_j))
    tmp2       = 4d0*diffus_sum/(rad_sum*sqrt(speedsq_i + speedsq_j))
    bckernel  = 4d0*const%pi*rad_sum*diffus_sum/(tmp1 + tmp2)

  end subroutine kernel_brown_helper

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_coag_kernel_brown
