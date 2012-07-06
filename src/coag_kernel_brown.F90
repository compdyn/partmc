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
  use pmc_fractal
  
contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Compute the Brownian coagulation kernel.
  !!
  !! Uses equation (16.28) of M. Z. Jacobson, Fundamentals of
  !! Atmospheric Modeling, Cambridge University Press, 1999.
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

    call kernel_brown_helper(aero_data, v1, d1, v2, d2, env_state%temp, &
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
          call kernel_brown_helper(aero_data, v1, d1, v2, d2, &
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
  !! Uses equation (16.28) of M. Z. Jacobson, Fundamentals of
  !! Atmospheric Modeling, Cambridge University Press, 1999.
  subroutine kernel_brown_helper(aero_data, v1, d1, v2, d2, tk, press, bckernel)

    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Volume of first particle (m^3).
    real(kind=dp), intent(in) :: v1
    !> Density of first particle (kg/m^3).
    real(kind=dp), intent(in) :: d1
    !> Volume of second particle (m^3).
    real(kind=dp), intent(in) :: v2
    !> Density of second particle (kg/m^3).
    real(kind=dp), intent(in) :: d2
    !> Temperature (K).
    real(kind=dp), intent(in) :: tk
    !> Pressure (Pa).
    real(kind=dp), intent(in) :: press
    !> Kernel k(a,b) (m^3/s).
    real(kind=dp), intent(out) :: bckernel

    integer, parameter :: nbin_maxd = 1000
    integer, save :: nbin = 0
    real(kind=dp), save :: rad_sv(nbin_maxd)
    real(kind=dp) :: avogad, bckernel1, boltz, cunning, deltasq_i, &
         deltasq_j, den_i, den_j, diffus_i, diffus_j, diffus_sum, &
         freepath, gasfreepath, gasspeed, knud, mwair, rad_i, rad_j, &
         rad_sum, rgas, rhoair, speedsq_i, speedsq_j, tmp1, tmp2, &
         viscosd, viscosk, vol_i, vol_j, Rme_i, Rme_j

    ! boltz   = boltzmann's constant (erg/K = g*cm^2/s/K)
    ! avogad  = avogadro's number (molecules/mol)
    ! mwair   = molecular weight of air (g/mol)
    ! rgas    = gas constant (atmos/(mol/liter)/K)
    ! rhoair  = air density (g/cm^3)
    ! viscosd = air dynamic viscosity (g/cm/s)
    ! viscosk = air kinematic viscosity (cm^2/s)
    ! gasspeed    = air molecule mean thermal velocity (cm/s)
    ! gasfreepath = air molecule mean free path (cm)

    boltz = const%boltzmann * 1d7 ! J/K to erg/K
    avogad = const%avagadro
    mwair = const%air_molec_weight * 1d3 ! kg/mole to g/mole
    rgas = const%univ_gas_const * 1d-2 ! J/mole/K to atmos/(mol/liter)/K
    
    rhoair = 0.001d0 * ((press/1.01325d5)*mwair/(rgas*tk))
    
    viscosd = (1.8325d-04*(296.16d0+120d0)/(tk+120d0)) * (tk/296.16d0)**1.5d0
    viscosk = viscosd/rhoair
    gasspeed = sqrt(8d0*boltz*tk*avogad/(const%pi*mwair))
    gasfreepath = 2d0*viscosk/gasspeed

    ! coagulation kernel from eqn 16.28 of jacobson (1999) text
    !
    ! diffus_i/j  = particle brownian diffusion coefficient  (cm^2/s)
    ! speedsq_i/j = square of particle mean thermal velocity (cm/s)
    ! freepath    = particle mean free path (cm)
    ! cunning     = cunningham slip-flow correction factor
    ! deltasq_i/j = square of "delta_i" in eqn 16.29d0
    !
    ! bckernel1   = brownian coagulation kernel (cm3/s)

    den_i     = d1 * 1.0d-3   ! particle wet density (g/cm3)
    vol_i     = v1 * 1.0d+6   ! particle wet volume (cm3)
    rad_i     = vol2rad(vol_i, aero_data%fractal)  ! particle wet radius (cm)
    Rme_i     = vol2Rme(vol_i/1d+6, tk, press, aero_data%fractal)*1d+2
    
    knud      = gasfreepath/Rme_i
    cunning   = 1d0 + knud*(1.249d0 + 0.42d0*exp(-0.87d0/knud))
    diffus_i  = boltz*tk*cunning/(6d0*const%pi*Rme_i*viscosd)
    speedsq_i = 8d0*boltz*tk/(const%pi*den_i*vol_i)
    freepath  = 8d0*diffus_i/(const%pi*sqrt(speedsq_i))
    tmp1      = (2d0*rad_i + freepath)**3
    tmp2      = (4d0*rad_i*rad_i + freepath*freepath)**1.5d0
    deltasq_i = ( (tmp1-tmp2)/(6d0*rad_i*freepath) - 2d0*rad_i )**2
    
    den_j     = d2 * 1.0d-3
    vol_j     = v2 * 1.0d+6
    rad_j     = vol2rad(vol_j, aero_data%fractal)
    Rme_j     = vol2Rme(vol_j/1d+6, tk, press, aero_data%fractal)*1d+2
    
    knud      = gasfreepath/Rme_j
    cunning   = 1d0 + knud*(1.249d0 + 0.42d0*exp(-0.87d0/knud))
    diffus_j  = boltz*tk*cunning/(6d0*const%pi*Rme_j*viscosd)
    speedsq_j = 8d0*boltz*tk/(const%pi*den_j*vol_j)
    freepath  = 8d0*diffus_j/(const%pi*sqrt(speedsq_j))
    tmp1      = (2d0*rad_j + freepath)**3
    tmp2      = (4d0*rad_j*rad_j + freepath*freepath)**1.5d0
    deltasq_j = ( (tmp1-tmp2)/(6d0*rad_j*freepath) - 2d0*rad_j )**2
    
    rad_sum    = rad_i + rad_j
    diffus_sum = diffus_i + diffus_j 
    tmp1       = rad_sum/(rad_sum + sqrt(deltasq_i + deltasq_j))
    tmp2       = 4d0*diffus_sum/(rad_sum*sqrt(speedsq_i + speedsq_j))
    bckernel1  = 4d0*const%pi*rad_sum*diffus_sum/(tmp1 + tmp2)
    
    bckernel   = bckernel1 * 1.0d-6
    
  end subroutine kernel_brown_helper
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
end module pmc_coag_kernel_brown
