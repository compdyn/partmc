! Copyright (C) 2005-2008 Nicole Riemer and Matthew West
! Copyright (C) 2007 Richard Easter
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.
    
!> \file
!> The pmc_kernel_brown module.

!> Brownian coagulation kernel.
module pmc_kernel_brown

  use pmc_env_state
  use pmc_constants
  use pmc_util
  
contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Compute the Brownian coagulation kernel.
  !!
  !! Uses equation (16.28) of M. Z. Jacobson, Fundamentals of
  !! Atmospheric Modeling, Cambridge University Press, 1999.
  subroutine kernel_brown(v1, v2, env_state, bckernel)

    !> Volume of first particle (m^3).
    real*8, intent(in) :: v1
    !> Volume of second particle (m^3).
    real*8, intent(in) :: v2
    !> Environment state.
    type(env_state_t), intent(in) :: env_state
    !> Kernel k(a,b) (m^3/s).
    real*8, intent(out) :: bckernel

    integer, parameter :: nbin_maxd = 1000
    integer, save :: nbin = 0
    integer :: i, j, k, m, n, lundiag1, lundiag2
    real*8, save :: rad_sv(nbin_maxd)

    real*8 :: avogad, bckernel1, boltz, cunning, deltasq_i, &
         deltasq_j, den_i, den_j, diffus_i, diffus_j, diffus_sum, &
         freepath, gasfreepath, gasspeed, knud, mwair, rad_i, rad_j, &
         rad_sum, rgas, rhoair, speedsq_i, speedsq_j, tmp1, tmp2, &
         viscosd, viscosk, vol_i, vol_j, tk, press, den_i_inp, den_j_inp

    ! boltz   = boltzmann's constant (erg/K = g*cm^2/s/K)
    ! avogad  = avogadro's number (molecules/mol)
    ! mwair   = molecular weight of air (g/mol)
    ! rgas    = gas constant (atmos/(mol/liter)/K)
    ! rhoair  = air density (g/cm^3)
    ! viscosd = air dynamic viscosity (g/cm/s)
    ! viscosk = air kinematic viscosity (cm^2/s)
    ! gasspeed    = air molecule mean thermal velocity (cm/s)
    ! gasfreepath = air molecule mean free path (cm)

    den_i_inp = 1.8d3
    den_j_inp = 1.8d3

    lundiag1 = -91
    lundiag2 = -92

    tk = env_state%temp
    press = env_state%pressure

    boltz = const%boltzmann * 1d7 ! J/K to erg/K
    avogad = const%avagadro
    mwair = const%air_molec_weight * 1d3 ! kg/mole to g/mole
    rgas = const%univ_gas_const * 1d-2 ! J/mole/K to ??? (FIXME)
    
    rhoair = 0.001d0 * ((press/1.01325d5)*mwair/(rgas*tk))
    
    viscosd = (1.8325d-04*(296.16d0+120d0)/(tk+120d0)) * (tk/296.16d0)**1.5d0
    viscosk = viscosd/rhoair
    gasspeed = sqrt(8d0*boltz*tk*avogad/(const%pi*mwair))
    gasfreepath = 2d0*viscosk/gasspeed

    ! following code attempts to construct the bin radius values
    !    by saving/organizing the input radius values 
    ! it is purely for diagnostic purposes
    i = -2
    j = -1
    if (lundiag2 > 0) then
       if (nbin == 0) rad_sv(:) = -1d0
       if (nbin == 0) write(*,*) '*** den_i,j =', den_i_inp, den_j_inp
       
       vol_i = v1 * 1.0d+6
       vol_j = v2 * 1.0d+6
       rad_i = vol2rad(vol_i)
       rad_j = vol2rad(vol_j)
       
       do k = 1, 2
          tmp1 = rad_i
          if (k == 2) tmp1 = rad_j
          m = -1
          do n = 1, nbin
             if (abs((tmp1/rad_sv(n))-1d0) <= 1.0d-5) then
                m = n
                exit
             end if
          end do
          if (m <= 0) then
             nbin = nbin + 1
             if (nbin > nbin_maxd) then
                write(*,*) '*** nbin > nbin_maxd'
                write(*,*) '    rad_i, rad_j =', rad_i, rad_j
                stop
             end if
             rad_sv(nbin) = tmp1
             m = nbin
          end if
          if (k == 1) i = m
          if (k == 2) j = m
       end do
    end if

    ! coagulation kernel from eqn 16.28 of jacobson (1999) text
    !
    ! diffus_i/j  = particle brownian diffusion coefficient  (cm^2/s)
    ! speedsq_i/j = square of particle mean thermal velocity (cm/s)
    ! freepath    = particle mean free path (cm)
    ! cunning     = cunningham slip-flow correction factor
    ! deltasq_i/j = square of "delta_i" in eqn 16.29d0
    !
    ! bckernel1   = brownian coagulation kernel (cm3/s)

    if (lundiag1 > 0) then
       write(lundiag1,'(/a,1p,2d12.4)') 'tk, patm', tk, (press/1.01325d5)
       write(lundiag1,'(a)')   &
            'i, j, coagcoef (cm3/s), dpwet_i,j (um), denswet_i,j (g/cm3)'
    end if
    
    den_i     = den_i_inp * 1.0d-3   ! particle wet density (g/cm3)
    vol_i     = v1 * 1.0d+6   ! particle wet volume (cm3)
    rad_i     = vol2rad(vol_i)       ! particle wet radius (cm)
    
    knud      = gasfreepath/rad_i  
    cunning   = 1d0 + knud*(1.249d0 + 0.42d0*exp(-0.87d0/knud))
    diffus_i  = boltz*tk*cunning/(6d0*const%pi*rad_i*viscosd)
    speedsq_i = 8d0*boltz*tk/(const%pi*den_i*vol_i)
    freepath  = 8d0*diffus_i/(const%pi*sqrt(speedsq_i))
    tmp1      = (2d0*rad_i + freepath)**3
    tmp2      = (4d0*rad_i*rad_i + freepath*freepath)**1.5d0
    deltasq_i = ( (tmp1-tmp2)/(6d0*rad_i*freepath) - 2d0*rad_i )**2
    
    den_j     = den_j_inp * 1.0d-3
    vol_j     = v2 * 1.0d+6
    rad_j     = vol2rad(vol_j)
    
    knud      = gasfreepath/rad_j  
    cunning   = 1d0 + knud*(1.249d0 + 0.42d0*exp(-0.87d0/knud))
    diffus_j  = boltz*tk*cunning/(6d0*const%pi*rad_j*viscosd)
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
    
    if ((lundiag1 > 0) .and. (i <= j)) then
       write(lundiag1,'(1p,2i4,5d12.4)')   &
            i, j, bckernel1, 2.0d4*rad_i, 2.0d4*rad_j, den_i, den_j
       write(lundiag1,'(1p,2i4,5d12.4)')   &
            i, j, bckernel, vol_i, v1, vol_j, v2
       if (lundiag2 > 0) then
          write(lundiag2,  '(1p,2i4,5d12.4)')   &
               i, j, bckernel1, 2.0d4*rad_i, 2.0d4*rad_j, den_i, den_j
       end if
    end if
    
  end subroutine kernel_brown
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
end module pmc_kernel_brown
