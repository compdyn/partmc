! Copyright (C) 2005-2007 Nicole Riemer and Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.
!
! Brownian coagulation kernel.
! See Seinfeld, Atmospheric chemistry and physics of air pollution,
! page 394 (equation 10.18)
! This expression is based on the assumption that the continuum regime applies.
    
module mod_kernel_brown
contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine kernel_brown(v1, v2, env, k)

    ! Compute the Brownian coagulation kernel.

    use mod_environ
    use mod_constants

    real*8, intent(in) :: v1            ! volume of first particle (m^3)
    real*8, intent(in) :: v2            ! volume of second particle (m^3)
    type(environ), intent(in) :: env    ! environment state
    real*8, intent(out) :: k            ! kernel k(a,b) (m^3/s)

!   real*8 c_1, a_third, b_third
    real*8 dens1, dens2                 ! particle densities (kg/m^3)
    integer lundiag1, lundiag2

!   c_1 = 2d0 * const%k_b * env%T / (3.d0 * const%mu)
!   a_third = v1**(1.d0/3.d0)
!   b_third = v2**(1.d0/3.d0)
    
!   k = c_1 * (a_third + b_third) *(1d0/a_third + 1d0/b_third) 
    

! calculate brownian kernel using coagsolv code
    dens1 = 1.8d3
    dens2 = 1.8d3
    lundiag1 = -91
    lundiag2 = -92
    call brownian_kernel( v1, v2, dens1, dens2, env%t, env%p,   &
                          lundiag1, lundiag2, k )


  end subroutine kernel_brown
      

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      subroutine brownian_kernel(   &
        vol_i_inp, vol_j_inp, den_i_inp, den_j_inp,   &
        tk, press, lundiag1, lundiag2, bckernel )
!
! this routine calculates brownian coagulation kernel
! using on eqn 16.28 of   
!    jacobson,  m. z. (1999) fundamentals of atmospheric modeling.
!       cambridge university press, new york, 656 pp.
!
! r. easter, june 2007
!
      implicit none

! arguments
      real*8,  intent(in) :: vol_i_inp, vol_j_inp  ! wet (ambient) particle volumes (m3)
      real*8,  intent(in) :: den_i_inp, den_j_inp  ! wet (ambient) particle densities (kg/m3)
      real*8,  intent(in) :: tk                    ! air temperature (K)
      real*8,  intent(in) :: press                 ! air pressure (Pa)
      integer, intent(in) :: lundiag1, lundiag2    ! logical units for diagnostic output

      real*8,  intent(out) :: bckernel             ! brownian coag kernel (m3/s)

! local variables
      integer, parameter :: nbin_maxd = 1000
      integer, save :: nbin = 0
      integer :: i, j
      integer :: k, m, n

      real*8, save :: rad_sv(nbin_maxd)
      real*8 ::   &
             avogad,   &
             bckernel1, boltz,   &
             cunning,   &
             deltasq_i, deltasq_j,   &
             den_i, den_j,   &
             diffus_i, diffus_j, diffus_sum,   &
             freepath,   &
             gasfreepath, gasspeed,   &
             knud,   &
             mwair,   &
             pi, pi43inv,   &
             rad_i, rad_j, rad_sum,   &
             rgas, rhoair,   &
             speedsq_i, speedsq_j,   &
             third, tmp1, tmp2,   &
             viscosd, viscosk,   &
             vol_i, vol_j

!
! boltz   = boltzmann's constant (erg/K = g*cm2/s/K)
! avogad  = avogadro's number (molecules/mol)
! mwair   = molecular weight of air (g/mol)
! rgas    = gas constant (atmos/(mol/liter)/K)
! rhoair  = air density (g/cm3)
! viscosd = air dynamic viscosity (g/cm/s)
! viscosk = air kinematic viscosity (cm2/s)
! gasspeed    = air molecule mean thermal velocity (cm/s)
! gasfreepath = air molecule mean free path (cm)
!
      third  = 1.0/3.0
      pi      = 3.1415926536
      pi43inv = 3.0/(4.0*pi)

      boltz  = 1.38054e-16
      avogad = 6.02252e+23
      mwair  = 28.966
      rgas = 0.08206

      rhoair    = 0.001 * ((press/1.01325d5)*mwair/(rgas*tk))

      viscosd = (1.8325e-04*(296.16+120.0)/(tk+120.0)) * (tk/296.16)**1.5
      viscosk = viscosd/rhoair
      gasspeed = sqrt(8.0*boltz*tk*avogad/(pi*mwair))
      gasfreepath = 2.0*viscosk/gasspeed

!
! following code attempts to construct the bin radius values
!    by saving/organizing the input radius values 
! it is purely for diagnostic purposes
!
      i = -2
      j = -1
      if (lundiag2 > 0) then
         if (nbin == 0) rad_sv(:) = -1.0
         if (nbin == 0) write(*,*) '*** den_i,j =', den_i_inp, den_j_inp

         vol_i = vol_i_inp * 1.0d+6
         vol_j = vol_j_inp * 1.0d+6
         rad_i = (vol_i*pi43inv)**third
         rad_j = (vol_j*pi43inv)**third

         do k = 1, 2
            tmp1 = rad_i
            if (k == 2) tmp1 = rad_j
            m = -1
            do n = 1, nbin
               if (abs((tmp1/rad_sv(n))-1.0) <= 1.0e-5) then
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

!
! coagulation kernel from eqn 16.28 of jacobson (1999) text
!
! diffus_i/j  = particle brownian diffusion coefficient  (cm2/s)
! speedsq_i/j = square of particle mean thermal velocity (cm/s)
! freepath    = particle mean free path (cm)
! cunning     = cunningham slip-flow correction factor
! deltasq_i/j = square of "delta_i" in eqn 16.29
!
! bckernel1   = brownian coagulation kernel (cm3/s)
!
       if (lundiag1 > 0) then
          write(lundiag1,'(/a,1p,2e12.4)') 'tk, patm', tk, (press/1.01325d5)
          write(lundiag1,'(a)')   &
          'i, j, coagcoef (cm3/s), dpwet_i,j (um), denswet_i,j (g/cm3)'
       end if

       den_i     = den_i_inp * 1.0d-3   ! particle wet density (g/cm3)
       vol_i     = vol_i_inp * 1.0d+6   ! particle wet volume (cm3)
       rad_i     = (vol_i*pi43inv)**third   ! particle wet radius (cm)

       knud      = gasfreepath/rad_i  
       cunning   = 1.0 + knud*(1.249 + 0.42*exp(-0.87/knud))
       diffus_i  = boltz*tk*cunning/(6.0*pi*rad_i*viscosd)
       speedsq_i = 8.0*boltz*tk/(pi*den_i*vol_i)
       freepath  = 8.0*diffus_i/(pi*sqrt(speedsq_i))
       tmp1      = (2.0*rad_i + freepath)**3
       tmp2      = (4.0*rad_i*rad_i + freepath*freepath)**1.5
       deltasq_i = ( (tmp1-tmp2)/(6.0*rad_i*freepath) - 2.0*rad_i )**2

       den_j     = den_j_inp * 1.0d-3
       vol_j     = vol_j_inp * 1.0d+6
       rad_j     = (vol_j*pi43inv)**third

       knud      = gasfreepath/rad_j  
       cunning   = 1.0 + knud*(1.249 + 0.42*exp(-0.87/knud))
       diffus_j  = boltz*tk*cunning/(6.0*pi*rad_j*viscosd)
       speedsq_j = 8.0*boltz*tk/(pi*den_j*vol_j)
       freepath  = 8.0*diffus_j/(pi*sqrt(speedsq_j))
       tmp1      = (2.0*rad_j + freepath)**3
       tmp2      = (4.0*rad_j*rad_j + freepath*freepath)**1.5
       deltasq_j = ( (tmp1-tmp2)/(6.0*rad_j*freepath) - 2.0*rad_j )**2

       rad_sum    = rad_i + rad_j
       diffus_sum = diffus_i + diffus_j 
       tmp1       = rad_sum/(rad_sum + sqrt(deltasq_i + deltasq_j))
       tmp2       = 4.0*diffus_sum/(rad_sum*sqrt(speedsq_i + speedsq_j))
       bckernel1  = 4.0*pi*rad_sum*diffus_sum/(tmp1 + tmp2)

       bckernel   = bckernel1 * 1.0d-6

       if ((lundiag1 > 0) .and. (i <= j)) then
          write(lundiag1,'(1p,2i4,5e12.4)')   &
             i, j, bckernel1, 2.0e4*rad_i, 2.0e4*rad_j, den_i, den_j
          write(lundiag1,'(1p,2i4,5e12.4)')   &
             i, j, bckernel, vol_i, vol_i_inp, vol_j, vol_j_inp
          if (lundiag2 > 0) then
              write(lundiag2,  '(1p,2i4,5e12.4)')   &
              i, j, bckernel1, 2.0e4*rad_i, 2.0e4*rad_j, den_i, den_j
          end if
       end if

      return
      end subroutine brownian_kernel


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


end module mod_kernel_brown

