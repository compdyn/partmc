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
    integer lundiag

!   c_1 = 2d0 * const%k_b * env%T / (3.d0 * const%mu)
!   a_third = v1**(1.d0/3.d0)
!   b_third = v2**(1.d0/3.d0)
    
!   k = c_1 * (a_third + b_third) *(1d0/a_third + 1d0/b_third) 
    

! calculate brownian kernel using coagsolv code
    dens1 = 1.8d3
    dens2 = 1.8d3
    lundiag = -91
    call coagsolv_kernel( v1, v2, dens1, dens2, env%t, env%p, lundiag, k )


  end subroutine kernel_brown
      

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine coagsolv(   &
        vol_i_inp, vol_j_inp, den_i_inp, den_j_inp,   &
        tk, press, lunout, coagcoef )
!
! *********************************************************************
! ***************** written by mark jacobson (1993) *******************
! ***           (c) copyright, 1993 by mark z. jacobson             ***
! ***              this version modified december, 2001             *** 
! ***                       (650) 723-6836                          ***
! *********************************************************************
!                                                                     *
!  cccccc   ooooo      a      gggggg   ssssss  ooooo   l     v      v *
! c        o     o    a a    g        s       o     o  l      v    v  *
! c        o     o   a   a   g   ggg   sssss  o     o  l       v  v   *
! c        o     o  aaaaaaa  g     g        s o     o  l        v v   *
!  cccccc   ooooo  a       a  gggggg  ssssss   ooooo   lllllll   v    *
!                                                                     *
! *********************************************************************
! the use of this module implies that the user agrees not to sell the *
! module or any portion of it, not to remove the copyright notice from*   
! it, and not to change the name of the module, "coagsolv". users may *
! modify the module as needed or incorporate it into a larger model.  *
! any special requests with regard to this module may be directed to  *
! jacobson@stanford.edu.                                              *
! *********************************************************************
!
! *********************************************************************
! * this routine calculates brownian coagulation kernels              *
! *                                                                   *
! * it was adapted (by r. easter, june-2007)                          *
! * from a box-model version of "coagsolv"                            *
! * that was obtained from mark jacobson in 2003                      *
! *                                                                   * 
! *********************************************************************
! *                         references                                *  
! *                                                                   * 
! * jacobson m. z., turco r. p., jensen, e. j. and toon o. b. (1994)  *
! *  modeling coagulation among particles of different compostion     *
! *  and size. atmos. environ. 28a, 1327 - 1338.                      *
! *                                                                   *
! * jacobson m. z. (1999) fundamentals of atmospheric modeling.       *
! *  cambridge university press, new york, 656 pp.                    *
! *                                                                   *
! *********************************************************************
!
      implicit none

!   IN arguments
      real*8, intent(in) :: vol_i_inp, vol_j_inp  ! wet (ambient) particle volumes (m3)
      real*8, intent(in) :: den_i_inp, den_j_inp  ! wet (ambient) particle densities (kg/m3)
      real*8, intent(in) :: tk                    ! air temperature (K)
      real*8, intent(in) :: press                 ! air pressure (Pa)
      integer, intent(in) :: lunout               ! logical unit for warning & diagnostic output

!   OUT arguments
      real*8, intent(out) :: coagcoef              ! coag kernel (m3/s)

!   local variables
      integer i, j, lunb

      real*8 aknud, amu, avg,   &
             boltg, bpm,   &
             cbr, consmu, cpi,   &
             deltp1, deltr,   &
             divis, dti1, dti2,   &
             fourpi,   &
             ggr, gmfp,   &
             onepi,   &
             patm, pmfp,   &
             rgas2, rho3, rsuma, rsumsq,   &
             sumdc,   &
             term1, term2, third, tworad,   &
             viscosk, vthermg,   &
             wtair

      real*8 den_i, den_j,         &  ! particle wet density (g/cm3)
             difcof_i, difcof_j,   &
             rad_i, rad_j,         &  ! particle wet radius (cm)
             sumdp_i, sumdp_j,     &
             sumvt_i, sumvt_j,     &
             vol_i, vol_j,         &  ! particle wet volume (cm3)
             vthermp_i, vthermp_j

      integer, parameter :: nbin_maxd = 1000
      integer, save :: nbin = 0
      integer :: k, m, n
      real*8, save :: rad_sv(nbin_maxd)
      real*8 :: duma


! *********************************************************************
! the "coagsolv" entry point cannot be used
! *********************************************************************
      write( *, '(/a/)' ) '*** the "coagsolv" entry cannot be used ***'
      stop
      return


! *********************************************************************
! this is the "active" entry point
! *********************************************************************
      entry coagsolv_kernel(   &
        vol_i_inp, vol_j_inp, den_i_inp, den_j_inp,   &
        tk, press, lunout, coagcoef )


! *********************************************************************
! set some parameters
! *********************************************************************
!
! boltg   = boltzmann's 1.38054e-16        (erg k-1 = g cm**2 sec-1 k-1)
! wtair   = molecular weight of air (g mol-1)
! avg     = avogadro's number (molec. mol-1)
! rgas2   = gas constant (l-atm mol-1 k-1)
! amu     = dynamic viscosity of air (g cm-1 sec-1)
!           est value at 20 c = 1.815e-04
! patm    = air pressure (atmospheres)
! rho3    = air density (g cm-3)
! viscosk = kinematic viscosity = amu / denair = (cm**2 sec-1)
! vthermg = mean thermal velocity of air molecules (cm sec-1)
! gmfp    = mean free path of an air molecule  = 2 x viscosk /
!           thermal velocity of air molecules (gmfp units of cm)
!
      third     = 1. / 3.
      onepi     = 3.14159265358979
      fourpi    = 4. * onepi
      cpi       = fourpi / 3.

      boltg     = 1.38054e-16
      wtair     = 28.966
      avg       = 6.02252e+23
      rgas2     = 0.08206
      consmu    = 1.8325e-04 * (296.16 + 120.0)

      patm = press / 1.01325d5
      amu       = (consmu / (tk + 120.)) * (tk / 296.16)**1.5
      rho3      = patm * wtair * 0.001 / (rgas2 * tk)
      viscosk   = amu / rho3
      vthermg   = sqrt(8. * boltg * tk * avg / (onepi * wtair))
      gmfp      = 2.0 * viscosk / vthermg


!
! following code attempts to construct the bin radius values
!    by saving/organizing the input radius values 
! it is purely for diagnostic purposes
!
      i = -2
      j = -1
      lunb = max( lunout+1, 92 )

      if (lunout > 0) then
         if (nbin == 0) rad_sv(:) = -1.0
         if (nbin == 0) write(*,*) '*** den_i,j =', den_i_inp, den_j_inp

         vol_i = vol_i_inp * 1.0d+6
         vol_j = vol_j_inp * 1.0d+6
         rad_i = (vol_i/cpi)**third
         rad_j = (vol_j/cpi)**third

         do k = 1, 2
            duma = rad_i
            if (k == 2) duma = rad_j
            m = -1
            do n = 1, nbin
               if (abs((duma/rad_sv(n))-1.0) <= 1.0e-5) then
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
               rad_sv(nbin) = duma
               m = nbin
            end if
            if (k == 1) i = m
            if (k == 2) j = m
         end do
      end if

!
! *********************************************************************
!              coagulation kernel from fuchs equations
! *********************************************************************
! difcof  = brownian particle diffus coef  (cm**2 sec-1)
!         = boltg * t * bpm / (6 * pi * mu * r(i))
!         = 5.25e-04 (diam = 0.01um); 6.23e-7 (diam = 0.5 um) seinfeld p.325.
! vthermp = avg thermal vel of particle    (cm sec-1)
!         = (8 * boltg * t / (pi * masmolec)**0.5
! pmfp    = mean free path of particle     (cm)
!         = 8 * di / (pi * vthermp)
! bpm     = correction for particle shape and for effects of low mean
!           free path of air. (toon et al., 1989, physical processes in
!           polar stratospheric ice clouds, jgr 94, 11,359. f1 and f2 are
!           included in the expression below (shape effects correction)
!         = 1 + kn[a + bexp(-c/kn)]
! deltp1  = if particles with mean free path pmfp leave the surface of
!           an absorbing sphere in all directions, each being probably
!           equal, then deltp1 is the mean distance from the surface of the
!           sphere reached by particles after covering a distance pmfp. (cm).
! cbr     = coag kernel (cm3 partic-1 s-1) due to brownian motion (fuchs, 1964)
! rrate   = coag kernal (cm3 partic-1 s-1) * time step (s)
!

       if (lunout > 0) then
          write(lunout,'(/a,1p,2e12.4)') 'tk, patm', tk, patm
          write(lunout,'(a)')   &
          'i, j, coagcoef (cm3/s), dpwet_i,j (um), denswet_i,j (g/cm3)'
       end if

         den_i      = den_i_inp * 1.0d-3
         vol_i      = vol_i_inp * 1.0d+6
         rad_i      = (vol_i/cpi)**third

         tworad     = rad_i + rad_i  
         aknud      = gmfp/rad_i  
         bpm        = 1. + aknud*(1.257 + 0.42*exp(-1.1/aknud))
         difcof_i   = boltg*tk*bpm/(6.*onepi*rad_i*amu)
         vthermp_i  = sqrt( 8.*boltg*tk/(onepi*den_i*vol_i) )
         sumvt_i    = vthermp_i*vthermp_i 
         pmfp       = 8.*difcof_i/(onepi*vthermp_i)
         dti1       = tworad + pmfp
         dti2       = (4.*rad_i*rad_i + pmfp*pmfp)**1.5
         divis      = 0.166667/(rad_i*pmfp)
         deltp1     = divis*(dti1*dti1*dti1 - dti2) - tworad
         sumdp_i    = deltp1   * deltp1

           den_j      = den_j_inp * 1.0d-3
           vol_j      = vol_j_inp * 1.0d+6
           rad_j      = (vol_j/cpi)**third

           tworad     = rad_j + rad_j  
           aknud      = gmfp/rad_j  
           bpm        = 1. + aknud*(1.257 + 0.42*exp(-1.1/aknud))
           difcof_j   = boltg*tk*bpm/(6.*onepi*rad_j*amu)
           vthermp_j  = sqrt( 8.*boltg*tk/(onepi*den_j*vol_j) )
           sumvt_j    = vthermp_j*vthermp_j 
           pmfp       = 8.*difcof_j/(onepi*vthermp_j)
           dti1       = tworad + pmfp
           dti2       = (4.*rad_j*rad_j + pmfp*pmfp)**1.5
           divis      = 0.166667/(rad_j*pmfp)
           deltp1     = divis*(dti1*dti1*dti1 - dti2) - tworad
           sumdp_j    = deltp1   * deltp1

           rsuma      = rad_i  + rad_j
           rsumsq     = rsuma*rsuma
           sumdc      = difcof_i + difcof_j 
           ggr        = sqrt(sumvt_i + sumvt_j)
           deltr      = sqrt(sumdp_i + sumdp_j)
           term1      = rsuma/(rsuma + deltr)
           term2      = 4.*sumdc/(rsuma*ggr)
           cbr        = fourpi*rsuma*sumdc/(term1 + term2)

           coagcoef   = cbr * 1.0d-6

           if ((lunout > 0) .and. (i <= j)) then
               write(lunout,'(1p,2i4,5e12.4)')   &
                  i, j, cbr, 2.0e4*rad_i, 2.0e4*rad_j, den_i, den_j
               write(lunout,'(1p,2i4,5e12.4)')   &
                  i, j, coagcoef, vol_i, vol_i_inp, vol_j, vol_j_inp
               write(lunb,  '(1p,2i4,5e12.4)')   &
                  i, j, cbr, 2.0e4*rad_i, 2.0e4*rad_j, den_i, den_j
           end if

      return
      end subroutine coagsolv


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


end module mod_kernel_brown

