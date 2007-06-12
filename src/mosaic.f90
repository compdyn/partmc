! Copyright (C) 2007 Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.
!
! Interface to MOSAIC

module mod_mosaic
contains
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine singlestep_mosaic(M, n_spec, n_bin, MH, VH, &
       bin_v, bin_g, bin_gs, bin_n, dlnr, t, del_t, env, mat, gas)

    use mod_constants
    use mod_util
    use mod_array
    use mod_bin 
    use mod_condensation
    use mod_environ
    use mod_material
    use mod_state
    use mod_gas

    use module_data_mosaic_aero, only: nbin_a, aer, mGAS_AER_XFER, &
         mDYNAMIC_SOLVER, alpha_ASTEM, rtol_eqb_ASTEM, ptol_mol_ASTEM, &
         jtotal, water_a

    use module_data_mosaic_main, only: tbeg_sec, tcur_sec, tmid_sec, &
         dt_sec, dt_min, dt_aeroptic_min, rlon, rlat, zalt_m, RH, te, &
         pr_atm, cnn, cair_mlc, cair_molm3, h2o, o2, h2, ppb, avogad, deg2rad, &
         mmode, mgas, maer, mcld, maeroptic, mshellcore, msolar, mphoto

    implicit none    

    integer, intent(inout) :: M         ! actual number of particles
    integer, intent(in) :: n_spec       ! number of species
    integer, intent(in) :: n_bin        ! number of bins
    integer, intent(inout) :: MH(n_bin) ! number of particles per bin
    type(bin_p), intent(inout) :: VH(n_bin) ! particle volumes (m^3)
    
    real*8, intent(in) :: bin_v(n_bin)  ! volume of particles in bins (m^3)
    real*8, intent(out) :: bin_g(n_bin) ! volume in bins  
    real*8, intent(out) :: bin_gs(n_bin,n_spec) ! species volume in bins
    integer, intent(out) :: bin_n(n_bin) ! number in bins
    real*8, intent(in) :: dlnr          ! bin scale factor
    
    real*8, intent(in) :: t             ! current time (s)
    real*8, intent(in) :: del_t         ! timestep for coagulation
    
    type(environ), intent(inout) :: env ! environment state
    type(material), intent(in) :: mat   ! material properties
    type(gas_chem), intent(inout) :: gas ! gas chemistry

    ! MOSAIC function interfaces
    interface
       subroutine SolarZenithAngle()
       end subroutine SolarZenithAngle
       subroutine LoadPeroxyParameters()
       end subroutine LoadPeroxyParameters
       subroutine init_data_modules()
       end subroutine init_data_modules
       subroutine IntegrateChemistry()
       end subroutine IntegrateChemistry
       function WaterVapor(RH, cair_mlc, te, pr_atm)
       end function WaterVapor
!       subroutine DoMassBalance()
!       end subroutine DoMassBalance
    end interface


    ! local variables
    real*8 :: time_UTC ! 24-hr UTC clock time (hr)
    real*8 :: tmar21_sec ! time at noon, march 21, UTC (s)
    real*8 :: conv_fac(20), dum_var
    integer :: i_bin, i_num, i_spec, i_mosaic, i_spec_mosaic

    logical first
    save first
    data first/.true./



!---------------------------------------------------------
! set MOSAIC data once at the beginning of the simulation
    if(first)then
      first=.false.

! parameters
    mmode = 1               ! 1 = time integration, 2 = parametric analysis
    mgas = 1                ! 1 = gas chem on, 0 = gas chem off
    maer = 1                ! 1 = aer chem on, 0 = aer chem off
    mcld = 0                ! 1 = cld chem on, 0 = cld chem off
    maeroptic = 1           ! 1 = aer_optical on, 0 = aer_optical off
    mshellcore = 1          ! 0 = no shellcore, 1 = core is BC only
                            ! 2 = core is BC and DUST
    msolar = 1              ! 1 = diurnally varying phot, 2 = fixed phot
    mphoto = 2              ! 1 = Rick's param, 2 = Yang's param
    mGAS_AER_XFER = 1       ! 1 = gas-aerosol partitioning, 0 = do not partition
    mDYNAMIC_SOLVER = 1     ! 1 = astem, 2 = lsodes
    alpha_ASTEM = 0.5d0     ! solver parameter. range: 0.01 - 1.0
    rtol_eqb_ASTEM = 0.01d0 ! relative eqb tolerance. range: 0.01 - 0.03
    ptol_mol_ASTEM = 0.01d0 ! percent mol tolerance.  range: 0.01 - 1.0

! time variables
      dt_sec = del_t				! time-step (s)
      tmar21_sec = (79*24 + 12)*3600		! noon, mar 21, UTC
      tbeg_sec = env%start_day*24*3600 + &	! time since the beg of year 00:00, UTC (s)
                 env%start_time 
      time_UTC = env%start_time/3600		! 24-hr UTC clock time (hr)

! geographic location
      rlon = env%longitude * deg2rad		! longitude
      rlat = env%latitude * deg2rad		! latitude
      zalt_m = env%altitude			! altitude (m)

! environmental parameters: map PartMC -> MOSAIC
      RH = env%RH
      te = env%T
      pr_atm = env%p / const%atm

      call init_data_modules			! initialize many indices and variables
      call LoadPeroxyParameters			! Aperox and Bperox only once

      cair_mlc = avogad*pr_atm/(82.056*te)	! air conc [molec/cc]
      cair_molm3 = 1.e6*pr_atm/(82.056*te)	! air conc [mol/m^3]
      ppb = 1.e+9

    endif
!----------------------------------------------------------

! update time variables
    tcur_sec = tbeg_sec + t			! current (old) time since the beg of year 00:00, UTC (s)
    time_UTC = time_UTC + dt_sec/3600.

    if(time_UTC .ge. 24.0)then
      time_UTC = time_UTC - 24.0
    endif

    tmid_sec = tcur_sec + 0.5*dt_sec
    if(tmid_sec .ge. tmar21_sec)then
      tmid_sec = tmid_sec - tmar21_sec		! seconds since noon, march 21
    else
      tmid_sec = tmid_sec + &
                 ((365-79)*24 - 12)*3600	! seconds since noon, march 21
    endif


    ! transport timestep (min)
    dt_min = dt_sec/60.
    ! aerosol optics timestep (min)
    dt_aeroptic_min = 0d0


!----------------------------------------------------------
! aerosol data: map PartMC -> MOSAIC

    do i_spec = 1, n_spec
      conv_fac(i_spec) = 1.D9*mat%rho(i_spec)/(mat%M_w(i_spec)*env%V_comp)  ! converts m^3(species) to nmol(species)/m^3(air)
    enddo


    nbin_a = M
    i_mosaic = 0 ! MOSAIC bin number
    aer = 0d0
    do i_bin = 1,n_bin
       do i_num = 1,MH(i_bin)
          i_mosaic = i_mosaic + 1
          do i_spec = 1,n_spec
             i_spec_mosaic = mat%mosaic_index(i_spec)
             if (i_spec_mosaic > 0) then
                aer(i_spec_mosaic, jtotal, i_mosaic) &		! convert m^3(species) to nmol(species)/m^3(air)
                     = VH(i_bin)%p(i_num, i_spec)*conv_fac(i_spec)
             end if
          end do
          ! handle water specially
          water_a(i_mosaic) = VH(i_bin)%p(i_num,mat%i_water)*mat%rho(mat%i_water)/env%V_comp  ! convert m^3(water) to kg(water)/m^3(air)
       end do
    end do

! gas chemistry: map PartMC -> MOSAIC
    cnn = 0d0
    do i_spec = 1,gas%n_spec
       i_spec_mosaic = gas%mosaic_index(i_spec)
       if (i_spec_mosaic > 0) then
          cnn(i_spec_mosaic) = gas%conc(i_spec)*cair_mlc/ppb	! convert ppbv to molec/cc
       end if
    end do




    if(msolar.eq.1)then
      call SolarZenithAngle
    end if

    call IntegrateChemistry

!    call DoMassBalance




! environmental parameters: map MOSAIC -> PartMC
    env%RH = RH
    env%T = te
    env%p = pr_atm * const%atm

! aerosol data: map MOSAIC -> PartMC
    nbin_a = M
    i_mosaic = 0 ! MOSAIC bin number
    aer = 0d0
    do i_bin = 1,n_bin
       do i_num = 1,MH(i_bin)
          i_mosaic = i_mosaic + 1
          do i_spec = 1,n_spec
             i_spec_mosaic = mat%mosaic_index(i_spec)
             if (i_spec_mosaic > 0) then
                VH(i_bin)%p(i_num, i_spec) = &
                     aer(i_spec_mosaic, jtotal, i_mosaic)/conv_fac(i_spec)	! convert nmol(species)/m^3(air) to m^3(species)
             end if
          end do
          ! handle water specially
          VH(i_bin)%p(i_num,mat%i_water) = water_a(i_mosaic)/mat%rho(mat%i_water)*env%V_comp  ! convert kg(water)/m^3(air) to m^3(water)
       end do
    end do

! gas chemistry: map MOSAIC -> PartMC
    do i_spec = 1,gas%n_spec
       i_spec_mosaic = gas%mosaic_index(i_spec)
       if (i_spec_mosaic > 0) then
          gas%conc(i_spec) = cnn(i_spec_mosaic)/cair_mlc*ppb  ! convert molec/cc to ppbv
       end if
    end do

  end subroutine singlestep_mosaic

end module mod_mosaic
