! Copyright (C) 2007 Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.
!
! Interface to MOSAIC aerosol and gas phase chemistry code.

module pmc_mosaic
contains
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  logical function mosaic_support()

    ! Whether MOSAIC support is compiled in.

#ifdef PMC_USE_MOSAIC
    mosaic_support = .true.
#else
    mosaic_support = .false.
#endif

  end function mosaic_support

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine mosaic_init(bin_grid, env, del_t)

    ! Initialize all MOSAIC data-structures.
    
    use pmc_constants
    use pmc_util
    use pmc_bin_grid 
    use pmc_env
    
#ifdef PMC_USE_MOSAIC
    use module_data_mosaic_aero, only: alpha_ASTEM, rtol_eqb_ASTEM, &
         ptol_mol_ASTEM, mGAS_AER_XFER, mDYNAMIC_SOLVER
    
    use module_data_mosaic_main, only: tbeg_sec, dt_sec, rlon, rlat, &
         zalt_m, RH, te, pr_atm, cair_mlc, cair_molm3, ppb, avogad, &
         deg2rad, mmode, mgas, maer, mcld, maeroptic, mshellcore, &
         msolar, mphoto, &
         lun_aeroptic
#endif
    
    type(bin_grid_t), intent(in) :: bin_grid ! bin grid
    type(env_t), intent(inout) :: env   ! environment state
    real*8, intent(in) :: del_t         ! timestep for coagulation

#ifdef PMC_USE_MOSAIC
    ! MOSAIC function interfaces
    interface
       subroutine LoadPeroxyParameters()
       end subroutine LoadPeroxyParameters
       subroutine init_data_modules()
       end subroutine init_data_modules
    end interface
    
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
    mGAS_AER_XFER = 1       ! 1 = gas-aerosol partitioning, 0 = no partition
    mDYNAMIC_SOLVER = 1     ! 1 = astem, 2 = lsodes
    alpha_ASTEM = 0.5d0     ! solver parameter. range: 0.01 - 1.0
    rtol_eqb_ASTEM = 0.01d0 ! relative eqb tolerance. range: 0.01 - 0.03
    ptol_mol_ASTEM = 0.01d0 ! percent mol tolerance.  range: 0.01 - 1.0
    
    ! time variables
    dt_sec = del_t                           ! time-step (s)
    tbeg_sec = env%start_day*24*3600 + &     ! time since the beg of
         nint(env%start_time)                ! year 00:00, UTC (s)
    
    ! geographic location
    rlon = env%longitude * deg2rad           ! longitude
    rlat = env%latitude * deg2rad            ! latitude
    zalt_m = env%altitude                    ! altitude (m)
    
    ! environmental parameters: map PartMC -> MOSAIC
    RH = env%rel_humid * 100.d0              ! relative humidity (%)
    te = env%temp                            ! temperature (K)
    pr_atm = env%pressure / const%atm        ! pressure (atm)
    
    call init_data_modules                   ! initialize indices and vars
    call LoadPeroxyParameters                ! Aperox and Bperox only once
    
    cair_mlc = avogad*pr_atm/(82.056d0*te)   ! air conc [molec/cc]
    cair_molm3 = 1d6*pr_atm/(82.056d0*te)    ! air conc [mol/m^3]
    ppb = 1d9

! matt -- this is for the separate aerosol optical output 
!         that jim wants.  i needed to get something
!         working for him initially. -- dick
    ! get unit for aerosol optical output
    if (lun_aeroptic <= 0 ) lun_aeroptic = get_unit()

#endif
    
  end subroutine mosaic_init

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine mosaic_from_partmc(bin_grid, env, aero_data, &
       aero_state, gas_data, gas_state, time)

    ! Map all data PartMC -> MOSAIC.
    
    use pmc_constants
    use pmc_util
    use pmc_aero_state
    use pmc_bin_grid 
    use pmc_env
    use pmc_aero_data
    use pmc_output_state
    use pmc_gas_data
    use pmc_gas_state
    
#ifdef PMC_USE_MOSAIC
    use module_data_mosaic_aero, only: nbin_a, aer, num_a, jhyst_leg, &
         jtotal, water_a
    
    use module_data_mosaic_main, only: tbeg_sec, tcur_sec, tmid_sec, &
         dt_sec, dt_min, dt_aeroptic_min, RH, te, pr_atm, cnn, cair_mlc, &
         ppb, msolar
#endif
    
    type(bin_grid_t), intent(in) :: bin_grid ! bin grid
    type(env_t), intent(in) :: env      ! environment state
    type(aero_data_t), intent(in) :: aero_data ! aerosol data
    type(aero_state_t), intent(in) :: aero_state ! aerosol state
    type(gas_data_t), intent(in) :: gas_data ! gas data
    type(gas_state_t), intent(in) :: gas_state ! gas state
    real*8, intent(in) :: time          ! current time (s)

#ifdef PMC_USE_MOSAIC
    ! local variables
    real*8 :: time_UTC ! 24-hr UTC clock time (hr)
    real*8 :: tmar21_sec ! time at noon, march 21, UTC (s)
    real*8 :: conv_fac(aero_data%n_spec), dum_var
    integer :: i_bin, i_part, i_spec, i_mosaic, i_spec_mosaic
    type(aero_particle_t), pointer :: particle

    ! update time variables
    tmar21_sec = dble((79*24 + 12)*3600)        ! noon, mar 21, UTC
    tcur_sec = dble(tbeg_sec) + time            ! current (old) time since
                                                ! the beg of year 00:00, UTC (s)

    time_UTC = env%start_time/3600d0         ! 24-hr UTC clock time (hr)
    time_UTC = time_UTC + dt_sec/3600d0

    do while (time_UTC >= 24d0)
       time_UTC = time_UTC - 24d0
    end do

    tmid_sec = tcur_sec + 0.5d0*dt_sec
    if(tmid_sec .ge. tmar21_sec)then
      tmid_sec = tmid_sec - tmar21_sec          ! seconds since noon, march 21
    else
      tmid_sec = tmid_sec + &
                 dble(((365-79)*24 - 12)*3600)  ! seconds since noon, march 21
    endif

    ! transport timestep (min)
    dt_min = dt_sec/60d0
    ! aerosol optics timestep (min)
    dt_aeroptic_min = 0d0

    ! compute aerosol conversion factors
    do i_spec = 1,aero_data%n_spec
       ! converts m^3(species) to nmol(species)/m^3(air)
       conv_fac(i_spec) = 1.D9 * aero_data%density(i_spec) &
            / (aero_data%molec_weight(i_spec) * aero_state%comp_vol)
    enddo

    ! environmental parameters: map PartMC -> MOSAIC
    RH = env%rel_humid * 100.d0              ! relative humidity (%)
    te = env%temp                            ! temperature (K)
    pr_atm = env%pressure / const%atm        ! pressure (atm)
    
    ! aerosol data: map PartMC -> MOSAIC
    nbin_a = total_particles(aero_state)
    i_mosaic = 0 ! MOSAIC bin number
    aer = 0d0    ! initialize to zero
    do i_bin = 1,bin_grid%n_bin
       do i_part = 1,aero_state%bins(i_bin)%n_part
          particle => aero_state%bins(i_bin)%particle(i_part)
          i_mosaic = i_mosaic + 1
          do i_spec = 1,aero_data%n_spec
             i_spec_mosaic = aero_data%mosaic_index(i_spec)
             if (i_spec_mosaic > 0) then
                ! convert m^3(species) to nmol(species)/m^3(air)
                aer(i_spec_mosaic, 3, i_mosaic) &   ! nmol/m^3(air)
                     = particle%vol(i_spec) * conv_fac(i_spec)
             end if
          end do
          ! handle water specially
          ! convert m^3(water) to kg(water)/m^3(air)
          water_a(i_mosaic) = particle%vol(aero_data%i_water) &
               * aero_data%density(aero_data%i_water) / aero_state%comp_vol
          num_a(i_mosaic) = 1d-6 / aero_state%comp_vol ! num conc (#/cc(air))
          jhyst_leg(i_mosaic) = 1
       end do
    end do

    ! gas chemistry: map PartMC -> MOSAIC
    cnn = 0d0
    do i_spec = 1,gas_data%n_spec
       i_spec_mosaic = gas_data%mosaic_index(i_spec)
       if (i_spec_mosaic > 0) then
          ! convert ppbv to molec/cc
          cnn(i_spec_mosaic) = gas_state%conc(i_spec) * cair_mlc / ppb
       end if
    end do
#endif

  end subroutine mosaic_from_partmc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine mosaic_to_partmc(bin_grid, env, aero_data, &
       aero_state, aero_binned, gas_data, gas_state)
    
    use pmc_constants
    use pmc_util
    use pmc_aero_state
    use pmc_bin_grid 
    use pmc_env
    use pmc_aero_data
    use pmc_output_state
    use pmc_gas_data
    use pmc_gas_state
    use pmc_aero_binned
    
#ifdef PMC_USE_MOSAIC
    use module_data_mosaic_aero, only: nbin_a, aer, num_a, jhyst_leg, &
         jtotal, water_a
    
    use module_data_mosaic_main, only: tbeg_sec, tcur_sec, tmid_sec, &
         dt_sec, dt_min, dt_aeroptic_min, RH, te, pr_atm, cnn, cair_mlc, &
         ppb, msolar
#endif
    
    type(bin_grid_t), intent(in) :: bin_grid ! bin grid
    type(env_t), intent(inout) :: env   ! environment state
    type(aero_data_t), intent(in) :: aero_data ! aerosol data
    type(aero_state_t), intent(inout) :: aero_state ! aerosol state
    type(aero_binned_t), intent(inout) :: aero_binned ! binned aerosol data
    type(gas_data_t), intent(in) :: gas_data ! gas data
    type(gas_state_t), intent(inout) :: gas_state ! gas state

#ifdef PMC_USE_MOSAIC
    ! local variables
    real*8 :: conv_fac(aero_data%n_spec), dum_var
    integer :: i_bin, i_part, i_spec, i_mosaic, i_spec_mosaic
    type(aero_particle_t), pointer :: particle

    ! compute aerosol conversion factors
    do i_spec = 1,aero_data%n_spec
       ! converts m^3(species) to nmol(species)/m^3(air)
       conv_fac(i_spec) = 1.D9 * aero_data%density(i_spec) &
            / (aero_data%molec_weight(i_spec) * aero_state%comp_vol)
    enddo

    ! environmental parameters: map MOSAIC -> PartMC
    env%rel_humid = RH / 100d0
    env%temp = te
    env%pressure = pr_atm * const%atm

    ! aerosol data: map MOSAIC -> PartMC
    i_mosaic = 0 ! MOSAIC bin number
    do i_bin = 1,bin_grid%n_bin
       do i_part = 1,aero_state%bins(i_bin)%n_part
          i_mosaic = i_mosaic + 1
          particle => aero_state%bins(i_bin)%particle(i_part)
          do i_spec = 1,aero_data%n_spec
             i_spec_mosaic = aero_data%mosaic_index(i_spec)
             if (i_spec_mosaic > 0) then
                particle%vol(i_spec) = &
                     ! convert nmol(species)/m^3(air) to m^3(species)
                     aer(i_spec_mosaic, 3, i_mosaic) &
                     / conv_fac(i_spec)
             end if
          end do
          ! handle water specially
          ! convert kg(water)/m^3(air) to m^3(water)
          particle%vol(aero_data%i_water) = water_a(i_mosaic) &
               / aero_data%density(aero_data%i_water) * aero_state%comp_vol
       end do
    end do
    call aero_state_resort(bin_grid, aero_state)
    call aero_state_to_binned(bin_grid, aero_data, aero_state, aero_binned)

    ! gas chemistry: map MOSAIC -> PartMC
    do i_spec = 1,gas_data%n_spec
       i_spec_mosaic = gas_data%mosaic_index(i_spec)
       if (i_spec_mosaic > 0) then
          ! convert molec/cc to ppbv
          gas_state%conc(i_spec) = cnn(i_spec_mosaic) / cair_mlc * ppb
       end if
    end do
#endif

  end subroutine mosaic_to_partmc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine mosaic_timestep(bin_grid, env, aero_data, &
       aero_state, aero_binned, gas_data, gas_state, time)

    ! Do one timestep with MOSAIC.
    
    use pmc_aero_state
    use pmc_bin_grid 
    use pmc_env
    use pmc_aero_data
    use pmc_gas_data
    use pmc_gas_state
    use pmc_aero_binned
    
#ifdef PMC_USE_MOSAIC
    use module_data_mosaic_main, only: msolar
#endif
    
    type(bin_grid_t), intent(in) :: bin_grid ! bin grid
    type(env_t), intent(inout) :: env   ! environment state
    type(aero_data_t), intent(in) :: aero_data ! aerosol data
    type(aero_state_t), intent(inout) :: aero_state ! aerosol state
    type(aero_binned_t), intent(inout) :: aero_binned ! binned aerosol data
    type(gas_data_t), intent(in) :: gas_data ! gas data
    type(gas_state_t), intent(inout) :: gas_state ! gas state
    real*8, intent(in) :: time          ! current time (s)

#ifdef PMC_USE_MOSAIC
    ! MOSAIC function interfaces
    interface
       subroutine SolarZenithAngle()
       end subroutine SolarZenithAngle
       subroutine IntegrateChemistry()
       end subroutine IntegrateChemistry
!DEBUG
       subroutine aerosol_optical()
       end subroutine aerosol_optical
!DEBUG
    end interface
    
    ! map PartMC -> MOSAIC
    call mosaic_from_partmc(bin_grid, env, aero_data, aero_state, &
         gas_data, gas_state, time)

    if (msolar == 1) then
      call SolarZenithAngle
    end if

    call IntegrateChemistry
    call aerosol_optical

    ! map MOSAIC -> PartMC
    call mosaic_to_partmc(bin_grid, env, aero_data, aero_state, &
         aero_binned, gas_data, gas_state)

    call mosaic_aero_optical(bin_grid, env, aero_data, &
         aero_state, gas_data, gas_state, time)
#endif

  end subroutine mosaic_timestep

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine mosaic_aero_optical(bin_grid, env, aero_data, &
       aero_state, gas_data, gas_state, time)

    ! Compute the optical properties of each aerosol particle.
    
    use pmc_aero_state
    use pmc_bin_grid 
    use pmc_env
    use pmc_aero_data
    use pmc_gas_data
    use pmc_gas_state
    use pmc_util
    
#ifdef PMC_USE_MOSAIC
    use module_data_mosaic_aero, only: ri_shell_a, ri_core_a, &
         ext_cross, scat_cross, asym_particle
#endif
    
    type(bin_grid_t), intent(in) :: bin_grid ! bin grid
    type(env_t), intent(in) :: env      ! environment state
    type(aero_data_t), intent(in) :: aero_data ! aerosol data
    type(aero_state_t), intent(inout) :: aero_state ! aerosol state
    type(gas_data_t), intent(in) :: gas_data ! gas data
    type(gas_state_t), intent(in) :: gas_state ! gas state
    real*8, intent(in) :: time          ! current time (s)

#ifdef PMC_USE_MOSAIC
    ! MOSAIC function interfaces
    interface
       subroutine aerosol_optical()
       end subroutine aerosol_optical
    end interface

    integer :: i_bin, i_part, i_mosaic
    type(aero_particle_t), pointer :: particle
    
    ! map PartMC -> MOSAIC
!    call mosaic_from_partmc(bin_grid, env, aero_data, aero_state, &
!         gas_data, gas_state, time)

!    call aerosol_optical

    ! map MOSAIC -> PartMC
    i_mosaic = 0 ! MOSAIC bin number
    do i_bin = 1,bin_grid%n_bin
       do i_part = 1,aero_state%bins(i_bin)%n_part
          i_mosaic = i_mosaic + 1
          particle => aero_state%bins(i_bin)%particle(i_part)
          particle%absorb_cross_sect = (ext_cross(i_mosaic) &
               - scat_cross(i_mosaic)) / 1d4                       ! (m^2)
          particle%scatter_cross_sect = scat_cross(i_mosaic) / 1d4 ! (m^2)
          particle%asymmetry = asym_particle(i_mosaic)             ! (1)
          particle%refract_shell = cmplx(ri_shell_a(i_mosaic), kind = 8) ! (1)
          particle%refract_core = cmplx(ri_core_a(i_mosaic), kind = 8)   ! (1)
          ! temporary debugging code follows
          !particle%core_vol = diam2vol(dp_core_a(i_mosaic))       ! (m^3)
          particle%core_vol = 0d0                                  ! (m^3)
       end do
    end do
#endif

  end subroutine mosaic_aero_optical

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_mosaic
