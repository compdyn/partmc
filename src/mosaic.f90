! Copyright (C) 2007 Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.
!
! Interface to MOSAIC aerosol and gas phase chemistry code.

module pmc_mosaic
contains
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine mosaic_init(bin_grid, env, del_t)
    
    use pmc_constants
    use pmc_util
    use pmc_bin_grid 
    use pmc_env
    
    use module_data_mosaic_aero, only: alpha_ASTEM, rtol_eqb_ASTEM, &
         ptol_mol_ASTEM, mGAS_AER_XFER, mDYNAMIC_SOLVER
    
    use module_data_mosaic_main, only: tbeg_sec, dt_sec, rlon, rlat, &
         zalt_m, RH, te, pr_atm, cair_mlc, cair_molm3, ppb, avogad, &
         deg2rad, mmode, mgas, maer, mcld, maeroptic, mshellcore, &
         msolar, mphoto
    
    type(bin_grid_t), intent(in) :: bin_grid ! bin grid
    type(env_t), intent(inout) :: env   ! environment state
    real*8, intent(in) :: del_t         ! timestep for coagulation

    ! MOSAIC function interfaces
    interface
       subroutine LoadPeroxyParameters()
       end subroutine LoadPeroxyParameters
       subroutine init_data_modules()
       end subroutine init_data_modules
    end interface
    
    ! local variables
    real*8 :: time_UTC ! 24-hr UTC clock time (hr)
    real*8 :: tmar21_sec ! time at noon, march 21, UTC (s)

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
    tmar21_sec = dble((79*24 + 12)*3600)     ! noon, mar 21, UTC
    tbeg_sec = env%start_day*24*3600 + &     ! time since the beg of
         nint(env%start_time)                ! year 00:00, UTC (s)
    time_UTC = env%start_time/3600d0         ! 24-hr UTC clock time (hr)
    
    ! geographic location
    rlon = env%longitude * deg2rad           ! longitude
    rlat = env%latitude * deg2rad            ! latitude
    zalt_m = env%altitude                    ! altitude (m)
    
    ! environmental parameters: map PartMC -> MOSAIC
    RH = env%rel_humid * 100.d0                     ! relative humidity (%)
    te = env%temp                               ! temperature (K)
    pr_atm = env%pressure / const%atm               ! pressure (atm)
    
    call init_data_modules                   ! initialize indices and vars
    call LoadPeroxyParameters                ! Aperox and Bperox only once
    
    cair_mlc = avogad*pr_atm/(82.056d0*te)   ! air conc [molec/cc]
    cair_molm3 = 1d6*pr_atm/(82.056d0*te)    ! air conc [mol/m^3]
    ppb = 1d9
    
  end subroutine mosaic_init

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine mosaic_timestep(bin_grid, env, aero_data, &
       aero_state, aero_binned, gas_data, gas_state, t, del_t)
    
    use pmc_constants
    use pmc_util
    use pmc_aero_state
    use pmc_bin_grid 
    use pmc_condensation
    use pmc_env
    use pmc_aero_data
    use pmc_output_state
    use pmc_gas_data
    use pmc_gas_state
    use pmc_aero_binned
    
    use module_data_mosaic_aero, only: nbin_a, aer, num_a, jhyst_leg, &
         jtotal, water_a
    
    use module_data_mosaic_main, only: tbeg_sec, tcur_sec, tmid_sec, &
         dt_sec, dt_min, dt_aeroptic_min, RH, te, pr_atm, cnn, cair_mlc, &
         ppb, msolar
    
    type(bin_grid_t), intent(in) :: bin_grid ! bin grid
    type(env_t), intent(inout) :: env   ! environment state
    type(aero_data_t), intent(in) :: aero_data ! aerosol data
    type(aero_state_t), intent(inout) :: aero_state ! aerosol state
    type(aero_binned_t), intent(inout) :: aero_binned ! binned aerosol data
    type(gas_data_t), intent(in) :: gas_data ! gas data
    type(gas_state_t), intent(inout) :: gas_state ! gas state
    real*8, intent(in) :: t             ! current time (s)
    real*8, intent(in) :: del_t         ! timestep for coagulation

    ! MOSAIC function interfaces
    interface
       subroutine SolarZenithAngle()
       end subroutine SolarZenithAngle
       subroutine IntegrateChemistry()
       end subroutine IntegrateChemistry
    end interface
    
    ! local variables
    real*8 :: time_UTC ! 24-hr UTC clock time (hr)
    real*8 :: tmar21_sec ! time at noon, march 21, UTC (s)
    real*8 :: conv_fac(aero_data%n_spec), dum_var
    integer :: i_bin, i_num, i_spec, i_mosaic, i_spec_mosaic
    type(aero_particle_t), pointer :: particle

    ! update time variables
    tmar21_sec = dble((79*24 + 12)*3600)        ! noon, mar 21, UTC
    tcur_sec = dble(tbeg_sec) + t               ! current (old) time since
                                                ! the beg of year 00:00, UTC (s)

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

    ! aerosol data: map PartMC -> MOSAIC
    nbin_a = total_particles(aero_state)
    i_mosaic = 0 ! MOSAIC bin number
    aer = 0d0    ! initialize to zero
    do i_bin = 1,bin_grid%n_bin
       do i_num = 1,aero_state%bins(i_bin)%n_part
          particle => aero_state%bins(i_bin)%particle(i_num)
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

    if (msolar == 1) then
      call SolarZenithAngle
    end if

    call IntegrateChemistry

    ! environmental parameters: map MOSAIC -> PartMC
    env%rel_humid = RH / 100d0
    env%temp = te
    env%pressure = pr_atm * const%atm

    ! aerosol data: map MOSAIC -> PartMC
    i_mosaic = 0 ! MOSAIC bin number
    do i_bin = 1,bin_grid%n_bin
       do i_num = 1,aero_state%bins(i_bin)%n_part
          i_mosaic = i_mosaic + 1
          particle => aero_state%bins(i_bin)%particle(i_num)
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
    call aero_state_to_binned(bin_grid, aero_data, aero_state, aero_binned)


     dum_var = particle%vol(2)
    

    ! gas chemistry: map MOSAIC -> PartMC
    do i_spec = 1,gas_data%n_spec
       i_spec_mosaic = gas_data%mosaic_index(i_spec)
       if (i_spec_mosaic > 0) then
          ! convert molec/cc to ppbv
          gas_state%conc(i_spec) = cnn(i_spec_mosaic) / cair_mlc * ppb
       end if
    end do

  end subroutine mosaic_timestep

end module pmc_mosaic
