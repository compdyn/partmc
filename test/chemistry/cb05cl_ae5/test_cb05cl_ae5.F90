! Copyright (C) 2017 Matt Dawson
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_test_cb05cl_ae5 program

!> Test for the cb05cl_ae5 mechanism from MONARCH. This program runs the
!! MONARCH CB5 code and the Phlex-chem version and compares the output.
program pmc_test_cb05cl_ae5

#define DEBUG
     
  use pmc_constants,                    only: const
  use pmc_util,                         only: i_kind, dp, assert, assert_msg, &
                                              almost_equal, string_t, &
                                              to_string, warn_assert_msg
  use pmc_phlex_core
  use pmc_phlex_state
  use pmc_phlex_solver_data
  use pmc_chem_spec_data
  use pmc_rxn_data
  use pmc_rxn_photolysis
  use pmc_property
#ifdef PMC_USE_JSON
  use json_module
#endif

  ! EBI Solver
  use module_bsc_chem_data

  implicit none

  ! New-line character
  character(len=*), parameter :: new_line = char(10)
  ! EBI solver output file unit
  integer(kind=i_kind), parameter :: EBI_FILE_UNIT = 10
  ! KPP solver output file unit
  integer(kind=i_kind), parameter :: KPP_FILE_UNIT = 11
  ! Phlex-chem output file unit
  integer(kind=i_kind), parameter :: PHLEX_FILE_UNIT = 12
  ! Number of timesteps to integrate over
  integer(kind=i_kind), parameter :: NUM_TIME_STEPS = 100
  ! Number of EBI-solver species
  integer(kind=i_kind), parameter :: NUM_EBI_SPEC = 72
  ! Number of EBI-solever photolysis reactions
  integer(kind=i_kind), parameter :: NUM_EBI_PHOTO_RXN = 23
  ! Small number for minimum concentrations
  real(kind=dp), parameter :: SMALL_NUM = 1.0d-30
  ! Used to check availability of a solver  
  type(phlex_solver_data_t), pointer :: phlex_solver_data

#ifdef DEBUG
  integer(kind=i_kind), parameter :: DEBUG_UNIT = 13
   
  open(unit=DEBUG_UNIT, file="out/debug_cb05cl_ae.txt", status="replace", action="write")
#endif

  phlex_solver_data => phlex_solver_data_t()

  if (.not.phlex_solver_data%is_solver_available()) then
    write(*,*) "CB5 mechanism test - no solver available - PASS"
  else if (run_cb05cl_ae5_tests()) then
    write(*,*) "CB5 mechanism tests - PASS"
  else
    write(*,*) "CB5 mechanism tests - FAIL"
  end if

#ifdef DEBUG
  close(DEBUG_UNIT)
#endif

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Run all CB5 tests
  logical function run_cb05cl_ae5_tests() result(passed)

    passed = run_standard_cb05cl_ae5_test()

  end function run_cb05cl_ae5_tests

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Run the cb05cl_ae5 mechanism under standard conditions using the original 
  !! MONARCH ebi-solver code, the KPP CB5 module and the Phlex-chem version
  logical function run_standard_cb05cl_ae5_test() result(passed)

    ! EBI Solver
    use EXT_HRDATA
    use EXT_RXCM,                               only : NRXNS, RXLABEL

    ! KPP Solver
    use cb05cl_ae5_Initialize,                  only : KPP_Initialize => Initialize
    use cb05cl_ae5_Model,                       only : KPP_NSPEC => NSPEC, &
                                                       KPP_STEPMIN => STEPMIN, &
                                                       KPP_STEPMAX => STEPMAX, &
                                                       KPP_RTOL => RTOL, &
                                                       KPP_ATOL => ATOL, &
                                                       KPP_TIME => TIME, &
                                                       KPP_C => C, &
                                                       KPP_RCONST => RCONST, &
                                                       KPP_Update_RCONST => Update_RCONST, &
                                                       KPP_INTEGRATE => INTEGRATE, &
                                                       KPP_SPC_NAMES => SPC_NAMES, &
                                                       KPP_PHOTO_RATES => PHOTO_RATES, &
                                                       KPP_TEMP => TEMP, &
                                                       KPP_PRESS => PRESS, &
                                                       KPP_SUN => SUN, &
                                                       KPP_M => M, &
                                                       KPP_N2 => N2, &
                                                       KPP_O2 => O2, &
                                                       KPP_H2 => H2, &
                                                       KPP_H2O => H2O, &
                                                       KPP_N2O => N2O, &
                                                       KPP_CH4 => CH4, &
                                                       KPP_NVAR => NVAR, &
                                                       KPP_NREACT => NREACT, &
                                                       KPP_DT => DT
    use cb05cl_ae5_Parameters,                  only : KPP_IND_O2 => IND_O2
    use cb05cl_ae5_Initialize, ONLY: Initialize

    ! EBI-solver species names
    type(string_t), dimension(NUM_EBI_SPEC) :: ebi_spec_names

    ! KPP reaction labels
    type(string_t), allocatable :: kpp_rxn_labels(:)
    ! KPP rstate
    real(kind=dp) :: KPP_RSTATE(20)
    ! #/cc -> ppm conversion factor
    real(kind=dp) :: conv

    ! Flag for sunlight
    logical :: is_sunny
    ! Photolysis rates (\min)
    real, allocatable :: photo_rates(:)
    ! Temperature (K)
    real :: temperature = 272.5
    ! Pressure (atm)
    real :: pressure = 0.8
    ! Water vapor concentration (ppmV)
    real :: water_conc = 0.0 ! (Set by Phlex-chem initial concentration)
    
    ! Phlex-chem core
    type(phlex_core_t), pointer :: phlex_core
    ! Phlex-chem state
    type(phlex_state_t), pointer :: phlex_state, phlex_state_comp
    ! Phlex-chem species names
    type(string_t), allocatable :: phlex_spec_names(:)
    ! EBI -> Phlex-chem species map
    integer(kind=i_kind), dimension(NUM_EBI_SPEC) :: spec_map

    ! Computation timer variables
    real(kind=dp) :: comp_start, comp_end, comp_ebi, comp_kpp, comp_phlex

    class(rxn_data_t), pointer :: rxn
    type(property_t), pointer :: prop_set
    character(len=:), allocatable :: key, spec_name, string_val, phlex_input_file
    real(kind=dp) :: real_val, phlex_rate, phlex_rate_const
    integer(kind=i_kind) :: i_spec, j_spec, i_mech, i_rxn, i_ebi_rxn, i_kpp_rxn, &
            i_time, i_repeat, n_gas_spec

    integer(kind=i_kind) :: i_M, i_O2, i_N2, i_H2O, i_CH4, i_H2
    integer(kind=i_kind), allocatable :: ebi_rxn_map(:), kpp_rxn_map(:)
    integer(kind=i_kind), allocatable :: ebi_spec_map(:), kpp_spec_map(:)
    type(string_t) :: str_temp

    ! Arrays to hold starting concentrations
    real(kind=dp), allocatable :: ebi_init(:), kpp_init(:), phlex_init(:)

    passed = .false.
   
    ! Set the #/cc -> ppm conversion factor
    conv = 1.0d0/ (const%avagadro /const%univ_gas_const * 10.0d0**(-12.0d0) * &
            (pressure*101325.d0) /temperature)
#ifdef DEBUG
    write(DEBUG_UNIT,*) "#/cc -> ppm conversion factor: ", conv
#endif

    ! Load the EBI solver species names
    call set_ebi_species(ebi_spec_names)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! Initialize the EBI solver !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call cpu_time(comp_start)
    ! Set the BSC chem parameters
    call init_bsc_chem_data()
    ! Set the output unit
    LOGDEV = 6
    ! Set the aerosol flag
    L_AE_VRSN = .false.
    ! Set the aq. chem flag
    L_AQ_VRSN = .false.
    ! Initialize the solver
    call EXT_HRINIT
    RKI(:) = 0.0
    RXRAT(:) = 0.0
    YC(:) = 0.0
    YC0(:) = 0.0
    YCP(:) = 0.0
    PROD(:) = 0.0
    LOSS(:) = 0.0
    PNEG(:) = 0.0
    ! Set the timestep (min)
    EBI_TMSTEP = 0.1
    ! Set the number of timesteps
    N_EBI_STEPS = 1
    ! Set the number of internal timesteps
    N_INR_STEPS = 1
    call cpu_time(comp_end)
    write(*,*) "EBI initialization time: ", comp_end-comp_start," s"

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! Initialize the KPP CB5 module !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    call cpu_time(comp_start)
    ! Set the step limits
    KPP_STEPMIN = 0.0d0
    KPP_STEPMAX = 0.0d0
    KPP_SUN = 1.0
    ! Set the tolerances
    do i_spec = 1, KPP_NVAR
      KPP_RTOL(i_spec) = 1.0d-4
      KPP_ATOL(i_spec) = 1.0d-3
    end do
    CALL KPP_Initialize()
    call cpu_time(comp_end)
    write(*,*) "KPP initialization time: ", comp_end-comp_start," s"


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! Initialize phlex-chem !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    call cpu_time(comp_start)
    phlex_input_file = "config_cb05cl_ae5.json"
    phlex_core => phlex_core_t(phlex_input_file)
    
    ! Initialize the model
    call phlex_core%initialize()
   
    ! Find the CB5 mechanism
    key = "cb05cl_ae5"
    call assert(418262750, phlex_core%find_mechanism(key, i_mech))

    ! Set the photolysis rate ids
    key = "rxn id"
    do i_rxn = 1, phlex_core%mechanism(i_mech)%size()
      rxn => phlex_core%mechanism(i_mech)%get_rxn(i_rxn)
      select type(rxn) 
        type is (rxn_photolysis_t)
          call assert(265614917, rxn%property_set%get_string(key, string_val))
          if (trim(string_val).eq."jo2") then
            ! Set O2 + hv rate constant to 0 (not present in ebi version)
            call rxn%set_photo_id(0)
          else
            call rxn%set_photo_id(1)
          end if
      end select
    end do

    ! Initialize the solver
    call phlex_core%solver_initialize()

    ! Get an new state variable
    phlex_state => phlex_core%new_state()
    
    ! Set the environmental conditions
    phlex_state%env_state%temp = temperature
    phlex_state%env_state%pressure = pressure * const%air_std_press
    call phlex_state%update_env_state()

    call cpu_time(comp_end)
    write(*,*) "Phlex-chem initialization time: ", comp_end-comp_start," s"

    ! Get a phlex-state for rate comparisons
    phlex_state_comp => phlex_core%new_state()
    phlex_state_comp%env_state%temp = phlex_state%env_state%temp
    phlex_state_comp%env_state%pressure = phlex_state%env_state%pressure
    call phlex_state_comp%update_env_state()
  
    ! Find the constant species in the CB5 mechanism
    spec_name = "M"
    i_M = phlex_core%chem_spec_data%gas_state_id(spec_name)
    spec_name = "O2"
    i_O2 = phlex_core%chem_spec_data%gas_state_id(spec_name)
    spec_name = "N2"
    i_N2 = phlex_core%chem_spec_data%gas_state_id(spec_name)
    spec_name = "H2O"
    i_H2O = phlex_core%chem_spec_data%gas_state_id(spec_name)
    spec_name = "CH4"
    i_CH4 = phlex_core%chem_spec_data%gas_state_id(spec_name)
    spec_name = "H2"
    i_H2 = phlex_core%chem_spec_data%gas_state_id(spec_name)

    ! Set the photolysis rates (dummy values for solver comparison)
    is_sunny = .true.
    allocate(photo_rates(NUM_EBI_PHOTO_RXN))
    photo_rates(:) = 0.0001 * 60.0 ! EBI solver wants rates in min^-1
    KPP_PHOTO_RATES(:) = 0.0001
    ! Set O2 + hv rate constant to 0 in KPP (not present in ebi version)
    KPP_PHOTO_RATES(1) = 0.0
    ! Set the phlex-chem photolysis rate constants
    call phlex_core%set_photo_rate(1, real(0.0001, kind=dp))

    ! Make sure the right number of reactions is present
    ! (KPP includes two Cl rxns with rate constants set to zero that are not
    !  present in phlex-chem)
    call assert_msg(396732632, phlex_core%mechanism(i_mech)%size().eq.186, &
            "Wrong number of phlex-chem reactions: "// &
            trim(to_string(phlex_core%mechanism(i_mech)%size())))

    ! Set the initial concentrations
    key = "init conc"
    YC(:) = 0.0
    KPP_C(:) = 0.0
    phlex_state%state_var(:) = 0.0

    ! Set the initial concentrations in each module
    do i_spec = 1, NUM_EBI_SPEC

      ! Get initial concentrations from phlex-chem input data
      call assert(787326679, phlex_core%chem_spec_data%get_property_set( &
              ebi_spec_names(i_spec)%string, prop_set))
      if (prop_set%get_real(key, real_val)) then
        
        ! Set the EBI solver concetration (ppm)
        YC(i_spec) = real_val

        ! Set the phlex-chem concetration (ppm)
        phlex_state%state_var( &
                phlex_core%chem_spec_data%gas_state_id( &
                ebi_spec_names(i_spec)%string)) = real_val
#ifdef DEBUG
        write(DEBUG_UNIT,*) "Species ", ebi_spec_names(i_spec)%string, &
                ", Phlex-chem id: ", phlex_core%chem_spec_data%gas_state_id( &
                ebi_spec_names(i_spec)%string), ", init conc: ", real_val
      else
        write(DEBUG_UNIT,*) "No initial concentration for species ", &
                ebi_spec_names(i_spec)%string
#endif
      end if
      
      ! Set KPP species concentrations (#/cc)
      do j_spec = 1, KPP_NSPEC
        if (trim(ebi_spec_names(i_spec)%string).eq.trim(KPP_SPC_NAMES(j_spec))) then
          KPP_C(j_spec) = YC(i_spec) / conv
        end if
      end do
    end do
    
    ! Set EBI solver constant species concentrations in Phlex-chem
    spec_name = "M"
    call assert(273497194, phlex_core%chem_spec_data%get_property_set(spec_name, prop_set))
    call assert(740666066, associated(prop_set))
    call assert(907464197, prop_set%get_real(key, real_val))
    phlex_state%state_var(i_M) = real_val
    KPP_M = real_val / conv
    spec_name = "O2"
    call assert(557877977, phlex_core%chem_spec_data%get_property_set(spec_name, prop_set))
    call assert(729136508, associated(prop_set))
    call assert(223930103, prop_set%get_real(key, real_val))
    phlex_state%state_var(i_O2) = real_val
    KPP_O2 = real_val / conv
    KPP_C(KPP_IND_O2) = real_val / conv
    ! KPP has variable O2 concentration  
    do j_spec = 1, KPP_NSPEC
      if (trim(KPP_SPC_NAMES(j_spec)).eq.'O2') then
        KPP_C(j_spec) = real_val / conv
      end if
    end do
    spec_name = "N2"
    call assert(329882514, phlex_core%chem_spec_data%get_property_set(spec_name, prop_set))
    call assert(553715297, associated(prop_set))
    call assert(666033642, prop_set%get_real(key, real_val))
    phlex_state%state_var(i_N2) = real_val
    KPP_N2 = real_val / conv
    spec_name = "H2O"
    call assert(101887051, phlex_core%chem_spec_data%get_property_set(spec_name, prop_set))
    call assert(160827237, associated(prop_set))
    call assert(273145582, prop_set%get_real(key, real_val))
    phlex_state%state_var(i_H2O) = real_val
    KPP_H2O = real_val / conv
    spec_name = "CH4"
    call assert(208941089, phlex_core%chem_spec_data%get_property_set(spec_name, prop_set))
    call assert(667939176, associated(prop_set))
    call assert(780257521, prop_set%get_real(key, real_val))
    phlex_state%state_var(i_CH4) = real_val
    KPP_CH4 = real_val / conv
    ! KPP has variable CH4 concentration  
    do j_spec = 1, KPP_NSPEC
      if (trim(KPP_SPC_NAMES(j_spec)).eq.'CH4') then
        KPP_C(j_spec) = real_val / conv
      end if
    end do
    spec_name = "H2"
    call assert(663478776, phlex_core%chem_spec_data%get_property_set(spec_name, prop_set))
    call assert(892575866, associated(prop_set))
    call assert(722418962, prop_set%get_real(key, real_val))
    phlex_state%state_var(i_H2) = real_val
    KPP_H2 = real_val / conv

    ! Set the water concentration for EBI solver (ppmV)
    water_conc = phlex_state%state_var(i_H2O)

#ifdef DEBUG
    call phlex_core%print(DEBUG_UNIT)
    write(DEBUG_UNIT,*) "*** Constant species concentrations ***"
    write(DEBUG_UNIT,*) "[M] = ", phlex_state%state_var(i_M), i_M
    write(DEBUG_UNIT,*) "[O2] = ", phlex_state%state_var(i_O2), i_O2
    write(DEBUG_UNIT,*) "[N2] = ", phlex_state%state_var(i_N2), i_N2
    write(DEBUG_UNIT,*) "[H2O] = ", phlex_state%state_var(i_H2O), i_H2O
    write(DEBUG_UNIT,*) "[CH4] = ", phlex_state%state_var(i_CH4), i_CH4
    write(DEBUG_UNIT,*) "[H2] = ", phlex_state%state_var(i_H2), i_H2
#endif

    ! Set up the output files
    open(EBI_FILE_UNIT, file="out/cb05cl_ae5_ebi_results.txt", status="replace", &
            action="write")
    open(KPP_FILE_UNIT, file="out/cb05cl_ae5_kpp_results.txt", status="replace", &
            action="write")
    open(PHLEX_FILE_UNIT, file="out/cb05cl_ae5_phlex_results.txt", status="replace", &
            action="write")
    n_gas_spec = phlex_core%chem_spec_data%size(spec_phase=CHEM_SPEC_GAS_PHASE)
    allocate(phlex_spec_names(n_gas_spec))
    do i_spec = 1, n_gas_spec
      phlex_spec_names(i_spec)%string = phlex_core%chem_spec_data%gas_state_name(i_spec)
    end do
    write(EBI_FILE_UNIT,*) "time ", (ebi_spec_names(i_spec)%string//" ", i_spec=1, &
            NUM_EBI_SPEC)
    write(KPP_FILE_UNIT,*) "time ", (trim(KPP_SPC_NAMES(i_spec))//" ", i_spec=1, &
            KPP_NSPEC)
    write(PHLEX_FILE_UNIT,*) "time ", (phlex_spec_names(i_spec)%string//" ", i_spec=1, &
            size(phlex_spec_names))

    ! Set up the reaction map between phlex-chem, kpp and ebi solvers
    key = "rxn id"
    allocate(ebi_rxn_map(phlex_core%mechanism(i_mech)%size()))
    ebi_rxn_map(:) = 0
    allocate(kpp_rxn_map(phlex_core%mechanism(i_mech)%size()))
    kpp_rxn_map(:) = 0
    call get_kpp_rxn_labels(kpp_rxn_labels)
    do i_rxn = 1, phlex_core%mechanism(i_mech)%size()
      rxn => phlex_core%mechanism(i_mech)%get_rxn(i_rxn)
      call assert_msg(917216189, associated(rxn), "Missing rxn "//to_string(i_rxn))
      call assert(656034097, rxn%property_set%get_string(key, string_val))
      do i_ebi_rxn = 1, NRXNS 
        if (trim(RXLABEL(i_ebi_rxn)).eq.trim(string_val)) then
          ebi_rxn_map(i_rxn) = i_ebi_rxn
          exit
        end if
      end do
      if (trim(string_val).ne.'jo2') then ! jo2 rxn O2 + hv is not in EBI solver
        call assert_msg(921715481, ebi_rxn_map(i_rxn).ne.0, "EBI missing rxn "//string_val)
      end if
      do i_kpp_rxn = 1, KPP_NREACT
        if (trim(kpp_rxn_labels(i_kpp_rxn)%string).eq.trim(string_val)) then
          kpp_rxn_map(i_rxn) = i_kpp_rxn
          exit
        end if
      end do
      call assert_msg(360769001, kpp_rxn_map(i_rxn).ne.0, "KPP missing rxn "//string_val)
    end do
  
    ! Set up the species map between phlex-chem, kpp and ebi solvers
    allocate(ebi_spec_map(phlex_core%chem_spec_data%size()))
    allocate(kpp_spec_map(phlex_core%chem_spec_data%size()))
    ebi_spec_map(:) = 0
    kpp_spec_map(:) = 0
    do i_spec = 1, NUM_EBI_SPEC
      j_spec = phlex_core%chem_spec_data%gas_state_id(ebi_spec_names(i_spec)%string)
      call assert_msg(194404050, j_spec.gt.0, "Missing EBI species: "//trim(ebi_spec_names(i_spec)%string))
      ebi_spec_map(j_spec) = i_spec
    end do
    do i_spec = 1, KPP_NSPEC
      spec_name = trim(KPP_SPC_NAMES(i_spec))
      j_spec = phlex_core%chem_spec_data%gas_state_id(spec_name)
      call assert_msg(194404050, j_spec.gt.0, "Missing KPP species: "//trim(KPP_SPC_NAMES(i_spec)))
      kpp_spec_map(j_spec) = i_spec
    end do

    ! Reset the computation timers
    comp_ebi = 0.0
    comp_kpp = 0.0
    comp_phlex = 0.0

    ! Compare the rates for the initial conditions
    call EXT_HRCALCKS( NUM_EBI_PHOTO_RXN,       & ! Number of EBI solver photolysis reactions
                       is_sunny,                & ! Flag for sunlight
                       photo_rates,             & ! Photolysis rates
                       temperature,             & ! Temperature (K)
                       pressure,                & ! Air pressure (atm)
                       water_conc,              & ! Water vapor concentration (ppmV)
                       RKI)                       ! Rate constants
    KPP_DT = EBI_TMSTEP * 60.0d0 ! KPP Time step in seconds
    KPP_TIME = 0.0
    KPP_TEMP = temperature
    KPP_PRESS = pressure * 1013.25 ! KPP pressure in hPa
    CALL KPP_Update_RCONST()
    phlex_state_comp%state_var(:) = phlex_state%state_var(:)
    call compare_rates(phlex_core, phlex_state_comp, ebi_spec_names, conv, &
            ebi_rxn_map, kpp_rxn_map)

    ! Save the initial states for repeat calls
    allocate(ebi_init(size(YC)))
    allocate(kpp_init(size(KPP_C)))
    allocate(phlex_init(size(phlex_state%state_var)))
    ebi_init(:) = YC(:)
    kpp_init(:) = KPP_C(:)
    phlex_init(:) = phlex_state%state_var(:)

    ! Repeatedly solve the mechanism
    do i_repeat = 1, 100

    YC(:) = ebi_init(:)
    KPP_C(:) = kpp_init(:)
    phlex_state%state_var(:) = phlex_init(:)

    ! Solve the mechanism
    do i_time = 1, NUM_TIME_STEPS

#ifdef DEBUG
      if (i_repeat.eq.1) then
      write(DEBUG_UNIT,*) "Time step: ", i_time
      write(DEBUG_UNIT,*) "EBI state: ", (ebi_spec_names(i_spec)%string, &
              YC(i_spec), i_spec=1,size(YC))
      write(DEBUG_UNIT,*) "KPP state: ", (TRIM(KPP_SPC_NAMES(i_spec)), &
              KPP_C(i_spec)*conv, i_spec=1,size(KPP_C))
      write(DEBUG_UNIT,*) "Phlex state: ", (phlex_spec_names(i_spec)%string, &
              phlex_state%state_var(i_spec), i_spec=1,size(phlex_state%state_var))
      end if
#endif

      ! Set minimum concentrations in all solvers
      YC(:) = MAX(YC(:), SMALL_NUM)
      KPP_C(:) = MAX(KPP_C(:), SMALL_NUM/conv)
      phlex_state%state_var(:) = max(phlex_state%state_var(:), SMALL_NUM)

      ! Output current time step
      if (i_repeat.eq.1) then
      write(EBI_FILE_UNIT,*) i_time*EBI_TMSTEP, YC(:)
      write(KPP_FILE_UNIT,*) i_time*KPP_DT/60.0d0, KPP_C(:)*conv
      write(PHLEX_FILE_UNIT,*) i_time*EBI_TMSTEP, phlex_state%state_var(:)
      end if

      ! EBI solver
      call cpu_time(comp_start)
      call EXT_HRCALCKS( NUM_EBI_PHOTO_RXN,       & ! Number of EBI solver photolysis reactions
                         is_sunny,                & ! Flag for sunlight
                         photo_rates,             & ! Photolysis rates
                         temperature,             & ! Temperature (K)
                         pressure,                & ! Air pressure (atm)
                         water_conc,              & ! Water vapor concentration (ppmV)
                         RKI)                       ! Rate constants
      call EXT_HRSOLVER( 2018012, 070000, 1, 1, 1)  ! These dummy variables are just for output
      call cpu_time(comp_end)
      comp_ebi = comp_ebi + (comp_end-comp_start)

      ! Set KPP and phlex-chem concentrations to those of EBI at first time step to match steady-state
      ! EBI species
      if (i_time.eq.1) then
        ! Set KPP species in #/cc
        do i_spec = 1, NUM_EBI_SPEC
          do j_spec = 1, KPP_NSPEC
            if (trim(ebi_spec_names(i_spec)%string).eq.trim(KPP_SPC_NAMES(j_spec))) then
              KPP_C(j_spec) = YC(i_spec) / conv
            end if
          end do
          phlex_state%state_var( &
                  phlex_core%chem_spec_data%gas_state_id( &
                  ebi_spec_names(i_spec)%string)) = YC(i_spec)
        end do
      end if

      ! KPP module
      call cpu_time(comp_start)
      KPP_TIME = (i_time-1)*KPP_DT
      KPP_TEMP = temperature
      KPP_PRESS = pressure * 1013.25 ! KPP pressure in hPa
      CALL KPP_Update_RCONST()
      CALL KPP_INTEGRATE( TIN = KPP_TIME, TOUT = (KPP_TIME+KPP_DT), &
              RSTATUS_U = KPP_RSTATE, ICNTRL_U = (/ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) )
      call cpu_time(comp_end)
      comp_kpp = comp_kpp + (comp_end-comp_start)

      ! Phlex-chem
      call cpu_time(comp_start)
      call phlex_core%solve(phlex_state, real(EBI_TMSTEP*60.0, kind=dp))
      call cpu_time(comp_end)
      comp_phlex = comp_phlex + (comp_end-comp_start)

    end do 
    end do

    ! Output final timestep
    write(EBI_FILE_UNIT,*) i_time*EBI_TMSTEP, YC(:)
    write(KPP_FILE_UNIT,*) i_time*EBI_TMSTEP, KPP_C(:)*conv
    write(PHLEX_FILE_UNIT,*) i_time*EBI_TMSTEP, phlex_state%state_var(:)

    ! Output the computational time
    write(*,*) "EBI calculation time: ", comp_ebi," s"
    write(*,*) "KPP calculation time: ", comp_kpp," s"
    write(*,*) "Phlex-chem calculation time: ", comp_phlex," s"

    ! Compare the results
    ! EBI <-> Phlex-chem
    do i_spec = 1, NUM_EBI_SPEC
      call warn_assert_msg(749090387, almost_equal(real(YC(i_spec), kind=dp), &
          phlex_state%state_var( &
                  phlex_core%chem_spec_data%gas_state_id( &
                  ebi_spec_names(i_spec)%string)), 5.0d-2), &
          "Species "//ebi_spec_names(i_spec)%string//" has different result. "// &
          "EBI solver: "//trim(to_string(real(YC(i_spec), kind=dp)))// &
          "; Phlex-chem: "// &
          trim(to_string( phlex_state%state_var( &
                  phlex_core%chem_spec_data%gas_state_id( &
                  ebi_spec_names(i_spec)%string)))))
    end do
    ! KPP <-> Phlex-chem
    do i_spec = 1, KPP_NSPEC
      str_temp%string = trim(KPP_SPC_NAMES(i_spec))
      call warn_assert_msg(749090436, almost_equal(real(KPP_C(i_spec)*conv, kind=dp), &
          phlex_state%state_var( &
                  phlex_core%chem_spec_data%gas_state_id( &
                  str_temp%string)), 5.0d-2), &
          "Species "//str_temp%string//" has different result. "// &
          "KPP solver: "//trim(to_string(real(KPP_C(i_spec)*conv, kind=dp)))// &
          "; Phlex-chem: "// &
          trim(to_string( phlex_state%state_var( &
                  phlex_core%chem_spec_data%gas_state_id( &
                  str_temp%string)))))
    end do

    ! Close the output files
    close(EBI_FILE_UNIT)
    close(KPP_FILE_UNIT)
    close(PHLEX_FILE_UNIT)

    ! Create the gnuplot script
    open(unit=PHLEX_FILE_UNIT, file="out/plot_cb05cl_ae.conf", status="replace", action="write")
    write(PHLEX_FILE_UNIT,*) "# plot_cb05cl_ae5.conf"
    write(PHLEX_FILE_UNIT,*) "# Run as: gnuplot plot_cb05cl_ae5.conf"
    write(PHLEX_FILE_UNIT,*) "set terminal png truecolor"
    write(PHLEX_FILE_UNIT,*) "set autoscale"
    do i_spec = 1, phlex_core%chem_spec_data%size()
      call assert(476528007, phlex_core%chem_spec_data%get_name(i_spec, spec_name))
      write(PHLEX_FILE_UNIT,*) "set output 'cb05cl_ae5_"//spec_name//".png'"
      write(PHLEX_FILE_UNIT,*) "plot\"
      if (ebi_spec_map(i_spec).gt.0) then
        write(PHLEX_FILE_UNIT,*) " 'cb05cl_ae5_ebi_results.txt'\"
        write(PHLEX_FILE_UNIT,*) " using 1:"//trim(to_string(ebi_spec_map(i_spec)+1))//" title '"// &
              spec_name//" (ebi)',\"
      end if
      if (kpp_spec_map(i_spec).gt.0) then
        write(PHLEX_FILE_UNIT,*) " 'cb05cl_ae5_kpp_results.txt'\"
        write(PHLEX_FILE_UNIT,*) " using 1:"//trim(to_string(kpp_spec_map(i_spec)+1))//" title '"// &
              spec_name//" (kpp)',\"
      end if
      write(PHLEX_FILE_UNIT,*) " 'cb05cl_ae5_phlex_results.txt'\"
      write(PHLEX_FILE_UNIT,*) " using 1:"//trim(to_string(i_spec+1))//" title '"// &
              spec_name//" (phlex)'"
    end do
    close(PHLEX_FILE_UNIT)
    
    passed = .true.

  end function run_standard_cb05cl_ae5_test

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Set the EBI-solver species names
  subroutine set_ebi_species(spec_names)

    !> EBI solver species names
    type(string_t), dimension(NUM_EBI_SPEC) :: spec_names

    spec_names(1)%string = "NO2"
    spec_names(2)%string = "NO"
    spec_names(3)%string = "O"
    spec_names(4)%string = "O3"
    spec_names(5)%string = "NO3"
    spec_names(6)%string = "O1D"
    spec_names(7)%string = "OH"
    spec_names(8)%string = "HO2"
    spec_names(9)%string = "N2O5"
    spec_names(10)%string = "HNO3"
    spec_names(11)%string = "HONO"
    spec_names(12)%string = "PNA"
    spec_names(13)%string = "H2O2"
    spec_names(14)%string = "XO2"
    spec_names(15)%string = "XO2N"
    spec_names(16)%string = "NTR"
    spec_names(17)%string = "ROOH"
    spec_names(18)%string = "FORM"
    spec_names(19)%string = "ALD2"
    spec_names(20)%string = "ALDX"
    spec_names(21)%string = "PAR"
    spec_names(22)%string = "CO"
    spec_names(23)%string = "MEO2"
    spec_names(24)%string = "MEPX"
    spec_names(25)%string = "MEOH"
    spec_names(26)%string = "HCO3"
    spec_names(27)%string = "FACD"
    spec_names(28)%string = "C2O3"
    spec_names(29)%string = "PAN"
    spec_names(30)%string = "PACD"
    spec_names(31)%string = "AACD"
    spec_names(32)%string = "CXO3"
    spec_names(33)%string = "PANX"
    spec_names(34)%string = "ROR"
    spec_names(35)%string = "OLE"
    spec_names(36)%string = "ETH"
    spec_names(37)%string = "IOLE"
    spec_names(38)%string = "TOL"
    spec_names(39)%string = "CRES"
    spec_names(40)%string = "TO2"
    spec_names(41)%string = "TOLRO2"
    spec_names(42)%string = "OPEN"
    spec_names(43)%string = "CRO"
    spec_names(44)%string = "MGLY"
    spec_names(45)%string = "XYL"
    spec_names(46)%string = "XYLRO2"
    spec_names(47)%string = "ISOP"
    spec_names(48)%string = "ISPD"
    spec_names(49)%string = "ISOPRXN"
    spec_names(50)%string = "TERP"
    spec_names(51)%string = "TRPRXN"
    spec_names(52)%string = "SO2"
    spec_names(53)%string = "SULF"
    spec_names(54)%string = "SULRXN"
    spec_names(55)%string = "ETOH"
    spec_names(56)%string = "ETHA"
    spec_names(57)%string = "CL2"
    spec_names(58)%string = "CL"
    spec_names(59)%string = "HOCL"
    spec_names(60)%string = "CLO"
    spec_names(61)%string = "FMCL"
    spec_names(62)%string = "HCL"
    spec_names(63)%string = "TOLNRXN"
    spec_names(64)%string = "TOLHRXN"
    spec_names(65)%string = "XYLNRXN"
    spec_names(66)%string = "XYLHRXN"
    spec_names(67)%string = "BENZENE"
    spec_names(68)%string = "BENZRO2"
    spec_names(69)%string = "BNZNRXN"
    spec_names(70)%string = "BNZHRXN"
    spec_names(71)%string = "SESQ"
    spec_names(72)%string = "SESQRXN"

  end subroutine set_ebi_species

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Load a string array with KPP reaction labels
  subroutine get_kpp_rxn_labels(kpp_rxn_labels)

    use cb05cl_ae5_Model,                       only : NREACT

    type(string_t), allocatable :: kpp_rxn_labels(:)
    integer(kind=i_kind) :: i_rxn = 1
 
    allocate(kpp_rxn_labels(NREACT))

     kpp_rxn_labels(i_rxn)%string = 'R1'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R2'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R3'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R4'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R5'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R6'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R7'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R8'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R9'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R10'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R11'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R12'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R13'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R14'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R15'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R16'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R17'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R18'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R19'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R20'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R21'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R22'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R23'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R24'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R25'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R26'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R27'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R28'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R29'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R30'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R31'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R32'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R33'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R34'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R35'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R36'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R37'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R38'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R39'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R40'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R41'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R42'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R43'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R44'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R45'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R46'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R47'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R48'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R49'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R50'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R51'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R52'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R53'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R54'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R55'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R56'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R57'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R58'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R59'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R60'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R61'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R62'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R63'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R64'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R65'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R66'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R67'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R68'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R69'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R70'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R71'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R72'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R73'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R74'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R75'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R76'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R77'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R78'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R79'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R80'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R81'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R82'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R83'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R84'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R85'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R86'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R87'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R88'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R89'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R90'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R91'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R92'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R93'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R94'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R95'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R96'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R97'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R98'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R99'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R100'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R101'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R102'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R103'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R104'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R105'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R106'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R107'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R108'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R109'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R110'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R111'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R112'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R113'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R114'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R115'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R116'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R117'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R118'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R119'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R120'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R121'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R122'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R123'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R124'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R125'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R126'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R127'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R128'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R129'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R130'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R131'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R132'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R133'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R134'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R135'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R136'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R137'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R138'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R139'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R140'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R141'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R142'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R143'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R144'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R145'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R146'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R147'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R148'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R149'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R150'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R151'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R152'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R153'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R154'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R155'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'R156'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'CL1'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'CL2'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'CL3'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'CL4'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'CL5'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'CL6'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'CL7'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'CL8'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'CL9'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'CL10'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'CL11'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'CL12'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'CL13'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'CL14'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'CL15'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'CL16'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'CL17'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'CL18'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'CL19'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'CL20'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'CL21'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'SA01'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'SA02'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'SA03'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'SA04'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'SA05'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'SA06'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'SA07'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'SA08'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'SA09'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'SA10'
     i_rxn = i_rxn + 1
     kpp_rxn_labels(i_rxn)%string = 'jo2'

     call assert_msg(313903826, i_rxn.eq.NREACT, "Labeled "//trim(to_string(i_rxn))// &
             " out of "//trim(to_string(NREACT))//" KPP reactions")

  end subroutine get_kpp_rxn_labels

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Compare calculated rates between the modules
  subroutine compare_rates(phlex_core, phlex_state, ebi_spec_names, conv, &
                  ebi_rxn_map, kpp_rxn_map)

    use EXT_RXCM,                               only : NRXNS, RXLABEL
    use EXT_HRDATA,                             only : EBI_PROD => PROD, &
                                                       EBI_LOSS => LOSS, &
                                                       EBI_YC => YC, &
                                                       EBI_RXRAT => RXRAT
    use cb05cl_ae5_Function,                    only : KPP_Fun => Fun, &
                                                       KPP_A => A
    use cb05cl_ae5_Model,                       only : KPP_NVAR => NVAR, &
                                                       KPP_SPC_NAMES => SPC_NAMES, &
                                                       KPP_VAR => VAR, &
                                                       KPP_FIX => FIX, &
                                                       KPP_RCONST => RCONST

    !> Phlex-chem core
    type(phlex_core_t), pointer :: phlex_core
    !> Phlex-chem state
    type(phlex_state_t), pointer :: phlex_state
    !> Species map Phlex-chem <-> EBI
    type(string_t), intent(in) :: ebi_spec_names(:)
    !> Coversion factor #/cc -> ppm
    real(kind=dp), intent(in) :: conv
    !> Reaction map EBI <-> phlex
    integer(kind=i_kind), allocatable :: ebi_rxn_map(:)
    !> Reaction map KPP <-> phlex
    integer(kind=i_kind), allocatable :: kpp_rxn_map(:)

    real(kind=dp) :: KPP_VDOT(KPP_NVAR)               
    integer(kind=i_kind) :: i_ebi_spec, i_kpp_spec, i_phlex_spec, i_mech, i_rxn
    real(kind=dp), allocatable :: phlex_rates(:)
    type(string_t) :: spec_name

    character(len=:), allocatable :: key

    write(*,*) "Comparing rates"

    call KPP_Fun ( KPP_VAR, KPP_FIX, KPP_RCONST, KPP_VDOT )

    call EXT_HRRATES()
    call EXT_HRPRODLOSS()

    ! Find the phlex-core mechanism
    key = "cb05cl_ae5"
    call assert(331333207, phlex_core%find_mechanism(key, i_mech))

    ! Get the current reaction rates
    phlex_rates = phlex_core%get_rxn_rates(phlex_state)    
    call assert_msg(743625443, size(phlex_rates).eq.phlex_core%mechanism(i_mech)%size(), &
            "Rate array is wrong size. Expected "// &
            to_string(phlex_core%mechanism(i_mech)%size())// &
            " but got "//to_string(size(phlex_rates)))

    ! Compare the calculated rates
    ! phlex <-> KPP
    do i_rxn = 1, phlex_core%mechanism(i_mech)%size()
      if ((.not.almost_equal(KPP_A(kpp_rxn_map(i_rxn))*conv, phlex_rates(i_rxn), 1.0d-2))) then
        write(*,*) "Rxn(", trim(RXLABEL(ebi_rxn_map(i_rxn))), ") rate: (kpp):", &
              KPP_A(kpp_rxn_map(i_rxn))*conv, " (ebi):", EBI_RXRAT(ebi_rxn_map(i_rxn))/60.0d0, &
              " (phlex):", phlex_rates(i_rxn)
      end if
    end do

    ! Compare the calculated rates
    ! EBI <-> KPP
    do i_ebi_spec = 1, NUM_EBI_SPEC
      do i_kpp_spec = 1, KPP_NVAR
        ! Skip EBI PSSA species
        if (trim(ebi_spec_names(i_ebi_spec)%string).eq.'NO2' .or. &
            trim(ebi_spec_names(i_ebi_spec)%string).eq.'NO' .or. &    
            trim(ebi_spec_names(i_ebi_spec)%string).eq.'O' .or. &    
            trim(ebi_spec_names(i_ebi_spec)%string).eq.'O3' .or. &    
            trim(ebi_spec_names(i_ebi_spec)%string).eq.'NO3' .or. &    
            trim(ebi_spec_names(i_ebi_spec)%string).eq.'O1D' .or. &    
            trim(ebi_spec_names(i_ebi_spec)%string).eq.'OH' .or. &    
            trim(ebi_spec_names(i_ebi_spec)%string).eq.'HO2' .or. &    
            trim(ebi_spec_names(i_ebi_spec)%string).eq.'N2O5' .or. &    
            trim(ebi_spec_names(i_ebi_spec)%string).eq.'HONO' .or. &    
            trim(ebi_spec_names(i_ebi_spec)%string).eq.'C2O3' .or. &    
            trim(ebi_spec_names(i_ebi_spec)%string).eq.'PAN' .or. &    
            trim(ebi_spec_names(i_ebi_spec)%string).eq.'PNA') cycle
        if (trim(ebi_spec_names(i_ebi_spec)%string).eq.trim(KPP_SPC_NAMES(i_kpp_spec))) then
          call warn_assert_msg(616862348, almost_equal(real(KPP_VDOT(i_kpp_spec)*60.0d0*conv, kind=dp), &
              real(EBI_PROD(i_ebi_spec)-EBI_LOSS(i_ebi_spec), kind=dp), 1.0d-2), &
              "Species "//ebi_spec_names(i_ebi_spec)%string//" has different time deriv. "// &
              "EBI solver: "//trim(to_string(real((EBI_PROD(i_ebi_spec)-EBI_LOSS(i_ebi_spec))/60.0d0, kind=dp)))// &
              "; KPP solver: "//trim(to_string(real(KPP_VDOT(i_kpp_spec)*conv, kind=dp))))
        end if
      end do
    end do

  end subroutine compare_rates

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program pmc_test_cb05cl_ae5
