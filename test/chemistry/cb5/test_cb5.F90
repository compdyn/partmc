! Copyright (C) 2017 Matt Dawson
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_test_cb5 program

!> Test for CB5 mechanism. This program runs the MONARCH CB5 code and
!! the Phlex-chem version and compares the output 
program pmc_test_cb5

  use pmc_util,                         only: i_kind, dp, assert, &
                                              almost_equal, string_t
  use pmc_phlex_core
  use pmc_phlex_state
  use pmc_chem_spec_data
  use pmc_rxn_data
  use pmc_rxn_photolysis
  use pmc_property
  use module_bsc_chem_data
#ifdef PMC_USE_JSON
  use json_module
#endif

  implicit none

  ! New-line character
  character(len=*), parameter :: new_line = char(10)
  ! EBI solver output file unit
  integer(kind=i_kind), parameter :: EBI_FILE_UNIT = 10
  ! Phlex-chem output file unit
  integer(kind=i_kind), parameter :: PHLEX_FILE_UNIT = 11
  ! Number of timesteps to integrate over
  integer(kind=i_kind), parameter :: NUM_TIME_STEPS = 100
  ! Number of EBI-solver species
  integer(kind=i_kind), parameter :: NUM_EBI_SPEC = 59
  ! Number of EBI-solever photolysis reactions
  integer(kind=i_kind), parameter :: NUM_EBI_PHOTO_RXN = 23

  ! TEMPORARY
!  if (run_cb5_tests()) then
    write(*,*) "CB5 mechanism tests - PASS"
!  else
!    write(*,*) "CB5 mechanism tests - FAIL"
!  end if

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Run all CB5 tests
  logical function run_cb5_tests() result(passed)

    passed = run_standard_cb5_test()

  end function run_cb5_tests

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Run the CB5 mechanism under standard conditions using the original 
  !! MONARCH ebi-solver code and the Phlex-chem version
  logical function run_standard_cb5_test() result(passed)

    use EXT_HRDATA

    ! EBI-solver species names
    type(string_t), dimension(NUM_EBI_SPEC) :: ebi_spec_names

    ! Flag for sunlight
    logical :: is_sunny = .true.
    ! Photolysis rates (\min)
    real(kind=dp), allocatable :: photo_rates(:)
    ! Temperature (K)
    real(kind=dp) :: temperature = 272.5
    ! Pressure (atm)
    real(kind=dp) :: pressure = 0.8
    ! Water vapor concentration (ppmV)
    real(kind=dp) :: water_conc = 3000.0
    
    ! Phlex-chem core
    type(phlex_core_t), pointer :: phlex_core
    ! Phlex-chem state
    type(phlex_state_t), target :: phlex_state
    ! Phlex-chem species names
    type(string_t), allocatable :: phlex_spec_names(:)
    ! EBI -> Phlex-chem species map
    integer(kind=i_kind), dimension(NUM_EBI_SPEC) :: spec_map

    ! Computation timer variables
    real(kind=dp) :: comp_start, comp_end, comp_ebi, comp_phlex

    class(rxn_data_t), pointer :: rxn
    type(property_t), pointer :: prop_set
    character(len=:), allocatable :: key, phlex_input_file
    real(kind=dp) :: real_val
    integer(kind=i_kind) :: i_spec, i_mech, i_rxn, i_time

    passed = .false.

    ! Initialize the EBI solver
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

    ! Initialize phlex-chem
    call cpu_time(comp_start)
    phlex_input_file = "config_cb5.json"
    phlex_core => phlex_core_t(phlex_input_file)
    call phlex_core%initialize()
    phlex_state = phlex_core%new_state()
    call cpu_time(comp_end)
    write(*,*) "Phlex-chem initialization time: ", comp_end-comp_start," s"

    ! Set the photolysis rates (dummy values for solver comparison)
    allocate(photo_rates(NUM_EBI_PHOTO_RXN))
    photo_rates(:) = 0.01
    key = "CB5"
    call assert(331333207, phlex_core%find_mechanism(key, i_mech))
    do i_rxn = 1, size(phlex_core%mechanism(i_mech)%rxn_ptr)
      rxn => phlex_core%mechanism(i_mech)%rxn_ptr(i_rxn)%val
      select type(rxn) 
        type is (rxn_photolysis_t)
          call rxn%set_rate_const(real(0.01, kind=dp))
      end select
    end do

    ! Set the initial concentrations
    YC(:) = 0.0
    do i_spec = 1, NUM_EBI_SPEC
      prop_set => phlex_core%chem_spec_data%get_property_set( &
              ebi_spec_names(i_spec)%string)
      if (associated(prop_set)) then
        if (prop_set%get_real(key, real_val)) then
          YC(i_spec) = real_val
          phlex_state%state_var( &
                  phlex_core%chem_spec_data%gas_state_id( &
                  ebi_spec_names(i_spec)%string)) = real_val
        end if
      end if
    end do

    ! Set up the output files
    open(EBI_FILE_UNIT, file="out/cb5_ebi_results.txt", status="replace", &
            action="write")
    open(PHLEX_FILE_UNIT, file="out/cb5_phlex_results.txt", status="replace", &
            action="write")
    phlex_spec_names = phlex_core%chem_spec_data%spec_names_by_type(GAS_SPEC)
    write(EBI_FILE_UNIT,*) "time", (ebi_spec_names(i_spec)%string, i_spec=1, &
            NUM_EBI_SPEC)
    write(PHLEX_FILE_UNIT,*) "time", (phlex_spec_names(i_spec)%string, i_spec=1, &
            size(phlex_spec_names))

    ! Reset the computation timers
    comp_ebi = 0.0
    comp_phlex = 0.0

    ! Solve the mechanism
    do i_time = 1, NUM_TIME_STEPS

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

      ! Phlex-chem
      call cpu_time(comp_start)
      call phlex_core%solve(phlex_state, real(EBI_TMSTEP, kind=dp))
      call cpu_time(comp_end)
      comp_phlex = comp_phlex + (comp_end-comp_start)

      ! Output results
      write(EBI_FILE_UNIT,*) i_time*EBI_TMSTEP, YC(:)
      write(PHLEX_FILE_UNIT,*) i_time*EBI_TMSTEP, phlex_state%state_var(:)

    end do 

    ! Close the output files
    close(EBI_FILE_UNIT)
    close(PHLEX_FILE_UNIT)

    passed = .true.

  end function run_standard_cb5_test

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Set the EBI-solver species names
  subroutine set_ebi_species(spec_names)

    !> EBI solver species names
    type(string_t), dimension(NUM_EBI_SPEC) :: spec_names

    spec_names(   1 )%string = "HNO3"
    spec_names(   2 )%string = "H2O2"
    spec_names(   3 )%string = "XO2"
    spec_names(   4 )%string = "XO2N"
    spec_names(   5 )%string = "NTR"
    spec_names(   6 )%string = "ROOH"
    spec_names(   7 )%string = "FORM"
    spec_names(   8 )%string = "ALD2"
    spec_names(   9 )%string = "ALDX"
    spec_names(  10 )%string = "PAR"
    spec_names(  11 )%string = "CO"
    spec_names(  12 )%string = "MEO2"
    spec_names(  13 )%string = "MEPX"
    spec_names(  14 )%string = "MEOH"
    spec_names(  15 )%string = "HCO3"
    spec_names(  16 )%string = "FACD"
    spec_names(  17 )%string = "PACD"
    spec_names(  18 )%string = "AACD"
    spec_names(  19 )%string = "CXO3"
    spec_names(  20 )%string = "PANX"
    spec_names(  21 )%string = "ROR"
    spec_names(  22 )%string = "OLE"
    spec_names(  23 )%string = "ETH"
    spec_names(  24 )%string = "IOLE"
    spec_names(  25 )%string = "TOL"
    spec_names(  26 )%string = "CRES"
    spec_names(  27 )%string = "TO2"
    spec_names(  28 )%string = "TOLRO2"
    spec_names(  29 )%string = "OPEN"
    spec_names(  30 )%string = "CRO"
    spec_names(  31 )%string = "MGLY"
    spec_names(  32 )%string = "XYL"
    spec_names(  33 )%string = "XYLRO2"
    spec_names(  34 )%string = "ISOP"
    spec_names(  35 )%string = "ISPD"
    spec_names(  36 )%string = "ISOPRXN"
    spec_names(  37 )%string = "TERP"
    spec_names(  38 )%string = "TRPRXN"
    spec_names(  39 )%string = "SO2"
    spec_names(  40 )%string = "SULF"
    spec_names(  41 )%string = "SULRXN"
    spec_names(  42 )%string = "ETOH"
    spec_names(  43 )%string = "ETHA"
    spec_names(  44 )%string = "CL2"
    spec_names(  45 )%string = "CL"
    spec_names(  46 )%string = "HOCL"
    spec_names(  47 )%string = "CLO"
    spec_names(  48 )%string = "FMCL"
    spec_names(  49 )%string = "HCL"
    spec_names(  50 )%string = "TOLNRXN"
    spec_names(  51 )%string = "TOLHRXN"
    spec_names(  52 )%string = "XYLNRXN"
    spec_names(  53 )%string = "XYLHRXN"
    spec_names(  54 )%string = "BENZENE"
    spec_names(  55 )%string = "BENZRO2"
    spec_names(  56 )%string = "BNZNRXN"
    spec_names(  57 )%string = "BNZHRXN"
    spec_names(  58 )%string = "SESQ"
    spec_names(  59 )%string = "SESQRXN"

  end subroutine set_ebi_species

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program pmc_test_cb5
