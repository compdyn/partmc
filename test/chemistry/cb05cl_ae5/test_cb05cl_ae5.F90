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
  integer(kind=i_kind), parameter :: NUM_EBI_SPEC = 72
  ! Number of EBI-solever photolysis reactions
  integer(kind=i_kind), parameter :: NUM_EBI_PHOTO_RXN = 23

#ifdef DEBUG
  integer(kind=i_kind), parameter :: DEBUG_UNIT = 12
   
  open(unit=DEBUG_UNIT, file="out/debug_cb05cl_ae.txt", status="replace", action="write")
#endif

  if (run_cb05cl_ae5_tests()) then
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
  !! MONARCH ebi-solver code and the Phlex-chem version
  logical function run_standard_cb05cl_ae5_test() result(passed)

    use EXT_HRDATA
    use EXT_RXCM,                               only : NRXNS, RXLABEL

    ! EBI-solver species names
    type(string_t), dimension(NUM_EBI_SPEC) :: ebi_spec_names

    ! Flag for sunlight
    logical :: is_sunny = .true.
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
    type(phlex_state_t), target :: phlex_state, phlex_state_comp
    ! Phlex-chem species names
    type(string_t), allocatable :: phlex_spec_names(:)
    ! EBI -> Phlex-chem species map
    integer(kind=i_kind), dimension(NUM_EBI_SPEC) :: spec_map

    ! Computation timer variables
    real(kind=dp) :: comp_start, comp_end, comp_ebi, comp_phlex

    class(rxn_data_t), pointer :: rxn
    type(property_t), pointer :: prop_set
    character(len=:), allocatable :: key, spec_name, string_val, phlex_input_file
    real(kind=dp) :: real_val, phlex_rate, phlex_rate_const
    integer(kind=i_kind) :: i_spec, i_mech, i_rxn, i_ebi_rxn, i_time

    integer(kind=i_kind) :: i_M, i_O2, i_N2, i_H2O, i_CH4, i_H2
    integer(kind=i_kind), allocatable :: rxn_map(:)

    passed = .false.
    
    ! Load the EBI solver species names
    call set_ebi_species(ebi_spec_names)

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
    phlex_input_file = "config_cb05cl_ae5.json"
    phlex_core => phlex_core_t(phlex_input_file)
    call phlex_core%initialize()
    phlex_state = phlex_core%new_state()
    phlex_state%env_state%temp = temperature
    phlex_state%env_state%pressure = pressure * const%air_std_press
    call cpu_time(comp_end)
    write(*,*) "Phlex-chem initialization time: ", comp_end-comp_start," s"

    ! Get a phlex-state for rate comparisons and find constant species
    phlex_state_comp = phlex_core%new_state()
    phlex_state_comp%env_state%temp = phlex_state%env_state%temp
    phlex_state_comp%env_state%pressure = phlex_state%env_state%pressure
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
    allocate(photo_rates(NUM_EBI_PHOTO_RXN))
    photo_rates(:) = 0.01 * 60.0 ! EBI solver wants rates in min^-1
    key = "cb05cl_ae5"
    call assert(331333207, phlex_core%find_mechanism(key, i_mech))
    do i_rxn = 1, phlex_core%mechanism(i_mech)%size()
      rxn => phlex_core%mechanism(i_mech)%rxn_ptr(i_rxn)%val
      select type(rxn) 
        type is (rxn_photolysis_t)
          call rxn%set_rate_const(real(0.01, kind=dp))
      end select
    end do

    ! Set the initial concentrations
    key = "init conc"
    YC(:) = 0.0
    phlex_state%state_var(:) = 0.0
    do i_spec = 1, NUM_EBI_SPEC
      prop_set => phlex_core%chem_spec_data%get_property_set( &
              ebi_spec_names(i_spec)%string)
      if (associated(prop_set)) then
        if (prop_set%get_real(key, real_val)) then
          YC(i_spec) = real_val
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
#ifdef DEBUG
      else
        write(DEBUG_UNIT,*) "No Properties for species ", ebi_spec_names(i_spec)%string
#endif
      end if
    end do
    
    ! Set EBI solver constant species concentrations in Phlex-chem
    spec_name = "M"
    prop_set => phlex_core%chem_spec_data%get_property_set(spec_name)
    call assert(740666066, associated(prop_set))
    call assert(907464197, prop_set%get_real(key, real_val))
    phlex_state%state_var(i_M) = real_val
    spec_name = "O2"
    prop_set => phlex_core%chem_spec_data%get_property_set(spec_name)
    call assert(729136508, associated(prop_set))
    call assert(223930103, prop_set%get_real(key, real_val))
    phlex_state%state_var(i_O2) = real_val
    spec_name = "N2"
    prop_set => phlex_core%chem_spec_data%get_property_set(spec_name)
    call assert(553715297, associated(prop_set))
    call assert(666033642, prop_set%get_real(key, real_val))
    phlex_state%state_var(i_N2) = real_val
    spec_name = "H2O"
    prop_set => phlex_core%chem_spec_data%get_property_set(spec_name)
    call assert(160827237, associated(prop_set))
    call assert(273145582, prop_set%get_real(key, real_val))
    phlex_state%state_var(i_H2O) = real_val
    spec_name = "CH4"
    prop_set => phlex_core%chem_spec_data%get_property_set(spec_name)
    call assert(667939176, associated(prop_set))
    call assert(780257521, prop_set%get_real(key, real_val))
    phlex_state%state_var(i_CH4) = real_val
    spec_name = "H2"
    prop_set => phlex_core%chem_spec_data%get_property_set(spec_name)
    call assert(892575866, associated(prop_set))
    call assert(722418962, prop_set%get_real(key, real_val))
    phlex_state%state_var(i_H2) = real_val

    ! Set the water concentration for EBI solver (ppmV)
    water_conc = phlex_state%state_var(i_H2O)

    ! EBI solver species have a minimum value set
    YC(:) = MAX(YC(:), 1e-30)
    phlex_state%state_var(:) = MAX(phlex_state%state_var(:), 1e-30)

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
    open(PHLEX_FILE_UNIT, file="out/cb05cl_ae5_phlex_results.txt", status="replace", &
            action="write")
    phlex_spec_names = phlex_core%chem_spec_data%spec_names_by_type(GAS_SPEC)
    write(EBI_FILE_UNIT,*) "time ", (ebi_spec_names(i_spec)%string//" ", i_spec=1, &
            NUM_EBI_SPEC)
    write(PHLEX_FILE_UNIT,*) "time ", (phlex_spec_names(i_spec)%string//" ", i_spec=1, &
            size(phlex_spec_names))

    ! Set up the reaction map between phlex-chem and ebi solver
    key = "rxn id"
    allocate(rxn_map(phlex_core%mechanism(i_mech)%size()))
    rxn_map(:) = 0
    do i_rxn = 1, phlex_core%mechanism(i_mech)%size()
      rxn => phlex_core%mechanism(i_mech)%rxn_ptr(i_rxn)%val
      call assert_msg(917216189, associated(rxn), "Missing rxn "//to_string(i_rxn))
      call rxn%get_test_info(phlex_state, phlex_rate, phlex_rate_const, prop_set)
      call assert(656034097, prop_set%get_string(key, string_val))
      do i_ebi_rxn = 1, NRXNS 
        if (trim(RXLABEL(i_ebi_rxn)).eq.trim(string_val)) then
          rxn_map(i_rxn) = i_ebi_rxn
          exit
        end if
      end do
      call assert_msg(921715481, rxn_map(i_rxn).ne.0, "Missing rxn "//string_val)
    end do

    ! Reset the computation timers
    comp_ebi = 0.0
    comp_phlex = 0.0

    ! Solve the mechanism
    do i_time = 1, NUM_TIME_STEPS

#ifdef DEBUG
      write(*,*) "Time step: ", i_time
      write(DEBUG_UNIT,*) "Phlex state: ", phlex_state%state_var(:)
#endif

      ! Output current time step
      write(EBI_FILE_UNIT,*) i_time*EBI_TMSTEP, YC(:)
      write(PHLEX_FILE_UNIT,*) i_time*EBI_TMSTEP, phlex_state%state_var(:)

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

#ifdef DEBUG
      ! Compare rate constants
      phlex_state_comp%state_var(:) = 1.0
      phlex_state_comp%state_var(i_M)   = phlex_state%state_var(i_M)
      phlex_state_comp%state_var(i_O2)  = phlex_state%state_var(i_O2)
      phlex_state_comp%state_var(i_N2)  = phlex_state%state_var(i_N2)
      phlex_state_comp%state_var(i_H2O) = phlex_state%state_var(i_H2O)
      phlex_state_comp%state_var(i_CH4) = phlex_state%state_var(i_CH4)
      phlex_state_comp%state_var(i_H2)  = phlex_state%state_var(i_H2)
      do i_rxn = 1, phlex_core%mechanism(i_mech)%size()
        rxn => phlex_core%mechanism(i_mech)%rxn_ptr(i_rxn)%val
        call rxn%get_test_info(phlex_state_comp, phlex_rate, phlex_rate_const, prop_set)
        ! EBI solver rates are in min^-1
        if (.not.almost_equal(phlex_rate*60.0, RKI(rxn_map(i_rxn)), 1.0d-4)) then
          write(DEBUG_UNIT,*) "Different rate constants for reaction "// &
             trim(RXLABEL(rxn_map(i_rxn)))//": ",phlex_rate*60.0, RKI(rxn_map(i_rxn))
        endif
        call rxn%get_test_info(phlex_state, phlex_rate, phlex_rate_const, prop_set)
        write(DEBUG_UNIT,*) trim(RXLABEL(rxn_map(i_rxn))), "Rate constant: ", &
                phlex_rate_const, "; rate: ", phlex_rate
      end do
#endif

      ! Phlex-chem
      call cpu_time(comp_start)
      call phlex_core%solve(phlex_state, real(EBI_TMSTEP, kind=dp))
      call cpu_time(comp_end)
      comp_phlex = comp_phlex + (comp_end-comp_start)

      ! Keep concentrations positive
      phlex_state%state_var(:) = MAX(phlex_state%state_var(:), 1.0d-30)

    end do 

    ! Output final timestep
    write(EBI_FILE_UNIT,*) i_time*EBI_TMSTEP, YC(:)
    write(PHLEX_FILE_UNIT,*) i_time*EBI_TMSTEP, phlex_state%state_var(:)

    ! Output the computational time
    write(*,*) "EBI calculation time: ", comp_ebi," s"
    write(*,*) "Phlex-chem calculation time: ", comp_phlex," s"

    ! Compare the results
    do i_spec = 1, NUM_EBI_SPEC
      call warn_assert_msg(749090387, almost_equal(real(YC(i_spec), kind=dp), &
          phlex_state%state_var( &
                  phlex_core%chem_spec_data%gas_state_id( &
                  ebi_spec_names(i_spec)%string)), 1.0d-2), &
          "Species "//ebi_spec_names(i_spec)%string//" has different result. "// &
          "EBI solver: "//trim(to_string(real(YC(i_spec), kind=dp)))// &
          "; Phlex-chem: "// &
          trim(to_string( phlex_state%state_var( &
                  phlex_core%chem_spec_data%gas_state_id( &
                  ebi_spec_names(i_spec)%string)))))
    end do

    ! Close the output files
    close(EBI_FILE_UNIT)
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

end program pmc_test_cb05cl_ae5
