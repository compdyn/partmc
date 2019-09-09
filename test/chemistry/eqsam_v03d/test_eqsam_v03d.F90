! Copyright (C) 2017 Matt Dawson
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_test_eqsam_v03d program

!> Test for the eqsam v03d inorganic module from MONARCH. This program runs the
!! MONARCH eqsam code and the CAMP-chem version and compares the output.
program pmc_test_eqsam_v03d

#define DEBUG

  use pmc_constants,                    only: const
  use pmc_util,                         only: i_kind, dp, assert, assert_msg, &
                                              almost_equal, string_t, &
                                              to_string, warn_assert_msg, &
                                              die_msg
  use pmc_camp_core
  use pmc_camp_state
  use pmc_camp_solver_data
  use pmc_aero_rep_data
  use pmc_chem_spec_data
  use pmc_rxn_data
  use pmc_rxn_photolysis
  use pmc_property
#ifdef PMC_USE_JSON
  use json_module
#endif

  ! EBI Solver

  implicit none

  ! New-line character
  character(len=*), parameter :: new_line = char(10)
  ! EQSAM output file unit
  integer(kind=i_kind), parameter :: EQSAM_FILE_UNIT = 10
  ! CAMP-chem output file unit
  integer(kind=i_kind), parameter :: CAMP_FILE_UNIT = 12
  ! Number of timesteps to integrate over
  integer(kind=i_kind), parameter :: NUM_TIME_STEPS = 100
  ! Small number for minimum concentrations
  real(kind=dp), parameter :: SMALL_NUM = 1.0d-30
  ! Used to check availability of a solver
  type(camp_solver_data_t), pointer :: camp_solver_data

#ifdef DEBUG
  integer(kind=i_kind), parameter :: DEBUG_UNIT = 13

  open(unit=DEBUG_UNIT, file="out/debug_eqsam.txt", status="replace", action="write")
#endif

  camp_solver_data => camp_solver_data_t()

  if (.not.camp_solver_data%is_solver_available()) then
    write(*,*) "EQSAM test - no solver available - PASS"
  else if (run_eqsam_tests()) then
    write(*,*) "EQSAM tests - PASS"
  else
    write(*,*) "EQSAM tests - FAIL"
  end if

  deallocate(camp_solver_data)

#ifdef DEBUG
  close(DEBUG_UNIT)
#endif

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Run all EQSAM tests
  logical function run_eqsam_tests() result(passed)

    passed = run_standard_eqsam_test()

  end function run_eqsam_tests

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Run the eqsam mechanism under standard conditions using the original
  !! MONARCH eqsam code and the CAMP-chem version
  logical function run_standard_eqsam_test() result(passed)

    ! CAMP-chem core
    type(camp_core_t), pointer :: camp_core
    ! CAMP-chem state
    type(camp_state_t), pointer :: camp_state

    ! Computation timer variables
    real(kind=dp) :: comp_start, comp_end, comp_eqsam, comp_camp

    ! Temperature (K)
    real :: temperature = 272.5
    ! Pressure (atm)
    real :: pressure = 0.8

    ! CAMP-chem configuration file
    character(len=:), allocatable :: camp_input_file

    type(chem_spec_data_t), pointer :: chem_spec_data
    character(len=:), allocatable :: key, str_val, spec_name, phase_name, rep_name
    integer(kind=i_kind) :: i_mech, i_spec, n_spec, int_val, i_state_elem
    real(kind=dp) :: real_val, init_conc
    type(property_t), pointer :: prop_set
    class(aero_rep_data_t), pointer :: aero_rep
    integer(kind=i_kind), allocatable :: spec_ids(:)
    type(string_t), allocatable :: spec_names(:), unique_names(:)

    !!!!!!!!!!!!!!!!!!!!!!!
    !!! EQSAM variables !!!
    !!!!!!!!!!!!!!!!!!!!!!!

    ! Number of input variables
    integer, parameter :: eqsam_nca = 11
    ! Number of output variables
    integer, parameter :: eqsam_nco = 36
    ! Phase selection (1=metastable liquid; 2=solid)
    integer, parameter :: eqsam_iopt = 1
    ! Number of sets to solve
    integer, parameter :: eqsam_loop = 1
    ! Maximum number of sets
    integer, parameter :: eqsam_imax = 1
    ! Input array
    real, dimension(eqsam_imax, eqsam_nca) :: eqsam_yi
    ! Output array
    real, dimension(eqsam_imax, eqsam_nco) :: eqsam_yo


    passed = .true.

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! Initialize camp-chem !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call cpu_time(comp_start)
    camp_input_file = "config_eqsam.json"
    camp_core => camp_core_t(camp_input_file)

    ! Initialize the model
    call camp_core%initialize()

    ! Initialize the solver
    call camp_core%solver_initialize()

    ! Get a state variable
    camp_state => camp_core%new_state()

    call cpu_time(comp_end)
    write(*,*) "CAMP-chem initialization time: ", (comp_end-comp_start)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! Set the initial conditions !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! CAMP-chem environmental conditions
    camp_state%env_state%temp = temperature
    camp_state%env_state%pressure = pressure * const%air_std_press
    call camp_state%update_env_state()

    ! Get the chemical species data
    call assert(410085820, camp_core%get_chem_spec_data(chem_spec_data))

    ! CAMP-chem species concentrations
    key = "init conc"
    phase_name = "aqueous aerosol"
    rep_name = "single particle"
    call assert(522998221, camp_core%get_aero_rep(rep_name, aero_rep))
    spec_names = chem_spec_data%get_spec_names()
    do i_spec = 1, size(spec_names)
      call assert(929748071, chem_spec_data%get_property_set( &
              spec_names(i_spec)%string, prop_set))
      call assert(368608641, prop_set%get_real(key, real_val))
      call assert(514801839, chem_spec_data%get_phase( &
              spec_names(i_spec)%string, int_val))
      select case (int_val)
        case (CHEM_SPEC_GAS_PHASE)
          int_val = chem_spec_data%gas_state_id(spec_names(i_spec)%string)
          call assert_msg(891348329, int_val.gt.0, "Cannot find gas-phase species "// &
                  spec_names(i_spec)%string)
          camp_state%state_var(int_val) = real_val
        case (CHEM_SPEC_AERO_PHASE)
          unique_names = aero_rep%unique_names( &
                  phase_name = phase_name, spec_name = spec_names(i_spec)%string)
          call assert_msg(284639316, size(unique_names).eq.1, &
                  "Cannot find aerosol-phase species "//spec_names(i_spec)%string)
          i_state_elem = aero_rep%spec_state_id(unique_names(1)%string)
          camp_state%state_var(i_state_elem) = real_val
          deallocate(unique_names)
        case default
          call die_msg(837607961, "Unknown type for species "//spec_names(i_spec)%string)
      end select
    end do

    ! EQSAM environmental conditions
    eqsam_yi(1:eqsam_imax, 1) = temperature                             ! Temperature (K)
    eqsam_yi(1:eqsam_imax, 11) = pressure * 1.0e-2                      ! Pressure (hPa)

    ! EQSAM species concentrations

    ! RH
    eqsam_yi(1:eqsam_imax, 2) = water_conc_to_RH( &
            real(temperature, kind=dp), real(pressure, kind=dp), &
            get_eqsam_gas_conc(camp_core, camp_state, "H2O"))
    ! Ammonium
    eqsam_yi(1:eqsam_imax, 3) = &
            get_eqsam_gas_conc(camp_core, camp_state, "NH3") &
            + get_eqsam_aero_conc(aero_rep, phase_name, camp_state, "NH3_aq") &
            + get_eqsam_aero_conc(aero_rep, phase_name, camp_state, "NH4_p") &
            + 2*get_eqsam_aero_conc(aero_rep, phase_name, camp_state, "(NH4)2SO4") &
            + get_eqsam_aero_conc(aero_rep, phase_name, camp_state, "NH4NO3") &
            + get_eqsam_aero_conc(aero_rep, phase_name, camp_state, "NH4Cl") &
            + get_eqsam_aero_conc(aero_rep, phase_name, camp_state, "NH4HSO4") &
            + 3*get_eqsam_aero_conc(aero_rep, phase_name, camp_state, "(NH4)3H(SO4)2")

    ! Sulfate
    eqsam_yi(1:eqsam_imax, 4) = &
            get_eqsam_gas_conc(camp_core, camp_state, "H2SO4") &
            + get_eqsam_aero_conc(aero_rep, phase_name, camp_state, "H2SO4_aq") &
            + get_eqsam_aero_conc(aero_rep, phase_name, camp_state, "HSO4_m") &
            + get_eqsam_aero_conc(aero_rep, phase_name, camp_state, "SO4_mm") &
            + get_eqsam_aero_conc(aero_rep, phase_name, camp_state, "(Na)2SO4") &
            + get_eqsam_aero_conc(aero_rep, phase_name, camp_state, "(NH4)2SO4") &
            + get_eqsam_aero_conc(aero_rep, phase_name, camp_state, "NH4HSO4") &
            + get_eqsam_aero_conc(aero_rep, phase_name, camp_state, "NaHSO4") &
            + 2*get_eqsam_aero_conc(aero_rep, phase_name, camp_state, "(NH4)3H(SO4)2")
    ! Nitrate
    eqsam_yi(1:eqsam_imax, 5) = &
            get_eqsam_gas_conc(camp_core, camp_state, "HNO3") &
            + get_eqsam_aero_conc(aero_rep, phase_name, camp_state, "HNO3_aq") &
            + get_eqsam_aero_conc(aero_rep, phase_name, camp_state, "NO3_m") &
            + get_eqsam_aero_conc(aero_rep, phase_name, camp_state, "NaNO3") &
            + get_eqsam_aero_conc(aero_rep, phase_name, camp_state, "NH4NO3")
    ! Sodium
    eqsam_yi(1:eqsam_imax, 6) = &
            get_eqsam_aero_conc(aero_rep, phase_name, camp_state, "Na_p") &
            + get_eqsam_aero_conc(aero_rep, phase_name, camp_state, "NaCl") &
            + 2*get_eqsam_aero_conc(aero_rep, phase_name, camp_state, "(Na)2SO4") &
            + get_eqsam_aero_conc(aero_rep, phase_name, camp_state, "NaNO3") &
            + get_eqsam_aero_conc(aero_rep, phase_name, camp_state, "NaHSO4")
    ! Chloride
    eqsam_yi(1:eqsam_imax, 7) = &
            get_eqsam_gas_conc(camp_core, camp_state, "HCl") &
            + get_eqsam_aero_conc(aero_rep, phase_name, camp_state, "HCl_aq") &
            + get_eqsam_aero_conc(aero_rep, phase_name, camp_state, "Cl_m") &
            + get_eqsam_aero_conc(aero_rep, phase_name, camp_state, "NaCl") &
            + get_eqsam_aero_conc(aero_rep, phase_name, camp_state, "NH4Cl")
    ! Potassium
    eqsam_yi(1:eqsam_imax, 8) = &
            get_eqsam_aero_conc(aero_rep, phase_name, camp_state, "K_p")
    ! Calcium
    eqsam_yi(1:eqsam_imax, 9) = &
            get_eqsam_aero_conc(aero_rep, phase_name, camp_state, "Ca_pp")
    ! Magnesium
    eqsam_yi(1:eqsam_imax, 10) = &
            get_eqsam_aero_conc(aero_rep, phase_name, camp_state, "Mg_pp")


    !!!!!!!!!!!!!!!!!!!!!!!!!
    !!! Verify model data !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!

    ! TODO finish building full eqsam mechanism for test

    ! Find the mechanism
!    key = "eqsam"
!    call assert(369786347, camp_core%find_mechanism(key, i_mech))

    ! Make sure the right number of reactions is present
!    call assert_msg(218328842, camp_core%mechanism(i_mech)%size().eq.19, &
!            "Wrong number of camp-chem reactions: "// &
!            trim(to_string(camp_core%mechanism(i_mech)%size())))

    deallocate(spec_names)
    deallocate(camp_state)
    deallocate(camp_core)

  end function run_standard_eqsam_test

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Convert [H2O] (ppm) to RH (0-1)
  ! from https://www.vaisala.com/sites/default/fiels/documents/Humidity_Conversion_Formulas_B210973EN-F.pdf
  real(kind=dp) function water_conc_to_RH(temperature, pressure, water_conc)

    !> Temperature (K)
    real(kind=dp), intent(in) :: temperature
    !> Pressure (Pa)
    real(kind=dp), intent(in) :: pressure
    !> Water concentration (ppm)
    real(kind=dp), intent(in) :: water_conc

    real(kind=dp) :: v

    v = 1.0d0 - temperature / 647.096d0

    water_conc_to_RH = 2206.4d0 / 101325.0d0 * exp( &
            647.096d0 / temperature * ( &
                -7.85951783d0 * v + &
                1.84408259 * v**(1.5) + &
                -11.7866497 * v**(3) + &
                22.6807411 * v**(3.5) + &
                -15.9618719 * v**(4) + &
                1.80122502 * v**(7.5) &
            ))

  end function water_conc_to_RH

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get gas-phase species concentration for eqsam from camp-chem
  !! EQSAM units : umol/m^3; camp-chem units : ppm
  real(kind=dp) function get_eqsam_gas_conc(camp_core, camp_state, spec_name)

    !> CAMP-chem core
    type(camp_core_t), pointer, intent(in) :: camp_core
    !> CAMP-chem state
    type(camp_state_t), pointer, intent(in) :: camp_state
    !> Species name
    character(len=*), intent(in) :: spec_name

    type(chem_spec_data_t), pointer :: chem_spec_data
    character(len=:), allocatable :: spec_name_def
    integer(kind=i_kind) :: spec_id

    spec_name_def = spec_name

    call assert(285785622, camp_core%get_chem_spec_data(chem_spec_data))

    spec_id = chem_spec_data%gas_state_id(spec_name_def)
    call assert_msg(326761745, spec_id.gt.0, &
            "Error getting gas-phase concentration for "//trim(spec_name))
    get_eqsam_gas_conc = camp_state%state_var(spec_id) * &
            camp_state%env_state%pressure / const%univ_gas_const / &
            camp_state%env_state%temp

  end function get_eqsam_gas_conc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get aerosol-phase species concentration for eqsam from camp-chem
  !! EQSAM units: umol/m^3; camp-chem units: umol/m^3
  real(kind=dp) function get_eqsam_aero_conc(aero_rep, phase_name, camp_state, spec_name)

    !> Aerosol representation
    class(aero_rep_data_t), pointer, intent(in) :: aero_rep
    !> Phase name
    character(len=:), allocatable :: phase_name
    !> CAMP-chem state
    type(camp_state_t), pointer, intent(in) :: camp_state
    !> Species name
    character(len=*), intent(in) :: spec_name

    character(len=:), allocatable :: spec_name_def
    type(string_t), allocatable :: unique_names(:)
    integer(kind=i_kind) :: i_state_elem

    spec_name_def = spec_name

    unique_names = aero_rep%unique_names( &
            phase_name = phase_name, spec_name = spec_name_def)
    call assert_msg(370434317, size(unique_names).eq.1, &
                  "Cannot find aerosol-phase species "//spec_name)
    i_state_elem = aero_rep%spec_state_id(unique_names(1)%string)
    get_eqsam_aero_conc = camp_state%state_var(i_state_elem)

  end function get_eqsam_aero_conc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program pmc_test_eqsam_v03d
