! Copyright (C) 2017 Matt Dawson
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_test_eqsam program

!> Test for the eqsam inorganic module from MONARCH. This program runs the
!! MONARCH eqsam code and the Phlex-chem version and compares the output.
program pmc_test_eqsam

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
  ! EQSAM output file unit
  integer(kind=i_kind), parameter :: EQSAM_FILE_UNIT = 10
  ! Phlex-chem output file unit
  integer(kind=i_kind), parameter :: PHLEX_FILE_UNIT = 12
  ! Number of timesteps to integrate over
  integer(kind=i_kind), parameter :: NUM_TIME_STEPS = 100
  ! Small number for minimum concentrations
  real(kind=dp), parameter :: SMALL_NUM = 1.0d-30
  ! Used to check availability of a solver  
  type(phlex_solver_data_t), pointer :: phlex_solver_data

#ifdef DEBUG
  integer(kind=i_kind), parameter :: DEBUG_UNIT = 13
   
  open(unit=DEBUG_UNIT, file="out/debug_eqsam.txt", status="replace", action="write")
#endif

  phlex_solver_data => phlex_solver_data_t()

  if (.not.phlex_solver_data%is_solver_available()) then
    write(*,*) "EQSAM test - no solver available - PASS"
  else if (run_eqsam_tests()) then
    write(*,*) "EQSAM tests - PASS"
  else
    write(*,*) "EQSAM tests - FAIL"
  end if

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
  !! MONARCH eqsam code and the Phlex-chem version
  logical function run_standard_eqsam_test() result(passed)

    ! Phlex-chem core
    type(phlex_core_t), pointer :: phlex_core
    ! Phlex-chem state
    type(phlex_state_t), pointer :: phlex_state, phlex_state_comp
    ! Phlex-chem species names
    type(string_t), allocatable :: phlex_spec_names(:)

    ! Computation timer variables
    real(kind=dp) :: comp_start, comp_end, comp_eqsam, comp_phlex

    ! Phlex-chem configuration file
    character(len=:), allocatable :: phlex_input_file

    run_standard_eqsam_test = .true.

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! Initialize phlex-chem !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    call cpu_time(comp_start)
    phlex_input_file = "config_eqsam.json"
    phlex_core => phlex_core_t(phlex_input_file)
    
    ! Initialize the model
    call phlex_core%initialize()

    ! Initialize the solver
    call phlex_core%solver_initialize()
    call cpu_time(comp_end)
    write(*,*) "Phlex-chem initialization time: ", (comp_end-comp_start)



  end function run standard_eqsam_test

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program pmc_test_eqsam
