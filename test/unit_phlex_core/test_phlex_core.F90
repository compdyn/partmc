! Copyright (C) 2017 Matt Dawson
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_test_phlex_core program

!> Test class for the phlex_core_t type
program pmc_test_phlex_core

  use pmc_util,                         only: i_kind, dp, assert, &
                                              almost_equal, string_t
  use pmc_phlex_core
#ifdef PMC_USE_JSON
  use json_module
#endif

  implicit none

  ! New-line character
  character(len=*), parameter :: new_line = char(10)

  if (run_pmc_phlex_core_tests()) then
    write(*,*) "Model data tests - PASS"
  else
    write(*,*) "Model data tests - FAIL"
  end if

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Run all pmc_phlex_core tests
  logical function run_pmc_phlex_core_tests() result(passed)

    passed = load_phlex_core_test()

  end function run_pmc_phlex_core_tests

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Load model data from a test input file
  logical function load_phlex_core_test()

    type(phlex_core_t), pointer :: phlex_core
    character(len=:), allocatable :: input_file_path
    integer :: i_mech
    character(len=:), allocatable :: key_name

    load_phlex_core_test = .false.

    input_file_path = 'test_run/unit_phlex_core/test_mech_config.json'

    phlex_core => phlex_core_t(input_file_path)

    ! Check the number of species in the model
    call assert(822520018, phlex_core%chem_spec_data%size().eq.7)

    ! Make sure one mechanism has been loaded
    call assert(533327223, size(phlex_core%mechanism).eq.1)

    ! Get the mechanism index
    key_name = "lunch mechanism"
    call assert(589233468, phlex_core%find_mechanism(key_name, i_mech))

    ! Check the mechanism name
    call assert(636308667, phlex_core%mechanism(i_mech)%name().eq."lunch mechanism")

    ! Make sure all three reactions were loaded
    call assert(360948482, phlex_core%mechanism(i_mech)%size().eq.3)

    load_phlex_core_test = .true.

  end function load_phlex_core_test

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program pmc_test_phlex_core
