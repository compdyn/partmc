! Copyright (C) 2017 Matt Dawson
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_test_model_data program

!> Test class for the model_data_t type
program pmc_test_model_data

  use pmc_util,                         only: i_kind, dp, assert, &
                                              almost_equal, string_t
  use pmc_model_data
  use pmc_model_state
#ifdef PMC_USE_JSON
  use json_module
#endif

  implicit none

  ! New-line character
  character(len=*), parameter :: new_line = char(10)

  if (run_pmc_model_data_tests()) then
    write(*,*) "Model data tests - PASS"
  else
    write(*,*) "Model data tests - FAIL"
  end if

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Run all pmc_model_data tests
  logical function run_pmc_model_data_tests() result(passed)

    passed = load_model_data_test()

  end function run_pmc_model_data_tests

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Load model data from a test input file
  logical function load_model_data_test()

    type(model_data_t), pointer :: model_data
    character(len=:), allocatable :: input_file_path
    integer :: i_mech
    character(len=:), allocatable :: key_name

    load_model_data_test = .false.

    input_file_path = 'test_run/unit_model_data/test_mech_config.json'

    model_data => model_data_t(input_file_path)

    ! Check the number of species in the model
    call assert(822520018, model_data%chem_spec_data%size().eq.7)

    ! Make sure one mechanism has been loaded
    call assert(533327223, size(model_data%mechanism).eq.1)

    ! Get the mechanism index
    key_name = "lunch mechanism"
    call assert(589233468, model_data%find_mechanism(key_name, i_mech))

    ! Check the mechanism name
    call assert(636308667, model_data%mechanism(i_mech)%name().eq."lunch mechanism")

    ! Make sure all three reactions were loaded
    call assert(360948482, model_data%mechanism(i_mech)%size().eq.3)

    load_model_data_test = .true.

  end function load_model_data_test

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program pmc_test_model_data
