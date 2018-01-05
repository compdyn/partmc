! Copyright (C) 2017 Matt Dawson
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_test_rxn_arrhenius program

!> Test class for the rxn_test_t type, which extends the abstract
!! rxn_arrhenius_t type.
program pmc_test_rxn_arrhenius

  use pmc_util,                         only: i_kind, dp, assert, &
                                              almost_equal
  use pmc_phlex_state
  use pmc_chem_spec_data
  use pmc_rxn_arrhenius
  use pmc_mpi
#ifdef PMC_USE_JSON
  use json_module
#endif

  implicit none

  ! New-line character
  character(len=*), parameter :: new_line = char(10)

  call pmc_mpi_init()

  if (run_pmc_rxn_arrhenius_tests()) then
    write(*,*) "Arrhenius reaction tests - PASS"
  else
    write(*,*) "Arrhenius reaction tests - FAIL"
  end if

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Run all pmc_rxn_arrhenius tests
  logical function run_pmc_rxn_arrhenius_tests() result(passed)

    passed = build_rxn_arrhenius_test()

  end function run_pmc_rxn_arrhenius_tests

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Build rxn_arrhenius
  logical function build_rxn_arrhenius_test()

    type(rxn_data_ptr), allocatable :: rxns(:)
    
    build_rxn_arrhenius_test = .false.

    ! TODO finish

    build_rxn_arrhenius_test = .true.

  end function build_rxn_arrhenius_test()

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program pmc_test_rxn_arrhenius
