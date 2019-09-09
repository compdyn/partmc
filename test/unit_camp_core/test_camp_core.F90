! Copyright (C) 2017 Matt Dawson
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_test_camp_core program

!> Test class for the camp_core_t type
program pmc_test_camp_core

#ifdef PMC_USE_JSON
  use json_module
#endif
  use pmc_chem_spec_data
  use pmc_mechanism_data
  use pmc_mpi
  use pmc_camp_core
  use pmc_util,                         only: i_kind, dp, assert, &
                                              almost_equal, string_t

  implicit none

  ! New-line character
  character(len=*), parameter :: new_line = char(10)

  !> initialize mpi
  call pmc_mpi_init()

  if (run_pmc_camp_core_tests()) then
    if (pmc_mpi_rank().eq.0) write(*,*) "CAMP-core tests - PASS"
  else
    if (pmc_mpi_rank().eq.0) write(*,*) "CAMP-core tests - FAIL"
  end if

  !> finalize mpi
  call pmc_mpi_finalize()

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Run all pmc_camp_core tests
  logical function run_pmc_camp_core_tests() result(passed)

    passed = load_camp_core_test()

  end function run_pmc_camp_core_tests

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Load model data from a test input file
  logical function load_camp_core_test()

    type(camp_core_t), pointer :: camp_core
    type(mechanism_data_t), pointer :: mechanism
    type(chem_spec_data_t), pointer :: chem_spec_data
    character(len=:), allocatable :: input_file_path
    character(len=:), allocatable :: key_name
#ifdef PMC_USE_MPI
    type(camp_core_t), pointer :: passed_core
    character, allocatable :: buffer(:), buffer_copy(:)
    integer(kind=i_kind) :: pos, pack_size, i_elem
#endif

    load_camp_core_test = .false.

    input_file_path = "test_run/unit_camp_core/test_mech_config.json"

    camp_core => camp_core_t(input_file_path)

    ! Get the chemical species data
    call assert(974584023, camp_core%get_chem_spec_data(chem_spec_data))

    ! Check the number of species in the model
    call assert(822520018, chem_spec_data%size().eq.7)

    ! Get the mechanism index
    key_name = "lunch mechanism"
    call assert(589233468, camp_core%get_mechanism(key_name, mechanism))

    ! Check the mechanism name
    call assert(636308667, mechanism%name().eq."lunch mechanism")

    ! Make sure all three reactions were loaded
    call assert(360948482, mechanism%size().eq.3)

#ifdef PMC_USE_MPI
    call camp_core%initialize()
    pack_size = camp_core%pack_size()
    allocate(buffer(pack_size))
    pos = 0
    call camp_core%bin_pack(buffer, pos)
    call assert(738553750, pos.eq.pack_size)
    passed_core => camp_core_t()
    pos = 0
    call passed_core%bin_unpack(buffer, pos)
    call assert(614253552, pos.eq.pack_size)
    call assert(893370028, passed_core%pack_size().eq.camp_core%pack_size())
    allocate(buffer_copy(pack_size))
    pos = 0
    call passed_core%bin_pack(buffer_copy, pos)
    do i_elem = 1, pack_size
      call assert_msg(303337336, buffer(i_elem).eq.buffer_copy(i_elem), &
              "Mismatch in element "//trim(to_string(i_elem)))
    end do
    deallocate(buffer)
    deallocate(buffer_copy)
    deallocate(passed_core)
#endif

    deallocate(camp_core)

    load_camp_core_test = .true.

  end function load_camp_core_test

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program pmc_test_camp_core
