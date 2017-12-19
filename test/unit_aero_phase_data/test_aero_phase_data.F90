! Copyright (C) 2017 Matt Dawson
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_test_aero_phase_data program

!> Test class for the aero_phase_data_t type
program pmc_test_aero_phase_data

  use pmc_util,                         only: i_kind, dp, assert, &
                                              almost_equal
  use pmc_property
  use pmc_aero_phase_data
  use pmc_chem_spec_data
#ifdef PMC_USE_JSON
  use json_module
#endif

  implicit none

  ! New-line character
  character(len=*), parameter :: new_line = char(10)

  if (run_pmc_aero_phase_data_tests()) then
    write(*,*) "Aerosol phase data tests - PASS"
  else
    write(*,*) "Aerosol phase data tests - FAIL"
  end if

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Run all pmc_aero_phase_data tests
  logical function run_pmc_aero_phase_data_tests() result(passed)

    passed = build_aero_phase_data_set_test()

  end function run_pmc_aero_phase_data_tests

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Build aero_phase_data set
  logical function build_aero_phase_data_set_test()

    type(aero_phase_data_t), pointer :: aero_phase_data
    type(aero_phase_data_ptr), allocatable :: aero_phase_data_set(:)
    type(chem_spec_data_t), pointer :: chem_spec_data
#ifdef PMC_USE_JSON
    type(json_file) :: j_file
    type(json_core), target :: json
    type(json_value), pointer :: j_obj, j_next

    integer(kind=i_kind) :: i_phase, i_spec
    type(property_t), pointer :: property_set
    character(len=:), allocatable :: key
    real(kind=dp) :: temp_real
    logical :: temp_logical

    call j_file%initialize()
    call j_file%get_core(json)
    call j_file%load_file(filename = &
            'test_run/unit_aero_phase_data/test_aero_phase_data.json')
    call j_file%get('pmc-data(1)',j_obj)

    build_aero_phase_data_set_test = .false.

    allocate(aero_phase_data_set(3))
    chem_spec_data => chem_spec_data_t()
    i_phase = 1
    i_spec = 1
    do while (associated(j_obj))
      if (i_phase.le.3) then
        aero_phase_data_set(i_phase)%val => aero_phase_data_t()
        call aero_phase_data_set(i_phase)%val%load(json, j_obj)
        i_phase = i_phase + 1
      else
        call chem_spec_data%load(json, j_obj)
        i_spec = i_spec + 1
      end if
      call json%get_next(j_obj, j_next)
      j_obj => j_next
    end do
 
    call assert(680635018, i_phase.eq.4)
    call assert(964927420, i_spec.eq.8)

    call assert(679181779, aero_phase_data_set(1)%val%name().eq."my test phase one")
    call assert(793291680, aero_phase_data_set(2)%val%name().eq."my test phase two")
    call assert(905610025, aero_phase_data_set(3)%val%name().eq."my last test phase")

    property_set => aero_phase_data_set(1)%val%get_property_set()
    key = "some property"
    call assert(313789343, property_set%get_real(key, temp_real))
    call assert(364910356, almost_equal(temp_real, real(12.2, kind=dp)))

    property_set => aero_phase_data_set(2)%val%get_property_set()
    key = "some other property"
    call assert(867210283, property_set%get_logical(key, temp_logical))
    call assert(979528628, .not.temp_logical)

    property_set => aero_phase_data_set(3)%val%get_property_set()
    key = "some property"
    call assert(191846974, property_set%get_real(key, temp_real))
    call assert(304165319, almost_equal(temp_real, real(13.75, kind=dp)))

    call aero_phase_data_set(1)%val%initialize(chem_spec_data)
    call aero_phase_data_set(2)%val%initialize(chem_spec_data)
    call aero_phase_data_set(3)%val%initialize(chem_spec_data)

    call assert(814209333, aero_phase_data_set(1)%val%state_size().eq.3)
    call assert(863424812, aero_phase_data_set(2)%val%state_size().eq.3)
    call assert(358218407, aero_phase_data_set(3)%val%state_size().eq.2)

    call assert(278773971, aero_phase_data_set(1)%val%size().eq.3)
    call assert(608559165, aero_phase_data_set(2)%val%size().eq.3)
    call assert(438402261, aero_phase_data_set(3)%val%size().eq.2)

    key = "species a"
    call assert(780621603, aero_phase_data_set(1)%val%state_id(key).ne.0)
    key = "species b"
    call assert(610464699, aero_phase_data_set(1)%val%state_id(key).ne.0)
    key = "species c"
    call assert(722783044, aero_phase_data_set(1)%val%state_id(key).ne.0)
    key = "species d"
    call assert(552626140, aero_phase_data_set(1)%val%state_id(key).eq.0)

    key = "species c"
    call assert(947419734, aero_phase_data_set(2)%val%state_id(key).ne.0)
    key = "species d"
    call assert(777262830, aero_phase_data_set(2)%val%state_id(key).ne.0)
    key = "species e"
    call assert(889581175, aero_phase_data_set(2)%val%state_id(key).ne.0)
    key = "species f"
    call assert(154473773, aero_phase_data_set(2)%val%state_id(key).eq.0)

    key = "species b"
    call assert(549267367, aero_phase_data_set(3)%val%state_id(key).ne.0)
    key = "species e"
    call assert(379110463, aero_phase_data_set(3)%val%state_id(key).ne.0)
    key = "species a"
    call assert(544003060, aero_phase_data_set(3)%val%state_id(key).eq.0)
#endif  
    build_aero_phase_data_set_test = .true.

  end function build_aero_phase_data_set_test

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program pmc_test_aero_phase_data
