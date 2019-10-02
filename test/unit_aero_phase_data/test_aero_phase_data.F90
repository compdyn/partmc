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
  use pmc_mpi
#ifdef PMC_USE_MPI
  use mpi
#endif

  implicit none

  ! New-line character
  character(len=*), parameter :: new_line = char(10)

  !> initialize mpi
  call pmc_mpi_init()

  if (run_pmc_aero_phase_data_tests()) then
    if (pmc_mpi_rank().eq.0) write(*,*) "Aerosol phase data tests - PASS"
  else
    if (pmc_mpi_rank().eq.0) write(*,*) "Aerosol phase data tests - FAIL"
    stop 3
  end if

  !> finalize mpi
  call pmc_mpi_finalize()

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
    type(json_core), pointer :: json
    type(json_value), pointer :: j_obj, j_next

    integer(kind=i_kind) :: i_phase, i_spec
    type(property_t), pointer :: property_set
    character(len=:), allocatable :: key
    real(kind=dp) :: temp_real
    logical :: temp_logical
#ifdef PMC_USE_MPI
    type(aero_phase_data_ptr), allocatable :: aero_phase_passed_data_set(:)
    character, allocatable :: buffer(:)
    integer(kind=i_kind) :: pos, pack_size, i_prop
#endif
    allocate(json)
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

    call assert(814209333, aero_phase_data_set(1)%val%size().eq.3)
    call assert(863424812, aero_phase_data_set(2)%val%size().eq.3)
    call assert(358218407, aero_phase_data_set(3)%val%size().eq.2)

    call assert(278773971, aero_phase_data_set(1)%val%size().eq.3)
    call assert(608559165, aero_phase_data_set(2)%val%size().eq.3)
    call assert(438402261, aero_phase_data_set(3)%val%size().eq.2)

#ifdef PMC_USE_MPI
    pack_size = 0
    do i_phase = 1, 3
      pack_size = pack_size + aero_phase_data_set(i_phase)%val%pack_size(MPI_COMM_WORLD)
    end do
    allocate(buffer(pack_size))
    pos = 0
    do i_phase = 1, 3
      call aero_phase_data_set(i_phase)%val%bin_pack(buffer, pos, MPI_COMM_WORLD)
    end do
    allocate(aero_phase_passed_data_set(3))
    pos = 0
    do i_phase = 1, 3
      aero_phase_passed_data_set(i_phase)%val => aero_phase_data_t()
      call aero_phase_passed_data_set(i_phase)%val%bin_unpack(buffer, pos, MPI_COMM_WORLD)
    end do
    do i_phase = 1, 3
      call assert(165060871, &
        size(aero_phase_data_set(i_phase)%val%condensed_data_real).eq. &
        size(aero_phase_passed_data_set(i_phase)%val%condensed_data_real))
      do i_prop = 1, size(aero_phase_data_set(i_phase)%val%condensed_data_real)
        call assert(513997759, &
          aero_phase_data_set(i_phase)%val%condensed_data_real(i_prop).eq. &
          aero_phase_passed_data_set(i_phase)%val%condensed_data_real(i_prop))
      end do
      call assert(104315834, &
        size(aero_phase_data_set(i_phase)%val%condensed_data_int).eq. &
        size(aero_phase_passed_data_set(i_phase)%val%condensed_data_int))
      do i_prop = 1, size(aero_phase_data_set(i_phase)%val%condensed_data_int)
        call assert(834158929, &
          aero_phase_data_set(i_phase)%val%condensed_data_int(i_prop).eq. &
          aero_phase_passed_data_set(i_phase)%val%condensed_data_int(i_prop))
      end do
    end do
    deallocate(aero_phase_passed_data_set)
    deallocate(buffer)
#endif


#endif

    call j_file%destroy()
    call json%destroy(j_obj)
    deallocate(json)
    deallocate(aero_phase_data_set)
    deallocate(chem_spec_data)

    build_aero_phase_data_set_test = .true.

  end function build_aero_phase_data_set_test

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program pmc_test_aero_phase_data
