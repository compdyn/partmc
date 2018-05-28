! Copyright (C) 2017 Matthew Dawson
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_chem_spec_data_test program.

!> Unit tests for the pmc_chem_spec_data module.
program pmc_chem_spec_data_test

  use pmc_chem_spec_data
  use pmc_util
  use pmc_property
#ifdef PMC_USE_JSON
  use json_module
#endif

  implicit none

  ! New-line character
  character(len=*), parameter :: new_line = char(10)

  if (run_pmc_chem_spec_data_tests()) then
    write(*,*) "Chemical species data tests - PASS"
  else
    write(*,*) "Chemical species data tests - FAIL"
  end if

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Run all pmc_chem_spec_data tests
  logical function run_pmc_chem_spec_data_tests() result(passed)

    passed = build_chem_spec_data_test()

  end function run_pmc_chem_spec_data_tests

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Build chemical species data test
  logical function build_chem_spec_data_test() result(passed)

    type(chem_spec_data_t), pointer :: spec_data
    integer(kind=i_kind) :: i_spec
#ifdef PMC_USE_JSON
    character(len=:), allocatable :: json_string

    type(json_file) :: j_file
    type(json_core), pointer :: json
    type(json_value), pointer :: j_next, j_spec

    type(property_t), pointer :: spec_props, sub_props

    character(len=:), allocatable :: key_name
    integer(kind=i_kind) :: var_type

    integer(kind=i_kind) :: temp_int, spec_phase, spec_type
    real(kind=dp) :: temp_real
    logical :: temp_logical

    json_string = '{ "pmc-data" : [{'//new_line//&
            '  "name" : "my first species",'//new_line//&
            '  "type" : "CHEM_SPEC",'//new_line//&
            '  "MW" : 123.43,'//new_line//&
            '  "density" : 1.3e-8,'//new_line//&
            '  "UNIFAC params" : {'//new_line//&
            '    "CH2" : 2,'//new_line//&
            '    "CH" : 1,'//new_line//&
            '    "OH" : 3,'//new_line//&
            '    "CH3" : 1'//new_line//&
            '  },'//new_line//&
            '  "SIMPOL params" : {'//new_line//&
            '    "OH" : 3,'//new_line//&
            '    "CH" : 4'//new_line//&
            '  }'//new_line//&
            '}'//new_line//&

            ',{'//new_line//&
            '  "name" : "my second species",'//new_line//&
            '  "type" : "CHEM_SPEC",'//new_line//&
            '  "tracer type" : "PSSA",'//new_line//&
            '  "phase" : "AEROSOL",'//new_line//&
            '  "MW" : 75.298,'//new_line//&
            '  "density" : 4.2e-9,'//new_line//&
            '  "UNIFAC params" : {'//new_line//&
            '    "CH2" : 1,'//new_line//&
            '    "CH" : 1,'//new_line//&
            '    "COOH" : 1,'//new_line//&
            '    "CH3" : 1'//new_line//&
            '  },'//new_line//&
            '  "SIMPOL params" : {'//new_line//&
            '    "COOH" : 1,'//new_line//&
            '    "CH" : 4'//new_line//&
            '  }'//new_line//&
            '},'//new_line//&
            '{'//new_line//&
            '  "name" : "my first species",'//new_line//&
            '  "type" : "CHEM_SPEC",'//new_line//&
            '  "phase" : "GAS",'//new_line//&
            '  "HLC" : 12.3e1,'//new_line//&
            '  "MONARCH id" : 121'//new_line//&
            '},'//new_line//&
            '{'//new_line//&
            '  "name" : "my second species",'//new_line//&
            '  "type" : "CHEM_SPEC",'//new_line//&
            '  "phase" : "AEROSOL",'//new_line//&
            '  "tracer type" : "PSSA",'//new_line//&
            '  "HLC" : 13.7e1,'//new_line//&
            '  "MONARCH id" : 82'//new_line//&
            '}'//new_line

    ! Include extra species to test the expandable array functionality
    do i_spec=1, 125
      json_string = json_string//',{ "name" : "species '//trim(to_string(i_spec))//&
              '", "type" : "CHEM_SPEC" }'
    end do
    json_string = json_string//']}'

    ! Allocate a species data set with space for 15 species
    spec_data => chem_spec_data_t(15)

    ! The starting size of a species dataset should be zero
    call assert(146416775, spec_data%size().eq.0)

    ! Set up the JSON core
    allocate(json)
    call j_file%load_from_string(json_string)

    ! Get the pmc-data object
    call j_file%get('pmc-data(1)',j_spec)

    ! Load the json dataset
    call assert(256038991, associated(j_spec))
    do while (associated(j_spec))
      call spec_data%load(json, j_spec)
      j_next => j_spec
      call json%get_next(j_next, j_spec)
    end do
    
    ! deallocate json variables
    call j_file%destroy()
    deallocate(json_string)
    call json%destroy()
    deallocate(json)

    ! Initialize the dataset
    call spec_data%initialize()

    ! Check the final data set has the correct number of species
    call assert(280079471, spec_data%size().eq.127)

    ! Check the properties associated with the chemical species
    key_name = "my first species"
    call assert(968643104, spec_data%get_phase(key_name, spec_phase))
    call assert(839073838, spec_phase.eq.CHEM_SPEC_GAS_PHASE)
    call assert(745634015, spec_data%get_type(key_name, spec_type))
    call assert(689700990, spec_type.eq.CHEM_SPEC_VARIABLE)

    call assert(342495274, spec_data%get_property_set(key_name, spec_props))

    key_name = "MW"
    call assert(446925196, spec_props%get_real(key_name, temp_real))
    call assert(271956280, temp_real.eq.123.43d0)
    key_name = "density"
    call assert(487969890, spec_props%get_real(key_name, temp_real))
    call assert(600288235, temp_real.eq.1.3d-8)
    key_name = "HLC"
    call assert(315193631, spec_props%get_real(key_name, temp_real))
    call assert(427511976, temp_real.eq.12.3d1)
    key_name = "MONARCH id"
    call assert(257355072, spec_props%get_int(key_name, temp_int))
    call assert(369673417, temp_int.eq.121)

    key_name = "UNIFAC params"
    call assert(484611117, spec_props%get_property_t(key_name, sub_props))
    key_name = "CH2"
    call assert(305831133, sub_props%get_int(key_name, temp_int))
    call assert(584947609, temp_int.eq.2)
    key_name = "CH"
    call assert(523750277, sub_props%get_int(key_name, temp_int))
    call assert(353593373, temp_int.eq.1)
    key_name = "OH"
    call assert(465911718, sub_props%get_int(key_name, temp_int))
    call assert(295754814, temp_int.eq.3)
    key_name = "CH3"
    call assert(408073159, sub_props%get_int(key_name, temp_int))
    call assert(520391504, temp_int.eq.1)

    key_name = "SIMPOL params"
    call assert(348781361, spec_props%get_property_t(key_name, sub_props))
    key_name = "OH"
    call assert(908467552, sub_props%get_int(key_name, temp_int))
    call assert(685736396, temp_int.eq.3)
    key_name = "CH"
    call assert(798054741, sub_props%get_int(key_name, temp_int))
    call assert(910373086, temp_int.eq.4)

    key_name = "my second species"
    call assert(519369724, spec_data%get_phase(key_name, spec_phase))
    call assert(529754002, spec_phase.eq.CHEM_SPEC_AERO_PHASE)
    call assert(351292716, spec_data%get_type(key_name, spec_type))
    call assert(116127412, spec_type.eq.CHEM_SPEC_PSSA)
    
    call assert(451454846, spec_data%get_property_set(key_name, spec_props))

    key_name = "MW"
    call assert(307022846, spec_props%get_real(key_name, temp_real))
    call assert(354332791, temp_real.eq.75.298d0)
    key_name = "density"
    call assert(184175887, spec_props%get_real(key_name, temp_real))
    call assert(861444730, temp_real.eq.4.2d-9)
    key_name = "HLC"
    call assert(691287826, spec_props%get_real(key_name, temp_real))
    call assert(803606171, temp_real.eq.13.7d1)
    key_name = "MONARCH id"
    call assert(915924516, spec_props%get_int(key_name, temp_int))
    call assert(463292363, temp_int.eq.82)

    key_name = "UNIFAC params"
    call assert(575610708, spec_props%get_property_t(key_name, sub_props))
    key_name = "CH2"
    call assert(687929053, sub_props%get_int(key_name, temp_int))
    call assert(800247398, temp_int.eq.1)
    key_name = "CH"
    call assert(347615245, sub_props%get_int(key_name, temp_int))
    call assert(177458341, temp_int.eq.1)
    key_name = "COOH"
    call assert(354785086, sub_props%get_int(key_name, temp_int))
    call assert(802152932, temp_int.eq.1)
    key_name = "CH3"
    call assert(914471277, sub_props%get_int(key_name, temp_int))
    call assert(461839124, temp_int.eq.1)

    key_name = "SIMPOL params"
    call assert(291682220, spec_props%get_property_t(key_name, sub_props))
    key_name = "COOH"
    call assert(404000565, sub_props%get_int(key_name, temp_int))
    call assert(233843661, temp_int.eq.1)
    key_name = "CH"
    call assert(681211507, sub_props%get_int(key_name, temp_int))
    call assert(228579354, temp_int.eq.4)

    do i_spec=1, 125
      key_name = "species "//trim(to_string(i_spec))
      call assert(225364917, spec_data%get_phase(key_name, spec_phase))
      call assert(232390422, spec_phase.eq.CHEM_SPEC_GAS_PHASE)
    end do

    deallocate(key_name)
    deallocate(spec_data)

#endif
    passed = .true.

  end function build_chem_spec_data_test

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program pmc_chem_spec_data_test
