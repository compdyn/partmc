! Copyright (C) 2017 Matthew Dawson
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_property_test program.

!> Unit tests for the pmc_property module.
program pmc_property_test 

#ifdef PMC_USE_JSON
  use json_module
#endif
  use pmc_mpi
  use pmc_property
  use pmc_util

  implicit none

  ! New-line character
  character(len=*), parameter :: new_line = char(10)

  !> initialize mpi
  call pmc_mpi_init()

  if (run_pmc_property_tests()) then
    if (pmc_mpi_rank().eq.0) write(*,*) "Property tests - PASS"
  else
    if (pmc_mpi_rank().eq.0) write(*,*) "Property tests - FAIL"
  end if

  !> finalize mpi
  call pmc_mpi_finalize()

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Run all pmc_property tests
  logical function run_pmc_property_tests() result(passed)

    passed = build_property_set_test()
    if (passed) passed = load_property_set_test()
    if (passed) passed = move_update_property_set_test()

  end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Build property links test
  logical function build_property_set_test()

    type(property_t), pointer :: prop_set, sub_set
    character(len=:), allocatable :: key_name

    integer(kind=phlex_int) :: temp_int
    real(kind=phlex_real) :: temp_real
    logical :: temp_logical
    character(len=:), allocatable :: temp_string

    build_property_set_test = .false.

    ! Build individual property links
    prop_set => property_t()
    key_name = "int_prop"
    call prop_set%put(key_name, 27)
    key_name = "bool_prop"
    call prop_set%put(key_name, .true.)
    key_name = "real_prop"
    call prop_set%put(key_name, real(12.32e0, kind=phlex_real))
    key_name = "string_prop"
    call prop_set%put(key_name, "kd ks8*2alf  s")

    ! Build a property subset
    sub_set => property_t()
    key_name = "sub_int"
    call sub_set%put(key_name, 832)
    key_name = "sub_bool"
    call sub_set%put(key_name, .false.)
    key_name = "sub_real"
    call sub_set%put(key_name, real(593.e12, kind=phlex_real))
    key_name = "sub_string"
    call sub_set%put(key_name, "nlsd98*)@ur soi87")
    key_name = "sub_prop"
    call prop_set%put(key_name, sub_set)
    deallocate(sub_set)

    ! Check the values of the individual links
    key_name = "int_prop"
    call assert(372080023, prop_set%get_int(key_name, temp_int))
    call assert(903024217, temp_int.eq.27)
    key_name = "bool_prop"
    call assert(760156071, prop_set%get_logical(key_name, temp_logical))
    call assert(375307967, temp_logical.eqv..true.)
    key_name = "real_prop"
    call assert(872474416, prop_set%get_real(key_name, temp_real))
    call assert(848067695, almost_equal(temp_real, real(12.32e0, kind=phlex_real), real(1.0e-6, kind=phlex_real)))
    key_name = "string_prop"
    call assert(984792761, prop_set%get_string(key_name, temp_string))
    call assert(137439884, temp_string.eq."kd ks8*2alf  s")

    ! Check the values in the link that contains a subset of properties
    key_name = "sub_prop"
    call assert(273946937, prop_set%get_property_t(key_name, sub_set))

    key_name = "sub_int"
    call assert(729187696, sub_set%get_int(key_name, temp_int))
    call assert(918083289, temp_int.eq.832)
    key_name = "sub_bool"
    call assert(780308709, sub_set%get_logical(key_name, temp_logical))
    call assert(888112914, temp_logical.eqv..false.)
    key_name = "sub_real"
    call assert(319053476, sub_set%get_real(key_name, temp_real))
    call assert(265323857, almost_equal(temp_real, real(593.e12, kind=phlex_real), real(1.0e-6, kind=phlex_real)))
    key_name = "sub_string"
    call assert(433277355, sub_set%get_string(key_name, temp_string))
    call assert(154911046, temp_string.eq."nlsd98*)@ur soi87")
    
    ! Make sure that requests for a key that is not present in the set
    ! are returned a null pointer
    key_name = "bad_key"
    call assert(997072486, .not.prop_set%get_int(key_name, temp_int))

    ! Make sure that requests for string values do not return pointers
    ! to the stored data
    key_name = "sub_string"
    call assert(743362206, sub_set%get_string(key_name, temp_string))
    temp_string = "some other string"
    call assert(682164874, sub_set%get_string(key_name, temp_string))
    call assert(154911046, temp_string.eq."nlsd98*)@ur soi87")

    ! Deallocate variables
    deallocate(prop_set)
    deallocate(key_name)
    deallocate(temp_string)

    build_property_set_test = .true.

  end function build_property_set_test

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Test loading properties from a json
  logical function load_property_set_test()
#ifdef PMC_USE_JSON
    ! Property set
    type(property_t), pointer :: props, sub_props, prop2
    ! JSON string
    character(kind=json_ck, len=*), parameter :: json_string = &
            ' {'//new_line//&
            '   "real_prop" : 124.3e7,'//new_line//&
            '   "int_prop" : 429,'//new_line//&
            '   "bool_prop" : true,'//new_line//&
            '   "string_prop" : "nd*oi 2l3j(",'//new_line//&
            '   "sub_props" : {'//new_line//&
            '     "sub_int" : 10294,'//new_line//&
            '     "sub_real" : 417e-7,'//new_line//&
            '     "sub_bool" : false,'//new_line//&
            '     "sub_str" : "kdm2i308vona aoid8$",'//new_line//&
            '     "other_real" : 123.4591'//new_line//&
            '   },'//new_line//&
            '   "other_real" : 129523.23e3'//new_line//&
            ' }'
    ! JSON core
    type(json_core), pointer :: json
    type(json_value), pointer :: j_obj, child, next

    character(:), allocatable :: key_name
    character(kind=json_ck, len=:), allocatable :: unicode_key_name

    integer(kind=phlex_int) :: temp_int
    real(kind=phlex_real) :: temp_real
    logical :: temp_logical
    character(len=:), allocatable :: temp_string

    ! Set up the JSON core
    allocate(json)
    call json%parse(j_obj, json_string)

    ! Initialize a property set
    props => property_t()

    ! Load the property set with data in the JSON string
    ! passed as a JSON object so that all the data is loaded
    call props%load(json, j_obj, .true.)

    ! Make sure the property set contains the right number of elements
    call assert(874625445, props%size().eq.6)

    ! Make sure requests for all the top-level elements return pointers
    ! to the expected data
    key_name = "int_prop"
    call assert(775492256, props%get_int(key_name, temp_int)) 
    call assert(152703199, temp_int.eq.429)
    
    key_name = "real_prop"
    call assert(320245189, props%get_real(key_name, temp_real))
    call assert(432563534, almost_equal(temp_real, real(124.3e7, kind=phlex_real), real(1.0e-6, kind=phlex_real)))
    
    key_name = "bool_prop"
    call assert(429914141, props%get_logical(key_name, temp_logical))
    call assert(207182985, temp_logical.eqv..true.)
    
    key_name = "string_prop"
    call assert(201918678, props%get_string(key_name, temp_string)) 
    call assert(596712272, temp_string.eq."nd*oi 2l3j(")
    
    key_name = "other_real"
    call assert(373981116, props%get_real(key_name, temp_real))
    call assert(768774710, almost_equal(temp_real, real(129523.23e3, kind=phlex_real), real(1.0e-6, kind=phlex_real)))
    
    ! Make sure the subset of properties is accessible and correct
    key_name = "sub_props"
    call assert(198559905, props%get_property_t(key_name, sub_props))
    call assert(823254496, associated(sub_props))

    call assert(486612181, sub_props%size().eq.5)

    key_name = "other_real"
    call assert(935572841, sub_props%get_real(key_name, temp_real))
    call assert(365358036, almost_equal(temp_real, real(123.4591e0, kind=phlex_real), real(1.0e-6, kind=phlex_real)))
    
    key_name = "sub_real"
    call assert(142626880, sub_props%get_real(key_name, temp_real))
    call assert(537420474, almost_equal(temp_real, real(417e-7, kind=phlex_real), real(1.0e-6, kind=phlex_real)))
    
    key_name = "sub_int"
    call assert(762057164, sub_props%get_int(key_name, temp_int))
    call assert(256850759, temp_int.eq.10294)
    
    key_name = "sub_bool"
    call assert(986693854, sub_props%get_logical(key_name, temp_logical))
    call assert(698954298, temp_logical.eqv..false.)
    
    key_name = "sub_str"
    call assert(476223142, sub_props%get_string(key_name, temp_string))
    call assert(135909334, temp_string.eq."kdm2i308vona aoid8$")
 
    ! Reload the JSON string to try passing individual key-value pairs to the
    ! property_t%load function
    call json%destroy(j_obj)
    call json%parse(j_obj, json_string)

    ! Initialize a new property set variable
    prop2 => property_t()

    ! Only send the two real elements to property_t%load()
    next => null()
    call json%get_child(j_obj, child)
    do while (associated(child))
      call json%info(child, name=unicode_key_name)
      if (unicode_key_name.eq."real_prop" .or. unicode_key_name.eq."other_real") then
        call prop2%load(json, child, .false.)
      end if
      call json%get_next(child, next)
      child => next
    end do

    ! Make sure the property set contains the right number of elements
    call assert(823504432, prop2%size().eq.2)

    ! Make sure requests for all the top-level elements return pointers
    ! to the expected data
    key_name = "real_prop"
    call assert(219188035, prop2%get_real(key_name, temp_real))
    call assert(778874226, almost_equal(temp_real, real(124.3e7, kind=phlex_real), real(1.0e-6, kind=phlex_real)))
    
    key_name = "other_real"
    call assert(838618319, prop2%get_real(key_name, temp_real))
    call assert(268403514, almost_equal(temp_real, real(129523.23e3, kind=phlex_real), real(1.0e-6, kind=phlex_real)))
    
    deallocate(props)
    deallocate(prop2)
    deallocate(key_name)
    deallocate(unicode_key_name)
    deallocate(temp_string)
    call json%destroy(j_obj)
    deallocate(json)

#endif
    load_property_set_test = .true.
  end function load_property_set_test

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Move property set test
  logical function move_update_property_set_test()

    type(property_t), pointer :: orig_set, dest_set, sub_set, update_set
    character(len=:), allocatable :: key_name

    integer(kind=phlex_int) :: temp_int
    real(kind=phlex_real) :: temp_real
    logical :: temp_logical
    character(len=:), allocatable :: temp_string

    move_update_property_set_test = .false.

    ! Create original property set
    orig_set => property_t()

    ! Build individual property links
    key_name = "int_prop"
    call orig_set%put(key_name, 27)
    key_name = "bool_prop"
    call orig_set%put(key_name, .true.)
    key_name = "real_prop"
    call orig_set%put(key_name, real(12.32e0, kind=phlex_real))
    key_name = "string_prop"
    call orig_set%put(key_name, "kd ks8*2alf  s")

    ! Build a property subset
    sub_set => property_t()
    key_name = "sub_int"
    call sub_set%put(key_name, 832)
    key_name = "sub_bool"
    call sub_set%put(key_name, .false.)
    key_name = "sub_real"
    call sub_set%put(key_name, real(593.e12, kind=phlex_real))
    key_name = "sub_string"
    call sub_set%put(key_name, "nlsd98*)@ur soi87")
    key_name = "sub_prop"
    call orig_set%put(key_name, sub_set)
    deallocate(sub_set)

    ! Create a property set to update with
    update_set => property_t()

    ! Build a few links to merge in
    key_name = "new_real_prop"
    call update_set%put(key_name, real(135.23e3, kind=phlex_real))
    sub_set => property_t()
    key_name = "new_sub_real_prop"
    call sub_set%put(key_name, real(1.6784e-14, kind=phlex_real))
    key_name = "sub_prop"
    call update_set%put(key_name, sub_set)
    deallocate(sub_set)
    
    sub_set => property_t()
    key_name = "new_sub_real_prop"
    call sub_set%put(key_name, real(5239.60e1, kind=phlex_real))
    key_name = "new_sub_prop"
    call update_set%put(key_name, sub_set)
    deallocate(sub_set)

    ! Update the original property data set
    call orig_set%update(update_set)
    deallocate(update_set)

    ! Move the property set to the destination variable
    dest_set => property_t()
    call orig_set%move(dest_set)

    ! Check that the properties have moved
    key_name = "int_prop"
    call assert(443699157, .not.orig_set%get_int(key_name, temp_int))
    call assert(103385349, dest_set%get_int(key_name, temp_int)) 
    call assert(498178943, temp_int.eq.27)

    key_name = "bool_prop"
    call assert(492914636, .not.orig_set%get_logical(key_name, temp_logical))
    call assert(887708230, dest_set%get_logical(key_name, temp_logical))
    call assert(382501825, temp_logical.eqv..true.)

    key_name = "real_prop"
    call assert(494820170, .not.orig_set%get_real(key_name, temp_real))
    call assert(889613764, dest_set%get_real(key_name, temp_real))
    call assert(101932110, almost_equal(temp_real, real(12.32e0, kind=phlex_real), real(1.0e-6, kind=phlex_real)))

    key_name = "string_prop"
    call assert(149242055, .not.orig_set%get_string(key_name, temp_string))
    call assert(326568800, dest_set%get_string(key_name, temp_string))
    call assert(438887145, temp_string.eq."kd ks8*2alf  s")

    key_name = "new_real_prop"
    call assert(609011286, dest_set%get_real(key_name, temp_real))
    call assert(609213628, almost_equal(temp_real, real(135.23e3, kind=phlex_real), real(1.0e-6, kind=phlex_real)))

    ! New and original elements in original sub-property set
    key_name = "sub_prop"
    call assert(433622838, .not.orig_set%get_property_t(key_name, sub_set))
    call assert(323210027, dest_set%get_property_t(key_name, sub_set))

    key_name = "sub_int"
    call assert(904921627, sub_set%get_int(key_name, temp_int))
    call assert(652995221, temp_int.eq.832)
    key_name = "sub_bool"
    call assert(335159117, sub_set%get_logical(key_name, temp_logical))
    call assert(430264065, temp_logical.eqv..false.)
    key_name = "sub_real"
    call assert(500051714, sub_set%get_real(key_name, temp_real))
    call assert(260107161, almost_equal(temp_real, real(593.e12, kind=phlex_real), real(1.0e-6, kind=phlex_real)))
    key_name = "sub_string"
    call assert(212312158, sub_set%get_string(key_name, temp_string))
    call assert(937376004, temp_string.eq."nlsd98*)@ur soi87")
    key_name = "new_sub_real_prop"
    call assert(646697207, sub_set%get_real(key_name, temp_real))
    call assert(888330104, almost_equal(temp_real, real(1.6784e-14, kind=phlex_real), real(1.0e-6, kind=phlex_real)))

    ! New sub-property set
    key_name = "new_sub_prop"
    call assert(308288933, dest_set%get_property_t(key_name, sub_set))
    key_name = "new_sub_real_prop"
    call assert(417248505, sub_set%get_real(key_name, temp_real))
    call assert(772652986, almost_equal(temp_real, real(5239.60e1, kind=phlex_real), real(1.0e-6, kind=phlex_real)))

    deallocate(orig_set)
    deallocate(dest_set)
    deallocate(key_name)
    deallocate(temp_string)

    move_update_property_set_test = .true.

 end function move_update_property_set_test 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program pmc_property_test
