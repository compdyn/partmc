! Copyright (C) 2017 Matthew Dawson
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_property_test program.

!> Unit tests for the pmc_property module.
program pmc_property_test 

  use pmc_util
  use pmc_property
#ifdef PMC_USE_JSON
  use json_module
#endif

  implicit none

  ! New-line character
  character(len=*), parameter :: new_line = char(10)

  if (run_pmc_property_tests()) then
    write(*,*) "Property tests - PASS"
  else
    write(*,*) "Property tests - FAIL"
  end if

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Run all pmc_property tests
  logical function run_pmc_property_tests() result(passed)

    passed = build_property_links_test()
    if (passed) passed = load_property_set_test()

  end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Build property links test
  logical function build_property_links_test()

    type(property_link_t), pointer :: a, b, c, d, e
    type(property_t), pointer :: prop_set, sub_set
    character(len=:), allocatable :: str_val, key_name

    build_property_links_test = .false.

    ! Build individual property links
    key_name = "int_prop"
    a => property_link_t(key_name, 27)
    key_name = "bool_prop"
    b => property_link_t(key_name, .true.)
    key_name = "real_prop"
    c => property_link_t(key_name, 12.32d0)
    key_name = "string_prop"
    d => property_link_t(key_name, "kd ks8*2alf  s")

    ! Build a property subset
    sub_set => property_t()
    key_name = "sub_int"
    call sub_set%put(key_name, 832)
    key_name = "sub_bool"
    call sub_set%put(key_name, .false.)
    key_name = "sub_real"
    call sub_set%put(key_name, 593.d12)
    key_name = "sub_string"
    call sub_set%put(key_name, "nlsd98*)@ur soi87")
    key_name = "sub_prop"
    e => property_link_t(key_name, sub_set)

    ! Check that the key names of the individual links are correct
    call assert(214833709, a%key().eq."int_prop")
    call assert(143948949, b%key().eq."bool_prop")
    call assert(841835324, c%key().eq."real_prop")
    call assert(420401760, d%key().eq."string_prop")
    call assert(903679690, e%key().eq."sub_prop")

    ! Check the values of the individual links
    call assert(903024217, a%value_int().eq.27)
    call assert(375307967, b%value_logical().eqv..true.)
    call assert(848067695, c%value_real().eq.12.32d0)
    call assert(137439884, d%value_string().eq."kd ks8*2alf  s")

    ! Check the values in the link that contains a subset of properties
    prop_set => e%value_property_t()
    call assert(273946937, associated(prop_set))
    deallocate(a)
    key_name = "sub_int"
    a => prop_set%get(key_name)
    call assert(918083289, a%value_int().eq.832)
    key_name = "sub_bool"
    a => prop_set%get(key_name)
    call assert(888112914, a%value_logical().eqv..false.)
    key_name = "sub_real"
    a => prop_set%get(key_name)
    call assert(265323857, a%value_real().eq.593.d12)
    key_name = "sub_string"
    a => prop_set%get(key_name)
    call assert(154911046, a%value_string().eq."nlsd98*)@ur soi87")
    
    ! Make sure that requests for a key that is not present in the set
    ! are returned a null pointer
    key_name = "bad_key"
    call assert(997072486, .not.associated(prop_set%get(key_name)))

    ! Make sure that requests for string values do not return pointers
    ! to the stored data
    str_val = a%value_string()
    str_val = "some other string"
    call assert(154911046, a%value_string().eq."nlsd98*)@ur soi87")

    deallocate(b)
    deallocate(c)
    deallocate(d)
    deallocate(e)
    deallocate(key_name)
    deallocate(str_val)

    build_property_links_test = .true.

  end function build_property_links_test

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Test loading properties from a json
  logical function load_property_set_test()
#ifdef PMC_USE_JSON
    ! Property set
    type(property_t), pointer :: props, sub_props, prop2
    ! Property value
    type(property_link_t), pointer :: prop_val
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

    ! Set up the JSON core
    allocate(json)
    call json%parse(j_obj, json_string)

    ! Initialize a property set
    allocate(props)
    props = property_t()

    ! Load the property set with data in the JSON string
    ! passed as a JSON object so that all the data is loaded
    call props%load(json, j_obj, .true.)

    ! Make sure the property set contains the right number of elements
    call assert(874625445, props%size().eq.6)

    ! Make sure requests for all the top-level elements return pointers
    ! to the expected data
    key_name = "real_prop"
    prop_val => props%get(key_name)
    call assert(521867099, associated(prop_val))
    call assert(894528147, prop_val%value_real().eq.124.3d7)
    
    key_name = "int_prop"
    prop_val => props%get(key_name)
    call assert(775492256, associated(prop_val))
    call assert(152703199, prop_val%value_int().eq.429)
    
    key_name = "bool_prop"
    prop_val => props%get(key_name)
    call assert(429914141, associated(prop_val))
    call assert(207182985, prop_val%value_logical().eqv..true.)
    
    key_name = "string_prop"
    prop_val => props%get(key_name)
    call assert(201918678, associated(prop_val))
    call assert(596712272, prop_val%value_string().eq."nd*oi 2l3j(")
    
    key_name = "other_real"
    prop_val => props%get(key_name)
    call assert(373981116, associated(prop_val))
    call assert(768774710, prop_val%value_real().eq.129523.23d3)
    
    ! Make sure the subset of properties is accessible and correct
    key_name = "sub_props"
    prop_val => props%get(key_name)
    call assert(198559905, associated(prop_val))
    sub_props => prop_val%value_property_t()
    call assert(823254496, associated(sub_props))

    call assert(486612181, sub_props%size().eq.5)

    key_name = "other_real"
    prop_val => sub_props%get(key_name)
    call assert(935572841, associated(prop_val))
    call assert(365358036, prop_val%value_real().eq.123.4591d0)
    
    key_name = "sub_real"
    prop_val => sub_props%get(key_name)
    call assert(142626880, associated(prop_val))
    call assert(537420474, prop_val%value_real().eq.417d-7)
    
    key_name = "sub_int"
    prop_val => sub_props%get(key_name)
    call assert(762057164, associated(prop_val))
    call assert(256850759, prop_val%value_int().eq.10294)
    
    key_name = "sub_bool"
    prop_val => sub_props%get(key_name)
    call assert(986693854, associated(prop_val))
    call assert(698954298, prop_val%value_logical().eqv..false.)
    
    key_name = "sub_str"
    prop_val => sub_props%get(key_name)
    call assert(476223142, associated(prop_val))
    call assert(135909334, prop_val%value_string().eq."kdm2i308vona aoid8$")
 
    ! Reload the JSON string to try passing individual key-value pairs to the
    ! property_t%load function
    call json%parse(j_obj, json_string)

    ! Initialize a new property set variable
    allocate(prop2)
    prop2 = property_t()

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
    prop_val => prop2%get(key_name)
    call assert(219188035, associated(prop_val))
    call assert(778874226, prop_val%value_real().eq.124.3d7)
    
    key_name = "other_real"
    prop_val => prop2%get(key_name)
    call assert(838618319, associated(prop_val))
    call assert(268403514, prop_val%value_real().eq.129523.23d3)
    
    deallocate(props)
    deallocate(prop2)
    deallocate(key_name)
    call json%destroy(j_obj)

#endif
    load_property_set_test = .true.
  end function load_property_set_test

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program pmc_property_test
