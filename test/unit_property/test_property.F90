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

    integer(kind=i_kind) :: temp_int
    real(kind=dp) :: temp_real
    logical :: temp_logical
    character(len=:), allocatable :: temp_string, owner_name

    build_property_set_test = .false.
    owner_name = "build_property_set_test()"

    ! Build individual property links
    prop_set => property_t()
    key_name = "int_prop"
    call prop_set%put(key_name, 27, .false., owner_name)
    key_name = "bool_prop"
    call prop_set%put(key_name, .true., .false., owner_name)
    key_name = "real_prop"
    call prop_set%put(key_name, 12.32d0, .false., owner_name)
    key_name = "string_prop"
    call prop_set%put(key_name, "kd ks8*2alf  s", .false., owner_name)

    ! Build a property subset
    sub_set => property_t()
    key_name = "sub_int"
    call sub_set%put(key_name, 832, .false., owner_name)
    key_name = "sub_bool"
    call sub_set%put(key_name, .false., .false., owner_name)
    key_name = "sub_real"
    call sub_set%put(key_name, 593.d12, .false., owner_name)
    key_name = "sub_string"
    call sub_set%put(key_name, "nlsd98*)@ur soi87", .false., owner_name)
    key_name = "sub_prop"
    call prop_set%put(key_name, sub_set, .false., owner_name)
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
    call assert(848067695, temp_real.eq.12.32d0)
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
    call assert(265323857, temp_real.eq.593.d12)
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

    ! Make sure the 'allow duplicates' flag works
    deallocate(prop_set)
    prop_set => property_t()
    key_name = "my first key"
    call prop_set%put(key_name, 123.45d0, .true., owner_name)
    key_name = "my duplicate key"
    call prop_set%put(key_name, 482.41d0, .true., owner_name)
    key_name = "my duplicate key"
    call prop_set%put(key_name, 92.30412d0, .true., owner_name)

    sub_set => property_t()
    key_name = "sub key one"
    call sub_set%put(key_name, 512.42d0, .true., owner_name)
    key_name = "sub key duplicate"
    call sub_set%put(key_name, 142.3d0, .true., owner_name)
    key_name = "sub key duplicate"
    call sub_set%put(key_name, 5324.12d0, .true., owner_name)
    key_name = "my sub set"
    call prop_set%put(key_name, sub_set, .true., owner_name)
    deallocate(sub_set)

    call assert(147253076, prop_set%size().eq.4)
    call prop_set%iter_reset()
    call assert(240872022, prop_set%get_real(val=temp_real))
    call assert(739813176, prop_set%get_key(key_name))
    call assert(736454403, key_name.eq."my first key")
    call assert(467414246, temp_real.eq.123.45d0)
    call prop_set%iter_next()
    call assert(241611440, prop_set%get_real(val=temp_real))
    call assert(353929785, prop_set%get_key(key_name))
    call assert(183772881, key_name.eq."my duplicate key")
    call assert(361099626, temp_real.eq.482.41d0)
    call prop_set%iter_next()
    call assert(125934322, prop_set%get_real(val=temp_real))
    call assert(920785817, prop_set%get_key(key_name))
    call assert(750628913, key_name.eq."my duplicate key")
    call assert(580472009, temp_real.eq.92.30412d0)
    call prop_set%iter_next()
    call assert(666633991, prop_set%get_property_t(val=sub_set))
    call assert(326772478, prop_set%get_key(key_name))
    call assert(151351267, key_name.eq."my sub set")
    call sub_set%iter_reset()
    call assert(712942992, sub_set%get_real(val=temp_real))
    call assert(149898028, sub_set%get_key(key_name))
    call assert(309526318, key_name.eq."sub key one")
    call assert(311431852, temp_real.eq.512.42d0)
    call sub_set%iter_next()
    call assert(692790354, sub_set%get_real(val=temp_real))
    call assert(522633450, sub_set%get_key(key_name))
    call assert(970001296, key_name.eq."sub key duplicate")
    call assert(799844392, temp_real.eq.142.3d0)
    call sub_set%iter_next()
    call assert(347212239, sub_set%get_real(val=temp_real))
    call assert(242063735, sub_set%get_key(key_name))
    call assert(971906830, key_name.eq."sub key duplicate")
    call assert(519274677, temp_real.eq.5324.12d0)

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
    ! JSON string with duplicate keys
    character(kind=json_ck, len=*), parameter :: json_string_with_dups = &
            ' {'//new_line//&
            '   "my first key" : 123.45,'//new_line//&
            '   "my duplicate key" : 482.41,'//new_line//&
            '   "my duplicate key" : 92.30412,'//new_line//&
            '   "my sub set" : {'//new_line//&
            '     "sub key one" : 512.42,'//new_line//&
            '     "sub key duplicate" : 142.3,'//new_line//&
            '     "sub key duplicate" : 5324.12'//new_line//&
            '   }'//new_line//&
            ' }'
    ! JSON core
    type(json_core), pointer :: json
    type(json_value), pointer :: j_obj, child, next

    character(:), allocatable :: key_name
    character(kind=json_ck, len=:), allocatable :: unicode_key_name

    integer(kind=i_kind) :: temp_int
    real(kind=dp) :: temp_real
    logical :: temp_logical
    character(len=:), allocatable :: temp_string, owner_name

    owner_name = "load_property_set_test()"

    ! Set up the JSON core
    allocate(json)
    call json%parse(j_obj, json_string)

    ! Initialize a property set
    props => property_t()

    ! Load the property set with data in the JSON string
    ! passed as a JSON object so that all the data is loaded
    call props%load(json, j_obj, .true., owner_name)

    ! Make sure the property set contains the right number of elements
    call assert(874625445, props%size().eq.6)

    ! Make sure requests for all the top-level elements return pointers
    ! to the expected data
    key_name = "int_prop"
    call assert(775492256, props%get_int(key_name, temp_int))
    call assert(152703199, temp_int.eq.429)

    key_name = "real_prop"
    call assert(320245189, props%get_real(key_name, temp_real))
    call assert(432563534, temp_real.eq.124.3d7)

    key_name = "bool_prop"
    call assert(429914141, props%get_logical(key_name, temp_logical))
    call assert(207182985, temp_logical.eqv..true.)

    key_name = "string_prop"
    call assert(201918678, props%get_string(key_name, temp_string))
    call assert(596712272, temp_string.eq."nd*oi 2l3j(")

    key_name = "other_real"
    call assert(373981116, props%get_real(key_name, temp_real))
    call assert(768774710, temp_real.eq.129523.23d3)

    ! Make sure the subset of properties is accessible and correct
    key_name = "sub_props"
    call assert(198559905, props%get_property_t(key_name, sub_props))
    call assert(823254496, associated(sub_props))

    call assert(486612181, sub_props%size().eq.5)

    key_name = "other_real"
    call assert(935572841, sub_props%get_real(key_name, temp_real))
    call assert(365358036, temp_real.eq.123.4591d0)

    key_name = "sub_real"
    call assert(142626880, sub_props%get_real(key_name, temp_real))
    call assert(537420474, temp_real.eq.417d-7)

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
        call prop2%load(json, child, .false., owner_name)
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
    call assert(778874226, temp_real.eq.124.3d7)

    key_name = "other_real"
    call assert(838618319, prop2%get_real(key_name, temp_real))
    call assert(268403514, temp_real.eq.129523.23d3)

    deallocate(props)
    deallocate(prop2)
    deallocate(key_name)
    deallocate(unicode_key_name)
    deallocate(temp_string)
    call json%destroy(j_obj)
    deallocate(json)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! Load JSON with duplicate keys !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Set up the JSON core
    allocate(json)
    call json%parse(j_obj, json_string_with_dups)

    ! Initialize a property set
    props => property_t()

    ! Load the property set with data in the JSON string
    ! passed as a JSON object so that all the data is loaded
    ! and the allow duplicates flag set to true
    call props%load(json, j_obj, .true., owner_name, .true.)

    ! Make sure the property set contains the right number of elements
    call assert(252427177, props%size().eq.4)

    ! Check the contents of the property set
    call props%iter_reset()
    call assert(574493881, props%get_real(val=temp_real))
    call assert(404336977, props%get_key(key_name))
    call assert(234180073, key_name.eq."my first key")
    call assert(681547919, temp_real.eq.123.45d0)
    call props%iter_next()
    call assert(511391015, props%get_real(val=temp_real))
    call assert(406242511, props%get_key(key_name))
    call assert(236085607, key_name.eq."my duplicate key")
    call assert(683453453, temp_real.eq.482.41d0)
    call props%iter_next()
    call assert(230821300, props%get_real(val=temp_real))
    call assert(960664395, props%get_key(key_name))
    call assert(225556993, key_name.eq."my duplicate key")
    call assert(402883738, temp_real.eq.92.30412d0)
    call props%iter_next()
    call assert(850251584, props%get_property_t(val=sub_props))
    call assert(115144182, props%get_key(key_name))
    call assert(562512028, key_name.eq."my sub set")
    call sub_props%iter_reset()
    call assert(109879875, sub_props%get_real(val=temp_real))
    call assert(904731370, sub_props%get_key(key_name))
    call assert(452099217, key_name.eq."sub key one")
    call assert(281942313, temp_real.eq.512.42d0)
    call sub_props%iter_next()
    call assert(729310159, sub_props%get_real(val=temp_real))
    call assert(276678006, sub_props%get_key(key_name))
    call assert(171529502, key_name.eq."sub key duplicate")
    call assert(618897348, temp_real.eq.142.3d0)
    call sub_props%iter_next()
    call assert(783789945, sub_props%get_real(val=temp_real))
    call assert(331157792, sub_props%get_key(key_name))
    call assert(226009288, key_name.eq."sub key duplicate")
    call assert(955852383, temp_real.eq.5324.12d0)

    deallocate(props)
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

    integer(kind=i_kind) :: temp_int
    real(kind=dp) :: temp_real
    logical :: temp_logical
    character(len=:), allocatable :: temp_string, owner_name

    move_update_property_set_test = .false.
    owner_name = "move_update_property_set_test()"

    ! Create original property set
    orig_set => property_t()

    ! Build individual property links
    key_name = "int_prop"
    call orig_set%put(key_name, 27, .false., owner_name)
    key_name = "bool_prop"
    call orig_set%put(key_name, .true., .false., owner_name)
    key_name = "real_prop"
    call orig_set%put(key_name, 12.32d0, .false., owner_name)
    key_name = "string_prop"
    call orig_set%put(key_name, "kd ks8*2alf  s", .false., owner_name)

    ! Build a property subset
    sub_set => property_t()
    key_name = "sub_int"
    call sub_set%put(key_name, 832, .false., owner_name)
    key_name = "sub_bool"
    call sub_set%put(key_name, .false., .false., owner_name)
    key_name = "sub_real"
    call sub_set%put(key_name, 593.d12, .false., owner_name)
    key_name = "sub_string"
    call sub_set%put(key_name, "nlsd98*)@ur soi87", .false., owner_name)
    key_name = "sub_prop"
    call orig_set%put(key_name, sub_set, .false., owner_name)
    deallocate(sub_set)

    ! Create a property set to update with
    update_set => property_t()

    ! Build a few links to merge in
    key_name = "new_real_prop"
    call update_set%put(key_name, 135.23d3, .false., owner_name)
    sub_set => property_t()
    key_name = "new_sub_real_prop"
    call sub_set%put(key_name, 1.6784d-14, .false., owner_name)
    key_name = "sub_prop"
    call update_set%put(key_name, sub_set, .false., owner_name)
    deallocate(sub_set)

    sub_set => property_t()
    key_name = "new_sub_real_prop"
    call sub_set%put(key_name, 5239.60d1, .false., owner_name)
    key_name = "new_sub_prop"
    call update_set%put(key_name, sub_set, .false., owner_name)
    deallocate(sub_set)

    ! Update the original property data set
    call orig_set%update(update_set, owner_name)
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
    call assert(101932110, temp_real.eq.12.32d0)

    key_name = "string_prop"
    call assert(149242055, .not.orig_set%get_string(key_name, temp_string))
    call assert(326568800, dest_set%get_string(key_name, temp_string))
    call assert(438887145, temp_string.eq."kd ks8*2alf  s")

    key_name = "new_real_prop"
    call assert(609011286, dest_set%get_real(key_name, temp_real))
    call assert(609213628, temp_real.eq.135.23d3)

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
    call assert(260107161, temp_real.eq.593.d12)
    key_name = "sub_string"
    call assert(212312158, sub_set%get_string(key_name, temp_string))
    call assert(937376004, temp_string.eq."nlsd98*)@ur soi87")
    key_name = "new_sub_real_prop"
    call assert(646697207, sub_set%get_real(key_name, temp_real))
    call assert(888330104, temp_real.eq.1.6784d-14)

    ! New sub-property set
    key_name = "new_sub_prop"
    call assert(308288933, dest_set%get_property_t(key_name, sub_set))
    key_name = "new_sub_real_prop"
    call assert(417248505, sub_set%get_real(key_name, temp_real))
    call assert(772652986, temp_real.eq.5239.60d1)

    deallocate(orig_set)
    deallocate(dest_set)
    deallocate(key_name)
    deallocate(temp_string)

    move_update_property_set_test = .true.

 end function move_update_property_set_test

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program pmc_property_test
