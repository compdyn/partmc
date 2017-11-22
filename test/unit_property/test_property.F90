! Copyright (C) 2017 Matthew Dawson
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_property_test program.

!> Unit tests for the pmc_property module.
program pmc_property_test 

  use pmc_util
  use pmc_property

  implicit none

  if (run_pmc_property_tests()) then
    write(*,*) "Property tests - PASS"
  else
    write(*,*) "Property tests - FAIL"
  end if

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Run all pmc_property tests
  logical function run_pmc_property_tests()

    run_pmc_property_tests = build_property_links_test()

  end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Build property links test
  logical function build_property_links_test()

    type(property_link_t), pointer :: a, b, c, d, e
    type(property_t), pointer :: prop_set, sub_set
    character(len=:), allocatable :: str_val

    build_property_links_test = .false.

    a => property_link_t("int_prop", 27)
    b => property_link_t("bool_prop", .true.)
    c => property_link_t("real_prop", 12.32d0)
    str_val = "string_prop"
    d => property_link_t(str_val, "kd ks8*2alf  s")
    deallocate(str_val)

    sub_set => property_t()
    call sub_set%put("sub_int", 832)
    call sub_set%put("sub_bool", .false.)
    call sub_set%put("sub_real", 593.d12)
    call sub_set%put("sub_string", "nlsd98*)@ur soi87")

    e => property_link_t("sub_prop", sub_set)

    call assert(214833709, a%key().eq."int_prop")
    call assert(143948949, b%key().eq."bool_prop")
    call assert(841835324, c%key().eq."real_prop")
    call assert(420401760, d%key().eq."string_prop")
    call assert(903679690, e%key().eq."sub_prop")

    call assert(903024217, a%value_int().eq.27)
    call assert(375307967, b%value_logical().eqv..true.)
    call assert(848067695, c%value_real().eq.12.32d0)
    call assert(137439884, d%value_string().eq."kd ks8*2alf  s")

    prop_set => e%value_property_t()
    call assert(273946937, associated(prop_set))
    deallocate(a)
    a => prop_set%get("sub_int")
    call assert(918083289, a%value_int().eq.832)
    a => prop_set%get("sub_bool")
    call assert(888112914, a%value_logical().eqv..false.)
    a => prop_set%get("sub_real")
    call assert(265323857, a%value_real().eq.593.d12)
    a => prop_set%get("sub_string")
    call assert(154911046, a%value_string().eq."nlsd98*)@ur soi87")
    call assert(997072486, .not.associated(prop_set%get("bad_key")))

    str_val = a%value_string()
    str_val = "some other string"
    call assert(154911046, a%value_string().eq."nlsd98*)@ur soi87")

    deallocate(b)
    deallocate(c)
    deallocate(d)
    deallocate(e)

    build_property_links_test = .true.

  end function build_property_links_test

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program pmc_property_test
