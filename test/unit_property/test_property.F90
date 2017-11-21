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

    run_pmc_property_tests = build_property_links_test

  end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Build property links test
  logical function build_property_links_test()

    class(property_link_t), pointer :: a, b, c, d

    a = property_link_t("int_prop", 27)
    b = property_link_t("bool_prop", .true.)
    c = property_link_t("real_prop", 12.32d0)
    d = property_link_t("string_prop", "hi there")

    assert(214833709, a%key.eq."int_prop")
    assert(143948949, b%key.eq."bool_prop")
    assert(841835324, c%key.eq."real_prop")
    assert(420401760, d%key.eq."string_prop")

    assert(328255399, a%value.eq.27)
    assert(375307967, b%value.eq..true.)
    assert(848067695, c%value.eq.12.32d0)
    assert(137439884, d%value.eq."hi there")

  end function build_property_links_test

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program pmc_property_test
