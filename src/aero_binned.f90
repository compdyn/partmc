! Copyright (C) 2005-2007 Nicole Riemer and Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.
!
! Functions that deal with the binned aerosol distributions.

module mod_aero_binned

  type aero_binned_t
     real*8, pointer :: v(:)            ! len n_bin, volume per bin (m^3)
     real*8, pointer :: vs(:,:)         ! n_bin x n_spec, vol per bin&spec (m^3)
     integer, pointer :: n(:)           ! len n_bin, number per bin (1)
  end type aero_binned_t

contains
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine alloc_aero_binned(n_bin, n_spec, aero_binned)

    ! Allocates a aero_binned.

    integer, intent(in) :: n_bin        ! number of bins
    integer, intent(in) :: n_spec       ! number of species
    type(aero_binned_t), intent(out) :: aero_binned ! bin distribution

    allocate(aero_binned%v(n_bin))
    allocate(aero_binned%vs(n_bin, n_spec))
    allocate(aero_binned%n(n_bin))

  end subroutine alloc_aero_binned

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine free_aero_binned(aero_binned)

    ! Frees all memory.

    type(aero_binned_t), intent(inout) :: aero_binned ! aero_binned to free

    deallocate(aero_binned%v)
    deallocate(aero_binned%vs)
    deallocate(aero_binned%n)

  end subroutine free_aero_binned

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
end module mod_aero_binned
