! Copyright (C) 2011 Nicole Riemer and Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_aero_sorted module.

!> The aero_sorted_t structure and assocated subroutines.
module pmc_aero_sorted

  use pmc_integer_array

  !> A sorted particle index set.
  type aero_sorted_t
     !> Array of integer arrays, one per bin.
     type(integer_array_t), allocatable, dimension(:) :: bin
  end type aero_sorted_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Allocate an empty structure.
  subroutine aero_sorted_allocate(aero_sorted)

    !> Structure to initialize.
    type(aero_sorted_t), intent(out) :: aero_sorted
    
    allocate(aero_sorted%bin(0))

  end subroutine aero_sorted_allocate
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Allocate a strcture with the given size.
  subroutine aero_sorted_allocate_size(aero_sorted, n_bin)

    !> Structure to initialize.
    type(aero_sorted_t), intent(out) :: aero_sorted
    !> Number of bins.
    integer, intent(in) :: n_bin
    
    allocate(aero_sorted%bin(n_bin))
    call integer_array_allocate(aero_sorted%bin)

  end subroutine aero_sorted_allocate_size
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Deallocates a previously allocated structure.
  subroutine aero_sorted_deallocate(aero_sorted)

    !> Structure to deallocate.
    type(aero_sorted_t), intent(inout) :: aero_sorted
    
    call integer_array_deallocate(aero_sorted%bin)
    deallocate(aero_sorted%bin)

  end subroutine aero_sorted_deallocate
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
end module pmc_aero_sorted
