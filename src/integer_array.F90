! Copyright (C) 2011 Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_integer_array module.

!> The integer_array_t structure and assocated subroutines.
module pmc_integer_array

  !> An array of integers packaged as a derived type.
  type integer_array_t
     !> Array of integer values.
     integer, allocatable, dimension(:) :: entry
  end type integer_array_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Allocates an empty structure.
  subroutine integer_array_allocate(integer_array)

    !> Structure to initialize.
    type(integer_array_t), intent(out) :: integer_array
    
    allocate(integer_array%entry(0))

  end subroutine integer_array_allocate
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Allocates a structure with the given size.
  subroutine integer_array_allocate_size(integer_array, n_entry)

    !> Structure to initialize.
    type(integer_array_t), intent(out) :: integer_array
    !> Number of entires.
    integer, intent(in) :: n_entry

    allocate(integer_array%entry(n_entry))

  end subroutine integer_array_allocate_size
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Deallocates a previously allocated structure.
  subroutine integer_array_deallocate(integer_array)

    !> Structure to deallocate.
    type(integer_array_t), intent(inout) :: integer_array
    
    deallocate(integer_array%entry)

  end subroutine integer_array_deallocate
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
end module pmc_integer_array
