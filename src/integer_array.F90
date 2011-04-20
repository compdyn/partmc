! Copyright (C) 2011 Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_integer_array module.

!> The integer_array_t structure and assocated subroutines.
module pmc_integer_array

  use pmc_util

  !> A variable-length 1D array of integers.
  !!
  !! The number of currently used entries in \c n_entry will generally
  !! be less than the allocated storage.
  type integer_array_t
     !> Number of currently used entries.
     integer :: n_entry
     !> Array of integer values.
     integer, allocatable, dimension(:) :: entry
  end type integer_array_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Allocates an empty structure.
  elemental subroutine integer_array_allocate(integer_array)

    !> Structure to initialize.
    type(integer_array_t), intent(out) :: integer_array

    integer_array%n_entry = 0
    allocate(integer_array%entry(0))

  end subroutine integer_array_allocate
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Allocates a structure with the given size.
  elemental subroutine integer_array_allocate_size(integer_array, n_entry)

    !> Structure to initialize.
    type(integer_array_t), intent(out) :: integer_array
    !> Number of entires.
    integer, intent(in) :: n_entry

    integer_array%n_entry = n_entry
    allocate(integer_array%entry(n_entry))
    integer_array%entry = 0

  end subroutine integer_array_allocate_size
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Deallocates a previously allocated structure.
  elemental subroutine integer_array_deallocate(integer_array)

    !> Structure to deallocate.
    type(integer_array_t), intent(inout) :: integer_array
    
    deallocate(integer_array%entry)

  end subroutine integer_array_deallocate
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Changes the given integer_array to exactly the given new_length.
  !!
  !! This function should not be called directly, but rather use
  !! integer_array_enlarge(), integer_array_enlarge_to() or
  !! integer_array_shrink().
  subroutine integer_array_reallocate(integer_array, new_length)

    !> Array to reallocate.
    type(integer_array_t), intent(inout) :: integer_array
    !> New length of the array.
    integer, intent(in) :: new_length

    integer, dimension(integer_array%n_entry) :: temp_array

    call assert(753399394, new_length >= integer_array%n_entry)
    if (new_length <= size(integer_array%entry)) return
    temp_array = integer_array%entry(1:integer_array%n_entry)
    deallocate(integer_array%entry)
    allocate(integer_array%entry(new_length))
    integer_array%entry(1:integer_array%n_entry) = temp_array
    
  end subroutine integer_array_reallocate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Enlarges the given integer_array by at least one element.
  !!
  !! Currently this at least doubles the length.
  subroutine integer_array_enlarge(integer_array)

    !> Array to enlarge.
    type(integer_array_t), intent(inout) :: integer_array

    integer :: length, new_length

    length = size(integer_array%entry)
    new_length = max(length * 2, length + 1)
    call integer_array_reallocate(integer_array, new_length)
    
  end subroutine integer_array_enlarge

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Enlarges the given array so that it is at least of size n.
  subroutine integer_array_enlarge_to(integer_array, n)

    !> Array to enlarge.
    type(integer_array_t), intent(inout) :: integer_array
    !> Minimum new size of array.
    integer, intent(in) :: n

    do while (size(integer_array%entry) < n)
       call integer_array_enlarge(integer_array)
    end do

  end subroutine integer_array_enlarge_to

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Possibly shrinks the storage of the given array, ensuring that
  !> it can still store the used entries.
  subroutine integer_array_shrink(integer_array)

    !> Array to shrink.
    type(integer_array_t), intent(inout) :: integer_array

    integer :: length, new_length

    length = size(integer_array%entry)
    new_length = length / 2
    do while ((integer_array%n_entry <= new_length) .and. (length > 0))
       call integer_array_reallocate(integer_array, new_length)
       length = size(integer_array%entry)
       new_length = length / 2
    end do

  end subroutine integer_array_shrink

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Adds the given number to the end of the array.
  subroutine integer_array_append(integer_array, val)

    !> Array to add to.
    type(integer_array_t), intent(inout) :: integer_array
    !> Value to add.
    integer, intent(in) :: val

    integer :: n

    n = integer_array%n_entry + 1
    call integer_array_enlarge_to(integer_array, n)
    integer_array%entry(n) = val
    integer_array%n_entry = n

  end subroutine integer_array_append

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Removes the entry at the given index, repacking values to
  !> maintain contiguous data.
  subroutine integer_array_remove_entry(integer_array, index)

    !> Array to remove from.
    type(integer_array_t), intent(inout) :: integer_array
    !> Index of entry to remove.
    integer, intent(in) :: index

    call assert(541032660, index >= 1)
    call assert(385739765, index <= integer_array%n_entry)
    if (index < integer_array%n_entry) then
       ! shift last entry into now-empty slot to preserve dense packing
       integer_array%entry(index) = integer_array%entry(integer_array%n_entry)
    end if
    ! clear now-unused last entry for safety
    integer_array%entry(integer_array%n_entry) = 0
    integer_array%n_entry = integer_array%n_entry - 1
    call integer_array_shrink(integer_array)

  end subroutine integer_array_remove_entry

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
end module pmc_integer_array
