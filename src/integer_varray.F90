! Copyright (C) 2011 Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_integer_varray module.

!> The integer_varray_t structure and assocated subroutines.
module pmc_integer_varray

  use pmc_util
  use pmc_mpi

  !> A variable-length 1D array of integers.
  !!
  !! The number of currently used entries in \c n_entry will generally
  !! be less than the allocated storage.
  type integer_varray_t
     !> Number of currently used entries.
     integer :: n_entry
     !> Array of integer values.
     integer, allocatable, dimension(:) :: entry
  end type integer_varray_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Allocates an empty structure.
  elemental subroutine integer_varray_allocate(integer_varray)

    !> Structure to initialize.
    type(integer_varray_t), intent(out) :: integer_varray

    integer_varray%n_entry = 0
    allocate(integer_varray%entry(0))

  end subroutine integer_varray_allocate
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Allocates a structure with the given size.
  elemental subroutine integer_varray_allocate_size(integer_varray, n_entry)

    !> Structure to initialize.
    type(integer_varray_t), intent(out) :: integer_varray
    !> Number of entires.
    integer, intent(in) :: n_entry

    integer_varray%n_entry = n_entry
    allocate(integer_varray%entry(n_entry))
    integer_varray%entry = 0

  end subroutine integer_varray_allocate_size
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Deallocates a previously allocated structure.
  elemental subroutine integer_varray_deallocate(integer_varray)

    !> Structure to deallocate.
    type(integer_varray_t), intent(inout) :: integer_varray
    
    deallocate(integer_varray%entry)

  end subroutine integer_varray_deallocate
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Changes the given integer_varray to exactly the given new_length.
  !!
  !! This function should not be called directly, but rather use
  !! integer_varray_enlarge(), integer_varray_enlarge_to() or
  !! integer_varray_shrink().
  subroutine integer_varray_reallocate(integer_varray, new_length)

    !> Array to reallocate.
    type(integer_varray_t), intent(inout) :: integer_varray
    !> New length of the array.
    integer, intent(in) :: new_length

    integer, dimension(integer_varray%n_entry) :: temp_array

    call assert(753399394, new_length >= integer_varray%n_entry)
    temp_array = integer_varray%entry(1:integer_varray%n_entry)
    deallocate(integer_varray%entry)
    allocate(integer_varray%entry(new_length))
    integer_varray%entry(1:integer_varray%n_entry) = temp_array
    
  end subroutine integer_varray_reallocate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Resets an integer_varray to have zero particles per bin.
  elemental subroutine integer_varray_zero(integer_varray)

    !> Structure to zero.
    type(integer_varray_t), intent(inout) :: integer_varray

    integer_varray%entry = 0
    integer_varray%n_entry = 0

  end subroutine integer_varray_zero
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Copies an integer_varray.
  subroutine integer_varray_copy(integer_varray_from, integer_varray_to)

    !> Structure to copy from.
    type(integer_varray_t), intent(in) :: integer_varray_from
    !> Structure to copy to.
    type(integer_varray_t), intent(inout) :: integer_varray_to
    
    call integer_varray_deallocate(integer_varray_to)
    call integer_varray_allocate_size(integer_varray_to, &
         integer_varray_from%n_entry)
    integer_varray_to%entry(1:integer_varray_from%n_entry) &
         = integer_varray_from%entry(1:integer_varray_from%n_entry)

  end subroutine integer_varray_copy
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Enlarges the given integer_varray by at least one element.
  !!
  !! Currently this at least doubles the length.
  subroutine integer_varray_enlarge(integer_varray)

    !> Array to enlarge.
    type(integer_varray_t), intent(inout) :: integer_varray

    integer :: length, new_length

    length = size(integer_varray%entry)
    new_length = max(length * 2, length + 1)
    call integer_varray_reallocate(integer_varray, new_length)
    
  end subroutine integer_varray_enlarge

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Enlarges the given array so that it is at least of size n.
  subroutine integer_varray_enlarge_to(integer_varray, n)

    !> Array to enlarge.
    type(integer_varray_t), intent(inout) :: integer_varray
    !> Minimum new size of array.
    integer, intent(in) :: n

    do while (size(integer_varray%entry) < n)
       call integer_varray_enlarge(integer_varray)
    end do

  end subroutine integer_varray_enlarge_to

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Possibly shrinks the storage of the given array, ensuring that
  !> it can still store the used entries.
  subroutine integer_varray_shrink(integer_varray)

    !> Array to shrink.
    type(integer_varray_t), intent(inout) :: integer_varray

    integer :: length, new_length

    length = size(integer_varray%entry)
    new_length = length / 2
    do while ((integer_varray%n_entry <= new_length) .and. (length > 0))
       call integer_varray_reallocate(integer_varray, new_length)
       length = size(integer_varray%entry)
       new_length = length / 2
    end do

  end subroutine integer_varray_shrink

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Adds the given number to the end of the array.
  subroutine integer_varray_append(integer_varray, val)

    !> Array to add to.
    type(integer_varray_t), intent(inout) :: integer_varray
    !> Value to add.
    integer, intent(in) :: val

    integer :: n

    n = integer_varray%n_entry + 1
    call integer_varray_enlarge_to(integer_varray, n)
    integer_varray%entry(n) = val
    integer_varray%n_entry = n

  end subroutine integer_varray_append

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Removes the entry at the given index, repacking values to
  !> maintain contiguous data.
  subroutine integer_varray_remove_entry(integer_varray, index)

    !> Array to remove from.
    type(integer_varray_t), intent(inout) :: integer_varray
    !> Index of entry to remove.
    integer, intent(in) :: index

    call assert(541032660, index >= 1)
    call assert(385739765, index <= integer_varray%n_entry)
    if (index < integer_varray%n_entry) then
       ! shift last entry into now-empty slot to preserve dense packing
       integer_varray%entry(index) = integer_varray%entry(integer_varray%n_entry)
    end if
    ! clear now-unused last entry for safety
    integer_varray%entry(integer_varray%n_entry) = 0
    integer_varray%n_entry = integer_varray%n_entry - 1
    call integer_varray_shrink(integer_varray)

  end subroutine integer_varray_remove_entry

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Determines the number of bytes required to pack the given value.
  integer function pmc_mpi_pack_size_integer_varray(val)

    !> Value to pack.
    type(integer_varray_t), intent(in) :: val

    integer :: i, total_size

    total_size = 0
    total_size = total_size &
         + pmc_mpi_pack_size_integer_array(val%entry(1:val%n_entry))
    pmc_mpi_pack_size_integer_varray = total_size

  end function pmc_mpi_pack_size_integer_varray

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Packs the given value into the buffer, advancing position.
  subroutine pmc_mpi_pack_integer_varray(buffer, position, val)

    !> Memory buffer.
    character, intent(inout) :: buffer(:)
    !> Current buffer position.
    integer, intent(inout) :: position
    !> Value to pack.
    type(integer_varray_t), intent(in) :: val

#ifdef PMC_USE_MPI
    integer :: prev_position, i

    prev_position = position
    call pmc_mpi_pack_integer_array(buffer, position, &
         val%entry(1:val%n_entry))
    call assert(230655880, &
         position - prev_position <= pmc_mpi_pack_size_integer_varray(val))
#endif

  end subroutine pmc_mpi_pack_integer_varray

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Unpacks the given value from the buffer, advancing position.
  subroutine pmc_mpi_unpack_integer_varray(buffer, position, val)

    !> Memory buffer.
    character, intent(inout) :: buffer(:)
    !> Current buffer position.
    integer, intent(inout) :: position
    !> Value to pack.
    type(integer_varray_t), intent(inout) :: val

#ifdef PMC_USE_MPI
    integer :: prev_position, i, n
    ! FIXME: should switch to allocatable arrays in pmc_mpi_unpack_*()
    integer, pointer, dimension(:) :: tmp_entry

    prev_position = position
    allocate(tmp_entry(0))
    call pmc_mpi_unpack_integer_array(buffer, position, tmp_entry)
    call integer_varray_deallocate(val)
    call integer_varray_allocate_size(val, size(tmp_entry))
    val%entry = tmp_entry
    call assert(355866103, &
         position - prev_position <= pmc_mpi_pack_size_integer_varray(val))
#endif

  end subroutine pmc_mpi_unpack_integer_varray

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
end module pmc_integer_varray
