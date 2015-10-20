! Copyright (C) 2011-2012 Matthew West
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
     integer, allocatable :: entry(:)
  end type integer_varray_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Return the current number of entries.
  elemental integer function integer_varray_n_entry(integer_varray)

    !> Array.
    type(integer_varray_t), intent(in) :: integer_varray

    if (allocated(integer_varray%entry)) then
       integer_varray_n_entry = integer_varray%n_entry
    else
       integer_varray_n_entry = 0
    end if

  end function integer_varray_n_entry

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Changes the given integer_varray to exactly the given new_length.
  !!
  !! This function should not be called directly, but rather use
  !! integer_varray_enlarge(), integer_varray_shrink().
  subroutine integer_varray_realloc(integer_varray, new_length)

    !> Array to reallocate.
    type(integer_varray_t), intent(inout) :: integer_varray
    !> New length of the array.
    integer, intent(in) :: new_length

    integer, allocatable :: new_entries(:)

    if (.not. allocated(integer_varray%entry)) then
       allocate(integer_varray%entry(new_length))
       integer_varray%n_entry = 0
       return
    end if

    call assert(479324776, new_length >= integer_varray%n_entry)
    allocate(new_entries(new_length))
    new_entries(:integer_varray%n_entry) &
         = integer_varray%entry(1:integer_varray%n_entry)
    call move_alloc(new_entries, integer_varray%entry)

  end subroutine integer_varray_realloc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Resets an integer_varray to have zero entries.
  elemental subroutine integer_varray_zero(integer_varray)

    !> Structure to zero.
    type(integer_varray_t), intent(inout) :: integer_varray

    integer_varray%n_entry = 0
    if (allocated(integer_varray%entry)) then
       deallocate(integer_varray%entry)
    end if

  end subroutine integer_varray_zero

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Enlarges the given array so that it is at least of size n.
  subroutine integer_varray_enlarge(integer_varray, n)

    !> Array to enlarge.
    type(integer_varray_t), intent(inout) :: integer_varray
    !> Minimum new size of array.
    integer, intent(in) :: n

    if (.not. allocated(integer_varray%entry)) then
       call integer_varray_realloc(integer_varray, pow2_above(n))
       return
    end if

    if (n <= size(integer_varray%entry)) return

    call integer_varray_realloc(integer_varray, pow2_above(n))

  end subroutine integer_varray_enlarge

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Possibly shrinks the storage of the given array, ensuring that
  !> it can still store the used entries.
  subroutine integer_varray_shrink(integer_varray)

    !> Array to shrink.
    type(integer_varray_t), intent(inout) :: integer_varray

    integer :: length, new_length

    if (.not. allocated(integer_varray%entry)) return

    new_length = pow2_above(integer_varray%n_entry)
    if (new_length < size(integer_varray%entry)) then
       call integer_varray_realloc(integer_varray, new_length)
    end if

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
    call integer_varray_enlarge(integer_varray, n)
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

    call assert(302759108, allocated(integer_varray%entry))
    call assert(541032660, index >= 1)
    call assert(385739765, index <= integer_varray%n_entry)
    if (index < integer_varray%n_entry) then
       ! shift last entry into now-empty slot to preserve dense packing
       integer_varray%entry(index) &
            = integer_varray%entry(integer_varray%n_entry)
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

    logical :: is_allocated
    integer, allocatable :: tmp_entry(:)
    integer :: total_size

    is_allocated = allocated(val%entry)
    total_size = pmc_mpi_pack_size_logical(is_allocated)
    if (is_allocated) then
       tmp_entry = val%entry(1:val%n_entry)
       total_size = total_size &
            + pmc_mpi_pack_size_integer_array(tmp_entry)
    end if
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
    logical :: is_allocated
    integer :: prev_position
    integer, allocatable :: tmp_entry(:)

    prev_position = position
    is_allocated = allocated(val%entry)
    call pmc_mpi_pack_logical(buffer, position, is_allocated)
    if (is_allocated) then
       tmp_entry = val%entry(1:val%n_entry)
       call pmc_mpi_pack_integer_array(buffer, position, tmp_entry)
    end if
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
    integer :: prev_position
    logical :: is_allocated
    integer, allocatable :: tmp_entry(:)

    prev_position = position
    call pmc_mpi_unpack_logical(buffer, position, is_allocated)
    if (is_allocated) then
       call pmc_mpi_unpack_integer_array(buffer, position, tmp_entry)
       call integer_varray_realloc(val, size(tmp_entry))
       val%entry(1:size(tmp_entry)) = tmp_entry
    else
       if (allocated(val%entry)) then
          deallocate(val%entry)
       end if
    end if
    call assert(355866103, &
         position - prev_position <= pmc_mpi_pack_size_integer_varray(val))
#endif

  end subroutine pmc_mpi_unpack_integer_varray

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_integer_varray
