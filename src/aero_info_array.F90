! Copyright (C) 2007-2012, 2017 Nicole Riemer and Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_aero_info_array module.

!> The aero_info_array_t structure and assoicated subroutines.
module pmc_aero_info_array

  use pmc_aero_info
  use pmc_util
  use pmc_spec_file
  use pmc_mpi
#ifdef PMC_USE_MPI
  use mpi
#endif

  !> 1-D arrays of aero_info_t structure.
  !!
  !! This type implements a variable-length array of aero_info_t
  !! structures. To give a reasonable tradeoff between frequent
  !! re-allocs and memory usage, the length of an aero_info_array is
  !! generally a bit longer than the number of particles stored in
  !! it. When the array is full then a larger array is allocated with
  !! new extra space. As a balance between memory usage and frequency
  !! of re-allocs the length of the array is currently doubled when
  !! necessary and halved when possible.
  !!
  !! The true allocated length of the aero_info_array can be obtained
  !! by size(aero_info_array%%aero_info), while the number of used
  !! particle slots in it is given by aero_info_array_n_item().
  type aero_info_array_t
     !> Number of items in the array (not the same as the length of
     !> the allocated memory).
     integer :: n_item
     !> Array of aero_info_t structures.
     type(aero_info_t), allocatable :: aero_info(:)
  end type aero_info_array_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Return the current number of items.
  elemental integer function aero_info_array_n_item(aero_info_array)

    !> Aero info array.
    type(aero_info_array_t), intent(in) :: aero_info_array

    if (allocated(aero_info_array%aero_info)) then
       aero_info_array_n_item = aero_info_array%n_item
    else
       aero_info_array_n_item = 0
    end if

  end function aero_info_array_n_item

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Sets an aero_info_array to contain zero data.
  subroutine aero_info_array_zero(aero_info_array)

    !> Structure to reset.
    type(aero_info_array_t), intent(inout) :: aero_info_array

    aero_info_array%n_item = 0
    if (allocated(aero_info_array%aero_info)) then
       deallocate(aero_info_array%aero_info)
    end if
    allocate(aero_info_array%aero_info(0))

  end subroutine aero_info_array_zero

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Changes the given aero_info_array to exactly the given
  !> new_length.
  !!
  !! This function should not be called directly, but rather use
  !! aero_info_array_enlarge(), aero_info_array_enlarge_to()
  !! or aero_info_array_shrink().
  subroutine aero_info_array_realloc(aero_info_array, new_length)

    !> Array to reallocate (must already be allocated on entry).
    type(aero_info_array_t), intent(inout) :: aero_info_array
    !> New length of the array.
    integer, intent(in) :: new_length

    integer :: i
    type(aero_info_t), allocatable :: new_items(:)

    if (.not. allocated(aero_info_array%aero_info)) then
       allocate(aero_info_array%aero_info(new_length))
       aero_info_array%n_item = 0
       return
    end if

    call assert(955874877, new_length >= aero_info_array%n_item)
    allocate(new_items(new_length))
    do i = 1,aero_info_array%n_item
       new_items(i) = aero_info_array%aero_info(i)
    end do
    call move_alloc(new_items, aero_info_array%aero_info)

  end subroutine aero_info_array_realloc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Possibly enlarges the given array, ensuring that it is at least of size n.
  subroutine aero_info_array_enlarge_to(aero_info_array, n)

    !> Array to enlarge.
    type(aero_info_array_t), intent(inout) :: aero_info_array
    !> Minimum new size of array.
    integer, intent(in) :: n

    if (.not. allocated(aero_info_array%aero_info)) then
       call aero_info_array_realloc(aero_info_array, pow2_above(n))
       return
    end if

    if (n <= size(aero_info_array%aero_info)) return

    call aero_info_array_realloc(aero_info_array, pow2_above(n))

  end subroutine aero_info_array_enlarge_to

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Possibly shrinks the storage of the given array, ensuring that
  !> it can still store the allocated particles.
  subroutine aero_info_array_shrink(aero_info_array)

    !> Array to shrink.
    type(aero_info_array_t), intent(inout) :: aero_info_array

    integer :: new_length

    if (.not. allocated(aero_info_array%aero_info)) return

    new_length = pow2_above(aero_info_array%n_item)
    if (new_length < size(aero_info_array%aero_info)) then
       call aero_info_array_realloc(aero_info_array, new_length)
    end if

  end subroutine aero_info_array_shrink

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Adds the given aero_info to the end of the array.
  subroutine aero_info_array_add_aero_info(aero_info_array, &
       aero_info)

    !> Array to add to.
    type(aero_info_array_t), intent(inout) :: aero_info_array
    !> Aero_info to add.
    type(aero_info_t), intent(in) :: aero_info

    integer :: n

    n = aero_info_array_n_item(aero_info_array) + 1
    call aero_info_array_enlarge_to(aero_info_array, n)
    aero_info_array%aero_info(n) = aero_info
    aero_info_array%n_item = n

  end subroutine aero_info_array_add_aero_info

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Removes the aero_info at the given index.
  subroutine aero_info_array_remove_aero_info(aero_info_array, &
       index)

    !> Array to remove from.
    type(aero_info_array_t), intent(inout) :: aero_info_array
    !> Index of aero_info to remove.
    integer, intent(in) :: index

    call assert(578870706, allocated(aero_info_array%aero_info))
    call assert(213892348, index >= 1)
    call assert(953927392, index <= aero_info_array%n_item)
    if (index < aero_info_array%n_item) then
       ! shift last aero_info into empty slot to preserve dense packing
       aero_info_array%aero_info(index) &
            = aero_info_array%aero_info(aero_info_array%n_item)
    end if
    aero_info_array%n_item = aero_info_array%n_item - 1
    call aero_info_array_shrink(aero_info_array)

  end subroutine aero_info_array_remove_aero_info

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Adds \c aero_info_array_delta to the end of \c aero_info_array.
  subroutine aero_info_array_add(aero_info_array, &
       aero_info_array_delta)

    !> Array to add to.
    type(aero_info_array_t), intent(inout) :: aero_info_array
    !> Aero_info to add.
    type(aero_info_array_t), intent(in) :: aero_info_array_delta

    integer :: i, n, n_delta, n_new

    n = aero_info_array%n_item
    n_delta = aero_info_array_delta%n_item
    n_new = n + n_delta
    call aero_info_array_enlarge_to(aero_info_array, n_new)
    do i = 1,n_delta
       aero_info_array%aero_info(n + i) = aero_info_array_delta%aero_info(i)
    end do
    aero_info_array%n_item = n_new

  end subroutine aero_info_array_add

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Determines the number of bytes required to pack the given value.
  integer function pmc_mpi_pack_size_aia(val)

    !> Value to pack.
    type(aero_info_array_t), intent(in) :: val

    integer :: i, total_size

    total_size = 0
    total_size = total_size &
         + pmc_mpi_pack_size_integer(aero_info_array_n_item(val))
    do i = 1,aero_info_array_n_item(val)
       total_size = total_size &
            + pmc_mpi_pack_size_aero_info(val%aero_info(i))
    end do
    pmc_mpi_pack_size_aia = total_size

  end function pmc_mpi_pack_size_aia

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Packs the given value into the buffer, advancing position.
  subroutine pmc_mpi_pack_aero_info_array(buffer, position, val)

    !> Memory buffer.
    character, intent(inout) :: buffer(:)
    !> Current buffer position.
    integer, intent(inout) :: position
    !> Value to pack.
    type(aero_info_array_t), intent(in) :: val

#ifdef PMC_USE_MPI
    integer :: prev_position, i

    prev_position = position
    call pmc_mpi_pack_integer(buffer, position, aero_info_array_n_item(val))
    do i = 1,aero_info_array_n_item(val)
       call pmc_mpi_pack_aero_info(buffer, position, val%aero_info(i))
    end do
    call assert(732927292, &
         position - prev_position <= pmc_mpi_pack_size_aia(val))
#endif

  end subroutine pmc_mpi_pack_aero_info_array

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Unpacks the given value from the buffer, advancing position.
  subroutine pmc_mpi_unpack_aero_info_array(buffer, position, val)

    !> Memory buffer.
    character, intent(inout) :: buffer(:)
    !> Current buffer position.
    integer, intent(inout) :: position
    !> Value to pack.
    type(aero_info_array_t), intent(inout) :: val

#ifdef PMC_USE_MPI
    integer :: prev_position, i, n

    prev_position = position
    call pmc_mpi_unpack_integer(buffer, position, n)
    call aero_info_array_realloc(val, n)
    val%n_item = n
    do i = 1,n
       call pmc_mpi_unpack_aero_info(buffer, position, val%aero_info(i))
    end do
    call assert(262838429, &
         position - prev_position <= pmc_mpi_pack_size_aia(val))
#endif

  end subroutine pmc_mpi_unpack_aero_info_array

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_aero_info_array
