! Copyright (C) 2007-2009 Nicole Riemer and Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_aero_info_array module.

!> The aero_info_array_t structure and assoicated subroutines.
module pmc_aero_info_array

  use pmc_aero_info
  use pmc_util
  use pmc_spec_read
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
  !! particle slots in it is given by aero_info_array%%n_item. It must
  !! be that aero_info_array%%n_item is less than or equal to
  !! size(aero_info_array%%aero_info).
  type aero_info_array_t
     !> Number of items in the array (not the same as the length of
     !> the allocated memory).
     integer :: n_item
     !> Array of aero_info_t structures.
     type(aero_info_t), pointer :: aero_info(:)
  end type aero_info_array_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Allocates and initializes.
  subroutine aero_info_array_alloc(aero_info_array, n_item)

    !> Result.
    type(aero_info_array_t), intent(inout) :: aero_info_array
    !> Number of items.
    integer, intent(in) :: n_item

    integer :: i

    aero_info_array%n_item = n_item
    allocate(aero_info_array%aero_info(n_item))
    do i = 1,n_item
       call aero_info_alloc(aero_info_array%aero_info(i))
    end do

  end subroutine aero_info_array_alloc
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Deallocates.
  subroutine aero_info_array_free(aero_info_array)

    !> Structure to deallocate.
    type(aero_info_array_t), intent(inout) :: aero_info_array

    integer :: i
    
    do i = 1,aero_info_array%n_item
       call aero_info_free(aero_info_array%aero_info(i))
    end do
    deallocate(aero_info_array%aero_info)

  end subroutine aero_info_array_free
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Copies aero_info_array_from to aero_info_array_to, both
  !> of which must already be allocated.
  subroutine aero_info_array_copy(aero_info_array_from, &
       aero_info_array_to)

    !> Origin structure.
    type(aero_info_array_t), intent(in) :: aero_info_array_from
    !> Destination structure.
    type(aero_info_array_t), intent(inout) :: aero_info_array_to

    integer :: i
    
    call aero_info_array_free(aero_info_array_to)
    call aero_info_array_alloc(aero_info_array_to, &
         aero_info_array_from%n_item)
    do i = 1,aero_info_array_from%n_item
       call aero_info_copy(aero_info_array_from%aero_info(i), &
            aero_info_array_to%aero_info(i))
    end do

  end subroutine aero_info_array_copy
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Resets an aero_info_array to contain zero particles.
  subroutine aero_info_array_zero(aero_info_array)

    !> Structure to reset.
    type(aero_info_array_t), intent(inout) :: aero_info_array

    call aero_info_array_free(aero_info_array)
    allocate(aero_info_array%aero_info(0))
    aero_info_array%n_item = 0

  end subroutine aero_info_array_zero
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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

    integer :: n_item, i
    type(aero_info_t), pointer :: new_particles(:)

    n_item = aero_info_array%n_item
    call assert(372938429, new_length >= n_item)
    allocate(new_particles(new_length))
    do i = 1,aero_info_array%n_item
       call aero_info_copy(aero_info_array%aero_info(i), &
            new_particles(i))
       call aero_info_free(aero_info_array%aero_info(i))
    end do
    deallocate(aero_info_array%aero_info)
    aero_info_array%aero_info => new_particles
    
  end subroutine aero_info_array_realloc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Enlarges the given aero_info_array by at least one element
  !!
  !! Currently this doubles the length.
  subroutine aero_info_array_enlarge(aero_info_array)

    !> Array to enlarge.
    type(aero_info_array_t), intent(inout) :: aero_info_array

    integer :: length, new_length

    length = size(aero_info_array%aero_info)
    new_length = max(length * 2, length + 1)
    call aero_info_array_realloc(aero_info_array, new_length)
    
  end subroutine aero_info_array_enlarge

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Enlarges the given array so that it is at least of size n.
  subroutine aero_info_array_enlarge_to(aero_info_array, n)

    !> Array to enlarge.
    type(aero_info_array_t), intent(inout) :: aero_info_array
    !> Minimum new size of array.
    integer, intent(in) :: n

    do while (size(aero_info_array%aero_info) < n)
       call aero_info_array_enlarge(aero_info_array)
    end do

  end subroutine aero_info_array_enlarge_to

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Possibly shrinks the storage of the given array, ensuring that
  !> it can still store the allocated particles.
  subroutine aero_info_array_shrink(aero_info_array)

    !> Array to shrink.
    type(aero_info_array_t), intent(inout) :: aero_info_array

    integer :: n_item, length, new_length

    n_item = aero_info_array%n_item
    length = size(aero_info_array%aero_info)
    new_length = length / 2
    do while ((n_item <= new_length) .and. (length > 0))
       call aero_info_array_realloc(aero_info_array, new_length)
       length = size(aero_info_array%aero_info)
       new_length = length / 2
    end do

  end subroutine aero_info_array_shrink

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Adds the given aero_info to the end of the array.
  subroutine aero_info_array_add_aero_info(aero_info_array, &
       aero_info)

    !> Array to add to.
    type(aero_info_array_t), intent(inout) :: aero_info_array
    !> Aero_info to add.
    type(aero_info_t), intent(in) :: aero_info

    integer :: n

    n = aero_info_array%n_item + 1
    call aero_info_array_enlarge_to(aero_info_array, n)
    call aero_info_alloc(aero_info_array%aero_info(n))
    call aero_info_copy(aero_info, aero_info_array%aero_info(n))
    aero_info_array%n_item = aero_info_array%n_item + 1

  end subroutine aero_info_array_add_aero_info

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Removes the aero_info at the given index.
  subroutine aero_info_array_remove_aero_info(aero_info_array, &
       index)

    !> Array to remove from.
    type(aero_info_array_t), intent(inout) :: aero_info_array
    !> Index of aero_info to remove.
    integer, intent(in) :: index

    call assert(213892348, index >= 1)
    call assert(953927392, index <= aero_info_array%n_item)
    call aero_info_free(aero_info_array%aero_info(index))
    if (index < aero_info_array%n_item) then
       ! shift last aero_info into empty slot to preserve dense packing
       call aero_info_copy( &
            aero_info_array%aero_info(aero_info_array%n_item), &
            aero_info_array%aero_info(index))
       call aero_info_free( &
            aero_info_array%aero_info(aero_info_array%n_item))
    end if
    aero_info_array%n_item = aero_info_array%n_item - 1
    call aero_info_array_shrink(aero_info_array)

  end subroutine aero_info_array_remove_aero_info

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Determines the number of bytes required to pack the given value.
  integer function pmc_mpi_pack_size_aia(val)

    !> Value to pack.
    type(aero_info_array_t), intent(in) :: val

    integer :: i, total_size

    total_size = 0
    total_size = total_size + pmc_mpi_pack_size_integer(val%n_item)
    do i = 1,val%n_item
       total_size = total_size &
            + pmc_mpi_pack_size_aero_info(val%aero_info(i))
    end do
    pmc_mpi_pack_size_aia = total_size

  end function pmc_mpi_pack_size_aia

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
    call pmc_mpi_pack_integer(buffer, position, val%n_item)
    do i = 1,val%n_item
       call pmc_mpi_pack_aero_info(buffer, position, val%aero_info(i))
    end do
    call assert(732927292, &
         position - prev_position == pmc_mpi_pack_size_apa(val))
#endif

  end subroutine pmc_mpi_pack_aero_info_array

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Unpacks the given value from the buffer, advancing position.
  subroutine pmc_mpi_unpack_aero_info_array(buffer, position, val)

    !> Memory buffer.
    character, intent(inout) :: buffer(:)
    !> Current buffer position.
    integer, intent(inout) :: position
    !> Value to pack.
    type(aero_info_array_t), intent(out) :: val

#ifdef PMC_USE_MPI
    integer :: prev_position, i

    prev_position = position
    call pmc_mpi_unpack_integer(buffer, position, val%n_item)
    allocate(val%aero_info(val%n_item))
    do i = 1,val%n_item
       call pmc_mpi_unpack_aero_info(buffer, position, val%aero_info(i))
    end do
    call assert(262838429, &
         position - prev_position == pmc_mpi_pack_size_apa(val))
#endif

  end subroutine pmc_mpi_unpack_aero_info_array

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
end module pmc_aero_info_array
