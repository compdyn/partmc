! Copyright (C) 2005-2011 Nicole Riemer and Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_aero_particle_array module.

!> The aero_particle_array_t structure and assoicated subroutines.
module pmc_aero_particle_array

  use pmc_aero_particle
  use pmc_util
  use pmc_spec_file
  use pmc_mpi
#ifdef PMC_USE_MPI
  use mpi
#endif

  !> 1-D arrays of particles, used by aero_state to build a ragged
  !> array.
  !!
  !! One aero_particle_array is generally a list of particles in a
  !! single size bin, but the basic type can be used for any list of
  !! particles.
  !!
  !! To give a reasonable tradeoff between frequent re-allocs and
  !! memory usage, the length of an aero_particle_array is generally a
  !! bit longer than the number of particles stored in it. When the
  !! array is full then a larger array is allocated with new extra
  !! space. As a balance between memory usage and frequency of
  !! re-allocs the length of the array is currently doubled when
  !! necessary and halved when possible.
  !!
  !! The true allocated length of the aero_particle_array can be
  !! obtained by size(aero_particle_array%%particle), while the number
  !! of used particle slots in it is given by
  !! aero_particle_array%%n_part. It must be that
  !! aero_particle_array%%n_part is less than or equal to
  !! size(aero_particle_array%%particle).
  type aero_particle_array_t
     !> Number of particles.
     integer :: n_part
     !> Particle array.
     type(aero_particle_t), pointer :: particle(:)
  end type aero_particle_array_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Allocates and initializes.
  subroutine aero_particle_array_allocate(aero_particle_array)

    !> Result.
    type(aero_particle_array_t), intent(out) :: aero_particle_array

    aero_particle_array%n_part = 0
    allocate(aero_particle_array%particle(0))

  end subroutine aero_particle_array_allocate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Allocates and initializes to the given size.
  subroutine aero_particle_array_allocate_size(aero_particle_array, n_part)

    !> Result.
    type(aero_particle_array_t), intent(out) :: aero_particle_array
    !> Number of particles.
    integer, intent(in) :: n_part

    integer :: i

    aero_particle_array%n_part = n_part
    allocate(aero_particle_array%particle(n_part))
    do i = 1,n_part
       call aero_particle_allocate(aero_particle_array%particle(i))
    end do

  end subroutine aero_particle_array_allocate_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Deallocates.
  subroutine aero_particle_array_deallocate(aero_particle_array)

    !> Structure to deallocate.
    type(aero_particle_array_t), intent(inout) :: aero_particle_array

    integer :: i

    do i = 1,aero_particle_array%n_part
       call aero_particle_deallocate(aero_particle_array%particle(i))
    end do
    deallocate(aero_particle_array%particle)

  end subroutine aero_particle_array_deallocate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Copies aero_particle_array_from to aero_particle_array_to, both
  !> of which must already be allocated.
  subroutine aero_particle_array_copy(aero_particle_array_from, &
       aero_particle_array_to)

    !> Origin structure.
    type(aero_particle_array_t), intent(in) :: aero_particle_array_from
    !> Destination structure.
    type(aero_particle_array_t), intent(inout) :: aero_particle_array_to

    integer :: i

    call aero_particle_array_deallocate(aero_particle_array_to)
    call aero_particle_array_allocate_size(aero_particle_array_to, &
         aero_particle_array_from%n_part)
    do i = 1,aero_particle_array_from%n_part
       call aero_particle_copy(aero_particle_array_from%particle(i), &
            aero_particle_array_to%particle(i))
    end do

  end subroutine aero_particle_array_copy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Resets an aero_particle_array to contain zero particles.
  subroutine aero_particle_array_zero(aero_particle_array)

    !> Structure to reset.
    type(aero_particle_array_t), intent(inout) :: aero_particle_array

    call aero_particle_array_deallocate(aero_particle_array)
    allocate(aero_particle_array%particle(0))
    aero_particle_array%n_part = 0

  end subroutine aero_particle_array_zero

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Changes the given aero_particle_array to exactly the given
  !> new_length.
  !!
  !! This function should not be called directly, but rather use
  !! aero_particle_array_enlarge(), aero_particle_array_enlarge_to()
  !! or aero_particle_array_shrink().
  subroutine aero_particle_array_realloc(aero_particle_array, new_length)

    !> Array to reallocate (must already be allocated on entry).
    type(aero_particle_array_t), intent(inout) :: aero_particle_array
    !> New length of the array.
    integer, intent(in) :: new_length

    integer :: n_part, i
    type(aero_particle_t), pointer :: new_particles(:)

    n_part = aero_particle_array%n_part
    call assert(867444847, new_length >= n_part)
    allocate(new_particles(new_length))
    do i = 1,aero_particle_array%n_part
       call aero_particle_shift(aero_particle_array%particle(i), &
            new_particles(i))
    end do
    deallocate(aero_particle_array%particle)
    aero_particle_array%particle => new_particles

  end subroutine aero_particle_array_realloc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Enlarges the given aero_particle_array by at least one element
  !!
  !! Currently this doubles the length.
  subroutine aero_particle_array_enlarge(aero_particle_array)

    !> Array to enlarge.
    type(aero_particle_array_t), intent(inout) :: aero_particle_array

    integer :: length, new_length

    length = size(aero_particle_array%particle)
    new_length = max(length * 2, length + 1)
    call aero_particle_array_realloc(aero_particle_array, new_length)

  end subroutine aero_particle_array_enlarge

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Enlarges the given array so that it is at least of size n.
  subroutine aero_particle_array_enlarge_to(aero_particle_array, n)

    !> Array to enlarge.
    type(aero_particle_array_t), intent(inout) :: aero_particle_array
    !> Minimum new size of array.
    integer, intent(in) :: n

    do while (size(aero_particle_array%particle) < n)
       call aero_particle_array_enlarge(aero_particle_array)
    end do

  end subroutine aero_particle_array_enlarge_to

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Possibly shrinks the storage of the given array, ensuring that
  !> it can still store the allocated particles.
  subroutine aero_particle_array_shrink(aero_particle_array)

    !> Array to shrink.
    type(aero_particle_array_t), intent(inout) :: aero_particle_array

    integer :: n_part, length, new_length

    n_part = aero_particle_array%n_part
    length = size(aero_particle_array%particle)
    new_length = length / 2
    do while ((n_part <= new_length) .and. (length > 0))
       call aero_particle_array_realloc(aero_particle_array, new_length)
       length = size(aero_particle_array%particle)
       new_length = length / 2
    end do

  end subroutine aero_particle_array_shrink

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Adds the given particle to the end of the array.
  subroutine aero_particle_array_add_particle(aero_particle_array, &
       aero_particle)

    !> Array to add to.
    type(aero_particle_array_t), intent(inout) :: aero_particle_array
    !> Particle to add.
    type(aero_particle_t), intent(in) :: aero_particle

    integer :: n

    n = aero_particle_array%n_part + 1
    call aero_particle_array_enlarge_to(aero_particle_array, n)
    call aero_particle_allocate(aero_particle_array%particle(n))
    call aero_particle_copy(aero_particle, &
         aero_particle_array%particle(n))
    aero_particle_array%n_part = aero_particle_array%n_part + 1

  end subroutine aero_particle_array_add_particle

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Removes the particle at the given index.
  subroutine aero_particle_array_remove_particle(aero_particle_array, &
       index)

    !> Array to remove from.
    type(aero_particle_array_t), intent(inout) :: aero_particle_array
    !> Index of particle to remove.
    integer, intent(in) :: index

    call assert(992946227, index >= 1)
    call assert(711246139, index <= aero_particle_array%n_part)
    call aero_particle_deallocate(aero_particle_array%particle(index))
    if (index < aero_particle_array%n_part) then
       ! shift last particle into empty slot to preserve dense packing
       call aero_particle_shift( &
            aero_particle_array%particle(aero_particle_array%n_part), &
            aero_particle_array%particle(index))
    end if
    aero_particle_array%n_part = aero_particle_array%n_part - 1
    call aero_particle_array_shrink(aero_particle_array)

  end subroutine aero_particle_array_remove_particle

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Determines the number of bytes required to pack the given value.
  integer function pmc_mpi_pack_size_apa(val)

    !> Value to pack.
    type(aero_particle_array_t), intent(in) :: val

    integer :: i, total_size

    total_size = 0
    total_size = total_size + pmc_mpi_pack_size_integer(val%n_part)
    do i = 1,val%n_part
       total_size = total_size &
            + pmc_mpi_pack_size_aero_particle(val%particle(i))
    end do
    pmc_mpi_pack_size_apa = total_size

  end function pmc_mpi_pack_size_apa

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Packs the given value into the buffer, advancing position.
  subroutine pmc_mpi_pack_aero_particle_array(buffer, position, val)

    !> Memory buffer.
    character, intent(inout) :: buffer(:)
    !> Current buffer position.
    integer, intent(inout) :: position
    !> Value to pack.
    type(aero_particle_array_t), intent(in) :: val

#ifdef PMC_USE_MPI
    integer :: prev_position, i

    prev_position = position
    call pmc_mpi_pack_integer(buffer, position, val%n_part)
    do i = 1,val%n_part
       call pmc_mpi_pack_aero_particle(buffer, position, val%particle(i))
    end do
    call assert(803856329, &
         position - prev_position <= pmc_mpi_pack_size_apa(val))
#endif

  end subroutine pmc_mpi_pack_aero_particle_array

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Unpacks the given value from the buffer, advancing position.
  subroutine pmc_mpi_unpack_aero_particle_array(buffer, position, val)

    !> Memory buffer.
    character, intent(inout) :: buffer(:)
    !> Current buffer position.
    integer, intent(inout) :: position
    !> Value to pack.
    type(aero_particle_array_t), intent(inout) :: val

#ifdef PMC_USE_MPI
    integer :: prev_position, i

    call aero_particle_array_deallocate(val)
    prev_position = position
    call pmc_mpi_unpack_integer(buffer, position, val%n_part)
    allocate(val%particle(val%n_part))
    do i = 1,val%n_part
       call aero_particle_allocate(val%particle(i))
       call pmc_mpi_unpack_aero_particle(buffer, position, val%particle(i))
    end do
    call assert(138783294, &
         position - prev_position <= pmc_mpi_pack_size_apa(val))
#endif

  end subroutine pmc_mpi_unpack_aero_particle_array

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_aero_particle_array
