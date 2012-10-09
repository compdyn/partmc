! Copyright (C) 2005-2012 Nicole Riemer and Matthew West
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

  !> 1-D array of particles, used by aero_state to store the
  !> particles.
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
  !! size(aero_particle_array%%particle). In user code, the \c
  !! aero_particle_array_n_part() getter function should be used.q
  !!
  !! For internal usage, if \c particle is not allocated then \c
  !! n_part is invalid. If \c particle is allocated then \c n_part
  !! must be valid.
  type aero_particle_array_t
     !> Number of particles.
     integer :: n_part
     !> Particle array.
     type(aero_particle_t), allocatable :: particle(:)
  end type aero_particle_array_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Return the current number of particles.
  elemental integer function aero_particle_array_n_part(aero_particle_array)

    !> Aerosol particle array.
    type(aero_particle_array_t), intent(in) :: aero_particle_array

    if (allocated(aero_particle_array%particle)) then
       aero_particle_array_n_part = aero_particle_array%n_part
    else
       aero_particle_array_n_part = 0
    end if

  end function aero_particle_array_n_part

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Resets an aero_particle_array to contain zero particles.
  subroutine aero_particle_array_zero(aero_particle_array)

    !> Structure to reset.
    type(aero_particle_array_t), intent(inout) :: aero_particle_array

    aero_particle_array%n_part = 0
    if (allocated(aero_particle_array%particle)) then
       deallocate(aero_particle_array%particle)
    end if

  end subroutine aero_particle_array_zero

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Changes the given aero_particle_array to exactly the given
  !> new_length.
  !!
  !! This function should not be called directly, but rather use
  !! aero_particle_array_enlarge() or aero_particle_array_shrink().
  subroutine aero_particle_array_realloc(aero_particle_array, new_length)

    !> Array to reallocate.
    type(aero_particle_array_t), intent(inout) :: aero_particle_array
    !> New length of the array.
    integer, intent(in) :: new_length

    integer :: i
    type(aero_particle_t), allocatable :: new_particles(:)

    if (.not. allocated(aero_particle_array%particle)) then
       allocate(aero_particle_array%particle(new_length))
       aero_particle_array%n_part = 0
       return
    end if

    call assert(867444847, new_length >= aero_particle_array%n_part)
    allocate(new_particles(new_length))
    do i = 1,aero_particle_array%n_part
       call aero_particle_shift(aero_particle_array%particle(i), &
            new_particles(i))
    end do
    call move_alloc(new_particles, aero_particle_array%particle)

  end subroutine aero_particle_array_realloc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Possibly enlarges the given array, ensuring that it is at least of size n.
  subroutine aero_particle_array_enlarge(aero_particle_array, n)

    !> Array to enlarge.
    type(aero_particle_array_t), intent(inout) :: aero_particle_array
    !> Minimum new size of array.
    integer, intent(in) :: n

    if (.not. allocated(aero_particle_array%particle)) then
       call aero_particle_array_realloc(aero_particle_array, pow2_above(n))
       return
    end if

    if (n <= size(aero_particle_array%particle)) return

    call aero_particle_array_realloc(aero_particle_array, pow2_above(n))

  end subroutine aero_particle_array_enlarge

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Possibly shrinks the storage of the given array, ensuring that
  !> it can still store the allocated particles.
  subroutine aero_particle_array_shrink(aero_particle_array)

    !> Array to shrink.
    type(aero_particle_array_t), intent(inout) :: aero_particle_array

    integer :: new_length

    if (.not. allocated(aero_particle_array%particle)) return

    new_length = pow2_above(aero_particle_array%n_part)
    if (new_length < size(aero_particle_array%particle)) then
       call aero_particle_array_realloc(aero_particle_array, new_length)
    end if

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

    n = aero_particle_array_n_part(aero_particle_array) + 1
    call aero_particle_array_enlarge(aero_particle_array, n)
    aero_particle_array%particle(n) = aero_particle
    aero_particle_array%n_part = n

  end subroutine aero_particle_array_add_particle

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Removes the particle at the given index.
  subroutine aero_particle_array_remove_particle(aero_particle_array, &
       index)

    !> Array to remove from.
    type(aero_particle_array_t), intent(inout) :: aero_particle_array
    !> Index of particle to remove.
    integer, intent(in) :: index

    call assert(883639923, allocated(aero_particle_array%particle))
    call assert(992946227, index >= 1)
    call assert(711246139, index <= aero_particle_array%n_part)
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
    total_size = total_size &
         + pmc_mpi_pack_size_integer(aero_particle_array_n_part(val))
    do i = 1,aero_particle_array_n_part(val)
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
    call pmc_mpi_pack_integer(buffer, position, &
         aero_particle_array_n_part(val))
    do i = 1,aero_particle_array_n_part(val)
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
    integer :: prev_position, i, n

    prev_position = position
    call pmc_mpi_unpack_integer(buffer, position, n)
    call aero_particle_array_realloc(val, n)
    val%n_part = n
    do i = 1,n
       call pmc_mpi_unpack_aero_particle(buffer, position, val%particle(i))
    end do
    call assert(138783294, &
         position - prev_position <= pmc_mpi_pack_size_apa(val))
#endif

  end subroutine pmc_mpi_unpack_aero_particle_array

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Check that the particle array data is consistent.
  subroutine aero_particle_array_check(aero_particle_array, aero_data, &
       continue_on_error)

    !> Aerosol particle array to check.
    type(aero_particle_array_t), intent(in) :: aero_particle_array
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Whether to continue despite error.
    logical, intent(in) :: continue_on_error

    integer :: i_part

    if (.not. allocated(aero_particle_array%particle)) return

    if (aero_particle_array%n_part < 0) then
       write(0, *) 'ERROR aero_particle_array A:'
       write(0, *) 'aero_particle_array%n_part', aero_particle_array%n_part
       call assert(250011397, continue_on_error)
    end if

    do i_part = 1,aero_particle_array%n_part
       call aero_particle_check(aero_particle_array%particle(i_part), &
            aero_data, continue_on_error)
    end do

  end subroutine aero_particle_array_check

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_aero_particle_array
