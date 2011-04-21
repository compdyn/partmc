! Copyright (C) 2011 Nicole Riemer and Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_aero_sorted module.

!> The aero_sorted_t structure and assocated subroutines.
module pmc_aero_sorted

  use pmc_integer_varray
  use pmc_aero_particle
  use pmc_aero_particle_array
  use pmc_mpi

  !> Bin index for particles sorted into bins.
  !!
  !! Both forward and reverse indexes are maintained. Particles are
  !! stored with both a linaer index \c i_part, and binned indexes
  !! <tt>(i_bin, i_entry)</tt>, indicating that the particle is number
  !! \c i_entry in bin number \c i_bin. The forward index satisfies
  !! \code
  !! i_part = aero_sorted%bin(i_bin)%entry(i_part)
  !! \endcode
  !! while the reverse index satisfies
  !! \code
  !! i_bin = aero_sorted%reverse_bin%entry(i_part)
  !! i_entry = aero_sorted%reverse_entry%entry(i_part)
  !! \endcode
  type aero_sorted_t
     !> Array of integer arrays, one per bin.
     type(integer_varray_t), allocatable, dimension(:) :: bin
     !> Reverse index array to bin numbers.
     type(integer_varray_t) :: reverse_bin
     !> Reverse index array to particle entry-in-bin numbers.
     type(integer_varray_t) :: reverse_entry
  end type aero_sorted_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Allocate an empty structure.
  subroutine aero_sorted_allocate(aero_sorted)

    !> Structure to initialize.
    type(aero_sorted_t), intent(out) :: aero_sorted
    
    allocate(aero_sorted%bin(0))
    call integer_varray_allocate(aero_sorted%reverse_bin)
    call integer_varray_allocate(aero_sorted%reverse_entry)

  end subroutine aero_sorted_allocate
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Allocate a strcture with the given size.
  subroutine aero_sorted_allocate_size(aero_sorted, n_bin)

    !> Structure to initialize.
    type(aero_sorted_t), intent(out) :: aero_sorted
    !> Number of bins.
    integer, intent(in) :: n_bin
    
    allocate(aero_sorted%bin(n_bin))
    call integer_varray_allocate(aero_sorted%bin)
    call integer_varray_allocate(aero_sorted%reverse_bin)
    call integer_varray_allocate(aero_sorted%reverse_entry)

  end subroutine aero_sorted_allocate_size
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Deallocates a previously allocated structure.
  subroutine aero_sorted_deallocate(aero_sorted)

    !> Structure to deallocate.
    type(aero_sorted_t), intent(inout) :: aero_sorted
    
    call integer_varray_deallocate(aero_sorted%bin)
    deallocate(aero_sorted%bin)
    call integer_varray_deallocate(aero_sorted%reverse_bin)
    call integer_varray_deallocate(aero_sorted%reverse_entry)

  end subroutine aero_sorted_deallocate
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Resets an aero_sorted to have zero particles per bin.
  subroutine aero_sorted_zero(aero_sorted)

    !> Structure to zero.
    type(aero_sorted_t), intent(inout) :: aero_sorted
    
    integer :: i_bin

    do i_bin = 1,size(aero_sorted%bin)
       call integer_varray_zero(aero_sorted%bin(i_bin))
    end do
    call integer_varray_zero(aero_sorted%reverse_bin)
    call integer_varray_zero(aero_sorted%reverse_entry)

  end subroutine aero_sorted_zero
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Add a new particle to both an aero_sorted and the corresponding
  !> aero_particle_array.
  subroutine aero_sorted_add_particle(aero_sorted, aero_particle_array, &
       aero_particle, i_bin)

    !> Sorted particle structure.
    type(aero_sorted_t), intent(inout) :: aero_sorted
    !> Aerosol particles.
    type(aero_particle_array_t), intent(inout) :: aero_particle_array
    !> Particle to add.
    type(aero_particle_t), intent(in) :: aero_particle
    !> Bin number of added particle.
    integer, intent(in) :: i_bin

    ! add the particle to the aero_particle_array
    call aero_particle_array_add_particle(aero_particle_array, aero_particle)

    ! update the forward index
    call integer_varray_append(aero_sorted%bin(i_bin), &
         aero_particle_array%n_part)

    ! update the reverse index
    call integer_varray_append(aero_sorted%reverse_bin, i_bin)
    call integer_varray_append(aero_sorted%reverse_entry, &
         aero_sorted%bin(i_bin)%n_entry)

  end subroutine aero_sorted_add_particle

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Remove a particle from both an aero_sorted and the corresponding
  !> aero_particle_array.
  subroutine aero_sorted_remove_particle(aero_sorted, aero_particle_array, &
       i_part)

    !> Sorted particle structure.
    type(aero_sorted_t), intent(inout) :: aero_sorted
    !> Aerosol particles.
    type(aero_particle_array_t), intent(inout) :: aero_particle_array
    !> Index of particle to remove.
    integer, intent(in) :: i_part

    integer :: i_bin, i_entry, i_part_shifted
    integer :: i_bin_fix, i_part_fix, i_entry_fix

    ! Deleting particles shifts the end particles into the empty slots
    ! in the aero_particle_array and the aero_sorted forward and
    ! reverse indexes. All must be fixed in the right order to
    ! maintain consistency.

    i_bin = aero_sorted%reverse_bin%entry(i_part)
    i_entry = aero_sorted%reverse_entry%entry(i_part)

    ! remove the particle from the aero_particle_array
    i_part_shifted = aero_particle_array%n_part ! old loc of shifted particle
    call aero_particle_array_remove_particle(aero_particle_array, i_part)

    if (i_part_shifted /= i_part) then
       ! fix up the forward index for the shifted particle
       i_bin_fix = aero_sorted%reverse_bin%entry(i_part_shifted)
       i_part_fix = aero_sorted%reverse_entry%entry(i_part_shifted)
       aero_sorted%bin(i_bin_fix)%entry(i_part_fix) = i_part
    end if

    ! remove the particle from the reverse index (with the side effect
    ! of fixing the reverse map for the shifted particle)
    call integer_varray_remove_entry(aero_sorted%reverse_bin, i_part)
    call integer_varray_remove_entry(aero_sorted%reverse_entry, i_part)

    ! remove the forward index entry
    i_entry_fix = aero_sorted%bin(i_bin)%n_entry
    i_part_fix = aero_sorted%bin(i_bin)%entry(i_entry_fix)
    call integer_varray_remove_entry(aero_sorted%bin(i_bin), i_entry)

    if (i_entry_fix /= i_entry) then
       ! fix reverse index
       aero_sorted%reverse_entry%entry(i_part_fix) = i_entry
    end if

  end subroutine aero_sorted_remove_particle

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Determines the number of bytes required to pack the given value.
  integer function pmc_mpi_pack_size_aero_sorted(val)

    !> Value to pack.
    type(aero_sorted_t), intent(in) :: val

    integer :: i, total_size

    total_size = 0
    total_size = total_size + pmc_mpi_pack_size_integer(size(val%bin))
    do i = 1,size(val%bin)
       total_size = total_size + pmc_mpi_pack_size_integer_varray(val%bin(i))
    end do
    total_size = total_size &
         + pmc_mpi_pack_size_integer_varray(val%reverse_bin)
    total_size = total_size &
         + pmc_mpi_pack_size_integer_varray(val%reverse_entry)
    pmc_mpi_pack_size_aero_sorted = total_size

  end function pmc_mpi_pack_size_aero_sorted

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Packs the given value into the buffer, advancing position.
  subroutine pmc_mpi_pack_aero_sorted(buffer, position, val)

    !> Memory buffer.
    character, intent(inout) :: buffer(:)
    !> Current buffer position.
    integer, intent(inout) :: position
    !> Value to pack.
    type(aero_sorted_t), intent(in) :: val

#ifdef PMC_USE_MPI
    integer :: prev_position, i

    prev_position = position
    call pmc_mpi_pack_integer(buffer, position, size(val%bin))
    do i = 1,size(val%bin)
       call pmc_mpi_pack_integer_varray(buffer, position, val%bin(i))
    end do
    call pmc_mpi_pack_integer_varray(buffer, position, val%reverse_bin)
    call pmc_mpi_pack_integer_varray(buffer, position, val%reverse_entry)
    call assert(178297816, &
         position - prev_position <= pmc_mpi_pack_size_aero_sorted(val))
#endif

  end subroutine pmc_mpi_pack_aero_sorted

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Unpacks the given value from the buffer, advancing position.
  subroutine pmc_mpi_unpack_aero_sorted(buffer, position, val)

    !> Memory buffer.
    character, intent(inout) :: buffer(:)
    !> Current buffer position.
    integer, intent(inout) :: position
    !> Value to pack.
    type(aero_sorted_t), intent(inout) :: val

#ifdef PMC_USE_MPI
    integer :: prev_position, i, n

    prev_position = position
    call pmc_mpi_unpack_integer(buffer, position, n)
    call aero_sorted_deallocate(val)
    call aero_sorted_allocate_size(val, n)
    do i = 1,size(val%bin)
       call pmc_mpi_unpack_integer_varray(buffer, position, val%bin(i))
    end do
    call pmc_mpi_unpack_integer_varray(buffer, position, val%reverse_bin)
    call pmc_mpi_unpack_integer_varray(buffer, position, val%reverse_entry)
    call assert(364064630, &
         position - prev_position <= pmc_mpi_pack_size_aero_sorted(val))
#endif

  end subroutine pmc_mpi_unpack_aero_sorted

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
end module pmc_aero_sorted
