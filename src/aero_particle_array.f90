! Copyright (C) 2005-2007 Nicole Riemer and Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.
!
! Basic particle structure.

module mod_aero_particle_array

  use mod_aero_particle

  type aero_particle_array_t
     integer :: n_part                  ! number of particles
     integer :: n_spec                  ! number of species
     type(aero_particle_t), pointer :: particles(:) ! particles
     ! NOTE: typically size(particles) > num as we allocate more meory
     ! than needed. num is the number of entries used in particles.
  end type aero_particle_array_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine aero_particle_array_alloc(n_part, n_spec, aero_particle_array)

    ! Allocates and initializes.

    integer, intent(in) :: n_part       ! number of particles
    integer, intent(in) :: n_spec       ! number of species
    type(aero_particle_array_t), intent(inout) :: aero_particle_array ! result

    integer :: i

    aero_particle_array%n_part = n_part
    aero_particle_array%n_spec = n_spec
    allocate(aero_particle_array%particles(n_part))
    do i = 1,n_part
       call aero_particle_alloc(aero_particle_array%particles(i), n_spec)
    end do

  end subroutine aero_particle_array_alloc
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine aero_particle_array_free(aero_particle_array)

    ! Deallocates.

    type(aero_particle_array_t), intent(inout) :: aero_particle_array

    integer :: i
    
    do i = 1,aero_particle_array%n_part
       call aero_particle_free(aero_particle_array%particles(i))
    end do
    deallocate(aero_particle_array%particles)

  end subroutine aero_particle_array_free
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine aero_particle_array_copy(aero_particle_array_from, &
       aero_particle_array_to)

    ! Copies aero_particle_array_from to aero_particle_array_to, both
    ! of which must already be allocated.

    type(aero_particle_array_t), intent(in) :: aero_particle_array_from
    type(aero_particle_array_t), intent(inout) :: aero_particle_array_to

    integer :: i
    
    call aero_particle_array_free(aero_particle_array_to)
    call aero_particle_array_alloc(aero_particle_array_from%n_part, &
         aero_particle_array_from%n_spec, aero_particle_array_to)
    do i = 1,aero_particle_array_from%n_part
       call aero_particle_copy(aero_particle_array_from%particles(i), &
            aero_particle_array_to%particles(i))
    end do

  end subroutine aero_particle_array_copy
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine aero_particle_array_zero(aero_particle_array)

    ! Resets an aero_particle_array to contain zero particles.

    type(aero_particle_array_t), intent(inout) :: aero_particle_array

    call aero_particle_array_free(aero_particle_array)
    allocate(aero_particle_array%particles(0))
    aero_particle_array%n_part = 0

  end subroutine aero_particle_array_zero
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine aero_particle_array_realloc(aero_particle_array, new_length)

    ! Changes the given aero_particle_array (which must be allocated)
    ! to exactly the given new_length. This function should not be
    ! called directly, but rather use aero_particle_array_enlarge(),
    ! aero_particle_array_enlarge_to() or aero_particle_array_shrink().

    use mod_util

    type(aero_particle_array_t), intent(inout) :: aero_particle_array
    integer, intent(in) :: new_length   ! new length of the array

    integer :: n_part, n_spec, i
    type(aero_particle_t), pointer :: new_particles(:)

    n_part = aero_particle_array%n_part
    n_spec = aero_particle_array%n_spec
    call assert(new_length >= n_part)
    allocate(new_particles(new_length))
    do i = 1,aero_particle_array%n_part
       call aero_particle_shift(aero_particle_array%particles(i), &
            new_particles(i))
    end do
    deallocate(aero_particle_array%particles)
    aero_particle_array%particles => new_particles
    
  end subroutine aero_particle_array_realloc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine aero_particle_array_enlarge(aero_particle_array)

    ! Enlarges the given aero_particle_array (which must be allocated)
    ! by at least one element (currently doubles the length).

    type(aero_particle_array_t), intent(inout) :: aero_particle_array

    integer :: length, new_length

    length = size(aero_particle_array%particles)
    new_length = max(length * 2, length + 1)
    call aero_particle_array_realloc(aero_particle_array, new_length)
    
  end subroutine aero_particle_array_enlarge

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine aero_particle_array_enlarge_to(aero_particle_array, n)

    ! Enlarges the given array so that it is at least of size n.

    type(aero_particle_array_t), intent(inout) :: aero_particle_array
    integer, intent(in) :: n            ! minimum new size of array

    do while (size(aero_particle_array%particles) < n)
       call aero_particle_array_enlarge(aero_particle_array)
    end do

  end subroutine aero_particle_array_enlarge_to

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine aero_particle_array_shrink(aero_particle_array)

    ! Possibly shrinks the storage of the given array, ensuring that
    ! it can still store the allocated particles.

    type(aero_particle_array_t), intent(inout) :: aero_particle_array

    integer :: n_part, n_spec, length, new_length

    n_part = aero_particle_array%n_part
    n_spec = aero_particle_array%n_spec
    length = size(aero_particle_array%particles)
    new_length = length / 2
    do while ((n_part <= new_length) .and. (length > 0))
       call aero_particle_array_realloc(aero_particle_array, new_length)
       length = size(aero_particle_array%particles)
       new_length = length / 2
    end do

  end subroutine aero_particle_array_shrink

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine aero_particle_array_add_particle(aero_particle_array, &
       aero_particle)

    ! Adds the given particle to the end of the array.

    type(aero_particle_array_t), intent(inout) :: aero_particle_array
    type(aero_particle_t), intent(in) :: aero_particle ! particle to add

    integer :: n

    n = aero_particle_array%n_part + 1
    call aero_particle_array_enlarge_to(aero_particle_array, n)
    call aero_particle_alloc(aero_particle_array%particles(n), 0)
    call aero_particle_copy(aero_particle, aero_particle_array%particles(n))
    aero_particle_array%n_part = aero_particle_array%n_part + 1

  end subroutine aero_particle_array_add_particle

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine aero_particle_array_remove_particle(aero_particle_array, &
       index)

    ! Removes the particle at the given index.

    use mod_util

    type(aero_particle_array_t), intent(inout) :: aero_particle_array
    integer, intent(in) :: index        ! index of particle to remove

    call assert(index >= 1)
    call assert(index <= aero_particle_array%n_part)
    call aero_particle_free(aero_particle_array%particles(index))
    if (index < aero_particle_array%n_part) then
       ! shift last particle into empty slot to preserve dense packing
       call aero_particle_shift( &
            aero_particle_array%particles(aero_particle_array%n_part), &
            aero_particle_array%particles(index))
    end if
    aero_particle_array%n_part = aero_particle_array%n_part - 1
    call aero_particle_array_shrink(aero_particle_array)

  end subroutine aero_particle_array_remove_particle

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine aero_particle_array_double(aero_particle_array)

    ! Doubles the number of particles by making a duplicate of each
    ! one.

    type(aero_particle_array_t), intent(inout) :: aero_particle_array

    integer :: n, i

    n = aero_particle_array%n_part
    call aero_particle_array_enlarge_to(aero_particle_array, 2 * n)
    do i = 1,n
       call aero_particle_alloc(aero_particle_array%particles(i + n), &
            aero_particle_array%n_spec)
       call aero_particle_copy(aero_particle_array%particles(i), &
            aero_particle_array%particles(i + n))
    end do
    aero_particle_array%n_part = 2 * n

  end subroutine aero_particle_array_double

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine inout_write_aero_particle_array(file, aero_particle_array)
    
    ! Write full state.
    
    use mod_inout
    
    type(inout_file_t), intent(inout) :: file ! file to write to
    type(aero_particle_array_t), intent(in) :: aero_particle_array

    integer :: i

    call inout_write_integer(file, "n_part", aero_particle_array%n_part)
    call inout_write_integer(file, "n_spec", aero_particle_array%n_spec)
    do i = 1,aero_particle_array%n_part
       call inout_write_integer(file, "particle_number", i)
       call inout_write_aero_particle(file, aero_particle_array%particles(i))
    end do
    
  end subroutine inout_write_aero_particle_array
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine inout_read_aero_particle_array(file, aero_particle_array)
    
    ! Read full state.
    
    use mod_inout
    
    type(inout_file_t), intent(inout) :: file ! file to write to
    type(aero_particle_array_t), intent(out) :: aero_particle_array

    integer :: i, check_i

    call inout_read_integer(file, "n_part", aero_particle_array%n_part)
    call inout_read_integer(file, "n_spec", aero_particle_array%n_spec)
    allocate(aero_particle_array%particles(aero_particle_array%n_part))
    do i = 1,aero_particle_array%n_part
       call inout_read_integer(file, "particle_number", check_i)
       call inout_check_index(file, i, check_i)
       call inout_read_aero_particle(file, aero_particle_array%particles(i))
    end do
    
  end subroutine inout_read_aero_particle_array
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
end module mod_aero_particle_array
