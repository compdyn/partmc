! Copyright (C) 2011 Nicole Riemer and Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_aero_sorted module.

!> The aero_sorted_t structure and assocated subroutines.
module pmc_aero_sorted

  use pmc_integer_array
  use pmc_aero_particle
  use pmc_aero_particle_array
  !>DEBUG
  use pmc_aero_state
  !<DEBUG

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
     type(integer_array_t), allocatable, dimension(:) :: bin
     !> Reverse index array to bin numbers.
     type(integer_array_t) :: reverse_bin
     !> Reverse index array to particle entry-in-bin numbers.
     type(integer_array_t) :: reverse_entry
  end type aero_sorted_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Allocate an empty structure.
  subroutine aero_sorted_allocate(aero_sorted)

    !> Structure to initialize.
    type(aero_sorted_t), intent(out) :: aero_sorted
    
    allocate(aero_sorted%bin(0))
    call integer_array_allocate(aero_sorted%reverse_bin)
    call integer_array_allocate(aero_sorted%reverse_entry)

  end subroutine aero_sorted_allocate
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Allocate a strcture with the given size.
  subroutine aero_sorted_allocate_size(aero_sorted, n_bin)

    !> Structure to initialize.
    type(aero_sorted_t), intent(out) :: aero_sorted
    !> Number of bins.
    integer, intent(in) :: n_bin
    
    allocate(aero_sorted%bin(n_bin))
    call integer_array_allocate(aero_sorted%bin)
    call integer_array_allocate(aero_sorted%reverse_bin)
    call integer_array_allocate(aero_sorted%reverse_entry)

  end subroutine aero_sorted_allocate_size
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Deallocates a previously allocated structure.
  subroutine aero_sorted_deallocate(aero_sorted)

    !> Structure to deallocate.
    type(aero_sorted_t), intent(inout) :: aero_sorted
    
    call integer_array_deallocate(aero_sorted%bin)
    deallocate(aero_sorted%bin)
    call integer_array_deallocate(aero_sorted%reverse_bin)
    call integer_array_deallocate(aero_sorted%reverse_entry)

  end subroutine aero_sorted_deallocate
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Fills in particle indexes from the given aero_state.
  subroutine aero_sorted_fill(aero_sorted, aero_particle_array, bin_grid)

    !> Sorted particle structure.
    type(aero_sorted_t), intent(inout) :: aero_sorted
    !> Aerosol particles.
    type(aero_particle_array_t), intent(in) :: aero_particle_array
    !> Bin grid.
    type(bin_grid_t), intent(in) :: bin_grid

    integer :: i_part, i_bin

    call aero_sorted_deallocate(aero_sorted)
    call aero_sorted_allocate_size(aero_sorted, bin_grid%n_bin)
    do i_part = 1,aero_particle_array%n_part
       i_bin = aero_particle_in_bin(aero_particle_array%particle(i_part), &
            bin_grid)

       ! fill in forward index
       call integer_array_append(aero_sorted%bin(i_bin), i_part)

       ! fill in reverse index
       call integer_array_append(aero_sorted%reverse_bin, i_bin)
       call integer_array_append(aero_sorted%reverse_entry, &
            aero_sorted%bin(i_bin)%n_entry)
    end do

  end subroutine aero_sorted_fill
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Check that the given aero_sorted is correct.
  subroutine aero_sorted_check(aero_sorted, aero_particle_array, bin_grid)

    !> Sorted particle structure.
    type(aero_sorted_t), intent(in) :: aero_sorted
    !> Aerosol particles.
    type(aero_particle_array_t), intent(in) :: aero_particle_array
    !> Bin grid.
    type(bin_grid_t), intent(in) :: bin_grid

    integer :: i_part, i_bin, i_entry

    if ((aero_particle_array%n_part /= aero_sorted%reverse_bin%n_entry) &
         .or. (aero_particle_array%n_part &
         /= aero_sorted%reverse_entry%n_entry)) then
       write(0,*) 'SORTED CHECK ERROR A'
       write(0,*) 'aero_particle_array%n_part', aero_particle_array%n_part
       write(0,*) 'reverse_bin%n_entry', aero_sorted%reverse_bin%n_entry
       write(0,*) 'reverse_entry%n_entry', aero_sorted%reverse_entry%n_entry
    end if

    do i_part = 1,aero_particle_array%n_part
       if ((i_part > aero_sorted%reverse_bin%n_entry) &
            .or. (i_part > aero_sorted%reverse_entry%n_entry)) then
          write(0,*) 'SORTED CHECK ERROR B'
          write(0,*) 'i_part', i_part
          write(0,*) 'ID', aero_particle_array%particle(i_part)%id
          write(0,*) 'reverse_bin%n_entry', aero_sorted%reverse_bin%n_entry
          write(0,*) 'reverse_entry%n_entry', &
               aero_sorted%reverse_entry%n_entry
       end if

       i_bin = aero_particle_in_bin(aero_particle_array%particle(i_part), &
            bin_grid)
       if (i_bin /= aero_sorted%reverse_bin%entry(i_part)) then
          write(0,*) 'SORTED CHECK ERROR C'
          write(0,*) 'i_part', i_part
          write(0,*) 'ID', aero_particle_array%particle(i_part)%id
          write(0,*) 'computed bin', i_bin
          write(0,*) 'reverse_bin%entry(i_part)', &
               aero_sorted%reverse_bin%entry(i_part)
          write(0,*) 'reverse_entry%entry(i_part)', &
               aero_sorted%reverse_entry%entry(i_part)
       end if

       i_entry = aero_sorted%reverse_entry%entry(i_part)
       if ((i_entry < 1) .or. (i_entry > aero_sorted%bin(i_bin)%n_entry)) then
          write(0,*) 'SORTED CHECK ERROR D'
          write(0,*) 'i_part', i_part
          write(0,*) 'ID', aero_particle_array%particle(i_part)%id
          write(0,*) 'computed bin', i_bin
          write(0,*) 'reverse_bin%entry(i_part)', &
               aero_sorted%reverse_bin%entry(i_part)
          write(0,*) 'reverse_entry%entry(i_part)', &
               aero_sorted%reverse_entry%entry(i_part)
          write(0,*) 'bin(i_bin)%n_entry', aero_sorted%bin(i_bin)%n_entry
       end if
       if (i_part /= aero_sorted%bin(i_bin)%entry(i_entry)) then
          write(0,*) 'SORTED CHECK ERROR E'
          write(0,*) 'i_part', i_part
          write(0,*) 'ID', aero_particle_array%particle(i_part)%id
          write(0,*) 'computed bin', i_bin
          write(0,*) 'reverse_bin%(i_part)', &
               aero_sorted%reverse_bin%entry(i_part)
          write(0,*) 'reverse_entry%entry(i_part)', &
               aero_sorted%reverse_entry%entry(i_part)
          write(0,*) 'forward i_part', &
               aero_sorted%bin(i_bin)%entry(i_entry)
       end if
    end do

    do i_bin = 1,bin_grid%n_bin
       do i_entry = 1,aero_sorted%bin(i_bin)%n_entry
          i_part = aero_sorted%bin(i_bin)%entry(i_entry)
          if ((i_part < 1) .or. (i_part > aero_particle_array%n_part)) then
             write(0,*) 'SORTED CHECK ERROR F'
             write(0,*) 'i_bin', i_bin
             write(0,*) 'i_entry', i_entry
             write(0,*) 'bin(i_bin)%entry(i_entry)', &
                  aero_sorted%bin(i_bin)%n_entry
             write(0,*) 'aero_particle_array%n_part', &
                  aero_particle_array%n_part
          end if
          if ((i_bin /= aero_sorted%reverse_bin%entry(i_part)) &
               .or. (i_entry /= aero_sorted%reverse_entry%entry(i_part))) then
             write(0,*) 'SORTED CHECK ERROR G'
             write(0,*) 'i_bin', i_bin
             write(0,*) 'i_entry', i_entry
             write(0,*) 'i_part', i_part
             write(0,*) 'reverse_bin%entry(i_part)', &
                  aero_sorted%reverse_bin%entry(i_part)
             write(0,*) 'reverse_entry%entry(i_part)', &
                  aero_sorted%reverse_entry%entry(i_part)
          end if
       end do
    end do

  end subroutine aero_sorted_check
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Temporary.
  subroutine aero_sorted_check_aero_state(aero_sorted, aero_particle_array, &
       bin_grid, aero_state)

    !> Sorted particle structure.
    type(aero_sorted_t), intent(in) :: aero_sorted
    !> Aerosol particles.
    type(aero_particle_array_t), intent(in) :: aero_particle_array
    !> Bin grid.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Aero state.
    type(aero_state_t), intent(in) :: aero_state

    integer :: i_bin, i_entry, i_entry_search, i_part
    integer :: i_entry_state, i_entry_sort, id_state, id_sort
    logical :: found
    integer, allocatable, dimension(:) :: data_state, data_sort
    integer, allocatable, dimension(:) :: perm_state, perm_sort

    do i_bin = 1,bin_grid%n_bin
       if (aero_state%bin(i_bin)%n_part &
            /= aero_sorted%bin(i_bin)%n_entry) then
          write(0,*) 'STATE CHECK ERROR A'
          write(0,*) 'i_bin', i_bin
          write(0,*) 'aero_state%bin(i_bin)%n_part', &
               aero_state%bin(i_bin)%n_part
          write(0,*) 'aero_sorted%bin(i_bin)%n_entry', &
               aero_sorted%bin(i_bin)%n_entry
       end if

       allocate(data_state(aero_state%bin(i_bin)%n_part))
       allocate(perm_state(aero_state%bin(i_bin)%n_part))
       allocate(data_sort(aero_state%bin(i_bin)%n_part))
       allocate(perm_sort(aero_state%bin(i_bin)%n_part))

       do i_entry = 1,aero_state%bin(i_bin)%n_part
          data_state(i_entry) = aero_state%bin(i_bin)%particle(i_entry)%id
          i_part = aero_sorted%bin(i_bin)%entry(i_entry)
          data_sort(i_entry) = aero_particle_array%particle(i_part)%id
       end do
       call integer_sort(data_state, perm_state)
       call integer_sort(data_sort, perm_sort)

       do i_entry = 1,aero_state%bin(i_bin)%n_part
          i_entry_state = perm_state(i_entry)
          i_entry_sort = perm_sort(i_entry)
          i_part = aero_sorted%bin(i_bin)%entry(i_entry_sort)
          id_state = aero_state%bin(i_bin)%particle(i_entry_state)%id
          id_sort = aero_particle_array%particle(i_part)%id
          if (id_state /= id_sort) then
             write(0,*) 'STATE CHECK ERROR B'
             write(0,*) 'i_bin', i_bin
             write(0,*) 'i_entry', i_entry
             write(0,*) 'i_entry_state', i_entry_state
             write(0,*) 'i_entry_sort', i_entry_sort
             write(0,*) 'i_part', i_part
             write(0,*) 'id_state', id_state
             write(0,*) 'id_sort', id_sort
          end if
       end do

       deallocate(data_state)
       deallocate(perm_state)
       deallocate(data_sort)
       deallocate(perm_sort)
    end do

  end subroutine aero_sorted_check_aero_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Add a new particle to both an aero_sorted and the corresponding
  !> aero_particle_array.
  subroutine aero_sorted_add_particle(aero_sorted, aero_particle_array, &
       bin_grid, aero_particle, i_bin)

    !> Sorted particle structure.
    type(aero_sorted_t), intent(inout) :: aero_sorted
    !> Aerosol particles.
    type(aero_particle_array_t), intent(inout) :: aero_particle_array
    !> Bin grid.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Particle to add.
    type(aero_particle_t), intent(in) :: aero_particle
    !> Bin number of added particle.
    integer, intent(in) :: i_bin

    ! add the particle to the aero_particle_array
    call aero_particle_array_add_particle(aero_particle_array, aero_particle)

    ! update the forward index
    call integer_array_append(aero_sorted%bin(i_bin), &
         aero_particle_array%n_part)

    ! update the reverse index
    call integer_array_append(aero_sorted%reverse_bin, i_bin)
    call integer_array_append(aero_sorted%reverse_entry, &
         aero_sorted%bin(i_bin)%n_entry)

  end subroutine aero_sorted_add_particle

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Remove a particle from both an aero_sorted and the corresponding
  !> aero_particle_array.
  subroutine aero_sorted_remove_particle(aero_sorted, aero_particle_array, &
       bin_grid, i_bin, i_entry)

    !> Sorted particle structure.
    type(aero_sorted_t), intent(inout) :: aero_sorted
    !> Aerosol particles.
    type(aero_particle_array_t), intent(inout) :: aero_particle_array
    !> Bin grid.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Bin number of particle to remove.
    integer, intent(in) :: i_bin
    !> Index in bin of particle to remove.
    integer, intent(in) :: i_entry

    integer :: i_part, i_part_shifted, i_bin_fix, i_part_fix, i_entry_fix

    ! Deleting particles shifts the end particles into the empty slots
    ! in the aero_particle_array and the aero_sorted forward and
    ! reverse indexes. All must be fixed in the right order to
    ! maintain consistency.

    ! remove the particle from the aero_particle_array
    i_part = aero_sorted%bin(i_bin)%entry(i_entry)
    i_part_shifted = aero_particle_array%n_part ! old loc of shifted particle
    !>DEBUG
    !write(*,*) 'asrp: remove i_part/id', i_part, &
    !     aero_particle_array%particle(i_part)%id
    !<DEBUG
    call aero_particle_array_remove_particle(aero_particle_array, i_part)

    if (i_part_shifted /= i_part) then
       ! fix up the forward index for the shifted particle
       i_bin_fix = aero_sorted%reverse_bin%entry(i_part_shifted)
       i_part_fix = aero_sorted%reverse_entry%entry(i_part_shifted)
       aero_sorted%bin(i_bin_fix)%entry(i_part_fix) = i_part
    end if

    ! remove the particle from the reverse index (with the side effect
    ! of fixing the reverse map for the shifted particle)
    call integer_array_remove_entry(aero_sorted%reverse_bin, i_part)
    call integer_array_remove_entry(aero_sorted%reverse_entry, i_part)

    ! remove the forward index entry
    i_entry_fix = aero_sorted%bin(i_bin)%n_entry
    i_part_fix = aero_sorted%bin(i_bin)%entry(i_entry_fix)
    call integer_array_remove_entry(aero_sorted%bin(i_bin), i_entry)

    if (i_entry_fix /= i_entry) then
       ! fix reverse index
       aero_sorted%reverse_entry%entry(i_part_fix) = i_entry
    end if

  end subroutine aero_sorted_remove_particle

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
end module pmc_aero_sorted
