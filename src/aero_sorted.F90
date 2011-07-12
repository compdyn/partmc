! Copyright (C) 2011 Nicole Riemer and Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_aero_sorted module.

!> The aero_sorted_t structure and assocated subroutines.
module pmc_aero_sorted

  use pmc_integer_varray
  use pmc_integer_rmap
  use pmc_aero_particle
  use pmc_aero_particle_array
  use pmc_bin_grid
  use pmc_mpi

  !> Bin index for particles sorted into bins.
  !!
  !! Both forward and reverse indexes are maintained. Particles are
  !! stored with both a linear index \c i_part, and binned indexes
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
     !> Bin grid for sorting.
     type(bin_grid_t) :: bin_grid
     !> Map of size bin numbers.
     type(integer_rmap_t) :: size
     !> Map of weight group numbers.
     type(integer_rmap_t) :: weight
     !> Reverse index array to bin numbers.
     type(integer_varray_t) :: reverse_bin
     !> Reverse index array to group numbers.
     type(integer_varray_t) :: reverse_group
     !> Array of integer arrays, one per size.
     type(integer_varray_t), allocatable, dimension(:) :: unif_bin
     !> Reverse index array to particle entry-in-unif_bin numbers.
     type(integer_varray_t) :: reverse_unif_entry
     !> Particle indices per weight group.
     type(integer_varray_t), allocatable, dimension(:) :: group
     !> Reverse index array to particle group numbers.
     type(integer_varray_t) :: reverse_group_entry
     !> Whether coagulation kernel bounds are valid.
     logical :: coag_kernel_bounds_valid
     !> Coagulation kernel lower bound.
     real(kind=dp), allocatable, dimension(:,:) :: coag_kernel_min
     !> Coagulation kernel upper bound.
     real(kind=dp), allocatable, dimension(:,:) :: coag_kernel_max
  end type aero_sorted_t

  !> How many size bins to use per decade of particle radius.
  real(kind=dp), parameter :: AERO_SORTED_BINS_PER_DECADE = 10d0
  !> Factor to extend size grid beyond largest/smallest particles.
  real(kind=dp), parameter :: AERO_SORTED_BIN_OVER_FACTOR = 10d0
  !> Size grid extension factor when we should regenerate grid.
  real(kind=dp), parameter :: AERO_SORTED_BIN_SAFETY_FACTOR = 3d0

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Allocate an empty structure.
  subroutine aero_sorted_allocate(aero_sorted)

    !> Structure to initialize.
    type(aero_sorted_t), intent(out) :: aero_sorted

    call bin_grid_allocate(aero_sorted%bin_grid)
    call integer_rmap_allocate(aero_sorted%size)
    call integer_rmap_allocate(aero_sorted%weight)
    call integer_varray_allocate(aero_sorted%reverse_bin)
    call integer_varray_allocate(aero_sorted%reverse_group)
    allocate(aero_sorted%unif_bin(0))
    call integer_varray_allocate(aero_sorted%reverse_unif_entry)
    allocate(aero_sorted%group(0))
    call integer_varray_allocate(aero_sorted%reverse_group_entry)
    aero_sorted%coag_kernel_bounds_valid = .false.
    allocate(aero_sorted%coag_kernel_min(0,0))
    allocate(aero_sorted%coag_kernel_max(0,0))

  end subroutine aero_sorted_allocate
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Allocate a strcture with the given size.
  subroutine aero_sorted_allocate_size(aero_sorted, n_bin, n_group)

    !> Structure to initialize.
    type(aero_sorted_t), intent(out) :: aero_sorted
    !> Number of bins.
    integer, intent(in) :: n_bin
    !> Number of weight groups.
    integer, intent(in) :: n_group

    call bin_grid_allocate_size(aero_sorted%bin_grid, n_bin)
    call integer_rmap_allocate_size(aero_sorted%size, n_bin)
    call integer_rmap_allocate_size(aero_sorted%weight, n_group)
    call integer_varray_allocate(aero_sorted%reverse_bin)
    call integer_varray_allocate(aero_sorted%reverse_group)
    allocate(aero_sorted%unif_bin(n_bin))
    call integer_varray_allocate(aero_sorted%unif_bin)
    call integer_varray_allocate(aero_sorted%reverse_unif_entry)
    allocate(aero_sorted%group(n_group))
    call integer_varray_allocate(aero_sorted%group)
    call integer_varray_allocate(aero_sorted%reverse_group_entry)
    aero_sorted%coag_kernel_bounds_valid = .false.
    allocate(aero_sorted%coag_kernel_min(n_bin, n_bin))
    allocate(aero_sorted%coag_kernel_max(n_bin, n_bin))

  end subroutine aero_sorted_allocate_size
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Deallocates a previously allocated structure.
  subroutine aero_sorted_deallocate(aero_sorted)

    !> Structure to deallocate.
    type(aero_sorted_t), intent(inout) :: aero_sorted

    call bin_grid_deallocate(aero_sorted%bin_grid)
    call integer_rmap_deallocate(aero_sorted%size)
    call integer_rmap_deallocate(aero_sorted%weight)
    call integer_varray_deallocate(aero_sorted%reverse_bin)
    call integer_varray_deallocate(aero_sorted%reverse_group)
    call integer_varray_deallocate(aero_sorted%unif_bin)
    deallocate(aero_sorted%unif_bin)
    call integer_varray_deallocate(aero_sorted%reverse_unif_entry)
    call integer_varray_deallocate(aero_sorted%group)
    deallocate(aero_sorted%group)
    call integer_varray_deallocate(aero_sorted%reverse_group_entry)
    aero_sorted%coag_kernel_bounds_valid = .false.
    deallocate(aero_sorted%coag_kernel_min)
    deallocate(aero_sorted%coag_kernel_max)

  end subroutine aero_sorted_deallocate
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Resets an aero_sorted to have zero particles per bin.
  subroutine aero_sorted_zero(aero_sorted)

    !> Structure to zero.
    type(aero_sorted_t), intent(inout) :: aero_sorted

    call integer_rmap_zero(aero_sorted%size)
    call integer_rmap_zero(aero_sorted%weight)
    call integer_varray_zero(aero_sorted%reverse_bin)
    call integer_varray_zero(aero_sorted%reverse_group)
    call integer_varray_zero(aero_sorted%unif_bin)
    call integer_varray_zero(aero_sorted%reverse_unif_entry)
    call integer_varray_zero(aero_sorted%group)
    call integer_varray_zero(aero_sorted%reverse_group_entry)
    aero_sorted%coag_kernel_bounds_valid = .false.
    aero_sorted%coag_kernel_min = 0d0
    aero_sorted%coag_kernel_max = 0d0

  end subroutine aero_sorted_zero
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Do a sorting of a set of aerosol particles.
  subroutine aero_sorted_set_bin_grid(aero_sorted, bin_grid, n_group)

    !> Aerosol sorted.
    type(aero_sorted_t), intent(inout) :: aero_sorted
    !> Bin grid.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Number of weight groups.
    integer, optional, intent(in) :: n_group

    integer :: use_n_group

    if (present(n_group)) then
       use_n_group = n_group
    else
       use_n_group = size(aero_sorted%group)
    end if
    call aero_sorted_deallocate(aero_sorted)
    call aero_sorted_allocate_size(aero_sorted, bin_grid%n_bin, use_n_group)
    call bin_grid_copy(bin_grid, aero_sorted%bin_grid)

  end subroutine aero_sorted_set_bin_grid

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Discard particles that don't fit the bin grid.
  subroutine aero_sorted_discard_outside_grid(aero_sorted, &
       aero_particle_array)

    !> Aerosol sorted.
    type(aero_sorted_t), intent(in) :: aero_sorted
    !> Aerosol particles to discard from.
    type(aero_particle_array_t), intent(inout) :: aero_particle_array
  
    integer :: i_part, i_bin

    ! Work backwards so we only shift particles that we've already
    ! tested.
    do i_part = aero_particle_array%n_part,1,-1
       i_bin = aero_sorted_particle_in_bin(aero_sorted, &
            aero_particle_array%particle(i_part))
       if ((i_bin < 1) .or. (i_bin > aero_sorted%bin_grid%n_bin)) then
          call warn_msg(954800836, "particle ID " &
               // trim(integer_to_string( &
               aero_particle_array%particle(i_part)%id)) &
               // " outside of bin_grid, discarding")
          call aero_particle_array_remove_particle(aero_particle_array, &
               i_part)
       end if
    end do

  end subroutine aero_sorted_discard_outside_grid

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !> Sort the particles.
    subroutine aero_sorted_sort_particles(aero_sorted, aero_particle_array)

    !> Aerosol sorted.
    type(aero_sorted_t), intent(inout) :: aero_sorted
    !> Aerosol particles to sort.
    type(aero_particle_array_t), intent(in) :: aero_particle_array

    integer :: i_part, i_bin, i_group

    call integer_rmap_zero(aero_sorted%size)
    call integer_rmap_zero(aero_sorted%weight)

    call integer_varray_zero(aero_sorted%reverse_bin)
    call integer_varray_zero(aero_sorted%reverse_group)
    call integer_varray_zero(aero_sorted%unif_bin)
    call integer_varray_zero(aero_sorted%reverse_unif_entry)
    call integer_varray_zero(aero_sorted%group)

    call assert(427582120, &
         size(aero_sorted%unif_bin) == aero_sorted%bin_grid%n_bin)

    do i_part = 1,aero_particle_array%n_part
       i_bin = aero_sorted_particle_in_bin(aero_sorted, &
            aero_particle_array%particle(i_part))
       call integer_rmap_append(aero_sorted%size, i_bin)

       call assert(754810952, i_bin >= 1)
       call assert(811054361, i_bin <= aero_sorted%bin_grid%n_bin)

       i_group = aero_particle_array%particle(i_part)%weight_group
       call integer_rmap_append(aero_sorted%weight, i_group)

       call assert(288643748, i_group >= 1)
       call assert(882759457, i_group <= size(aero_sorted%group))

       ! fill in reverse index for size/group bins
       call integer_varray_append(aero_sorted%reverse_bin, i_bin)
       call integer_varray_append(aero_sorted%reverse_group, i_group)

       ! fill in forward index for size-only bins
       call integer_varray_append(aero_sorted%unif_bin(i_bin), i_part)
       
       ! fill in reverse index for size-only bins
       call integer_varray_append(aero_sorted%reverse_unif_entry, &
            aero_sorted%unif_bin(i_bin)%n_entry)

       ! fill in forward index for group
       call integer_varray_append(aero_sorted%group(i_group), i_part)

       ! fill in reverse index for group entry
       call integer_varray_append(aero_sorted%reverse_group_entry, &
            aero_sorted%group(i_group)%n_entry)
    end do

  end subroutine aero_sorted_sort_particles

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Remake a sorting if particles are getting too close to the edges.
  subroutine aero_sorted_remake_if_needed(aero_sorted, aero_particle_array, &
       valid_sort, n_group, bin_grid, all_procs_same)

    !> Aerosol sorted to (possibly) remake.
    type(aero_sorted_t), intent(inout) :: aero_sorted
    !> Aerosol particles to sort.
    type(aero_particle_array_t), intent(inout) :: aero_particle_array
    !> Whether the given aero_sorted is valid.
    logical, intent(in) :: valid_sort
    !> Number of weight groups.
    integer, optional, intent(in) :: n_group
    !> An optional bin_grid to use for the sort.
    type(bin_grid_t), optional, intent(in) :: bin_grid
    !> Whether all processors should use the same bin grid.
    logical, optional, intent(in) :: all_procs_same

    integer :: i_bin, i_bin_min, i_bin_max, i_part, n_bin
    real(kind=dp) :: r, r_min, r_max, grid_r_min, grid_r_max
    real(kind=dp) :: local_r_min, local_r_max
    logical :: need_new_bin_grid
    type(bin_grid_t) :: new_bin_grid

    call assert(886415045, present(n_group) .or. valid_sort)

    if (present(bin_grid)) then
       call aero_sorted_set_bin_grid(aero_sorted, bin_grid, n_group)
       call aero_sorted_discard_outside_grid(aero_sorted, aero_particle_array)
       call aero_sorted_sort_particles(aero_sorted, aero_particle_array)
       return
    end if

    if (aero_particle_array%n_part == 0) then
       call assert(274242189, pmc_mpi_size() == 1)
       ! FIXME: this breaks on MPI: what if some procs have no
       ! particles and some do?
       call bin_grid_allocate(new_bin_grid)
       call bin_grid_make(new_bin_grid, n_bin=0, r_min=0d0, r_max=0d0)
       call aero_sorted_set_bin_grid(aero_sorted, new_bin_grid, n_group)
       call bin_grid_deallocate(new_bin_grid)
       return
    end if

    need_new_bin_grid = .false.

    ! determine r_min and r_max
    r_min = 0d0
    r_max = 0d0
    if (valid_sort) then
       ! use bin data to avoid looping over all particles
       i_bin_min = 0
       i_bin_max = 0
       do i_bin = 1,aero_sorted%bin_grid%n_bin
          if (aero_sorted%unif_bin(i_bin)%n_entry > 0) then
             if (i_bin_min == 0) then
                i_bin_min = i_bin
             end if
             i_bin_max = i_bin
          end if
       end do

       if (i_bin_min == 0) then
          ! there are't any particles, take r_min = upper edge, etc.
          call assert(333430891, i_bin_max == 0)
          r_min = aero_sorted%bin_grid%edge_radius( &
               aero_sorted%bin_grid%n_bin + 1)
          r_max = aero_sorted%bin_grid%edge_radius(1)
       else
          r_min = aero_sorted%bin_grid%edge_radius(i_bin_min)
          r_max = aero_sorted%bin_grid%edge_radius(i_bin_max + 1)
       end if
    else
       ! no bin data, need to loop over all particles
       do i_part = 1,aero_particle_array%n_part
          r = aero_particle_radius(aero_particle_array%particle(i_part))
          if (i_part == 1) then
             r_min = r
             r_max = r
          else
             r_min = min(r_min, r)
             r_max = max(r_max, r)
          end if
       end do
    end if

    if (present(all_procs_same)) then
       if (all_procs_same) then
          ! take global min/max
          local_r_min = r_min
          local_r_max = r_max
          call assert(539388373, r_min > 0d0) ! FIXME: not true if some
          call assert(539388373, r_max > 0d0) ! procs have no particles
          call pmc_mpi_allreduce_min_real(local_r_min, r_min)
          call pmc_mpi_allreduce_max_real(local_r_max, r_max)
          
          ! check that all the bin grids are really the same
          if (.not. pmc_mpi_allequal_bin_grid(aero_sorted%bin_grid)) then
             need_new_bin_grid = .true.
          end if
       end if
    end if

    if (aero_sorted%bin_grid%n_bin < 1) then
       need_new_bin_grid = .true.
    else
       grid_r_min = aero_sorted%bin_grid%edge_radius(1)
       grid_r_max &
            = aero_sorted%bin_grid%edge_radius(aero_sorted%bin_grid%n_bin + 1)
       
       ! We don't check to see whether we could make the bin grid
       ! smaller, as there doesn't seem much point. It would be easy
       ! to add if desired.
       if ((r_min / grid_r_min < AERO_SORTED_BIN_SAFETY_FACTOR) &
            .or. (grid_r_max / r_max < AERO_SORTED_BIN_SAFETY_FACTOR)) then
          need_new_bin_grid = .true.
       end if
    end if
       
    if (need_new_bin_grid) then
       grid_r_min = r_min / AERO_SORTED_BIN_OVER_FACTOR
       grid_r_max = r_max * AERO_SORTED_BIN_OVER_FACTOR
       n_bin = ceiling((log10(grid_r_max) - log10(grid_r_min)) &
            * AERO_SORTED_BINS_PER_DECADE)
       call bin_grid_allocate(new_bin_grid)
       call bin_grid_make(new_bin_grid, n_bin, grid_r_min, grid_r_max)
       call aero_sorted_set_bin_grid(aero_sorted, new_bin_grid, n_group)
       call bin_grid_deallocate(new_bin_grid)
       call aero_sorted_sort_particles(aero_sorted, aero_particle_array)
    else
       if (.not. valid_sort) then
          call aero_sorted_sort_particles(aero_sorted, aero_particle_array)
       end if
    end if

  end subroutine aero_sorted_remake_if_needed

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Find the bin number that contains a given particle.
  integer function aero_sorted_particle_in_bin(aero_sorted, aero_particle)

    !> Aerosol sort.
    type(aero_sorted_t), intent(in) :: aero_sorted
    !> Particle.
    type(aero_particle_t), intent(in) :: aero_particle
    
    aero_sorted_particle_in_bin &
         = bin_grid_particle_in_bin(aero_sorted%bin_grid, &
         aero_particle_radius(aero_particle))
    
  end function aero_sorted_particle_in_bin
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Add a new particle to both an aero_sorted and the corresponding
  !> aero_particle_array.
  subroutine aero_sorted_add_particle(aero_sorted, aero_particle_array, &
       aero_particle, n_group, allow_resort)

    !> Sorted particle structure.
    type(aero_sorted_t), intent(inout) :: aero_sorted
    !> Aerosol particles.
    type(aero_particle_array_t), intent(inout) :: aero_particle_array
    !> Particle to add.
    type(aero_particle_t), intent(in) :: aero_particle
    !> Number of weight groups.
    integer, intent(in) :: n_group
    !> Whether to allow a resort due to the add.
    logical, optional, intent(in) :: allow_resort

    integer :: i_bin, i_group

    i_bin = aero_sorted_particle_in_bin(aero_sorted, aero_particle)
    i_group = aero_particle%weight_group

    call assert(894889664, i_group >= 1)
    call assert(517084587, i_group <= size(aero_sorted%group))

    ! add the particle to the aero_particle_array
    call aero_particle_array_add_particle(aero_particle_array, aero_particle)

    if ((i_bin < 1) .or. (i_bin > aero_sorted%bin_grid%n_bin)) then
       ! particle doesn't fit in the current bin_grid, so remake the
       ! bin_grid if we are allowed
       if (present(allow_resort)) then
          if (.not. allow_resort) then
             ! FIXME: this could be avoided if the new bin_grid was an
             ! extension of the old one (only added bins, first bins
             ! are the same)
             call die_msg(134572570, "particle outside of bin_grid: " &
                  // "try reducing the timestep del_t")
          end if
       end if
       call aero_sorted_remake_if_needed(aero_sorted, aero_particle_array, &
            valid_sort=.false., n_group=n_group)
    else
       ! particle fits in the current bin_grid

       call integer_rmap_append(aero_sorted%size, i_bin)
       call integer_rmap_append(aero_sorted%weight, i_group)

       ! update the reverse index for size/group bins
       call integer_varray_append(aero_sorted%reverse_bin, i_bin)
       call integer_varray_append(aero_sorted%reverse_group, i_group)

       ! update the forward index for size-only bins
       call integer_varray_append(aero_sorted%unif_bin(i_bin), &
            aero_particle_array%n_part)
       
       ! update the reverse index for size-only bins
       call integer_varray_append(aero_sorted%reverse_unif_entry, &
            aero_sorted%unif_bin(i_bin)%n_entry)

       ! update the forward index for group
       call integer_varray_append(aero_sorted%group(i_group), &
            aero_particle_array%n_part)

       ! update the reverse index for group entry
       call integer_varray_append(aero_sorted%reverse_group_entry, &
            aero_sorted%group(i_group)%n_entry)
    end if

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

    integer :: i_bin, i_group, i_part_shifted
    integer :: i_bin_fix, i_group_fix, i_part_fix, i_entry_fix
    integer :: i_unif_entry, i_unif_entry_fix, i_unif_part_fix
    integer :: i_group_entry, i_group_entry_fix, i_group_part_fix

    ! Deleting particles shifts the end particles into the empty slots
    ! in the aero_particle_array and the aero_sorted forward and
    ! reverse indexes. All must be fixed in the right order to
    ! maintain consistency.

    i_bin = aero_sorted%reverse_bin%entry(i_part)
    i_group = aero_sorted%reverse_group%entry(i_part)
    i_unif_entry = aero_sorted%reverse_unif_entry%entry(i_part)
    i_group_entry = aero_sorted%reverse_group_entry%entry(i_part)

    ! remove the particle from the aero_particle_array
    i_part_shifted = aero_particle_array%n_part ! old loc of shifted particle
    call aero_particle_array_remove_particle(aero_particle_array, i_part)
    call integer_rmap_remove(aero_sorted%size, i_part)
    call integer_rmap_remove(aero_sorted%weight, i_part)

    if (i_part_shifted /= i_part) then
       ! fix up the forward index for the shifted particle
       i_bin_fix = aero_sorted%reverse_bin%entry(i_part_shifted)
       i_group_fix = aero_sorted%reverse_group%entry(i_part_shifted)

       ! fix up indices for size-only bins
       i_unif_entry_fix = aero_sorted%reverse_unif_entry%entry(i_part_shifted)
       aero_sorted%unif_bin(i_bin_fix)%entry(i_unif_entry_fix) = i_part

       ! fix up indices for group
       i_group_entry_fix &
            = aero_sorted%reverse_group_entry%entry(i_part_shifted)
       aero_sorted%group(i_group_fix)%entry(i_group_entry_fix) = i_part
    end if

    ! remove the particle from the reverse index (with the side effect
    ! of fixing the reverse map for the shifted particle)
    call integer_varray_remove_entry(aero_sorted%reverse_bin, i_part)
    call integer_varray_remove_entry(aero_sorted%reverse_group, i_part)
    call integer_varray_remove_entry(aero_sorted%reverse_unif_entry, i_part)
    call integer_varray_remove_entry(aero_sorted%reverse_group_entry, i_part)

    ! remove the forward index entry for size-only bins
    i_unif_entry_fix = aero_sorted%unif_bin(i_bin)%n_entry
    i_unif_part_fix = aero_sorted%unif_bin(i_bin)%entry(i_unif_entry_fix)
    call integer_varray_remove_entry(aero_sorted%unif_bin(i_bin), i_unif_entry)

    if (i_unif_entry_fix /= i_unif_entry) then
       ! fix reverse index
       aero_sorted%reverse_unif_entry%entry(i_unif_part_fix) = i_unif_entry
    end if

    ! remove the forward index entry for group
    i_group_entry_fix = aero_sorted%group(i_group)%n_entry
    i_group_part_fix = aero_sorted%group(i_group)%entry(i_group_entry_fix)
    call integer_varray_remove_entry(aero_sorted%group(i_group), i_group_entry)

    if (i_group_entry_fix /= i_group_entry) then
       ! fix reverse index
       aero_sorted%reverse_group_entry%entry(i_group_part_fix) = i_group_entry
    end if

  end subroutine aero_sorted_remove_particle

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Move a particle to a different bin and group.
  subroutine aero_sorted_move_particle(aero_sorted, i_part, new_bin, new_group)

    !> Aerosol sorted.
    type(aero_sorted_t), intent(inout) :: aero_sorted
    !> Particle number to move.
    integer, intent(in) :: i_part
    !> New bin to move particle to.
    integer, intent(in) :: new_bin
    !> New group to move particle to.
    integer, intent(in) :: new_group

    integer :: i_bin, i_group
    integer :: i_unif_entry, i_unif_part_shifted, new_unif_entry
    integer :: i_group_entry, i_group_part_shifted, new_group_entry

    call integer_rmap_change(aero_sorted%size, i_part, new_bin)
    call integer_rmap_change(aero_sorted%weight, i_part, new_group)

    i_bin = aero_sorted%reverse_bin%entry(i_part)
    i_group = aero_sorted%reverse_group%entry(i_part)
    i_unif_entry = aero_sorted%reverse_unif_entry%entry(i_part)
    i_group_entry = aero_sorted%reverse_group_entry%entry(i_part)
    !if ((i_bin == new_bin) .and. (i_group == new_group)) return

    if (i_bin /= new_bin) then
       ! remove the old forward maps
       call integer_varray_remove_entry(aero_sorted%unif_bin(i_bin), i_unif_entry)

       ! fix the reverse entry map for the last entry moved into the new slot
       if (i_unif_entry <= aero_sorted%unif_bin(i_bin)%n_entry) then
          i_unif_part_shifted = aero_sorted%unif_bin(i_bin)%entry(i_unif_entry)
          aero_sorted%reverse_unif_entry%entry(i_unif_part_shifted) = i_unif_entry
       end if

       ! add the new forward map
       call integer_varray_append(aero_sorted%unif_bin(new_bin), i_part)
       new_unif_entry = aero_sorted%unif_bin(new_bin)%n_entry

       ! fix the reverse maps
       aero_sorted%reverse_bin%entry(i_part) = new_bin
       aero_sorted%reverse_unif_entry%entry(i_part) = new_unif_entry
    end if


    if (i_group /= new_group) then
       ! remove the old forward maps
       call integer_varray_remove_entry(aero_sorted%group(i_group), i_group_entry)
       
       ! fix the reverse entry map for the last entry moved into the new slot
       if (i_group_entry <= aero_sorted%group(i_group)%n_entry) then
          i_group_part_shifted = aero_sorted%group(i_group)%entry(i_group_entry)
          aero_sorted%reverse_group_entry%entry(i_group_part_shifted) &
               = i_group_entry
       end if
       
       ! add the new forward map
       call integer_varray_append(aero_sorted%group(new_group), i_part)
       new_group_entry = aero_sorted%group(new_group)%n_entry
       
       ! fix the reverse maps
       aero_sorted%reverse_group%entry(i_part) = new_group
       aero_sorted%reverse_group_entry%entry(i_part) = new_group_entry
    end if

  end subroutine aero_sorted_move_particle

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Check that the sorted data is consistent.
  subroutine aero_sorted_check_base(name, n_domain, n_range, rmap, map, &
       index, continue_on_error)

    !> Check name.
    character(len=*), intent(in) :: name
    !> Number of domain items.
    integer, intent(in) :: n_domain
    !> Number of image items.
    integer, intent(in) :: n_range
    !> Reverse map from image to domain (multi-valued).
    type(integer_varray_t), intent(in) :: rmap(:)
    !> Forward map from domain to image (single valued).
    type(integer_varray_t), intent(in) :: map
    !> Forward map from domain to indexes (single valued).
    type(integer_varray_t), intent(in) :: index
    !> Whether to continue despite error.
    logical, intent(in) :: continue_on_error

    integer :: i_domain, i_range, i_index

    if ((n_domain /= map%n_entry) &
         .or. (n_domain /= index%n_entry) &
         .or. (n_range /= size(rmap))) then
       write(0,*) 'SORT CHECK ERROR A:', name
       write(0,*) 'n_domain', n_domain
       write(0,*) 'n_range', n_range
       write(0,*) 'map%n_entry', map%n_entry
       write(0,*) 'index%n_entry', index%n_entry
       write(0,*) 'size(rmap)', size(rmap)
       call assert(973643016, continue_on_error)
    end if

    do i_domain = 1,n_domain
       i_range = map%entry(i_domain)
       if ((i_range < 1) .or. (i_range > n_range)) then
          write(0,*) 'SORT CHECK ERROR B:', name
          write(0,*) 'i_domain', i_domain
          write(0,*) 'i_range', i_range
          write(0,*) 'n_range', n_range
          call assert(798857945, continue_on_error)
       end if

       i_index = index%entry(i_domain)
       if ((i_index < 1) .or. (i_index > rmap(i_range)%n_entry)) then
          write(0,*) 'SORT CHECK ERROR C:', name
          write(0,*) 'i_domain', i_domain
          write(0,*) 'i_range', i_range
          write(0,*) 'i_index', i_index
          write(0,*) 'rmap(i_range)%n_entry', rmap(i_range)%n_entry
          call assert(823748734, continue_on_error)
       end if
       if (i_domain /= rmap(i_range)%entry(i_index)) then
          write(0,*) 'SORT CHECK ERROR D:', name
          write(0,*) 'i_domain', i_domain
          write(0,*) 'i_range', i_range
          write(0,*) 'i_index', i_index
          write(0,*) 'rmap(i_range)%entry(i_index)', &
               rmap(i_range)%entry(i_index)
          call assert(735205557, continue_on_error)
       end if
    end do

    do i_range = 1,n_range
       do i_index = 1,rmap(i_range)%n_entry
          i_domain = rmap(i_range)%entry(i_index)
          if ((i_domain < 1) .or. (i_domain > n_domain)) then
             write(0,*) 'SORT CHECK ERROR E:', name
             write(0,*) 'i_range', i_range
             write(0,*) 'i_index', i_index
             write(0,*) 'i_domain', i_domain
             write(0,*) 'n_domain', n_domain
             call assert(502643520, continue_on_error)
          end if
          if ((i_range /= map%entry(i_domain)) &
               .or. (i_index /= index%entry(i_domain))) then
             write(0,*) 'SORT CHECK ERROR F:', name
             write(0,*) 'i_domain', i_domain
             write(0,*) 'i_range', i_range
             write(0,*) 'map%entry(i_domain)', map%entry(i_domain)
             write(0,*) 'i_index', i_index
             write(0,*) 'index%entry(i_domain)', index%entry(i_domain)
             call assert(544747928, continue_on_error)
          end if
       end do
    end do

  end subroutine aero_sorted_check_base

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Check that the sorted data is consistent.
  subroutine aero_sorted_check_rmap(integer_rmap, name, &
       rmap, map, index)

    !> Integer rmap to check.
    type(integer_rmap_t), intent(in) :: integer_rmap
    !> Check name.
    character(len=*), intent(in) :: name
    !> Reverse map from image to domain (multi-valued).
    type(integer_varray_t), intent(in) :: rmap(:)
    !> Forward map from domain to image (single valued).
    type(integer_varray_t), intent(in) :: map
    !> Forward map from domain to indexes (single valued).
    type(integer_varray_t), intent(in) :: index

    logical, parameter :: writeout = .false.

    integer :: i_domain, i_range, i_index

    if (writeout) then
       write(*,*) '##############################################################'
       write(*,*) name
       write(*,*) 'integer_rmap%forward%n_entry', integer_rmap%forward%n_entry
       write(*,*) 'map%n_entry', map%n_entry
    end if
    call assert(721760363, integer_rmap%forward%n_entry == map%n_entry)
    if (writeout) then
       write(*,*) 'integer_rmap%index%n_entry', integer_rmap%index%n_entry
       write(*,*) 'index%n_entry', index%n_entry
    end if
    call assert(499883188, integer_rmap%index%n_entry == index%n_entry)
    do i_domain = 1,map%n_entry
       if (writeout) then
          write(*,*) 'i_domain', i_domain
          write(*,*) 'integer_rmap%forward%entry(i_domain)', integer_rmap%forward%entry(i_domain)
          write(*,*) 'map%entry(i_domain)', map%entry(i_domain)
       end if
       call assert(678255436, &
            integer_rmap%forward%entry(i_domain) == map%entry(i_domain))
       if (writeout .or. (integer_rmap%index%entry(i_domain) /= index%entry(i_domain))) then
          write(*,*) 'integer_rmap%index%entry(i_domain)', integer_rmap%index%entry(i_domain)
          write(*,*) 'index%entry(i_domain)', index%entry(i_domain)
       end if
       call assert(313108040, &
            integer_rmap%index%entry(i_domain) == index%entry(i_domain))
    end do
    if (writeout) then
       write(*,*) 'size(integer_rmap%inverse)', size(integer_rmap%inverse)
       write(*,*) 'size(rmap)', size(rmap)
    end if
    call assert(130528257, size(integer_rmap%inverse) == size(rmap))
    do i_range = 1,size(rmap)
       if (writeout) then
          write(*,*) 'i_range', i_range
          write(*,*) 'integer_rmap%inverse(i_range)%n_entry', integer_rmap%inverse(i_range)%n_entry
          write(*,*) 'rmap(i_range)%n_entry', rmap(i_range)%n_entry
       end if
       call assert(378702184, &
            integer_rmap%inverse(i_range)%n_entry == rmap(i_range)%n_entry)
       do i_index = 1,rmap(i_range)%n_entry
          if (writeout) then
             write(*,*) 'i_index', i_index
             write(*,*) 'integer_rmap%inverse(i_range)%entry(i_index)', integer_rmap%inverse(i_range)%entry(i_index)
             write(*,*) 'rmap(i_range)%entry(i_index)', rmap(i_range)%entry(i_index)
          end if
          call assert(182642405, integer_rmap%inverse(i_range)%entry(i_index) &
               == rmap(i_range)%entry(i_index))
       end do
    end do

  end subroutine aero_sorted_check_rmap

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Check sorting.
  subroutine aero_sorted_check(aero_sorted, aero_particle_array, &
       n_group, continue_on_error)

    !> Aerosol sorted to check.
    type(aero_sorted_t), intent(in) :: aero_sorted
    !> Aerosol particles.
    type(aero_particle_array_t), intent(in) :: aero_particle_array
    !> Number of weight groups.
    integer, optional, intent(in) :: n_group
    !> Whether to continue despite error.
    logical, intent(in) :: continue_on_error

    integer :: i_part, i_bin

    call integer_rmap_check(aero_sorted%size, "size", &
         n_domain=aero_particle_array%n_part, &
         n_range=aero_sorted%bin_grid%n_bin, &
         continue_on_error=continue_on_error)
    do i_part = 1,aero_particle_array%n_part
       i_bin = aero_sorted_particle_in_bin(aero_sorted, &
            aero_particle_array%particle(i_part))
       if (i_bin /= aero_sorted%size%forward%entry(i_part)) then
          write(0,*) 'ERROR aero_sorted A: ', "size"
          write(0,*) 'i_part', i_part
          write(0,*) 'i_bin', i_bin
          write(0,*) 'aero_sorted%size%forward%entry(i_part)', &
               aero_sorted%size%forward%entry(i_part)
          call assert(553067208, continue_on_error)
       end if
    end do

    call integer_rmap_check(aero_sorted%weight, "weight", &
         n_domain=aero_particle_array%n_part, &
         n_range=n_group, &
         continue_on_error=continue_on_error)
    do i_part = 1,aero_particle_array%n_part
       if (aero_particle_array%particle(i_part)%weight_group &
            /= aero_sorted%weight%forward%entry(i_part)) then
          write(0,*) 'ERROR aero_sorted B: ', "group"
          write(0,*) 'i_part', i_part
          write(0,*) 'aero_particle_array%particle(i_part)%weight_group', &
               aero_particle_array%particle(i_part)%weight_group
          write(0,*) 'aero_sorted%weight%forward%entry(i_part)', &
               aero_sorted%weight%forward%entry(i_part)
          call assert(389482223, continue_on_error)
       end if
    end do

  end subroutine aero_sorted_check

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Determines the number of bytes required to pack the given value.
  integer function pmc_mpi_pack_size_aero_sorted(val)

    !> Value to pack.
    type(aero_sorted_t), intent(in) :: val

    integer :: i_bin, i_group, total_size

    total_size = 0
    total_size = total_size + pmc_mpi_pack_size_integer(size(val%unif_bin))
    total_size = total_size + pmc_mpi_pack_size_integer(size(val%group))
    total_size = total_size + pmc_mpi_pack_size_bin_grid(val%bin_grid)
    do i_bin = 1,size(val%unif_bin)
       total_size = total_size &
            + pmc_mpi_pack_size_integer_varray(val%unif_bin(i_bin))
    end do
    do i_group = 1,size(val%group)
       total_size = total_size &
            + pmc_mpi_pack_size_integer_varray(val%group(i_group))
    end do
    total_size = total_size &
         + pmc_mpi_pack_size_integer_varray(val%reverse_bin)
    total_size = total_size &
         + pmc_mpi_pack_size_integer_varray(val%reverse_group)
    total_size = total_size &
         + pmc_mpi_pack_size_integer_varray(val%reverse_unif_entry)
    total_size = total_size &
         + pmc_mpi_pack_size_integer_varray(val%reverse_group_entry)
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
    integer :: prev_position, i_bin, i_group

    prev_position = position
    call pmc_mpi_pack_integer(buffer, position, size(val%unif_bin))
    call pmc_mpi_pack_integer(buffer, position, size(val%group))
    call pmc_mpi_pack_bin_grid(buffer, position, val%bin_grid)
    do i_bin = 1,size(val%unif_bin)
       call pmc_mpi_pack_integer_varray(buffer, position, &
            val%unif_bin(i_bin))
    end do
    do i_group = 1,size(val%group)
       call pmc_mpi_pack_integer_varray(buffer, position, &
            val%group(i_group))
    end do
    call pmc_mpi_pack_integer_varray(buffer, position, val%reverse_bin)
    call pmc_mpi_pack_integer_varray(buffer, position, val%reverse_group)
    call pmc_mpi_pack_integer_varray(buffer, position, val%reverse_unif_entry)
    call pmc_mpi_pack_integer_varray(buffer, position, val%reverse_group_entry)
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
    integer :: prev_position, i_bin, i_group, n_bin, n_group

    prev_position = position
    call pmc_mpi_unpack_integer(buffer, position, n_bin)
    call pmc_mpi_unpack_integer(buffer, position, n_group)
    call aero_sorted_deallocate(val)
    call aero_sorted_allocate_size(val, n_bin, n_group)
    call pmc_mpi_unpack_bin_grid(buffer, position, val%bin_grid)
    do i_bin = 1,size(val%unif_bin)
       call pmc_mpi_unpack_integer_varray(buffer, position, &
            val%unif_bin(i_bin))
    end do
    do i_group = 1,size(val%group)
       call pmc_mpi_unpack_integer_varray(buffer, position, &
            val%group(i_group))
    end do
    call pmc_mpi_unpack_integer_varray(buffer, position, val%reverse_bin)
    call pmc_mpi_unpack_integer_varray(buffer, position, val%reverse_group)
    call pmc_mpi_unpack_integer_varray(buffer, position, &
         val%reverse_unif_entry)
    call pmc_mpi_unpack_integer_varray(buffer, position, &
         val%reverse_group_entry)
    call assert(364064630, &
         position - prev_position <= pmc_mpi_pack_size_aero_sorted(val))
#endif

  end subroutine pmc_mpi_unpack_aero_sorted

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
end module pmc_aero_sorted
