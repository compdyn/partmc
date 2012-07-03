! Copyright (C) 2011, 2012 Nicole Riemer and Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_aero_sorted module.

!> The aero_sorted_t structure and assocated subroutines.
module pmc_aero_sorted

  use pmc_integer_varray
  use pmc_integer_rmap
  use pmc_integer_rmap2
  use pmc_aero_particle
  use pmc_aero_data
  use pmc_aero_particle_array
  use pmc_bin_grid
  use pmc_mpi

  !> Sorting of particles into bins.
  !!
  !! Two different bin-sortings are maintained, one per size bin and
  !! weight class, and the other per weight group and weight class.
  !!
  !! A particle can thus be identified by its position \c i_part in an
  !! \c aero_particle_array_t, or by an entry in one of the two
  !! sortings.
  !!
  !! For example, for size bin \c i_bin and weight class \c i_class,
  !! the number of particles with this size and class are
  !! <pre>
  !! n = integer_varray_n_entry(aero_sorted%size_class%inverse(i_bin, i_class))
  !! </pre>
  !! For particle number \c i_entry in this size/class bin, the
  !! particle number is
  !! <pre>
  !! i_part = aero_sorted%%size_class%%inverse(i_bin, i_class)%%entry(i_entry)
  !! </pre>
  !! For particle number \c i_part, the size bin and weight class are
  !! <pre>
  !! i_bin = aero_sorted%%size_class%%forward1%%entry(i_part)
  !! i_class = aero_sorted%%size_class%%forward2%%entry(i_part)
  !! </pre>
  !!
  !! Similar relationships hold for \c aero_sorted%%group_class which
  !! sorts particles per weight group/class.
  type aero_sorted_t
     !> Bin grid for sorting.
     type(bin_grid_t) :: bin_grid
     !> Map of size bin and weight class numbers.
     type(integer_rmap2_t) :: size_class
     !> Map of weight group and weight class numbers.
     type(integer_rmap2_t) :: group_class
     !> Whether coagulation kernel bounds are valid.
     logical :: coag_kernel_bounds_valid
     !> Coagulation kernel lower bound [<tt>bin_grid_size(bin_grid) x
     !> bin_grid_size(bin_grid)</tt>].
     real(kind=dp), allocatable, dimension(:,:) :: coag_kernel_min
     !> Coagulation kernel upper bound [<tt>bin_grid_size(bin_grid) x
     !> bin_grid_size(bin_grid)</tt>].
     real(kind=dp), allocatable, dimension(:,:) :: coag_kernel_max
     !> Whether particle removal rate bounds are valid.
     logical :: removal_rate_bounds_valid
     !> Particle removal rate upper bound [<tt>bin_grid_size(bin_grid)</tt>].
     real(kind=dp), allocatable, dimension(:) :: removal_rate_max
  end type aero_sorted_t

  !> How many size bins to use per decade of particle radius.
  real(kind=dp), parameter :: AERO_SORTED_BINS_PER_DECADE = 10d0
  !> Factor to extend size grid beyond largest/smallest particles.
  real(kind=dp), parameter :: AERO_SORTED_BIN_OVER_FACTOR = 10d0
  !> Size grid extension factor when we should regenerate grid.
  real(kind=dp), parameter :: AERO_SORTED_BIN_SAFETY_FACTOR = 3d0

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the number of size bins.
  integer function aero_sorted_n_bin(aero_sorted)

    !> Aerosol sorting to use.
    type(aero_sorted_t), intent(in) :: aero_sorted

    aero_sorted_n_bin = size(aero_sorted%size_class%inverse, 1)

  end function aero_sorted_n_bin

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the number of weight groups.
  integer function aero_sorted_n_group(aero_sorted)

    !> Aerosol sorting to use.
    type(aero_sorted_t), intent(in) :: aero_sorted

    aero_sorted_n_group = size(aero_sorted%group_class%inverse, 1)

  end function aero_sorted_n_group

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the number of weight classes.
  integer function aero_sorted_n_class(aero_sorted)

    !> Aerosol sorting to use.
    type(aero_sorted_t), intent(in) :: aero_sorted

    aero_sorted_n_class = size(aero_sorted%size_class%inverse, 2)

  end function aero_sorted_n_class

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Set the bin grid to be used for sorting.
  subroutine aero_sorted_set_bin_grid(aero_sorted, bin_grid, n_group, n_class)

    !> Aerosol sorted.
    type(aero_sorted_t), intent(inout) :: aero_sorted
    !> Bin grid.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Number of weight groups.
    integer, intent(in) :: n_group
    !> Number of weight classes.
    integer, intent(in) :: n_class

    integer :: n_bin

    n_bin = bin_grid_size(bin_grid)
    call integer_rmap2_set_ranges(aero_sorted%size_class, n_bin, n_class)
    call integer_rmap2_set_ranges(aero_sorted%group_class, n_group, n_class)
    aero_sorted%coag_kernel_bounds_valid = .false.
    if (allocated(aero_sorted%coag_kernel_min)) then
       deallocate(aero_sorted%coag_kernel_min)
    end if
    allocate(aero_sorted%coag_kernel_min(n_bin, n_bin))
    if (allocated(aero_sorted%coag_kernel_max)) then
       deallocate(aero_sorted%coag_kernel_max)
    end if
    allocate(aero_sorted%coag_kernel_max(n_bin, n_bin))
    aero_sorted%removal_rate_bounds_valid = .false.
    if (allocated(aero_sorted%removal_rate_max)) then
       deallocate(aero_sorted%removal_rate_max)
    end if
    allocate(aero_sorted%removal_rate_max(n_bin))
    aero_sorted%bin_grid = bin_grid

  end subroutine aero_sorted_set_bin_grid

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Discard particles that don't fit the bin grid.
  subroutine aero_sorted_discard_outside_grid(aero_sorted, &
       aero_particle_array, aero_data)

    !> Aerosol sorted.
    type(aero_sorted_t), intent(in) :: aero_sorted
    !> Aerosol particles to discard from.
    type(aero_particle_array_t), intent(inout) :: aero_particle_array
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data

    integer :: i_part, i_bin

    ! Work backwards so we only shift particles that we've already
    ! tested.
    do i_part = aero_particle_array%n_part,1,-1
       i_bin = aero_sorted_particle_in_bin(aero_sorted, &
            aero_particle_array%particle(i_part), aero_data)
       if ((i_bin < 1) .or. (i_bin > bin_grid_size(aero_sorted%bin_grid))) then
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
    subroutine aero_sorted_sort_particles(aero_sorted, aero_particle_array, &
         aero_data)

    !> Aerosol sorted.
    type(aero_sorted_t), intent(inout) :: aero_sorted
    !> Aerosol particles to sort.
    type(aero_particle_array_t), intent(in) :: aero_particle_array
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data

    integer :: i_part, i_bin, i_group, i_class

    call integer_rmap2_zero(aero_sorted%size_class)
    call integer_rmap2_zero(aero_sorted%group_class)

    do i_part = 1,aero_particle_array%n_part
       i_bin = aero_sorted_particle_in_bin(aero_sorted, &
            aero_particle_array%particle(i_part), aero_data)
       i_group = aero_particle_array%particle(i_part)%weight_group
       i_class = aero_particle_array%particle(i_part)%weight_class
       call integer_rmap2_append(aero_sorted%size_class, i_bin, i_class)
       call integer_rmap2_append(aero_sorted%group_class, i_group, i_class)
    end do

  end subroutine aero_sorted_sort_particles

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Remake a sorting if particles are getting too close to the edges.
  subroutine aero_sorted_remake_if_needed(aero_sorted, aero_particle_array, &
       aero_data, valid_sort, n_group, n_class, bin_grid, all_procs_same)

    !> Aerosol sorted to (possibly) remake.
    type(aero_sorted_t), intent(inout) :: aero_sorted
    !> Aerosol particles to sort.
    type(aero_particle_array_t), intent(inout) :: aero_particle_array
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Whether the given aero_sorted is valid.
    logical, intent(in) :: valid_sort
    !> Number of weight groups.
    integer, optional, intent(in) :: n_group
    !> Number of weight classes.
    integer, optional, intent(in) :: n_class
    !> An optional bin_grid to use for the sort.
    type(bin_grid_t), optional, intent(in) :: bin_grid
    !> Whether all processors should use the same bin grid.
    logical, optional, intent(in) :: all_procs_same

    integer :: i_bin, i_bin_min, i_bin_max, i_part, n_bin, use_n_group
    integer :: use_n_class
    real(kind=dp) :: r, r_min, r_max, grid_r_min, grid_r_max
    real(kind=dp) :: local_r_min, local_r_max
    logical :: need_new_bin_grid
    type(bin_grid_t) :: new_bin_grid

    if (present(n_group)) then
       call assert(267881270, present(n_class))
       use_n_group = n_group
       use_n_class = n_class
    else
       call assert(352582858, valid_sort)
       use_n_group = aero_sorted_n_group(aero_sorted)
       use_n_class = aero_sorted_n_class(aero_sorted)
    end if

    if (present(bin_grid)) then
       call aero_sorted_set_bin_grid(aero_sorted, bin_grid, use_n_group, &
            use_n_class)
       call aero_sorted_discard_outside_grid(aero_sorted, &
            aero_particle_array, aero_data)
       call aero_sorted_sort_particles(aero_sorted, aero_particle_array, &
            aero_data)
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
       do i_bin = 1,bin_grid_size(aero_sorted%bin_grid)
          if (sum(integer_varray_n_entry( &
               aero_sorted%size_class%inverse(i_bin, :))) > 0) then
             if (i_bin_min == 0) then
                i_bin_min = i_bin
             end if
             i_bin_max = i_bin
          end if
       end do

       if (i_bin_min == 0) then
          ! there aren't any particles
          call assert(333430891, i_bin_max == 0)
          if (bin_grid_size(aero_sorted%bin_grid) > 0) then
             ! take r_min = upper edge, etc.
             r_min = aero_sorted%bin_grid%edges( &
                  bin_grid_size(aero_sorted%bin_grid) + 1)
             r_max = aero_sorted%bin_grid%edges(1)
          end if
       else
          r_min = aero_sorted%bin_grid%edges(i_bin_min)
          r_max = aero_sorted%bin_grid%edges(i_bin_max + 1)
       end if
    else
       ! no bin data, need to loop over all particles
       do i_part = 1,aero_particle_array%n_part
          r = aero_particle_radius(aero_particle_array%particle(i_part), &
               aero_data)
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
          local_r_max = r_max
          call pmc_mpi_allreduce_max_real(local_r_max, r_max)
          ! don't contaminate global min with zeros
          if (r_min == 0d0) then
             local_r_min = r_max
          else
             local_r_min = r_min
          end if
          call pmc_mpi_allreduce_min_real(local_r_min, r_min)

          ! check that all the bin grids are really the same
          if (.not. pmc_mpi_allequal_bin_grid(aero_sorted%bin_grid)) then
             need_new_bin_grid = .true.
          end if
       end if
    end if

    ! no particles and no existing useful bin_grid
    if (r_max == 0d0) then
       if (valid_sort) return
       call bin_grid_make(new_bin_grid, BIN_GRID_TYPE_LOG, n_bin=0, min=0d0, &
            max=0d0)
       call aero_sorted_set_bin_grid(aero_sorted, new_bin_grid, use_n_group, &
            use_n_class)
       return
    end if

    if (bin_grid_size(aero_sorted%bin_grid) < 1) then
       need_new_bin_grid = .true.
    else
       grid_r_min = aero_sorted%bin_grid%edges(1)
       grid_r_max = aero_sorted%bin_grid%edges( &
            bin_grid_size(aero_sorted%bin_grid) + 1)

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
       call bin_grid_make(new_bin_grid, BIN_GRID_TYPE_LOG, n_bin, grid_r_min, &
            grid_r_max)
       call aero_sorted_set_bin_grid(aero_sorted, new_bin_grid, use_n_group, &

            use_n_class)
       call aero_sorted_sort_particles(aero_sorted, aero_particle_array, &
            aero_data)
    else
       if (.not. valid_sort) then
          call aero_sorted_sort_particles(aero_sorted, aero_particle_array, &
               aero_data)
       end if
    end if

  end subroutine aero_sorted_remake_if_needed

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Find the bin number that contains a given particle.
  integer function aero_sorted_particle_in_bin(aero_sorted, aero_particle, &
       aero_data)

    !> Aerosol sort.
    type(aero_sorted_t), intent(in) :: aero_sorted
    !> Particle.
    type(aero_particle_t), intent(in) :: aero_particle
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data

    aero_sorted_particle_in_bin = bin_grid_find(aero_sorted%bin_grid, &
         aero_particle_radius(aero_particle, aero_data))

  end function aero_sorted_particle_in_bin

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Add a new particle to both an aero_sorted and the corresponding
  !> aero_particle_array.
  subroutine aero_sorted_add_particle(aero_sorted, aero_particle_array, &
       aero_particle, aero_data, allow_resort)

    !> Sorted particle structure.
    type(aero_sorted_t), intent(inout) :: aero_sorted
    !> Aerosol particles.
    type(aero_particle_array_t), intent(inout) :: aero_particle_array
    !> Particle to add.
    type(aero_particle_t), intent(in) :: aero_particle
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Whether to allow a resort due to the add.
    logical, optional, intent(in) :: allow_resort

    integer :: i_bin, i_group, i_class, n_bin, n_group, n_class

    i_bin = aero_sorted_particle_in_bin(aero_sorted, aero_particle, &
         aero_data)
    i_group = aero_particle%weight_group
    i_class = aero_particle%weight_class

    n_bin = bin_grid_size(aero_sorted%bin_grid)
    n_group = aero_sorted_n_group(aero_sorted)
    n_class = aero_sorted_n_class(aero_sorted)
    call assert(417177855, (i_group >= 1) .and. (i_group <= n_group))
    call assert(233133947, (i_class >= 1) .and. (i_class <= n_class))

    ! add the particle to the aero_particle_array
    call aero_particle_array_add_particle(aero_particle_array, aero_particle)

    if ((i_bin < 1) .or. (i_bin > n_bin)) then
       ! particle doesn't fit in the current bin_grid, so remake the
       ! bin_grid if we are allowed
       ! if bin_grid is unallocated, then i_bin will be -1 thus will remake
       ! the bin_grid.
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
            aero_data, valid_sort=.false., n_group=n_group, n_class=n_class)
    else
       ! particle fits in the current bin_grid
       call integer_rmap2_append(aero_sorted%size_class, i_bin, i_class)
       call integer_rmap2_append(aero_sorted%group_class, i_group, i_class)
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

    ! all of these shift the last item into the newly-empty slot
    call aero_particle_array_remove_particle(aero_particle_array, i_part)
    call integer_rmap2_remove(aero_sorted%size_class, i_part)
    call integer_rmap2_remove(aero_sorted%group_class, i_part)

  end subroutine aero_sorted_remove_particle

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Move a particle to a different bin and group.
  subroutine aero_sorted_move_particle(aero_sorted, i_part, new_bin, &
       new_group, new_class)

    !> Aerosol sorted.
    type(aero_sorted_t), intent(inout) :: aero_sorted
    !> Particle number to move.
    integer, intent(in) :: i_part
    !> New bin to move particle to.
    integer, intent(in) :: new_bin
    !> New weight group to move particle to.
    integer, intent(in) :: new_group
    !> New weight class to move particle to.
    integer, intent(in) :: new_class

    call integer_rmap2_change(aero_sorted%size_class, i_part, new_bin, &
         new_class)
    call integer_rmap2_change(aero_sorted%group_class, i_part, new_group, &
         new_class)

  end subroutine aero_sorted_move_particle

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Check sorting.
  subroutine aero_sorted_check(aero_sorted, aero_particle_array, &
       aero_data, n_group, n_class, continue_on_error)

    !> Aerosol sorted to check.
    type(aero_sorted_t), intent(in) :: aero_sorted
    !> Aerosol particles.
    type(aero_particle_array_t), intent(in) :: aero_particle_array
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Number of weight groups.
    integer, optional, intent(in) :: n_group
    !> Number of weight classes.
    integer, optional, intent(in) :: n_class
    !> Whether to continue despite error.
    logical, intent(in) :: continue_on_error

    integer :: i_part, i_bin

    call integer_rmap2_check(aero_sorted%size_class, "size_class", &
         n_domain=aero_particle_array%n_part, &
         n_range_1=bin_grid_size(aero_sorted%bin_grid), n_range_2=n_class, &
         continue_on_error=continue_on_error)
    do i_part = 1,aero_particle_array%n_part
       i_bin = aero_sorted_particle_in_bin(aero_sorted, &
            aero_particle_array%particle(i_part), aero_data)
       if ((i_bin /= aero_sorted%size_class%forward1%entry(i_part)) &
            .or. (i_bin /= aero_sorted%size_class%forward1%entry(i_part))) then
          write(0,*) 'ERROR aero_sorted A: ', "size_class"
          write(0,*) 'i_part', i_part
          write(0,*) 'i_bin', i_bin
          write(0,*) 'aero_sorted%size_class%forward1%entry(i_part)', &
               aero_sorted%size_class%forward1%entry(i_part)
          write(0,*) 'aero_sorted%size_class%forward2%entry(i_part)', &
               aero_sorted%size_class%forward2%entry(i_part)
          call assert(565030916, continue_on_error)
       end if
    end do

    call integer_rmap2_check(aero_sorted%group_class, "group_class", &
         n_domain=aero_particle_array%n_part, &
         n_range_1=n_group, n_range_2=n_class, &
         continue_on_error=continue_on_error)
    do i_part = 1,aero_particle_array%n_part
       if ((aero_particle_array%particle(i_part)%weight_group &
            /= aero_sorted%group_class%forward1%entry(i_part)) &
            .or. (aero_particle_array%particle(i_part)%weight_class &
            /= aero_sorted%group_class%forward2%entry(i_part))) then
          write(0,*) 'ERROR aero_sorted B: ', "group_class"
          write(0,*) 'i_part', i_part
          write(0,*) 'aero_particle_array%particle(i_part)%weight_group', &
               aero_particle_array%particle(i_part)%weight_group
          write(0,*) 'aero_particle_array%particle(i_part)%weight_class', &
               aero_particle_array%particle(i_part)%weight_class
          write(0,*) 'aero_sorted%group_class%forward1%entry(i_part)', &
               aero_sorted%group_class%forward1%entry(i_part)
          write(0,*) 'aero_sorted%group_class%forward2%entry(i_part)', &
               aero_sorted%group_class%forward2%entry(i_part)
          call assert(803595412, continue_on_error)
       end if
    end do

  end subroutine aero_sorted_check

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Determines the number of bytes required to pack the given value.
  integer function pmc_mpi_pack_size_aero_sorted(val)

    !> Value to pack.
    type(aero_sorted_t), intent(in) :: val

    integer :: total_size

    total_size = 0
    total_size = total_size + pmc_mpi_pack_size_bin_grid(val%bin_grid)
    total_size = total_size + pmc_mpi_pack_size_integer_rmap2(val%size_class)
    total_size = total_size + pmc_mpi_pack_size_integer_rmap2(val%group_class)
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
    integer :: prev_position

    prev_position = position
    call pmc_mpi_pack_bin_grid(buffer, position, val%bin_grid)
    call pmc_mpi_pack_integer_rmap2(buffer, position, val%size_class)
    call pmc_mpi_pack_integer_rmap2(buffer, position, val%group_class)
    call assert(786981367, &
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
    integer :: prev_position, n_bin, n_group, n_class

    prev_position = position
    call pmc_mpi_unpack_bin_grid(buffer, position, val%bin_grid)
    call pmc_mpi_unpack_integer_rmap2(buffer, position, val%size_class)
    call pmc_mpi_unpack_integer_rmap2(buffer, position, val%group_class)
    call assert(703866072, &
         position - prev_position <= pmc_mpi_pack_size_aero_sorted(val))
#endif

  end subroutine pmc_mpi_unpack_aero_sorted

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_aero_sorted
