! Copyright (C) 2005-2012 Nicole Riemer and Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_aero_state module.

!> The aero_state_t structure and assocated subroutines.
module pmc_aero_state

  use pmc_aero_particle_array
  use pmc_aero_sorted
  use pmc_integer_varray
  use pmc_bin_grid
  use pmc_aero_data
  use pmc_aero_particle
  use pmc_aero_dist
  use pmc_util
  use pmc_rand
  use pmc_aero_binned
  use pmc_mpi
  use pmc_spec_file
  use pmc_aero_info
  use pmc_aero_info_array
  use pmc_aero_weight
  use pmc_aero_weight_array
#ifdef PMC_USE_MPI
  use mpi
#endif

  !> MPI tag for mixing particles between processes.
  integer, parameter :: AERO_STATE_TAG_MIX     = 4987
  !> MPI tag for gathering between processes.
  integer, parameter :: AERO_STATE_TAG_GATHER  = 4988
  !> MPI tag for scattering between processes.
  integer, parameter :: AERO_STATE_TAG_SCATTER = 4989

  !> Single flat weighting scheme.
  integer, parameter :: AERO_STATE_WEIGHT_NONE = 1
  !> Single flat weighting scheme.
  integer, parameter :: AERO_STATE_WEIGHT_FLAT = 2
  !> Power-law weighting scheme.
  integer, parameter :: AERO_STATE_WEIGHT_POWER = 3
  !> Coupled number/mass weighting scheme.
  integer, parameter :: AERO_STATE_WEIGHT_NUMMASS = 4
  !> Flat weighting by source.
  integer, parameter :: AERO_STATE_WEIGHT_FLAT_SOURCE = 5
  !> Power-law weighting by source.
  integer, parameter :: AERO_STATE_WEIGHT_POWER_SOURCE = 6
  !> Coupled number/mass weighting by source.
  integer, parameter :: AERO_STATE_WEIGHT_NUMMASS_SOURCE = 7

  !> The current collection of aerosol particles.
  !!
  !! The particles in \c aero_state_t are stored in a single flat
  !! array (the \c apa data member), with a sorting into bins possibly
  !! stored in the \c aero_sorted data member.
  !!
  !! Every time we remove particles we keep track of the particle ID
  !! and the action performed in the aero_info_array_t structure. This
  !! is typically cleared each time we output data to disk.
  type aero_state_t
     type(aero_particle_array_t) :: apa
     type(aero_sorted_t) :: aero_sorted
     logical :: valid_sort
     !> Weighting functions.
     type(aero_weight_array_t) :: awa
     !> Ideal number of computational particles.
     real(kind=dp) :: n_part_ideal
     !> Information on removed particles.
     type(aero_info_array_t) :: aero_info_array
  end type aero_state_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Allocates aerosol arrays.
  subroutine aero_state_allocate(aero_state)

    !> Aerosol to initialize.
    type(aero_state_t), intent(out) :: aero_state
    
    call aero_particle_array_allocate(aero_state%apa)
    call aero_sorted_allocate(aero_state%aero_sorted)
    aero_state%valid_sort = .false.
    call aero_weight_array_allocate(aero_state%awa)
    aero_state%n_part_ideal = 0d0
    call aero_info_array_allocate(aero_state%aero_info_array)

  end subroutine aero_state_allocate
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Deallocates a previously allocated aerosol.
  subroutine aero_state_deallocate(aero_state)

    !> Aerosol to deallocate.
    type(aero_state_t), intent(inout) :: aero_state
    
    call aero_particle_array_deallocate(aero_state%apa)
    call aero_sorted_deallocate(aero_state%aero_sorted)
    aero_state%valid_sort = .false.
    call aero_weight_array_deallocate(aero_state%awa)
    call aero_info_array_deallocate(aero_state%aero_info_array)

  end subroutine aero_state_deallocate
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Resets an \c aero_state to an empty state.
  subroutine aero_state_reset(aero_state)

    !> Aerosol to reset.
    type(aero_state_t), intent(inout) :: aero_state

    call aero_state_deallocate(aero_state)
    call aero_state_allocate(aero_state)

  end subroutine aero_state_reset

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Copies aerosol to a destination that has already had
  !> aero_state_allocate() called on it.
  subroutine aero_state_copy(aero_state_from, aero_state_to)

    !> Reference aerosol.
    type(aero_state_t), intent(in) :: aero_state_from
    !> Already allocated.
    type(aero_state_t), intent(inout) :: aero_state_to

    call aero_particle_array_copy(aero_state_from%apa, aero_state_to%apa)
    aero_state_to%valid_sort = .false.
    call aero_state_copy_weight(aero_state_from, aero_state_to)
    aero_state_to%n_part_ideal = aero_state_from%n_part_ideal
    call aero_info_array_copy(aero_state_from%aero_info_array, &
         aero_state_to%aero_info_array)

  end subroutine aero_state_copy
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Copies weighting information for an \c aero_state.
  subroutine aero_state_copy_weight(aero_state_from, aero_state_to)

    !> Reference aerosol.
    type(aero_state_t), intent(in) :: aero_state_from
    !> Already allocated.
    type(aero_state_t), intent(inout) :: aero_state_to

    call aero_weight_array_copy(aero_state_from%awa, aero_state_to%awa)

  end subroutine aero_state_copy_weight

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Sets the weighting functions for an \c aero_state.
  subroutine aero_state_set_weight(aero_state, aero_data, weight_type, &
       exponent)

    !> Aerosol to set the weights on.
    type(aero_state_t), intent(inout) :: aero_state
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Type of weighting scheme to use.
    integer, intent(in) :: weight_type
    !> Exponent for power-law weighting (only used if \c weight_type
    !> is \c AERO_STATE_WEIGHT_POWER).
    real(kind=dp), intent(in), optional :: exponent

    aero_state%valid_sort = .false.
    call aero_weight_array_deallocate(aero_state%awa)
    select case(weight_type)
    case(AERO_STATE_WEIGHT_NONE)
       call aero_weight_array_allocate(aero_state%awa)
    case(AERO_STATE_WEIGHT_FLAT)
       call aero_weight_array_allocate_flat(aero_state%awa, 1)
    case(AERO_STATE_WEIGHT_POWER)
       call assert_msg(656670336, present(exponent), &
            "exponent parameter required for AERO_STATE_WEIGHT_POWER")
       call aero_weight_array_allocate_power(aero_state%awa, 1, exponent)
    case(AERO_STATE_WEIGHT_NUMMASS)
       call aero_weight_array_allocate_nummass(aero_state%awa, 1)
    case(AERO_STATE_WEIGHT_FLAT_SOURCE)
       call aero_weight_array_allocate_flat(aero_state%awa, aero_data%n_source)
    case(AERO_STATE_WEIGHT_POWER_SOURCE)
       call assert_msg(656670336, present(exponent), &
            "exponent parameter required for AERO_STATE_WEIGHT_POWER")
       call aero_weight_array_allocate_power(aero_state%awa, &
            aero_data%n_source, exponent)
    case(AERO_STATE_WEIGHT_NUMMASS_SOURCE)
       call aero_weight_array_allocate_nummass(aero_state%awa, &
            aero_data%n_source)
    case default
       call die_msg(969076992, "unknown weight_type: " &
            // trim(integer_to_string(weight_type)))
    end select

  end subroutine aero_state_set_weight
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the total number of particles in an aerosol distribution.
  integer function aero_state_total_particles(aero_state, i_group, i_set)

    !> Aerosol state.
    type(aero_state_t), intent(in) :: aero_state
    !> Weight group.
    integer, optional, intent(in) :: i_group
    !> Weight set.
    integer, optional, intent(in) :: i_set

    integer :: i_part

    if (present(i_group)) then
       call assert(908743823, present(i_set))
       if (aero_state%valid_sort) then
          aero_state_total_particles &
               = aero_state%aero_sorted%weight%inverse(i_group)%n_entry
       else
          ! FIXME: should we just sort?
          aero_state_total_particles = 0
          do i_part = 1,aero_state%apa%n_part
             if ((aero_state%apa%particle(i_part)%weight_group == i_group) &
                  .and. &
                  (aero_state%apa%particle(i_part)%weight_set == i_set)) then
                aero_state_total_particles = aero_state_total_particles + 1
             end if
          end do
       end if
    else
       aero_state_total_particles = aero_state%apa%n_part
    end if

  end function aero_state_total_particles

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the total number of particles across all processes.
  integer function aero_state_total_particles_all_procs(aero_state, i_group, &
       i_set)

    !> Aerosol state.
    type(aero_state_t), intent(in) :: aero_state
    !> Weight group.
    integer, optional, intent(in) :: i_group
    !> Weight set.
    integer, optional, intent(in) :: i_set

    call pmc_mpi_allreduce_sum_integer(&
         aero_state_total_particles(aero_state, i_group, i_set), &
         aero_state_total_particles_all_procs)

  end function aero_state_total_particles_all_procs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Resets an aero_state to have zero particles per bin. This must
  !> already have had aero_state_allocate() called on it. This
  !> function can be called more than once on the same state.
  subroutine aero_state_zero(aero_state)

    !> State to zero.
    type(aero_state_t), intent(inout) :: aero_state
    
    integer :: i, n_bin

    call aero_particle_array_zero(aero_state%apa)
    aero_state%valid_sort = .false.
    call aero_info_array_zero(aero_state%aero_info_array)

  end subroutine aero_state_zero
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Add the given particle.
  subroutine aero_state_add_particle(aero_state, aero_particle, allow_resort)

    !> Aerosol state.
    type(aero_state_t), intent(inout) :: aero_state
    !> Particle to add.
    type(aero_particle_t), intent(in) :: aero_particle
    !> Whether to allow a resort due to the add.
    logical, optional, intent(in) :: allow_resort

    if (aero_state%valid_sort) then
       call aero_sorted_add_particle(aero_state%aero_sorted, aero_state%apa, &
            aero_particle, size(aero_state%awa%weight), allow_resort)
    else
       call aero_particle_array_add_particle(aero_state%apa, aero_particle)
    end if

  end subroutine aero_state_add_particle

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Remove the given particle without recording it.
  subroutine aero_state_remove_particle_no_info(aero_state, i_part)

    !> Aerosol state.
    type(aero_state_t), intent(inout) :: aero_state
    !> Index of particle to remove.
    integer, intent(in) :: i_part

    if (aero_state%valid_sort) then
       call aero_sorted_remove_particle(aero_state%aero_sorted, &
            aero_state%apa, i_part)
    else
       call aero_particle_array_remove_particle(aero_state%apa, i_part)
    end if

  end subroutine aero_state_remove_particle_no_info

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Remove the given particle and record the removal.
  subroutine aero_state_remove_particle_with_info(aero_state, i_part, &
       aero_info)

    !> Aerosol state.
    type(aero_state_t), intent(inout) :: aero_state
    !> Index of particle to remove.
    integer, intent(in) :: i_part
    !> Removal info.
    type(aero_info_t), intent(in) :: aero_info

    call aero_state_remove_particle_no_info(aero_state, i_part)
    call aero_info_array_add_aero_info(aero_state%aero_info_array, &
         aero_info)

  end subroutine aero_state_remove_particle_with_info

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Remove the given particle and possibly record the removal.
  subroutine aero_state_remove_particle(aero_state, i_part, &
       record_removal, aero_info)

    !> Aerosol state.
    type(aero_state_t), intent(inout) :: aero_state
    !> Index of particle to remove.
    integer, intent(in) :: i_part
    !> Whether to record the removal in the aero_info_array.
    logical, intent(in) :: record_removal
    !> Removal info.
    type(aero_info_t), intent(in) :: aero_info

    if (record_removal) then
       call aero_state_remove_particle_with_info(aero_state, i_part, &
            aero_info)
    else
       call aero_state_remove_particle_no_info(aero_state, i_part)
    end if

  end subroutine aero_state_remove_particle

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Remove a randomly chosen particle from the given bin and return
  !> it.
  subroutine aero_state_remove_rand_particle_from_bin(aero_state, &
       i_bin, aero_particle)

    !> Aerosol state.
    type(aero_state_t), intent(inout) :: aero_state
    !> Bin number to remove particle from.
    integer, intent(in) :: i_bin
    !> Removed particle.
    type(aero_particle_t), intent(inout) :: aero_particle

    integer :: i_entry, i_part

    call assert(742996300, aero_state%valid_sort)
    call assert(392182617, &
         aero_state%aero_sorted%size%inverse(i_bin)%n_entry > 0)
    i_entry = pmc_rand_int(aero_state%aero_sorted%size%inverse(i_bin)%n_entry)
    i_part = aero_state%aero_sorted%size%inverse(i_bin)%entry(i_entry)
    call aero_particle_copy(aero_state%apa%particle(i_part), aero_particle)
    call aero_state_remove_particle_no_info(aero_state, i_part)

  end subroutine aero_state_remove_rand_particle_from_bin

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Add copies or remove a particle, with a given mean number of
  !> resulting particles.
  !!
  !! The particle number \c i_part is either removed, or zero or more
  !! copies are added, with a random number of copies with the given
  !! mean \c n_part_mean.  The final number of particles is either
  !! <tt>floor(n_part_mean)</tt> or <tt>ceiling(n_part_mean)</tt>,
  !! chosen randomly so the mean is \c n_part_mean.
  subroutine aero_state_dup_particle(aero_state, i_part, n_part_mean, &
       random_weight_group)

    !> Aerosol state.
    type(aero_state_t), intent(inout) :: aero_state
    !> Particle number.
    integer, intent(in) :: i_part
    !> Mean number of resulting particles.
    real(kind=dp), intent(in) :: n_part_mean
    !> Whether particle copies should be placed in a randomly chosen
    !> weight group.
    logical, optional, intent(in) :: random_weight_group

    integer :: n_copies, i_dup, new_group
    type(aero_particle_t), pointer :: aero_particle
    type(aero_particle_t) :: new_aero_particle
    type(aero_info_t) :: aero_info

    aero_particle => aero_state%apa%particle(i_part)
    n_copies = prob_round(n_part_mean)
    if (n_copies == 0) then
       call aero_info_allocate(aero_info)
       aero_info%id = aero_particle%id
       aero_info%action = AERO_INFO_WEIGHT
       aero_info%other_id = 0
       call aero_state_remove_particle_with_info(aero_state, &
            i_part, aero_info)
       call aero_info_deallocate(aero_info)
    elseif (n_copies > 1) then
       call aero_particle_allocate(new_aero_particle)
       do i_dup = 1,(n_copies - 1)
          call aero_particle_copy(aero_particle, new_aero_particle)
          call aero_particle_new_id(new_aero_particle)
          if (present(random_weight_group)) then
             if (random_weight_group) then
                new_group &
                     = aero_weight_array_rand_group(aero_state%awa, &
                     aero_particle%weight_set, &
                     aero_particle_radius(aero_particle))
                call aero_particle_set_group(new_aero_particle, new_group)
             end if
          end if
          call aero_state_add_particle(aero_state, new_aero_particle)
          ! re-get the particle pointer, which may have
          ! changed due to reallocations caused by adding
          aero_particle => aero_state%apa%particle(i_part)
       end do
       call aero_particle_deallocate(new_aero_particle)
    end if

  end subroutine aero_state_dup_particle

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> The number concentration of a single particle (m^{-3}).
  real(kind=dp) function aero_state_particle_num_conc(aero_state, &
       aero_particle)

    !> Aerosol state containing the particle.
    type(aero_state_t), intent(in) :: aero_state
    !> Aerosol particle.
    type(aero_particle_t), intent(in) :: aero_particle

    aero_state_particle_num_conc &
         = aero_weight_array_num_conc(aero_state%awa, aero_particle)

  end function aero_state_particle_num_conc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Save the correct number concentrations for later use by
  !> aero_state_reweight().
  subroutine aero_state_num_conc_for_reweight(aero_state, reweight_num_conc)

    !> Aerosol state.
    type(aero_state_t), intent(in) :: aero_state
    !> Number concentrations for later use by aero_state_reweight().
    real(kind=dp), intent(out) :: reweight_num_conc(aero_state%apa%n_part)

    integer :: i_part

    do i_part = 1,aero_state%apa%n_part
       reweight_num_conc(i_part) &
            = aero_weight_array_single_num_conc(aero_state%awa, &
            aero_state%apa%particle(i_part))
    end do

  end subroutine aero_state_num_conc_for_reweight

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Reweight all particles after their constituent volumes have been
  !> altered.
  !!
  !! The pattern for use should be like:
  !! <pre>
  !! call aero_state_num_conc_for_reweight(aero_state, reweight_num_conc)
  !! ... alter particle species volumes in aero_state ...
  !! call aero_state_reweight(aero_state, reweight_num_conc)
  !! </pre>
  subroutine aero_state_reweight(aero_state, reweight_num_conc)

    !> Aerosol state.
    type(aero_state_t), intent(inout) :: aero_state
    !> Number concentrations previously computed by
    !> aero_state_num_conc_for_reweight().
    real(kind=dp), intent(in) :: reweight_num_conc(aero_state%apa%n_part)

    integer :: i_part, i_group, i_set
    real(kind=dp) :: n_part_old(size(aero_state%awa%weight, 1), &
         size(aero_state%awa%weight, 2))
    real(kind=dp) :: n_part_new(size(aero_state%awa%weight, 1), &
         size(aero_state%awa%weight, 2))
    real(kind=dp) :: old_num_conc, new_num_conc, n_part_mean
    type(aero_particle_t), pointer :: aero_particle

    ! find average number of new particles in each weight group, if
    ! comp_vol is not changed
    n_part_old = 0d0
    n_part_new = 0d0
    do i_part = 1,aero_state%apa%n_part
       aero_particle => aero_state%apa%particle(i_part)
       old_num_conc = reweight_num_conc(i_part)
       new_num_conc = aero_weight_array_single_num_conc(aero_state%awa, &
            aero_particle)
       n_part_mean = old_num_conc / new_num_conc
       i_group = aero_particle%weight_group
       i_set = aero_particle%weight_set
       n_part_new(i_group, i_set) = n_part_new(i_group, i_set) + n_part_mean
       n_part_old(i_group, i_set) = n_part_old(i_group, i_set) + 1d0
    end do

    ! alter comp_vol to leave the number of computational particles
    ! per weight bin unchanged
    do i_group = 1,size(aero_state%awa%weight, 1)
       do i_set = 1,size(aero_state%awa%weight, 2)
          if (n_part_old(i_group, i_set) == 0d0) cycle
          call aero_weight_scale_comp_vol( &
               aero_state%awa%weight(i_group, i_set), &
               n_part_old(i_group, i_set) / n_part_new(i_group, i_set))
       end do
    end do

    ! work backwards so any additions and removals will only affect
    ! particles that we've already dealt with
    do i_part = aero_state%apa%n_part,1,-1
       aero_particle => aero_state%apa%particle(i_part)
       old_num_conc = reweight_num_conc(i_part)
       new_num_conc &
            = aero_weight_array_single_num_conc(aero_state%awa, aero_particle)
       n_part_mean = old_num_conc / new_num_conc
       call aero_state_dup_particle(aero_state, i_part, n_part_mean)
    end do

  end subroutine aero_state_reweight

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> <tt>aero_state += aero_state_delta</tt>, with adding the computational
  !> volumes, so the new concentration is the (volume-weighted)
  !> average of the two concentration.
  subroutine aero_state_add(aero_state, aero_state_delta)

    !> Aerosol state.
    type(aero_state_t), intent(inout) :: aero_state
    !> Increment.
    type(aero_state_t), intent(in) :: aero_state_delta

    call aero_state_add_particles(aero_state, aero_state_delta)
    call aero_weight_array_add_comp_vol(aero_state%awa, aero_state_delta%awa)

  end subroutine aero_state_add

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> <tt>aero_state += aero_state_delta</tt>, with the computational
  !> volume of \c aero_state left unchanged, so the new concentration is the
  !> sum of the two concentrations, computed with \c aero_state%comp_vol.
  subroutine aero_state_add_particles(aero_state, aero_state_delta)

    !> Aerosol state.
    type(aero_state_t), intent(inout) :: aero_state
    !> Increment.
    type(aero_state_t), intent(in) :: aero_state_delta

    integer :: i_part, i_bin

    do i_part = 1,aero_state_delta%apa%n_part
       call aero_state_add_particle(aero_state, &
            aero_state_delta%apa%particle(i_part))
    end do
    call aero_info_array_add(aero_state%aero_info_array, &
         aero_state_delta%aero_info_array)
    
  end subroutine aero_state_add_particles

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Generates a Poisson sample of an \c aero_dist, adding to \c
  !> aero_state. The sampled amount is <tt>sample_prop *
  !> aero_state%comp_vol</tt>.
  subroutine aero_state_add_aero_dist_sample(aero_state, aero_data, &
       aero_dist, sample_prop, create_time, n_part_add)

    !> Aero state to add to.
    type(aero_state_t), intent(inout) :: aero_state
    !> Aero data values.
    type(aero_data_t), intent(in) :: aero_data
    !> Distribution to sample.
    type(aero_dist_t), intent(in) :: aero_dist
    !> Volume fraction to sample (1).
    real(kind=dp), intent(in) :: sample_prop
    !> Creation time for new particles (s).
    real(kind=dp), intent(in) :: create_time
    !> Number of particles added.
    integer, intent(out), optional :: n_part_add

    real(kind=dp) :: n_samp_avg, radius, vol, mean_n_part, n_part_new
    real(kind=dp) :: comp_vol_ratio, n_part_ideal_local_group
    integer :: n_samp, i_mode, i_samp, i_group, i_set, n_group, n_set
    integer :: global_n_part
    type(aero_mode_t), pointer :: aero_mode
    type(aero_particle_t) :: aero_particle

    call aero_particle_allocate_size(aero_particle, aero_data%n_spec, &
         aero_data%n_source)

    n_group = size(aero_state%awa%weight, 1)
    n_set = size(aero_state%awa%weight, 2)
    if (present(n_part_add)) then
       n_part_add = 0
    end if
    do i_group = 1,n_group
       do i_mode = 1,aero_dist%n_mode
          aero_mode => aero_dist%mode(i_mode)
          i_set = aero_mode%source

          ! adjust comp_vol if necessary
          global_n_part = aero_state_total_particles_all_procs(aero_state, &
               i_group, i_set)
          mean_n_part = real(global_n_part, kind=dp) &
               / real(pmc_mpi_size(), kind=dp)
          if (aero_state%awa%weight(i_group, i_set)%comp_vol == 0d0) then
             ! FIXME: assert that n_part in this weight group is zero
             aero_state%awa%weight(i_group, i_set)%comp_vol = 1d0
          end if
          n_samp_avg = sample_prop * aero_dist_number(aero_dist, &
               aero_state%awa%weight(i_group, i_set))
          n_part_new = mean_n_part + n_samp_avg
          if (n_part_new == 0d0) cycle
          n_part_ideal_local_group = aero_state%n_part_ideal &
               / real(pmc_mpi_size(), kind=dp) / real(n_group * n_set, kind=dp)
          if ((n_part_new <  n_part_ideal_local_group / 2d0) &
               .or. (n_part_new > n_part_ideal_local_group * 2d0)) &
               then
             comp_vol_ratio = n_part_ideal_local_group / n_part_new
             call aero_state_scale_comp_vol(aero_state, i_group, i_set, &
                  comp_vol_ratio)
          end if

          ! sample particles
          n_samp_avg = sample_prop * aero_mode_number(aero_mode, &
               aero_state%awa%weight(i_group, i_set))
          n_samp = rand_poisson(n_samp_avg)
          if (present(n_part_add)) then
             n_part_add = n_part_add + n_samp
          end if
          do i_samp = 1,n_samp
             call aero_particle_zero(aero_particle)
             call aero_mode_sample_radius(aero_mode, &
                  aero_state%awa%weight(i_group, i_set), radius)
             vol = rad2vol(radius)
             call aero_particle_set_vols(aero_particle, &
                  aero_mode%vol_frac * vol)
             call aero_particle_new_id(aero_particle)
             call aero_particle_set_group(aero_particle, i_group)
             call aero_particle_set_create_time(aero_particle, create_time)
             call aero_particle_set_source(aero_particle, aero_mode%source)
             call aero_state_add_particle(aero_state, aero_particle)
          end do
       end do
    end do
    call aero_particle_deallocate(aero_particle)

  end subroutine aero_state_add_aero_dist_sample
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Choose a random particle from the aero_state.
  subroutine aero_state_rand_particle(aero_state, i_part)

    !> Original state.
    type(aero_state_t), intent(in) :: aero_state
    !> Chosen random particle number.
    integer, intent(out) :: i_part

    call assert(950725003, aero_state%apa%n_part > 0)
    i_part = pmc_rand_int(aero_state%apa%n_part)

  end subroutine aero_state_rand_particle

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Generates a random sample by removing particles from
  !> aero_state_from and adding them to aero_state_to, which must
  !> be already allocated (and should have its comp_vol set).
  !!
  !! None of the computational volumes are altered by this sampling,
  !! making this the equivalent of aero_state_add_particles().
  subroutine aero_state_sample_particles(aero_state_from, aero_state_to, &
       sample_prob, removal_action)

    !> Original state.
    type(aero_state_t), intent(inout) :: aero_state_from
    !> Destination state.
    type(aero_state_t), intent(inout) :: aero_state_to
    !> Probability of sampling each particle.
    real(kind=dp), intent(in) :: sample_prob
    !> Action for removal (see pmc_aero_info module for action
    !> parameters). Set to AERO_INFO_NONE to not log removal.
    integer, intent(in) :: removal_action
    
    integer :: n_transfer, i_transfer, i_part
    logical :: do_add, do_remove
    real(kind=dp) :: num_conc_from, num_conc_to
    type(aero_info_t) :: aero_info

    call assert(721006962, (sample_prob >= 0d0) .and. (sample_prob <= 1d0))
    call aero_state_reset(aero_state_to)
    call aero_state_copy_weight(aero_state_from, aero_state_to)
    n_transfer = rand_binomial(aero_state_total_particles(aero_state_from), &
         sample_prob)
    i_transfer = 0
    do while (i_transfer < n_transfer)
       if (aero_state_total_particles(aero_state_from) <= 0) exit
       call aero_state_rand_particle(aero_state_from, i_part)
       num_conc_from = aero_weight_array_num_conc(aero_state_from%awa, &
            aero_state_from%apa%particle(i_part))
       num_conc_to = aero_weight_array_num_conc(aero_state_to%awa, &
            aero_state_from%apa%particle(i_part))

       if (num_conc_to == num_conc_from) then ! add and remove
          do_add = .true.
          do_remove = .true.
       elseif (num_conc_to < num_conc_from) then ! always add, maybe remove
          do_add = .true.
          do_remove = .false.
          if (pmc_random() < num_conc_to / num_conc_from) then
             do_remove = .true.
          end if
       else ! always remove, maybe add
          do_add = .false.
          if (pmc_random() < num_conc_from / num_conc_to) then
             do_add = .true.
          end if
          do_remove = .true.
       end if
       if (do_add) then
          call aero_state_add_particle(aero_state_to, &
               aero_state_from%apa%particle(i_part))
       end if
       if (do_remove) then
          if (removal_action /= AERO_INFO_NONE) then
             call aero_info_allocate(aero_info)
             aero_info%id = aero_state_from%apa%particle(i_part)%id
             aero_info%action = removal_action
             call aero_state_remove_particle_with_info(aero_state_from, &
                  i_part, aero_info)
             call aero_info_deallocate(aero_info)
          else
             call aero_state_remove_particle_no_info(aero_state_from, &
                  i_part)
          end if
          i_transfer = i_transfer + 1
       end if
    end do
    
  end subroutine aero_state_sample_particles
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Generates a random sample by removing particles from
  !> aero_state_from and adding them to aero_state_to, transfering
  !> computational volume as well. This is the equivalent of
  !> aero_state_add().
  subroutine aero_state_sample(aero_state_from, aero_state_to, &
       sample_prob, removal_action)

    !> Original state.
    type(aero_state_t), intent(inout) :: aero_state_from
    !> Destination state (previous contents will be lost).
    type(aero_state_t), intent(inout) :: aero_state_to
    !> Probability of sampling each particle.
    real(kind=dp), intent(in) :: sample_prob
    !> Action for removal (see pmc_aero_info module for action
    !> parameters). Set to AERO_INFO_NONE to not log removal.
    integer, intent(in) :: removal_action

    integer :: n_transfer, i_transfer, i_part
    logical :: do_add, do_remove
    real(kind=dp) :: num_conc_from, num_conc_to
    type(aero_info_t) :: aero_info

    call assert(393205561, (sample_prob >= 0d0) .and. (sample_prob <= 1d0))
    call aero_state_reset(aero_state_to)
    call aero_state_copy_weight(aero_state_from, aero_state_to)
    call aero_weight_array_zero_comp_vol(aero_state_to%awa)
    n_transfer = rand_binomial(aero_state_total_particles(aero_state_from), &
         sample_prob)
    do i_transfer = 1,n_transfer
       if (aero_state_total_particles(aero_state_from) <= 0) exit
       call aero_state_rand_particle(aero_state_from, i_part)

       call aero_state_add_particle(aero_state_to, &
            aero_state_from%apa%particle(i_part))
       if (removal_action /= AERO_INFO_NONE) then
          call aero_info_allocate(aero_info)
          aero_info%id = aero_state_from%apa%particle(i_part)%id
          aero_info%action = removal_action
          call aero_state_remove_particle_with_info(aero_state_from, &
               i_part, aero_info)
          call aero_info_deallocate(aero_info)
       else
          call aero_state_remove_particle_no_info(aero_state_from, &
               i_part)
       end if
    end do
    call aero_weight_array_transfer_comp_vol(aero_state_from%awa, &
         aero_state_to%awa, sample_prob)

  end subroutine aero_state_sample
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Create the bin number and mass arrays from aero_state%v.
  subroutine aero_state_to_binned(bin_grid, aero_data, aero_state, &
       aero_binned)
    
    !> Bin grid.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Aerosol state.
    type(aero_state_t), intent(in) :: aero_state
    !> Binned distributions.
    type(aero_binned_t), intent(inout) :: aero_binned
    
    integer :: i_part, i_bin
    type(aero_particle_t), pointer :: aero_particle

    aero_binned%num_conc = 0d0
    aero_binned%vol_conc = 0d0
    do i_part = 1,aero_state%apa%n_part
       aero_particle => aero_state%apa%particle(i_part)
       i_bin = bin_grid_particle_in_bin(bin_grid, &
            aero_particle_radius(aero_particle))
       if ((i_bin < 1) .or. (i_bin > bin_grid%n_bin)) then
          call warn_msg(980232449, "particle ID " &
               // trim(integer_to_string(aero_particle%id)) &
               // " outside of bin_grid, discarding")
       else
          aero_binned%vol_conc(i_bin,:) = aero_binned%vol_conc(i_bin,:) &
               + aero_particle%vol &
               * aero_weight_array_num_conc(aero_state%awa, &
               aero_particle) / bin_grid%log_width
          aero_binned%num_conc(i_bin) = aero_binned%num_conc(i_bin) &
               + aero_weight_array_num_conc(aero_state%awa, &
               aero_particle) / bin_grid%log_width
       end if
    end do
    
  end subroutine aero_state_to_binned
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the diameters of all particles. The \c diameters array
  !> will be reallocated if necessary.
  subroutine aero_state_diameters(aero_state, diameters)

    !> Aerosol state.
    type(aero_state_t), intent(in) :: aero_state
    !> Diameters array (m).
    real(kind=dp), intent(inout), allocatable :: diameters(:)

    call ensure_real_array_size(diameters, aero_state%apa%n_part)
    diameters = aero_particle_diameter( &
         aero_state%apa%particle(1:aero_state%apa%n_part))

  end subroutine aero_state_diameters

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the masses of all particles. The \c masses array
  !> will be reallocated if necessary.
  !!
  !! If \c include is specified then only those species are included
  !! in computing the masses. If \c exclude is specified then all
  !! species except those species are included. If both \c include and
  !! \c exclude arguments are specified then only those species in \c
  !! include but not in \c exclude are included.
  subroutine aero_state_masses(aero_state, aero_data, masses, include, exclude)

    !> Aerosol state.
    type(aero_state_t), intent(in) :: aero_state
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Masses array (kg).
    real(kind=dp), intent(inout), allocatable :: masses(:)
    !> Species names to include in the mass.
    character(len=*), optional :: include(:)
    !> Species names to exclude in the mass.
    character(len=*), optional :: exclude(:)

    logical :: use_species(aero_data%n_spec)
    integer :: i_name, i_spec

    call ensure_real_array_size(masses, aero_state%apa%n_part)
    if ((.not. present(include)) .and. (.not. present(exclude))) then
       masses = aero_particle_mass( &
            aero_state%apa%particle(1:aero_state%apa%n_part), aero_data)
    else
       if (present(include)) then
          use_species = .false.
          do i_name = 1, size(include)
             i_spec = aero_data_spec_by_name(aero_data, include(i_name))
             call assert_msg(963163690, i_spec > 0, &
                  "unknown species: " // trim(include(i_name)))
             use_species(i_spec) = .true.
          end do
       else
          use_species = .true.
       end if
       if (present(exclude)) then
          do i_name = 1, size(exclude)
             i_spec = aero_data_spec_by_name(aero_data, exclude(i_name))
             call assert_msg(950847713, i_spec > 0, &
                  "unknown species: " // trim(exclude(i_name)))
             use_species(i_spec) = .false.
          end do
       end if
       masses = 0d0
       do i_spec = 1,aero_data%n_spec
          if (use_species(i_spec)) then
             masses = masses + aero_particle_species_mass( &
                  aero_state%apa%particle(1:aero_state%apa%n_part), &
                  i_spec, aero_data)
          end if
       end do
    end if

  end subroutine aero_state_masses

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the number concentrations of all particles. The \c
  !> num_concs array will be reallocated if necessary.
  subroutine aero_state_num_concs(aero_state, num_concs)

    !> Aerosol state.
    type(aero_state_t), intent(in) :: aero_state
    !> Number concentrations array (m^{-3}).
    real(kind=dp), intent(inout), allocatable :: num_concs(:)

    integer :: i_part

    call ensure_real_array_size(num_concs, aero_state%apa%n_part)
    do i_part = 1,aero_state%apa%n_part
       num_concs(i_part) = aero_state_particle_num_conc(aero_state, &
            aero_state%apa%particle(i_part))
    end do

  end subroutine aero_state_num_concs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the total number concentration.
  real(kind=dp) function aero_state_total_num_conc(aero_state)

    !> Aerosol state.
    type(aero_state_t), intent(in) :: aero_state

    integer :: i_part

    aero_state_total_num_conc = 0d0
    do i_part = 1,aero_state%apa%n_part
       aero_state_total_num_conc = aero_state_total_num_conc &
            + aero_state_particle_num_conc(aero_state, &
            aero_state%apa%particle(i_part))
    end do

  end function aero_state_total_num_conc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Does the same thing as aero_state_to_bin() but based on dry radius.
  subroutine aero_state_to_binned_dry(bin_grid, aero_data, aero_state, &
       aero_binned)
    
    !> Bin grid.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Aerosol state.
    type(aero_state_t), intent(in) :: aero_state
    !> Binned distributions.
    type(aero_binned_t), intent(inout) :: aero_binned
    
    integer :: i_part, i_bin
    type(aero_particle_t), pointer :: aero_particle
    
    aero_binned%num_conc = 0d0
    aero_binned%vol_conc = 0d0
    do i_part = 1,aero_state%apa%n_part
       aero_particle => aero_state%apa%particle(i_part)
       i_bin = bin_grid_particle_in_bin(bin_grid, &
            aero_particle_solute_radius(aero_particle, aero_data))
       if ((i_bin < 1) .or. (i_bin > bin_grid%n_bin)) then
          call warn_msg(503871022, "particle ID " &
               // trim(integer_to_string(aero_particle%id)) &
               // " outside of bin_grid, discarding")
       else
          aero_binned%vol_conc(i_bin,:) = aero_binned%vol_conc(i_bin,:) &
               + aero_particle%vol &
               * aero_weight_array_num_conc(aero_state%awa, &
               aero_particle) / bin_grid%log_width
          aero_binned%num_conc(i_bin) = aero_binned%num_conc(i_bin) &
               + aero_weight_array_num_conc(aero_state%awa, &
               aero_particle) / bin_grid%log_width
       end if
    end do
    
  end subroutine aero_state_to_binned_dry
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Doubles number of particles in the given weight group.
  subroutine aero_state_double(aero_state, i_group, i_set)
    
    !> Aerosol state.
    type(aero_state_t), intent(inout) :: aero_state
    !> Weight group to double.
    integer, intent(in) :: i_group
    !> Weight set to double.
    integer, intent(in) :: i_set

    integer :: i_part
    type(aero_particle_t) :: aero_particle

    call aero_particle_allocate(aero_particle)
    do i_part = 1,aero_state%apa%n_part
       if ((aero_state%apa%particle(i_part)%weight_group == i_group) &
            .and. (aero_state%apa%particle(i_part)%weight_set == i_set)) then
          call aero_particle_copy(aero_state%apa%particle(i_part), &
               aero_particle)
          call aero_particle_new_id(aero_particle)
          call aero_state_add_particle(aero_state, aero_particle)
       end if
    end do
    call aero_particle_deallocate(aero_particle)
    aero_state%valid_sort = .false.
    call aero_weight_scale_comp_vol(aero_state%awa%weight(i_group, i_set), 2d0)

  end subroutine aero_state_double
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Remove approximately half of the particles in the given weight group.
  subroutine aero_state_halve(aero_state, i_group, i_set)
    
    !> Aerosol state.
    type(aero_state_t), intent(inout) :: aero_state
    !> Weight group to halve.
    integer, intent(in) :: i_group
    !> Weight set to halve.
    integer, intent(in) :: i_set
    
    integer :: i_part
    type(aero_info_t) :: aero_info

    call aero_info_allocate(aero_info)
    do i_part = aero_state%apa%n_part,1,-1
       if ((aero_state%apa%particle(i_part)%weight_group == i_group) &
            .and. (aero_state%apa%particle(i_part)%weight_set == i_set)) then
          if (pmc_random() < 0.5d0) then
             aero_info%id = aero_state%apa%particle(i_part)%id
             aero_info%action = AERO_INFO_HALVED
             call aero_state_remove_particle_with_info(aero_state, i_part, &
                  aero_info)
          end if
       end if
    end do
    call aero_info_deallocate(aero_info)
    call aero_weight_scale_comp_vol(aero_state%awa%weight(i_group, i_set), &
         0.5d0)

  end subroutine aero_state_halve
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Double or halve the particle population in each weight group to
  !> maintain close to \c n_part_ideal particles per process,
  !> allocated equally amongst the weight groups.
  subroutine aero_state_rebalance(aero_state, allow_doubling, allow_halving, &
       initial_state_warning)

    !> Aerosol state.
    type(aero_state_t), intent(inout) :: aero_state
    !> Whether to allow doubling of the population.
    logical, intent(in) :: allow_doubling
    !> Whether to allow halving of the population.
    logical, intent(in) :: allow_halving
    !> Whether to warn due to initial state doubling/halving.
    logical, intent(in) :: initial_state_warning

    integer :: i_group, i_set, n_group, n_set, global_n_part

    n_group = size(aero_state%awa%weight, 1)
    n_set = size(aero_state%awa%weight, 2)

    ! if we have less than half the maximum number of particles then
    ! double until we fill up the array
    if (allow_doubling) then
       do i_group = 1,n_group
          do i_set = 1,n_set
             global_n_part &
                  = aero_state_total_particles_all_procs(aero_state, i_group, &
                  i_set)
             do while ((real(global_n_part, kind=dp) &
                  < aero_state%n_part_ideal &
                  / real(n_group * n_set, kind=dp) / 2d0) &
                  .and. (global_n_part > 0))
                if (initial_state_warning) then
                   call warn_msg(716882783, "doubling particles in initial " &
                        // "condition")
                end if
                call aero_state_double(aero_state, i_group, i_set)
                global_n_part &
                     = aero_state_total_particles_all_procs(aero_state, &
                     i_group, i_set)
             end do
          end do
       end do
    end if

    ! same for halving if we have too many particles
    if (allow_halving) then
       do i_group = 1,n_group
          do i_set = 1,n_set
             do while (real(aero_state_total_particles_all_procs(aero_state, &
                  i_group, i_set), kind=dp) > aero_state%n_part_ideal &
                  / real(n_group * n_set, kind=dp) * 2d0)
                if (initial_state_warning) then
                   call warn_msg(661936373, &
                        "halving particles in initial condition")
                end if
                call aero_state_halve(aero_state, i_group, i_set)
             end do
          end do
       end do
    end if

  end subroutine aero_state_rebalance

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Scale the computational volume of the given group by the given
  !> ratio, altering particle number as necessary to preserve the
  !> number concentration.
  subroutine aero_state_scale_comp_vol(aero_state, i_group, i_set, &
       comp_vol_ratio)

    !> Aerosol state.
    type(aero_state_t), intent(inout) :: aero_state
    !> Weight group number.
    integer, intent(in) :: i_group
    !> Weight set number.
    integer, intent(in) :: i_set
    !> Ratio of <tt>new_comp_vol / old_comp_vol</tt>.
    real(kind=dp), intent(in) :: comp_vol_ratio

    real(kind=dp) :: ratio
    integer :: i_part, i_remove, n_remove, i_entry
    type(aero_info_t) :: aero_info

    ! We could use the ratio > 1 case unconditionally, but that would
    ! have higher variance for the ratio < 1 case than the current
    ! scheme.

    call aero_weight_scale_comp_vol(aero_state%awa%weight(i_group, i_set), &
         comp_vol_ratio)

    if (aero_state%apa%n_part == 0) return

    call aero_state_sort(aero_state)

    if (comp_vol_ratio < 1d0) then
       n_remove = prob_round(comp_vol_ratio &
            * real(aero_state%aero_sorted%weight%inverse(i_group)%n_entry, &
            kind=dp))
       do i_remove = 1,n_remove
          i_entry = pmc_rand_int(aero_state%aero_sorted%weight%inverse( &
               i_group)%n_entry)
          i_part = aero_state%aero_sorted%weight%inverse(i_group)%entry( &
               i_entry)
          call aero_info_allocate(aero_info)
          aero_info%id = aero_state%apa%particle(i_part)%id
          aero_info%action = AERO_INFO_HALVED
          call aero_state_remove_particle(aero_state, i_part, .true., &
               aero_info)
          call aero_info_deallocate(aero_info)
       end do
    elseif (comp_vol_ratio > 1d0) then
       do i_entry = aero_state%aero_sorted%weight%inverse(i_group)%n_entry,1,-1
          i_part = aero_state%aero_sorted%weight%inverse(i_group)%entry( &
               i_entry)
          call aero_state_dup_particle(aero_state, i_part, comp_vol_ratio)
       end do
    end if

  end subroutine aero_state_scale_comp_vol

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Mix the aero_states between all processes. Currently uses a
  !> simple all-to-all diffusion.
  subroutine aero_state_mix(aero_state, del_t, mix_timescale, &
       aero_data, specify_prob_transfer)

    !> Aerosol state.
    type(aero_state_t), intent(inout) :: aero_state
    !> Timestep (s).
    real(kind=dp), intent(in) :: del_t
    !> Mixing timescale (s).
    real(kind=dp), intent(in) :: mix_timescale
    !> Aero data values.
    type(aero_data_t), intent(in) :: aero_data
    !> Transfer probability of each particle (0 means no mixing, 1
    !> means total mixing).
    real(kind=dp), optional, intent(in) :: specify_prob_transfer

#ifdef PMC_USE_MPI
    integer :: rank, n_proc, i_proc, ierr
    integer :: buffer_size, buffer_size_check
    character, allocatable :: buffer(:)
    type(aero_state_t), allocatable :: aero_state_sends(:)
    type(aero_state_t), allocatable :: aero_state_recvs(:)
    real(kind=dp) :: prob_transfer, prob_not_transferred
    real(kind=dp) :: prob_transfer_given_not_transferred

    ! process information
    rank = pmc_mpi_rank()
    n_proc = pmc_mpi_size()
    if (n_proc == 1) then
       ! buffer allocation below fails if n_proc == 1
       ! so bail out early (nothing to mix anyway)
       return
    end if

    ! allocate aero_state arrays
    allocate(aero_state_sends(n_proc))
    allocate(aero_state_recvs(n_proc))
    do i_proc = 0,(n_proc - 1)
       call aero_state_allocate(aero_state_sends(i_proc + 1))
       call aero_state_allocate(aero_state_recvs(i_proc + 1))
    end do

    ! compute the transfer probability
    if (present(specify_prob_transfer)) then
       prob_transfer = specify_prob_transfer / real(n_proc, kind=dp)
    else
       prob_transfer = (1d0 - exp(- del_t / mix_timescale)) &
            / real(n_proc, kind=dp)
    end if

    ! extract particles to send
    prob_not_transferred = 1d0
    do i_proc = 0,(n_proc - 1)
       if (i_proc /= rank) then
          ! because we are doing sequential sampling from the aero_state
          ! we need to scale up the later transfer probabilities, because
          ! the later particles are being transferred conditioned on the
          ! fact that they did not transfer already
          prob_transfer_given_not_transferred = prob_transfer &
               / prob_not_transferred
          call aero_state_sample(aero_state, &
               aero_state_sends(i_proc + 1), &
               prob_transfer_given_not_transferred, AERO_INFO_NONE)
          prob_not_transferred = prob_not_transferred - prob_transfer
       end if
    end do

    ! exchange the particles
    call aero_state_mpi_alltoall(aero_state_sends, aero_state_recvs)

    ! process the received particles
    do i_proc = 0,(n_proc - 1)
       if (i_proc /= rank) then
          call aero_state_add(aero_state, aero_state_recvs(i_proc + 1))
       end if
    end do

    ! cleanup
    do i_proc = 0,(n_proc - 1)
          call aero_state_deallocate(aero_state_sends(i_proc + 1))
          call aero_state_deallocate(aero_state_recvs(i_proc + 1))
    end do
    deallocate(aero_state_sends)
    deallocate(aero_state_recvs)
#endif

  end subroutine aero_state_mix

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Do an MPI all-to-all transfer of aerosol states.
  !!
  !! States without particles are not sent.
  subroutine aero_state_mpi_alltoall(send, recv)

    !> Array of aero_states to send (one per process).
    type(aero_state_t), intent(in) :: send(:)
    !> Array of aero_states to receives (one per process).
    type(aero_state_t), intent(inout) :: recv(size(send))

#ifdef PMC_USE_MPI
    character, allocatable :: sendbuf(:), recvbuf(:)
    integer :: sendcounts(size(send)), sdispls(size(send))
    integer :: recvcounts(size(send)), rdispls(size(send))
    integer :: i_proc, position, old_position, max_sendbuf_size, ierr

    call assert(978709842, size(send) == pmc_mpi_size())

    max_sendbuf_size = 0
    do i_proc = 1,pmc_mpi_size()
       if (aero_state_total_particles(send(i_proc)) > 0) then
          max_sendbuf_size = max_sendbuf_size &
               + pmc_mpi_pack_size_aero_state(send(i_proc))
       end if
    end do

    allocate(sendbuf(max_sendbuf_size))

    position = 0
    do i_proc = 1,pmc_mpi_size()
       old_position = position
       if (aero_state_total_particles(send(i_proc)) > 0) then
          call pmc_mpi_pack_aero_state(sendbuf, position, send(i_proc))
       end if
       sendcounts(i_proc) = position - old_position
    end do
    call assert(393267406, position <= max_sendbuf_size)

    call pmc_mpi_alltoall_integer(sendcounts, recvcounts)
    allocate(recvbuf(sum(recvcounts)))

    sdispls(1) = 0
    rdispls(1) = 0
    do i_proc = 2,pmc_mpi_size()
       sdispls(i_proc) = sdispls(i_proc - 1) + sendcounts(i_proc - 1)
       rdispls(i_proc) = rdispls(i_proc - 1) + recvcounts(i_proc - 1)
    end do

    call mpi_alltoallv(sendbuf, sendcounts, sdispls, MPI_CHARACTER, recvbuf, &
         recvcounts, rdispls, MPI_CHARACTER, MPI_COMM_WORLD, ierr)
    call pmc_mpi_check_ierr(ierr)

    position = 0
    do i_proc = 1,pmc_mpi_size()
       call assert(189739257, position == rdispls(i_proc))
       if (recvcounts(i_proc) > 0) then
          call pmc_mpi_unpack_aero_state(recvbuf, position, recv(i_proc))
       end if
    end do

    deallocate(sendbuf)
    deallocate(recvbuf)
#endif

  end subroutine aero_state_mpi_alltoall

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Set each aerosol particle to have its original total volume, but
  !> species volume ratios given by the total species volume ratio
  !> within each bin. This preserves the (weighted) total species
  !> volume per bin as well as per-particle total volumes.
  subroutine aero_state_bin_average_comp(aero_state, bin_grid, aero_data, &
       dry_volume)

    !> Aerosol state to average.
    type(aero_state_t), intent(inout) :: aero_state
    !> Bin grid to average within.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Whether to use dry volume (rather than wet).
    logical, intent(in) :: dry_volume

    real(kind=dp) :: species_volume_conc(aero_data%n_spec)
    real(kind=dp) :: total_volume_conc, particle_volume, num_conc
    integer :: i_bin, i_entry, i_part, i_spec
    type(aero_particle_t), pointer :: aero_particle

    call aero_state_sort(aero_state, bin_grid)

    do i_bin = 1,bin_grid%n_bin
       species_volume_conc = 0d0
       total_volume_conc = 0d0
       do i_entry = 1,aero_state%aero_sorted%size%inverse(i_bin)%n_entry
          i_part = aero_state%aero_sorted%size%inverse(i_bin)%entry(i_entry)
          aero_particle => aero_state%apa%particle(i_part)
          num_conc = aero_weight_array_num_conc(aero_state%awa, aero_particle)
          particle_volume = aero_particle_volume_maybe_dry(aero_particle, &
               aero_data, dry_volume)
          species_volume_conc = species_volume_conc &
               + num_conc * aero_particle%vol
          total_volume_conc = total_volume_conc + num_conc * particle_volume
       end do
       do i_entry = 1,aero_state%aero_sorted%size%inverse(i_bin)%n_entry
          i_part = aero_state%aero_sorted%size%inverse(i_bin)%entry(i_entry)
          aero_particle => aero_state%apa%particle(i_part)
          particle_volume = aero_particle_volume_maybe_dry(aero_particle, &
               aero_data, dry_volume)
          aero_particle%vol = particle_volume * species_volume_conc &
               / total_volume_conc
          if (dry_volume .and. (aero_data%i_water > 0)) then
             ! set water to zero if we are doing dry volume averaging
             aero_particle%vol(aero_data%i_water) = 0d0
          end if
       end do
    end do

  end subroutine aero_state_bin_average_comp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Set each aerosol particle to have its original species ratios,
  !> but total volume given by the average volume of all particles
  !> within each bin.
  !!
  !! This does not preserve the total species volume
  !! per bin. If the \c bin_center parameter is \c .true. then the
  !! particles in each bin are set to have the bin center volume,
  !! rather than the average volume of the particles in that bin.
  !!
  !! If the weighting function is not constant (AERO_WEIGHT_TYPE_NONE)
  !! then the averaging can be performed in either a number-preserving
  !! way or in a volume-preserving way. The volume-preserving way does
  !! not preserve species volume ratios in gernal, but will do so if
  !! the particle population has already been composition-averaged.
  subroutine aero_state_bin_average_size(aero_state, bin_grid, aero_data, &
       dry_volume, bin_center, preserve_number)

    !> Aerosol state to average.
    type(aero_state_t), intent(inout) :: aero_state
    !> Bin grid to average within.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Whether to use dry volume (rather than wet).
    logical, intent(in) :: dry_volume
    !> Whether to assign the bin center volume (rather than the average
    !> volume).
    logical, intent(in) :: bin_center
    !> Whether to use the number-preserving scheme (otherwise will use
    !> the volume-preserving scheme). This parameter has no effect if
    !> \c bin_center is \c .true.
    logical, intent(in) :: preserve_number

    real(kind=dp) :: total_volume_conc, particle_volume
    real(kind=dp) :: new_particle_volume, num_conc, total_num_conc
    real(kind=dp) :: lower_volume, upper_volume, center_volume
    real(kind=dp) :: lower_function, upper_function, center_function
    integer :: i_bin, i_entry, i_part, i_bisect, n_part
    logical :: monotone_increasing, monotone_decreasing
    type(aero_particle_t), pointer :: aero_particle

    call aero_state_sort(aero_state, bin_grid)

    call die_msg(324563039, "FIXME: i_set is not used correctly here")
    do i_bin = 1,bin_grid%n_bin
       if (aero_state%aero_sorted%size%inverse(i_bin)%n_entry == 0) then
          cycle
       end if

       n_part = aero_state%aero_sorted%size%inverse(i_bin)%n_entry
       total_num_conc = 0d0
       total_volume_conc = 0d0
       do i_entry = 1,aero_state%aero_sorted%size%inverse(i_bin)%n_entry
          i_part = aero_state%aero_sorted%size%inverse(i_bin)%entry(i_entry)
          aero_particle => aero_state%apa%particle(i_part)
          num_conc = aero_weight_array_num_conc(aero_state%awa, aero_particle)
          total_num_conc = total_num_conc + num_conc
          particle_volume = aero_particle_volume_maybe_dry(aero_particle, &
               aero_data, dry_volume)
          total_volume_conc = total_volume_conc &
               + num_conc * particle_volume
       end do

       ! determine the new_particle_volume for all particles in this bin
       if (bin_center) then
          new_particle_volume = rad2vol(bin_grid%center_radius(i_bin))
       elseif (aero_weight_array_check_flat(aero_state%awa)) then
          num_conc & ! any radius will have the same num_conc
               = aero_weight_array_num_conc_at_radius(aero_state%awa, 1, 1d0)
          new_particle_volume = total_volume_conc / num_conc &
               / real(aero_state%aero_sorted%size%inverse(i_bin)%n_entry, &
               kind=dp)
       elseif (preserve_number) then
          ! number-preserving scheme: Solve the implicit equation:
          ! n_part * W(new_vol) = total_num_conc
          !
          ! We assume that the weighting function is strictly monotone
          ! so this equation has a unique solution and the solution
          ! lies between the min and max particle volumes. We use
          ! bisection as this doesn't really need to be fast, just
          ! robust.

          call aero_weight_array_check_monotonicity(aero_state%awa, &
               monotone_increasing, monotone_decreasing)
          call assert_msg(214077200, &
               monotone_increasing .or. monotone_decreasing, &
               "monotone weight function required for averaging")

          ! initialize to min and max particle volumes
          do i_entry = 1,n_part
             i_part = aero_state%aero_sorted%size%inverse(i_bin)%entry(i_entry)
             aero_particle => aero_state%apa%particle(i_part)
             particle_volume = aero_particle_volume_maybe_dry(aero_particle, &
                  aero_data, dry_volume)
             if (i_part == 1) then
                lower_volume = particle_volume
                upper_volume = particle_volume
             end if
             lower_volume = min(lower_volume, particle_volume)
             upper_volume = max(upper_volume, particle_volume)
          end do

          lower_function = real(n_part, kind=dp) &
               * aero_weight_array_num_conc_at_radius( &
               aero_state%awa, 1, vol2rad(lower_volume)) - total_num_conc
          upper_function = real(n_part, kind=dp) &
               * aero_weight_array_num_conc_at_radius(&
               aero_state%awa, 1, vol2rad(upper_volume)) - total_num_conc

          ! do 50 rounds of bisection (2^50 = 10^15)
          do i_bisect = 1,50
             center_volume = (lower_volume + upper_volume) / 2d0
             center_function = real(n_part, kind=dp) &
                  * aero_weight_array_num_conc_at_radius(aero_state%awa, 1, &
                  vol2rad(center_volume)) - total_num_conc
             if ((lower_function > 0d0 .and. center_function > 0d0) &
                  .or. (lower_function < 0d0 .and. center_function < 0d0)) &
                  then
                lower_volume = center_volume
                lower_function = center_function
             else
                upper_volume = center_volume
                upper_function = center_function
             end if
          end do

          new_particle_volume = center_volume
       else
          ! volume-preserving scheme: Solve the implicit equation:
          ! n_part * W(new_vol) * new_vol = total_volume_conc
          !
          ! We assume that the weighting function is strictly monotone
          ! so this equation has a unique solution and the solution
          ! lies between the min and max particle volumes. We use
          ! bisection as this doesn't really need to be fast, just
          ! robust.

          call aero_weight_array_check_monotonicity(aero_state%awa, &
               monotone_increasing, monotone_decreasing)
          call assert_msg(483078128, &
               monotone_increasing .or. monotone_decreasing, &
               "monotone weight function required for averaging")

          ! initialize to min and max particle volumes
          do i_entry = 1,n_part
             i_part = aero_state%aero_sorted%size%inverse(i_bin)%entry(i_entry)
             aero_particle => aero_state%apa%particle(i_part)
             particle_volume = aero_particle_volume_maybe_dry(aero_particle, &
                  aero_data, dry_volume)
             if (i_part == 1) then
                lower_volume = particle_volume
                upper_volume = particle_volume
             end if
             lower_volume = min(lower_volume, particle_volume)
             upper_volume = max(upper_volume, particle_volume)
          end do

          lower_function = real(n_part, kind=dp) &
               * aero_weight_array_num_conc_at_radius( &
               aero_state%awa, 1, vol2rad(lower_volume)) &
               * lower_volume - total_volume_conc
          upper_function = real(n_part, kind=dp) &
               * aero_weight_array_num_conc_at_radius( &
               aero_state%awa, 1, vol2rad(upper_volume)) &
               * upper_volume - total_volume_conc

          ! do 50 rounds of bisection (2^50 = 10^15)
          do i_bisect = 1,50
             center_volume = (lower_volume + upper_volume) / 2d0
             center_function = real(n_part, kind=dp) &
                  * aero_weight_array_num_conc_at_radius( &
                  aero_state%awa, 1, vol2rad(center_volume)) &
                  * center_volume - total_volume_conc
             if ((lower_function > 0d0 .and. center_function > 0d0) &
                  .or. (lower_function < 0d0 .and. center_function < 0d0)) &
                  then
                lower_volume = center_volume
                lower_function = center_function
             else
                upper_volume = center_volume
                upper_function = center_function
             end if
          end do

          new_particle_volume = center_volume
       end if

       do i_entry = 1,n_part
          i_part = aero_state%aero_sorted%size%inverse(i_bin)%entry(i_entry)
          aero_particle => aero_state%apa%particle(i_part)
          particle_volume = aero_particle_volume_maybe_dry(aero_particle, &
               aero_data, dry_volume)
          aero_particle%vol = aero_particle%vol / particle_volume &
               * new_particle_volume
          if (dry_volume .and. (aero_data%i_water > 0)) then
             ! set water to zero if we are doing dry volume averaging
             aero_particle%vol(aero_data%i_water) = 0d0
          end if
       end do
    end do

  end subroutine aero_state_bin_average_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Make all particles dry (water set to zero).
  subroutine aero_state_make_dry(aero_state, aero_data)

    !> Aerosol state to make dry.
    type(aero_state_t), intent(inout) :: aero_state
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data

    integer :: i_part

    if (aero_data%i_water > 0) then
       do i_part = 1,aero_state%apa%n_part
          aero_state%apa%particle(i_part)%vol(aero_data%i_water) = 0d0
       end do
       aero_state%valid_sort = .false.
    end if

   end subroutine aero_state_make_dry

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Determines the number of bytes required to pack the given value.
  integer function pmc_mpi_pack_size_aero_state(val)

    !> Value to pack.
    type(aero_state_t), intent(in) :: val

    integer :: total_size, i_group

    total_size = 0
    total_size = total_size + pmc_mpi_pack_size_apa(val%apa)
    total_size = total_size + pmc_mpi_pack_size_aero_weight_array(val%awa)
    total_size = total_size + pmc_mpi_pack_size_real(val%n_part_ideal)
    total_size = total_size + pmc_mpi_pack_size_aia(val%aero_info_array)
    pmc_mpi_pack_size_aero_state = total_size

  end function pmc_mpi_pack_size_aero_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Packs the given value into the buffer, advancing position.
  subroutine pmc_mpi_pack_aero_state(buffer, position, val)

    !> Memory buffer.
    character, intent(inout) :: buffer(:)
    !> Current buffer position.
    integer, intent(inout) :: position
    !> Value to pack.
    type(aero_state_t), intent(in) :: val

#ifdef PMC_USE_MPI
    integer :: prev_position, i_group

    prev_position = position
    call pmc_mpi_pack_aero_particle_array(buffer, position, val%apa)
    call pmc_mpi_pack_aero_weight_array(buffer,position,val%awa)
    call pmc_mpi_pack_real(buffer, position, val%n_part_ideal)
    call pmc_mpi_pack_aero_info_array(buffer, position, val%aero_info_array)
    call assert(850997402, &
         position - prev_position <= pmc_mpi_pack_size_aero_state(val))
#endif

  end subroutine pmc_mpi_pack_aero_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Unpacks the given value from the buffer, advancing position.
  subroutine pmc_mpi_unpack_aero_state(buffer, position, val)

    !> Memory buffer.
    character, intent(inout) :: buffer(:)
    !> Current buffer position.
    integer, intent(inout) :: position
    !> Value to pack.
    type(aero_state_t), intent(inout) :: val

#ifdef PMC_USE_MPI
    integer :: prev_position, i_group, n_group

    val%valid_sort = .false.
    prev_position = position
    call pmc_mpi_unpack_aero_particle_array(buffer, position, val%apa)
    call pmc_mpi_unpack_aero_weight_array(buffer,position,val%awa)
    call pmc_mpi_unpack_real(buffer, position, val%n_part_ideal)
    call pmc_mpi_unpack_aero_info_array(buffer, position, val%aero_info_array)
    call assert(132104747, &
         position - prev_position <= pmc_mpi_pack_size_aero_state(val))
#endif

  end subroutine pmc_mpi_unpack_aero_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Gathers data from all processes into one aero_state on process 0.
  subroutine aero_state_mpi_gather(aero_state, aero_state_total)

    !> Local aero_state.
    type(aero_state_t), intent(in) :: aero_state
    !> Centralized aero_state (only on process 0).
    type(aero_state_t), intent(inout) :: aero_state_total

#ifdef PMC_USE_MPI
    type(aero_state_t) :: aero_state_transfer
    integer :: n_proc, ierr, status(MPI_STATUS_SIZE)
    integer :: buffer_size, max_buffer_size, i_proc, position
    character, allocatable :: buffer(:)
#endif

    if (pmc_mpi_rank() == 0) then
       call aero_state_copy(aero_state, aero_state_total)
    end if

#ifdef PMC_USE_MPI
    
    if (pmc_mpi_rank() /= 0) then
       ! send data from remote processes
       max_buffer_size = 0
       max_buffer_size = max_buffer_size &
            + pmc_mpi_pack_size_aero_state(aero_state)
       allocate(buffer(max_buffer_size))
       position = 0
       call pmc_mpi_pack_aero_state(buffer, position, aero_state)
       call assert(542772170, position <= max_buffer_size)
       buffer_size = position
       call mpi_send(buffer, buffer_size, MPI_CHARACTER, 0, &
            AERO_STATE_TAG_GATHER, MPI_COMM_WORLD, ierr)
       call pmc_mpi_check_ierr(ierr)
       deallocate(buffer)
    else
       ! root process receives data from each remote proc
       n_proc = pmc_mpi_size()
       do i_proc = 1,(n_proc - 1)
          ! determine buffer size at root process
          call mpi_probe(i_proc, AERO_STATE_TAG_GATHER, MPI_COMM_WORLD, &
               status, ierr)
          call pmc_mpi_check_ierr(ierr)
          call mpi_get_count(status, MPI_CHARACTER, buffer_size, ierr)
          call pmc_mpi_check_ierr(ierr)

          ! get buffer at root process
          allocate(buffer(buffer_size))
          call mpi_recv(buffer, buffer_size, MPI_CHARACTER, i_proc, &
               AERO_STATE_TAG_GATHER, MPI_COMM_WORLD, status, ierr)

          ! unpack it
          position = 0
          call aero_state_allocate(aero_state_transfer)
          call pmc_mpi_unpack_aero_state(buffer, position, &
               aero_state_transfer)
          call assert(518174881, position == buffer_size)
          deallocate(buffer)

          call aero_state_add(aero_state_total, aero_state_transfer)
          
          call aero_state_deallocate(aero_state_transfer)
       end do
    end if

#endif

  end subroutine aero_state_mpi_gather

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Write the aero particle dimension to the given NetCDF file if it
  !> is not already present and in any case return the associated
  !> dimid.
  subroutine aero_state_netcdf_dim_aero_particle(aero_state, ncid, &
       dimid_aero_particle)

    !> aero_state structure.
    type(aero_state_t), intent(in) :: aero_state
    !> NetCDF file ID, in data mode.
    integer, intent(in) :: ncid
    !> Dimid of the aero particle dimension.
    integer, intent(out) :: dimid_aero_particle

    integer :: status, i_part
    integer :: varid_aero_particle
    integer :: aero_particle_centers(aero_state%apa%n_part)

    ! try to get the dimension ID
    status = nf90_inq_dimid(ncid, "aero_particle", dimid_aero_particle)
    if (status == NF90_NOERR) return
    if (status /= NF90_EBADDIM) call pmc_nc_check(status)

    ! dimension not defined, so define now define it
    call pmc_nc_check(nf90_redef(ncid))

    call pmc_nc_check(nf90_def_dim(ncid, "aero_particle", &
         aero_state%apa%n_part, dimid_aero_particle))

    call pmc_nc_check(nf90_enddef(ncid))

    do i_part = 1,aero_state%apa%n_part
       aero_particle_centers(i_part) = i_part
    end do
    call pmc_nc_write_integer_1d(ncid, aero_particle_centers, &
         "aero_particle", (/ dimid_aero_particle /), &
         description="dummy dimension variable (no useful value)")

  end subroutine aero_state_netcdf_dim_aero_particle

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Write the aero removed dimension to the given NetCDF file if it
  !> is not already present and in any case return the associated
  !> dimid.
  subroutine aero_state_netcdf_dim_aero_removed(aero_state, ncid, &
       dimid_aero_removed)

    !> aero_state structure.
    type(aero_state_t), intent(in) :: aero_state
    !> NetCDF file ID, in data mode.
    integer, intent(in) :: ncid
    !> Dimid of the aero removed dimension.
    integer, intent(out) :: dimid_aero_removed

    integer :: status, i_remove, dim_size
    integer :: varid_aero_removed
    integer :: aero_removed_centers(max(aero_state%aero_info_array%n_item,1))

    ! try to get the dimension ID
    status = nf90_inq_dimid(ncid, "aero_removed", dimid_aero_removed)
    if (status == NF90_NOERR) return
    if (status /= NF90_EBADDIM) call pmc_nc_check(status)

    ! dimension not defined, so define now define it
    call pmc_nc_check(nf90_redef(ncid))

    dim_size = max(aero_state%aero_info_array%n_item, 1)
    call pmc_nc_check(nf90_def_dim(ncid, "aero_removed", &
         dim_size, dimid_aero_removed))

    call pmc_nc_check(nf90_enddef(ncid))

    do i_remove = 1,dim_size
       aero_removed_centers(i_remove) = i_remove
    end do
    call pmc_nc_write_integer_1d(ncid, aero_removed_centers, &
         "aero_removed", (/ dimid_aero_removed /), &
         description="dummy dimension variable (no useful value)")

  end subroutine aero_state_netcdf_dim_aero_removed

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Write full state.
  subroutine aero_state_output_netcdf(aero_state, ncid, aero_data, &
       record_removals, record_optical)
    
    !> aero_state to write.
    type(aero_state_t), intent(in) :: aero_state
    !> NetCDF file ID, in data mode.
    integer, intent(in) :: ncid
    !> aero_data structure.
    type(aero_data_t), intent(in) :: aero_data
    !> Whether to output particle removal info.
    logical, intent(in) :: record_removals
    !> Whether to output aerosol optical properties.
    logical, intent(in) :: record_optical

    integer :: dimid_aero_particle, dimid_aero_species, dimid_aero_source
    integer :: dimid_aero_removed
    integer :: i_part, i_remove
    type(aero_particle_t), pointer :: particle
    real(kind=dp) :: aero_particle_mass(aero_state%apa%n_part, &
         aero_data%n_spec)
    integer :: aero_n_orig_part(aero_state%apa%n_part, aero_data%n_source)
    integer :: aero_weight_group(aero_state%apa%n_part)
    real(kind=dp) :: aero_absorb_cross_sect(aero_state%apa%n_part)
    real(kind=dp) :: aero_scatter_cross_sect(aero_state%apa%n_part)
    real(kind=dp) :: aero_asymmetry(aero_state%apa%n_part)
    real(kind=dp) :: aero_refract_shell_real(aero_state%apa%n_part)
    real(kind=dp) :: aero_refract_shell_imag(aero_state%apa%n_part)
    real(kind=dp) :: aero_refract_core_real(aero_state%apa%n_part)
    real(kind=dp) :: aero_refract_core_imag(aero_state%apa%n_part)
    real(kind=dp) :: aero_core_vol(aero_state%apa%n_part)
    integer :: aero_water_hyst_leg(aero_state%apa%n_part)
    real(kind=dp) :: aero_comp_vol(aero_state%apa%n_part)
    integer :: aero_id(aero_state%apa%n_part)
    real(kind=dp) :: aero_least_create_time(aero_state%apa%n_part)
    real(kind=dp) :: aero_greatest_create_time(aero_state%apa%n_part)
    integer :: aero_removed_id(max(aero_state%aero_info_array%n_item,1))
    integer :: aero_removed_action(max(aero_state%aero_info_array%n_item,1))
    integer :: aero_removed_other_id(max(aero_state%aero_info_array%n_item,1))

    !> \page output_format_aero_state Output File Format: Aerosol Particle State
    !!
    !! The aerosol state consists of a set of individual aerosol
    !! particles, each with its own individual properties. The
    !! properties of all particles are stored in arrays, one per
    !! property. For example, <tt>aero_absorb_cross_sect(i)</tt> gives
    !! the absorption cross section of particle number \c i, while
    !! <tt>aero_particle_mass(i,s)</tt> gives the mass of species \c s
    !! in particle \c i. The aerosol species are described in \ref
    !! output_format_aero_data.
    !!
    !! Each aerosol particle \c i represents a volume of space known
    !! as the computational volume and given by
    !! <tt>aero_comp_vol(i)</tt>. Dividing a per-particle quantity by
    !! the respective computational volume gives the concentration of
    !! that quantity contributed by the particle. For example, summing
    !! <tt>aero_particle_mass(i,s) / aero_comp_vol(i)</tt> over all \c
    !! i gives the total concentration of species \c s in
    !! kg/m^3. Similarly, summing <tt>aero_absorb_cross_sect(i) /
    !! aero_comp_vol(i)</tt> over all \c i will give the concentration
    !! of scattering cross section in m^2/m^3.
    !!
    !! FIXME: the aero_weight is also output
    !!
    !! The aerosol state uses the \c aero_species NetCDF dimension as
    !! specified in the \ref output_format_aero_data section, as well
    !! as the NetCDF dimension:
    !!   - \b aero_particle: number of aerosol particles
    !!
    !! The aerosol state NetCDF variables are:
    !!   - \b aero_particle (dim \c aero_particle): dummy dimension variable
    !!     (no useful value)
    !!   - \b aero_particle_mass (unit kg,
    !!     dim <tt>aero_particle x aero_species</tt>): constituent masses of
    !!     each aerosol particle - <tt>aero_particle_mass(i,s)</tt> gives the
    !!     mass of species \c s in particle \c i
    !!   - \b aero_n_orig_part (dim <tt>aero_particle x
    !!     aero_source</tt>): number of original particles from each
    !!     source that formed each aerosol particle -
    !!     <tt>aero_n_orig_part(i,s)</tt> is the number of particles
    !!     from source \c s that contributed to particle \c i - when
    !!     particle \c i first enters the simulation (by emissions,
    !!     dilution, etc.) it has <tt>aero_n_orig_part(i,s) = 1</tt>
    !!     for the source number \c s it came from (otherwise zero)
    !!     and when two particles coagulate, their values of \c
    !!     aero_n_orig_part are added (the number of coagulation
    !!     events that formed each particle is thus
    !!     <tt>sum(aero_n_orig_part(i,:)) - 1</tt>)
    !!   - \b aero_weight_group (dim <tt>aero_particle</tt>): weight
    !!     group number (1 to <tt>aero_weight</tt>) of each aerosol
    !!     particle
    !!   - \b aero_absorb_cross_sect (unit m^2, dim \c aero_particle):
    !!     optical absorption cross sections of each aerosol particle
    !!   - \b aero_scatter_cross_sect (unit m^2, dim \c aero_particle):
    !!     optical scattering cross sections of each aerosol particle
    !!   - \b aero_asymmetry (dimensionless, dim \c aero_particle): optical
    !!     asymmetry parameters of each aerosol particle
    !!   - \b aero_refract_shell_real (dimensionless, dim \c aero_particle):
    !!     real part of the refractive indices of the shell of each
    !!     aerosol particle
    !!   - \b aero_refract_shell_imag (dimensionless, dim \c aero_particle):
    !!     imaginary part of the refractive indices of the shell of each
    !!     aerosol particle
    !!   - \b aero_refract_core_real (dimensionless, dim \c aero_particle):
    !!     real part of the refractive indices of the core of each
    !!     aerosol particle
    !!   - \b aero_refract_core_imag (dimensionless, dim \c aero_particle):
    !!     imaginary part of the refractive indices of the core of each
    !!     aerosol particle
    !!   - \b aero_core_vol (unit m^3, dim \c aero_particle): volume of the
    !!     optical cores of each aerosol particle
    !!   - \b aero_water_hyst_leg (dim \c aero_particle): integers
    !!     specifying which leg of the water hysteresis curve each
    !!     particle is on, using the MOSAIC numbering convention
    !!   - \b aero_comp_vol (unit m^3, dim \c aero_particle): computational
    !!     volume containing each particle
    !!   - \b aero_id (dim \c aero_particle): unique ID number of each
    !!     aerosol particle
    !!   - \b aero_least_create_time (unit s, dim \c aero_particle): least
    !!     (earliest) creation time of any original constituent particles
    !!     that coagulated to form each particle, measured from the start
    !!     of the simulation - a particle is said to be created when it
    !!     first enters the simulation (by emissions, dilution, etc.)
    !!   - \b aero_greatest_create_time (unit s, dim \c
    !!     aero_particle): greatest (latest) creation time of any
    !!     original constituent particles that coagulated to form each
    !!     particle, measured from the start of the simulation - a
    !!     particle is said to be created when it first enters the
    !!     simulation (by emissions, dilution, etc.)

    call aero_weight_array_output_netcdf(aero_state%awa, ncid)
    
    call aero_data_netcdf_dim_aero_species(aero_data, ncid, &
         dimid_aero_species)
    call aero_data_netcdf_dim_aero_source(aero_data, ncid, &
         dimid_aero_source)

    if (aero_state%apa%n_part > 0) then
       call aero_state_netcdf_dim_aero_particle(aero_state, ncid, &
            dimid_aero_particle)

       ! FIXME: replace this loop with statements like
       ! aero_n_orig_part = aero_state%apa%particle%n_orig_part
       do i_part = 1,aero_state%apa%n_part
          particle => aero_state%apa%particle(i_part)
          aero_particle_mass(i_part, :) = particle%vol * aero_data%density
          aero_n_orig_part(i_part, :) = particle%n_orig_part
          aero_weight_group(i_part) = particle%weight_group
          aero_water_hyst_leg(i_part) = particle%water_hyst_leg
          aero_comp_vol(i_part) &
               = 1d0 / aero_state_particle_num_conc(aero_state, particle)
          aero_id(i_part) = particle%id
          aero_least_create_time(i_part) = particle%least_create_time
          aero_greatest_create_time(i_part) = particle%greatest_create_time
          if (record_optical) then
             aero_absorb_cross_sect(i_part) = particle%absorb_cross_sect
             aero_scatter_cross_sect(i_part) = particle%scatter_cross_sect
             aero_asymmetry(i_part) = particle%asymmetry
             aero_refract_shell_real(i_part) = real(particle%refract_shell)
             aero_refract_shell_imag(i_part) = aimag(particle%refract_shell)
             aero_refract_core_real(i_part) = real(particle%refract_core)
             aero_refract_core_imag(i_part) = aimag(particle%refract_core)
             aero_core_vol(i_part) = particle%core_vol
          end if
       end do
       call pmc_nc_write_real_2d(ncid, aero_particle_mass, &
            "aero_particle_mass", (/ dimid_aero_particle, &
            dimid_aero_species /), unit="kg", &
            long_name="constituent masses of each aerosol particle")
       call pmc_nc_write_integer_2d(ncid, aero_n_orig_part, &
            "aero_n_orig_part", (/ dimid_aero_particle, &
            dimid_aero_source /), &
            long_name="number of original constituent particles from " &
            // "each source that coagulated to form each aerosol particle")
       call pmc_nc_write_integer_1d(ncid, aero_weight_group, &
            "aero_weight_group", (/ dimid_aero_particle /), &
            long_name="weight group number of each aerosol particle")
       call pmc_nc_write_integer_1d(ncid, aero_water_hyst_leg, &
            "aero_water_hyst_leg", (/ dimid_aero_particle /), &
            long_name="leg of the water hysteresis curve leg of each "&
            // "aerosol particle")
       call pmc_nc_write_real_1d(ncid, aero_comp_vol, &
            "aero_comp_vol", (/ dimid_aero_particle /), unit="m^3", &
            long_name="computational volume containing each particle")
       call pmc_nc_write_integer_1d(ncid, aero_id, &
            "aero_id", (/ dimid_aero_particle /), &
            long_name="unique ID number of each aerosol particle")
       call pmc_nc_write_real_1d(ncid, aero_least_create_time, &
            "aero_least_create_time", (/ dimid_aero_particle /), unit="s", &
            long_name="least creation time of each aerosol particle", &
            description="least (earliest) creation time of any original " &
            // "constituent particles that coagulated to form each " &
            // "particle, measured from the start of the simulation")
       call pmc_nc_write_real_1d(ncid, aero_greatest_create_time, &
            "aero_greatest_create_time", (/ dimid_aero_particle /), &
            unit="s", &
            long_name="greatest creation time of each aerosol particle", &
            description="greatest (latest) creation time of any original " &
            // "constituent particles that coagulated to form each " &
            // "particle, measured from the start of the simulation")
       if (record_optical) then
          call pmc_nc_write_real_1d(ncid, aero_absorb_cross_sect, &
               "aero_absorb_cross_sect", (/ dimid_aero_particle /), &
               unit="m^2", &
               long_name="optical absorption cross sections of each " &
               // "aerosol particle")
          call pmc_nc_write_real_1d(ncid, aero_scatter_cross_sect, &
               "aero_scatter_cross_sect", (/ dimid_aero_particle /), &
               unit="m^2", &
               long_name="optical scattering cross sections of each " &
               // "aerosol particle")
          call pmc_nc_write_real_1d(ncid, aero_asymmetry, &
               "aero_asymmetry", (/ dimid_aero_particle /), unit="1", &
               long_name="optical asymmetry parameters of each " &
               // "aerosol particle")
          call pmc_nc_write_real_1d(ncid, aero_refract_shell_real, &
               "aero_refract_shell_real", (/ dimid_aero_particle /), &
               unit="1", &
               long_name="real part of the refractive indices of the " &
               // "shell of each aerosol particle")
          call pmc_nc_write_real_1d(ncid, aero_refract_shell_imag, &
               "aero_refract_shell_imag", (/ dimid_aero_particle /), &
               unit="1", &
               long_name="imaginary part of the refractive indices of " &
               // "the shell of each aerosol particle")
          call pmc_nc_write_real_1d(ncid, aero_refract_core_real, &
               "aero_refract_core_real", (/ dimid_aero_particle /), &
               unit="1", &
               long_name="real part of the refractive indices of the core " &
               // "of each aerosol particle")
          call pmc_nc_write_real_1d(ncid, aero_refract_core_imag, &
               "aero_refract_core_imag", (/ dimid_aero_particle /), &
               unit="1", &
               long_name="imaginary part of the refractive indices of " &
               // "the core of each aerosol particle")
          call pmc_nc_write_real_1d(ncid, aero_core_vol, &
               "aero_core_vol", (/ dimid_aero_particle /), unit="m^3", &
               long_name="volume of the optical cores of each " &
               // "aerosol particle")
       end if
    end if

    !> \page output_format_aero_removed Output File Format: Aerosol Particle Removal Information
    !!
    !! When an aerosol particle is introduced into the simulation it
    !! is assigned a unique ID number. This ID number will persist
    !! over time, allowing tracking of a paticular particle's
    !! evolution. If the \c record_removals variable in the input spec
    !! file is \c yes, then the every time a particle is removed from
    !! the simulation its removal will be recorded in the removal
    !! information.
    !!
    !! The removal information written at timestep \c n contains
    !! information about every particle ID that is present at time
    !! <tt>(n - 1)</tt> but not present at time \c n.
    !!
    !! The removal information is always written in the output files,
    !! even if no particles were removed in the previous
    !! timestep. Unfortunately, NetCDF files cannot contain arrays of
    !! length 0. In the case of no particles being removed, the \c
    !! aero_removed dimension will be set to 1 and
    !! <tt>aero_removed_action(1)</tt> will be 0 (\c AERO_INFO_NONE).
    !!
    !! When two particles coagulate, the ID number of the combined
    !! particle will be the ID particle of the largest constituent, if
    !! possible (weighting functions can make this impossible to
    !! achieve). A given particle ID may thus be lost due to
    !! coagulation (if the resulting combined particle has a different
    !! ID), or the ID may be preserved (as the ID of the combined
    !! particle). Only if the ID is lost will the particle be recorded
    !! in the removal information, and in this case
    !! <tt>aero_removed_action(i)</tt> will be 2 (\c AERO_INFO_COAG)
    !! and <tt>aero_removed_other_id(i)</tt> will be the ID number of
    !! the combined particle.
    !!
    !! The aerosol removal information NetCDF dimensions are:
    !!   - \b aero_removed: number of aerosol particles removed from the
    !!     simulation during the previous timestep (or 1, as described
    !!     above)
    !!
    !! The aerosol removal information NetCDF variables are:
    !!   - \b aero_removed (dim \c aero_removed): dummy dimension variable
    !!     (no useful value)
    !!   - \b aero_removed_id (dim \c aero_removed): the ID number of each
    !!     removed particle
    !!   - \b aero_removed_action (dim \c aero_removed): the reasons for
    !!     removal for each particle, with values:
    !!     - 0 (\c AERO_INFO_NONE): no information (invalid entry)
    !!     - 1 (\c AERO_INFO_DILUTION): particle was removed due to dilution
    !!       with outside air
    !!     - 2 (\c AERO_INFO_COAG): particle was removed due to coagulation
    !!     - 3 (\c AERO_INFO_HALVED): particle was removed due to halving of
    !!       the aerosol population
    !!     - 4 (\c AERO_INFO_WEIGHT): particle was removed due to adjustments
    !!       in the particle's weighting function
    !!   - \b aero_removed_other_id (dim \c aero_removed): the ID number of
    !!     the combined particle formed by coagulation, if the removal reason
    !!     was coagulation (2, \c AERO_INFO_COAG). May be 0, if the new
    !!     coagulated particle was not created due to weighting.

    ! FIXME: move this to aero_info_array.F90, together with
    ! aero_state_netcdf_dim_aero_removed() ?
    if (record_removals) then
       call aero_state_netcdf_dim_aero_removed(aero_state, ncid, &
            dimid_aero_removed)
       if (aero_state%aero_info_array%n_item >= 1) then
          do i_remove = 1,aero_state%aero_info_array%n_item
             aero_removed_id(i_remove) = &
                  aero_state%aero_info_array%aero_info(i_remove)%id
             aero_removed_action(i_remove) = &
                  aero_state%aero_info_array%aero_info(i_remove)%action
             aero_removed_other_id(i_remove) = &
                  aero_state%aero_info_array%aero_info(i_remove)%other_id
          end do
       else
          aero_removed_id(1) = 0
          aero_removed_action(1) = AERO_INFO_NONE
          aero_removed_other_id(1) = 0
       end if
       call pmc_nc_write_integer_1d(ncid, aero_removed_id, &
            "aero_removed_id", (/ dimid_aero_removed /), &
            long_name="ID of removed particles")
       call pmc_nc_write_integer_1d(ncid, aero_removed_action, &
            "aero_removed_action", (/ dimid_aero_removed /), &
            long_name="reason for particle removal", &
            description="valid is 0 (invalid entry), 1 (removed due to " &
            // "dilution), 2 (removed due to coagulation -- combined " &
            // "particle ID is in \c aero_removed_other_id), 3 (removed " &
            // "due to populating halving), or 4 (removed due to " &
            // "weighting changes")
       call pmc_nc_write_integer_1d(ncid, aero_removed_other_id, &
            "aero_removed_other_id", (/ dimid_aero_removed /), &
            long_name="ID of other particle involved in removal", &
            description="if <tt>aero_removed_action(i)</tt> is 2 " &
            // "(due to coagulation), then " &
            // "<tt>aero_removed_other_id(i)</tt> is the ID of the " &
            // "resulting combined particle, or 0 if the new particle " &
            // "was not created")
    end if

  end subroutine aero_state_output_netcdf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read full state.
  subroutine aero_state_input_netcdf(aero_state, ncid, aero_data)
    
    !> aero_state to read.
    type(aero_state_t), intent(inout) :: aero_state
    !> NetCDF file ID, in data mode.
    integer, intent(in) :: ncid
    !> aero_data structure.
    type(aero_data_t), intent(in) :: aero_data

    integer :: dimid_aero_particle, dimid_aero_removed, n_info_item, n_part
    integer :: i_bin, i_part_in_bin, i_part, i_remove, status
    type(aero_particle_t) :: aero_particle
    character(len=1000) :: name

    real(kind=dp), allocatable :: aero_particle_mass(:,:)
    integer, allocatable :: aero_n_orig_part(:,:)
    integer, allocatable :: aero_weight_group(:)
    real(kind=dp), allocatable :: aero_absorb_cross_sect(:)
    real(kind=dp), allocatable :: aero_scatter_cross_sect(:)
    real(kind=dp), allocatable :: aero_asymmetry(:)
    real(kind=dp), allocatable :: aero_refract_shell_real(:)
    real(kind=dp), allocatable :: aero_refract_shell_imag(:)
    real(kind=dp), allocatable :: aero_refract_core_real(:)
    real(kind=dp), allocatable :: aero_refract_core_imag(:)
    real(kind=dp), allocatable :: aero_core_vol(:)
    integer, allocatable :: aero_water_hyst_leg(:)
    real(kind=dp), allocatable :: aero_comp_vol(:)
    integer, allocatable :: aero_id(:)
    real(kind=dp), allocatable :: aero_least_create_time(:)
    real(kind=dp), allocatable :: aero_greatest_create_time(:)
    integer, allocatable :: aero_removed_id(:)
    integer, allocatable :: aero_removed_action(:)
    integer, allocatable :: aero_removed_other_id(:)

    status = nf90_inq_dimid(ncid, "aero_particle", dimid_aero_particle)
    if (status == NF90_EBADDIM) then
       ! no aero_particle dimension means no particles present
       call aero_state_deallocate(aero_state)
       call aero_state_allocate(aero_state)
       call aero_weight_array_input_netcdf(aero_state%awa, ncid)
       return
    end if
    call pmc_nc_check(status)
    call pmc_nc_check(nf90_Inquire_Dimension(ncid, dimid_aero_particle, &
         name, n_part))

    allocate(aero_particle_mass(n_part, aero_data%n_spec))
    allocate(aero_n_orig_part(n_part, aero_data%n_source))
    allocate(aero_weight_group(n_part))
    allocate(aero_absorb_cross_sect(n_part))
    allocate(aero_scatter_cross_sect(n_part))
    allocate(aero_asymmetry(n_part))
    allocate(aero_refract_shell_real(n_part))
    allocate(aero_refract_shell_imag(n_part))
    allocate(aero_refract_core_real(n_part))
    allocate(aero_refract_core_imag(n_part))
    allocate(aero_core_vol(n_part))
    allocate(aero_water_hyst_leg(n_part))
    allocate(aero_comp_vol(n_part))
    allocate(aero_id(n_part))
    allocate(aero_least_create_time(n_part))
    allocate(aero_greatest_create_time(n_part))

    call pmc_nc_read_real_2d(ncid, aero_particle_mass, &
         "aero_particle_mass")
    call pmc_nc_read_integer_2d(ncid, aero_n_orig_part, &
         "aero_n_orig_part")
    call pmc_nc_read_integer_1d(ncid, aero_weight_group, &
         "aero_weight_group")
    call pmc_nc_read_real_1d(ncid, aero_absorb_cross_sect, &
         "aero_absorb_cross_sect", must_be_present=.false.)
    call pmc_nc_read_real_1d(ncid, aero_scatter_cross_sect, &
         "aero_scatter_cross_sect", must_be_present=.false.)
    call pmc_nc_read_real_1d(ncid, aero_asymmetry, &
         "aero_asymmetry", must_be_present=.false.)
    call pmc_nc_read_real_1d(ncid, aero_refract_shell_real, &
         "aero_refract_shell_real", must_be_present=.false.)
    call pmc_nc_read_real_1d(ncid, aero_refract_shell_imag, &
         "aero_refract_shell_imag", must_be_present=.false.)
    call pmc_nc_read_real_1d(ncid, aero_refract_core_real, &
         "aero_refract_core_real", must_be_present=.false.)
    call pmc_nc_read_real_1d(ncid, aero_refract_core_imag, &
         "aero_refract_core_imag", must_be_present=.false.)
    call pmc_nc_read_real_1d(ncid, aero_core_vol, &
         "aero_core_vol", must_be_present=.false.)
    call pmc_nc_read_integer_1d(ncid, aero_water_hyst_leg, &
         "aero_water_hyst_leg")
    call pmc_nc_read_real_1d(ncid, aero_comp_vol, &
         "aero_comp_vol")
    call pmc_nc_read_integer_1d(ncid, aero_id, &
         "aero_id")
    call pmc_nc_read_real_1d(ncid, aero_least_create_time, &
         "aero_least_create_time")
    call pmc_nc_read_real_1d(ncid, aero_greatest_create_time, &
         "aero_greatest_create_time")

    call aero_state_deallocate(aero_state)
    call aero_state_allocate(aero_state)
    
    call aero_weight_array_input_netcdf(aero_state%awa, ncid)
    aero_state%n_part_ideal = 0d0

    call aero_particle_allocate_size(aero_particle, aero_data%n_spec, &
         aero_data%n_source)
    do i_part = 1,n_part
       aero_particle%vol = aero_particle_mass(i_part, :) / aero_data%density
       aero_particle%n_orig_part = aero_n_orig_part(i_part, :)
       aero_particle%weight_group = aero_weight_group(i_part)
       aero_particle%absorb_cross_sect = aero_absorb_cross_sect(i_part)
       aero_particle%scatter_cross_sect = aero_scatter_cross_sect(i_part)
       aero_particle%asymmetry = aero_asymmetry(i_part)
       aero_particle%refract_shell = &
            cmplx(aero_refract_shell_real(i_part), &
            aero_refract_shell_imag(i_part), kind=dc)
       aero_particle%refract_core = cmplx(aero_refract_core_real(i_part), &
            aero_refract_core_imag(i_part), kind=dc)
       aero_particle%core_vol = aero_core_vol(i_part)
       aero_particle%water_hyst_leg = aero_water_hyst_leg(i_part)
       aero_particle%id = aero_id(i_part)
       aero_particle%least_create_time = aero_least_create_time(i_part)
       aero_particle%greatest_create_time = aero_greatest_create_time(i_part)

       call assert(314368871, almost_equal(aero_comp_vol(i_part), &
            1d0 / aero_weight_array_num_conc(aero_state%awa, aero_particle)))

       call aero_state_add_particle(aero_state, aero_particle)
    end do
    call aero_particle_deallocate(aero_particle)

    deallocate(aero_particle_mass)
    deallocate(aero_n_orig_part)
    deallocate(aero_weight_group)
    deallocate(aero_absorb_cross_sect)
    deallocate(aero_scatter_cross_sect)
    deallocate(aero_asymmetry)
    deallocate(aero_refract_shell_real)
    deallocate(aero_refract_shell_imag)
    deallocate(aero_refract_core_real)
    deallocate(aero_refract_core_imag)
    deallocate(aero_core_vol)
    deallocate(aero_water_hyst_leg)
    deallocate(aero_comp_vol)
    deallocate(aero_id)
    deallocate(aero_least_create_time)
    deallocate(aero_greatest_create_time)

    status = nf90_inq_dimid(ncid, "aero_removed", dimid_aero_removed)
    if ((status /= NF90_NOERR) .and. (status /= NF90_EBADDIM)) then
       call pmc_nc_check(status)
    end if
    if (status == NF90_NOERR) then
       call pmc_nc_check(nf90_Inquire_Dimension(ncid, dimid_aero_removed, &
            name, n_info_item))

       allocate(aero_removed_id(max(n_info_item,1)))
       allocate(aero_removed_action(max(n_info_item,1)))
       allocate(aero_removed_other_id(max(n_info_item,1)))

       call pmc_nc_read_integer_1d(ncid, aero_removed_id, &
            "aero_removed_id")
       call pmc_nc_read_integer_1d(ncid, aero_removed_action, &
            "aero_removed_action")
       call pmc_nc_read_integer_1d(ncid, aero_removed_other_id, &
            "aero_removed_other_id")

       if ((n_info_item > 1) .or. (aero_removed_id(1) /= 0)) then
          call aero_info_array_enlarge_to(aero_state%aero_info_array, &
               n_info_item)
          do i_remove = 1,n_info_item
             aero_state%aero_info_array%aero_info(i_remove)%id &
                  = aero_removed_id(i_remove)
             aero_state%aero_info_array%aero_info(i_remove)%action &
                  = aero_removed_action(i_remove)
             aero_state%aero_info_array%aero_info(i_remove)%other_id &
                  = aero_removed_other_id(i_remove)
          end do
       end if

       deallocate(aero_removed_id)
       deallocate(aero_removed_action)
       deallocate(aero_removed_other_id)
    end if

  end subroutine aero_state_input_netcdf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Sorts the particles if necessary.
  subroutine aero_state_sort(aero_state, bin_grid, all_procs_same)

    !> Aerosol state to sort.
    type(aero_state_t), intent(inout) :: aero_state
    !> Bin grid.
    type(bin_grid_t), optional, intent(in) :: bin_grid
    !> Whether all processors should use the same bin grid.
    logical, optional, intent(in) :: all_procs_same

    call aero_sorted_remake_if_needed(aero_state%aero_sorted, aero_state%apa, &
         aero_state%valid_sort, size(aero_state%awa%weight), bin_grid, &
         all_procs_same)
    aero_state%valid_sort = .true.

  end subroutine aero_state_sort
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Check that the sorted data is consistent.
  subroutine aero_state_check_sort(aero_state)

    !> Aerosol state to check.
    type(aero_state_t), intent(in) :: aero_state

    logical, parameter :: continue_on_error = .false.

    integer :: i_part, i_bin

    if (.not. aero_state%valid_sort) then
       write(0,*) 'SORTED CHECK ERROR: SORT NOT VALID'
       return
    end if

    call aero_sorted_check(aero_state%aero_sorted, aero_state%apa, &
         size(aero_state%awa%weight), continue_on_error)

  end subroutine aero_state_check_sort
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
end module pmc_aero_state
