! Copyright (C) 2005-2021 Nicole Riemer and Matthew West
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
  !! array (the \c apa data member), with a sorting into size bins and
  !! weight groups/classes possibly stored in the \c aero_sorted data
  !! member (if \c valid_sort is true).
  !!
  !! Every time we remove particles we keep track of the particle ID
  !! and the action performed in the aero_info_array_t structure. This
  !! is typically cleared each time we output data to disk.
  type aero_state_t
     !> Linear array of particles.
     type(aero_particle_array_t) :: apa
     !> Sorting of particles into size bins and weight groups/classes.
     type(aero_sorted_t) :: aero_sorted
     !> Whether the \c aero_sorted is a correct sorting.
     logical :: valid_sort
     !> Weighting functions.
     type(aero_weight_array_t) :: awa
     !> Ideal number of computational particles, per weight group and class.
     real(kind=dp), allocatable :: n_part_ideal(:, :)
     !> Information on removed particles.
     type(aero_info_array_t) :: aero_info_array
#ifdef PMC_USE_CAMP
     !> CAMP update number conc cookie
     type(aero_rep_update_data_single_particle_number_t) :: update_number
#endif
  end type aero_state_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Return the current number of particles.
  elemental integer function aero_state_n_part(aero_state)

    !> Aerosol state.
    type(aero_state_t), intent(in) :: aero_state

    aero_state_n_part = aero_particle_array_n_part(aero_state%apa)

  end function aero_state_n_part

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read the specification for a weighting type from a spec file.
  subroutine spec_file_read_aero_state_weighting_type(file, weighting_type, &
       exponent)

    !> Spec file.
    type(spec_file_t), intent(inout) :: file
    !> Aerosol weighting scheme.
    integer, intent(out) :: weighting_type
    !> Exponent for power-law weighting (only used if \c weight_type
    !> is \c AERO_STATE_WEIGHT_POWER).
    real(kind=dp), intent(out) :: exponent

    character(len=SPEC_LINE_MAX_VAR_LEN) :: weighting_name

    !> \page input_format_weight_type Input File Format: Type of aerosol size distribution weighting functions.
    !!
    !! The weighting function is specified by the parameters:
    !!   - \b weight_type (string): the type of weighting function ---
    !!     must be one of: \c flat for flat weighting, \c flat_source for
    !!     flat weighting by source, \c power for power weighting,
    !!     \c power_source for power source weighting, \c nummass for number
    !!     and mass weighting, and \c nummass_source for number and mass
    !!     weighting by source. If \c weight_type is \c power or \c
    !!     power_source then the next parameter must also be provided:
    !!     - \b weighting_exponent (real): the exponent for \c power or
    !!       \c power_source. Setting the \c exponent to 0 is equivalent to no
    !!       weighting, while setting the exponent negative uses more
    !!       computational particles at larger diameters and setting the
    !!       exponent positive uses more computaitonal partilces at smaller
    !!       diameters; in practice exponents between 0 and -3 are most useful.
    !!
    !! See also:
    !!   - \ref spec_file_format --- the input file text format

    call spec_file_read_string(file, 'weight_type', weighting_name)
    if (trim(weighting_name) == 'flat') then
       weighting_type = AERO_STATE_WEIGHT_FLAT
    elseif (trim(weighting_name) == 'power') then
       weighting_type = AERO_STATE_WEIGHT_POWER
       call spec_file_read_real(file, 'weighting_exponent', exponent)
    elseif (trim(weighting_name) == 'nummass') then
       weighting_type = AERO_STATE_WEIGHT_NUMMASS
    elseif (trim(weighting_name) == 'flat_source') then
       weighting_type = AERO_STATE_WEIGHT_FLAT_SOURCE
    elseif (trim(weighting_name) == 'power_source') then
       weighting_type = AERO_STATE_WEIGHT_POWER_SOURCE
       call spec_file_read_real(file, 'weighting_exponent', exponent)
    elseif (trim(weighting_name) == 'nummass_source') then
       weighting_type = AERO_STATE_WEIGHT_NUMMASS_SOURCE
    else
       call spec_file_die_msg(920321729, file, &
            "Unknown weighting type: " // trim(weighting_name))
    end if

  end subroutine spec_file_read_aero_state_weighting_type

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Copies weighting information for an \c aero_state.
  subroutine aero_state_copy_weight(aero_state_from, aero_state_to)

    !> Reference aerosol.
    type(aero_state_t), intent(in) :: aero_state_from
    !> Already allocated.
    type(aero_state_t), intent(inout) :: aero_state_to

    aero_state_to%awa = aero_state_from%awa

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
    select case(weight_type)
    case(AERO_STATE_WEIGHT_NONE)
    case(AERO_STATE_WEIGHT_FLAT)
       call aero_weight_array_set_flat(aero_state%awa, 1)
    case(AERO_STATE_WEIGHT_POWER)
       call assert_msg(656670336, present(exponent), &
            "exponent parameter required for AERO_STATE_WEIGHT_POWER")
       call aero_weight_array_set_power(aero_state%awa, 1, exponent)
    case(AERO_STATE_WEIGHT_NUMMASS)
       call aero_weight_array_set_nummass(aero_state%awa, 1)
    case(AERO_STATE_WEIGHT_FLAT_SOURCE)
       call aero_weight_array_set_flat(aero_state%awa, &
            aero_data_n_source(aero_data))
    case(AERO_STATE_WEIGHT_POWER_SOURCE)
       call assert_msg(102143848, present(exponent), &
            "exponent parameter required for AERO_STATE_WEIGHT_POWER")
       call aero_weight_array_set_power(aero_state%awa, &
            aero_data_n_source(aero_data), exponent)
    case(AERO_STATE_WEIGHT_NUMMASS_SOURCE)
       call aero_weight_array_set_nummass(aero_state%awa, &
            aero_data_n_source(aero_data))
    case default
       call die_msg(969076992, "unknown weight_type: " &
            // trim(integer_to_string(weight_type)))
    end select

  end subroutine aero_state_set_weight

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Set the ideal number of particles to the given value. The \c
  !> aero_state%%awa must be already set correctly.
  subroutine aero_state_set_n_part_ideal(aero_state, n_part)

    !> Aerosol state (with \c aero_state%%awa set).
    type(aero_state_t), intent(inout) :: aero_state
    !> Ideal total number of particles.
    real(kind=dp), intent(in) :: n_part

    integer :: n_group, n_class

    n_group = aero_weight_array_n_group(aero_state%awa)
    n_class = aero_weight_array_n_class(aero_state%awa)
    if (allocated(aero_state%n_part_ideal)) then
       deallocate(aero_state%n_part_ideal)
    end if
    allocate(aero_state%n_part_ideal(n_group, n_class))
    aero_state%n_part_ideal = n_part / real(n_group * n_class, kind=dp)

  end subroutine aero_state_set_n_part_ideal

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Determine the appropriate weight class for a source.
  integer function aero_state_weight_class_for_source(aero_state, source)

    !> Aerosol state.
    type(aero_state_t), intent(in) :: aero_state
    !> Source to find the class for.
    integer, intent(in) :: source

    integer :: n_class

    call assert(932390238, source >= 1)
    n_class = aero_weight_array_n_class(aero_state%awa)
    ! we are either using i_class = i_source or always i_class = n_class = 1
    if (n_class > 1) then
       call assert(765048788, source <= n_class)
       aero_state_weight_class_for_source = source
    else
       aero_state_weight_class_for_source = 1
    end if

  end function aero_state_weight_class_for_source

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the total number of particles in an aerosol distribution.
  integer function aero_state_total_particles(aero_state, i_group, i_class)

    !> Aerosol state.
    type(aero_state_t), intent(in) :: aero_state
    !> Weight group.
    integer, optional, intent(in) :: i_group
    !> Weight class.
    integer, optional, intent(in) :: i_class

    integer :: i_part

    if (present(i_group)) then
       call assert(908743823, present(i_class))
       if (aero_state%valid_sort) then
          aero_state_total_particles &
               = integer_varray_n_entry( &
               aero_state%aero_sorted%group_class%inverse(i_group, i_class))
       else
          ! FIXME: should we just sort?
          aero_state_total_particles = 0
          do i_part = 1,aero_state_n_part(aero_state)
             if ((aero_state%apa%particle(i_part)%weight_group == i_group) &
                  .and. &
                  (aero_state%apa%particle(i_part)%weight_class == i_class)) &
                  then
                aero_state_total_particles = aero_state_total_particles + 1
             end if
          end do
       end if
    else
       aero_state_total_particles = aero_state_n_part(aero_state)
    end if

  end function aero_state_total_particles

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the total number of particles across all processes.
  integer function aero_state_total_particles_all_procs(aero_state, i_group, &
       i_class)

    !> Aerosol state.
    type(aero_state_t), intent(in) :: aero_state
    !> Weight group.
    integer, optional, intent(in) :: i_group
    !> Weight class.
    integer, optional, intent(in) :: i_class

    call pmc_mpi_allreduce_sum_integer(&
         aero_state_total_particles(aero_state, i_group, i_class), &
         aero_state_total_particles_all_procs)

  end function aero_state_total_particles_all_procs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Resets an aero_state to have zero particles per bin.
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
  subroutine aero_state_add_particle(aero_state, aero_particle, aero_data, &
       allow_resort)

    !> Aerosol state.
    type(aero_state_t), intent(inout) :: aero_state
    !> Particle to add.
    type(aero_particle_t), intent(in) :: aero_particle
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Whether to allow a resort due to the add.
    logical, optional, intent(in) :: allow_resort

    if (aero_state%valid_sort) then
       call aero_sorted_add_particle(aero_state%aero_sorted, aero_state%apa, &
            aero_particle, aero_data, allow_resort)
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
       i_bin, i_class, aero_particle)

    !> Aerosol state.
    type(aero_state_t), intent(inout) :: aero_state
    !> Bin number to remove particle from.
    integer, intent(in) :: i_bin
    !> Weight class to remove particle from.
    integer, intent(in) :: i_class
    !> Removed particle.
    type(aero_particle_t), intent(inout) :: aero_particle

    integer :: i_entry, i_part

    call assert(742996300, aero_state%valid_sort)
    call assert(392182617, integer_varray_n_entry( &
         aero_state%aero_sorted%size_class%inverse(i_bin, i_class)) > 0)
    i_entry = pmc_rand_int(integer_varray_n_entry( &
         aero_state%aero_sorted%size_class%inverse(i_bin, i_class)))
    i_part = aero_state%aero_sorted%size_class%inverse(i_bin, &
         i_class)%entry(i_entry)
    aero_particle = aero_state%apa%particle(i_part)
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
  subroutine aero_state_dup_particle(aero_state, aero_data, i_part, &
       n_part_mean, random_weight_group)

    !> Aerosol state.
    type(aero_state_t), intent(inout) :: aero_state
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Particle number.
    integer, intent(in) :: i_part
    !> Mean number of resulting particles.
    real(kind=dp), intent(in) :: n_part_mean
    !> Whether particle copies should be placed in a randomly chosen
    !> weight group.
    logical, optional, intent(in) :: random_weight_group

    integer :: n_copies, i_dup, new_group
    type(aero_particle_t) :: new_aero_particle
    type(aero_info_t) :: aero_info

    n_copies = prob_round(n_part_mean)
    if (n_copies == 0) then
       aero_info%id = aero_state%apa%particle(i_part)%id
       aero_info%action = AERO_INFO_WEIGHT
       aero_info%other_id = 0
       call aero_state_remove_particle_with_info(aero_state, &
            i_part, aero_info)
    elseif (n_copies > 1) then
       do i_dup = 1,(n_copies - 1)
          new_aero_particle = aero_state%apa%particle(i_part)
          call aero_particle_new_id(new_aero_particle)
          if (present(random_weight_group)) then
             if (random_weight_group) then
                new_group &
                     = aero_weight_array_rand_group(aero_state%awa, &
                     aero_state%apa%particle(i_part)%weight_class, &
                     aero_particle_radius(aero_state%apa%particle(i_part), &
                          aero_data))
                call aero_particle_set_weight(new_aero_particle, new_group)
             end if
          end if
          call aero_state_add_particle(aero_state, new_aero_particle, &
               aero_data)
       end do
    end if

  end subroutine aero_state_dup_particle

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> The number concentration of a single particle (m^{-3}).
  real(kind=dp) function aero_state_particle_num_conc(aero_state, &
       aero_particle, aero_data)

    !> Aerosol state containing the particle.
    type(aero_state_t), intent(in) :: aero_state
    !> Aerosol particle.
    type(aero_particle_t), intent(in) :: aero_particle
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data

    aero_state_particle_num_conc &
         = aero_weight_array_num_conc(aero_state%awa, aero_particle, &
         aero_data)

  end function aero_state_particle_num_conc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Save the correct number concentrations for later use by
  !> aero_state_reweight().
  subroutine aero_state_num_conc_for_reweight(aero_state, aero_data, &
       reweight_num_conc)

    !> Aerosol state.
    type(aero_state_t), intent(in) :: aero_state
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Number concentrations for later use by aero_state_reweight().
    real(kind=dp), intent(out) &
         :: reweight_num_conc(aero_state_n_part(aero_state))

    integer :: i_part

    do i_part = 1,aero_state_n_part(aero_state)
       reweight_num_conc(i_part) &
            = aero_weight_array_single_num_conc(aero_state%awa, &
            aero_state%apa%particle(i_part), aero_data)
    end do

  end subroutine aero_state_num_conc_for_reweight

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Reweight all particles after their constituent volumes have been
  !> altered.
  !!
  !! The pattern for use should be like:
  !! <pre>
  !! call aero_state_num_conc_for_reweight(aero_state, aero_data,
  !!      reweight_num_conc)
  !! ... alter particle species volumes in aero_state ...
  !! call aero_state_reweight(aero_state, aero_data, reweight_num_conc)
  !! </pre>
  subroutine aero_state_reweight(aero_state, aero_data, reweight_num_conc)

    !> Aerosol state.
    type(aero_state_t), intent(inout) :: aero_state
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Number concentrations previously computed by
    !> aero_state_num_conc_for_reweight().
    real(kind=dp), intent(in) &
         :: reweight_num_conc(aero_state_n_part(aero_state))

    integer :: i_part, i_group, i_class
    real(kind=dp) :: n_part_old(size(aero_state%awa%weight, 1), &
         size(aero_state%awa%weight, 2))
    real(kind=dp) :: n_part_new(size(aero_state%awa%weight, 1), &
         size(aero_state%awa%weight, 2))
    real(kind=dp) :: old_num_conc, new_num_conc, n_part_mean

    ! find average number of new particles in each weight group, if
    ! weight is not changed
    n_part_old = 0d0
    n_part_new = 0d0
    do i_part = 1,aero_state_n_part(aero_state)
       old_num_conc = reweight_num_conc(i_part)
       new_num_conc = aero_weight_array_single_num_conc(aero_state%awa, &
            aero_state%apa%particle(i_part), aero_data)
       n_part_mean = old_num_conc / new_num_conc
       i_group = aero_state%apa%particle(i_part)%weight_group
       i_class = aero_state%apa%particle(i_part)%weight_class
       n_part_new(i_group, i_class) = n_part_new(i_group, i_class) &
            + n_part_mean
       n_part_old(i_group, i_class) = n_part_old(i_group, i_class) + 1d0
    end do

    ! alter weight to leave the number of computational particles
    ! per weight bin unchanged
    do i_group = 1,size(aero_state%awa%weight, 1)
       do i_class = 1,size(aero_state%awa%weight, 2)
          if (n_part_old(i_group, i_class) == 0d0) cycle
          call aero_weight_scale(aero_state%awa%weight(i_group, i_class), &
               n_part_new(i_group, i_class) / n_part_old(i_group, i_class))
       end do
    end do

    ! work backwards so any additions and removals will only affect
    ! particles that we've already dealt with
    do i_part = aero_state_n_part(aero_state),1,-1
       old_num_conc = reweight_num_conc(i_part)
       new_num_conc &
            = aero_weight_array_single_num_conc(aero_state%awa, &
            aero_state%apa%particle(i_part), aero_data)
       n_part_mean = old_num_conc / new_num_conc
       call aero_state_dup_particle(aero_state, aero_data, i_part, &
            n_part_mean)
    end do

  end subroutine aero_state_reweight

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> <tt>aero_state += aero_state_delta</tt>, including combining the
  !> weights, so the new concentration is the weighted average of the
  !> two concentrations.
  subroutine aero_state_add(aero_state, aero_state_delta, aero_data)

    !> Aerosol state.
    type(aero_state_t), intent(inout) :: aero_state
    !> Increment.
    type(aero_state_t), intent(in) :: aero_state_delta
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data

    call aero_state_add_particles(aero_state, aero_state_delta, aero_data)
    call aero_weight_array_combine(aero_state%awa, aero_state_delta%awa)

  end subroutine aero_state_add

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> <tt>aero_state += aero_state_delta</tt>, with the weight
  !> of \c aero_state left unchanged, so the new concentration is the
  !> sum of the two concentrations, computed with \c aero_state%%awa.
  subroutine aero_state_add_particles(aero_state, aero_state_delta, aero_data)

    !> Aerosol state.
    type(aero_state_t), intent(inout) :: aero_state
    !> Increment.
    type(aero_state_t), intent(in) :: aero_state_delta
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data

    integer :: i_part, i_bin

    do i_part = 1,aero_state_delta%apa%n_part
       call aero_state_add_particle(aero_state, &
            aero_state_delta%apa%particle(i_part), aero_data)
    end do
    call aero_info_array_add(aero_state%aero_info_array, &
         aero_state_delta%aero_info_array)

  end subroutine aero_state_add_particles

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Change the weight if necessary to ensure that the addition of
  !> about \c n_add computational particles will give the correct
  !> final particle number.
  subroutine aero_state_prepare_weight_for_add(aero_state, aero_data, &
       i_group, i_class, n_add, allow_doubling, allow_halving)

    !> Aero state to add to.
    type(aero_state_t), intent(inout) :: aero_state
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Weight group number to add to.
    integer, intent(in) :: i_group
    !> Weight class number to add to.
    integer, intent(in) :: i_class
    !> Approximate number of particles to be added at current weighting.
    real(kind=dp), intent(in) :: n_add
    !> Whether to allow doubling of the population.
    logical, intent(in) :: allow_doubling
    !> Whether to allow halving of the population.
    logical, intent(in) :: allow_halving

    integer :: global_n_part, n_group, n_class
    real(kind=dp) :: mean_n_part, n_part_new, weight_ratio
    real(kind=dp) :: n_part_ideal_local_group

    n_group = aero_weight_array_n_group(aero_state%awa)
    n_class = aero_weight_array_n_class(aero_state%awa)
    global_n_part = aero_state_total_particles_all_procs(aero_state, &
         i_group, i_class)
    mean_n_part = real(global_n_part, kind=dp) / real(pmc_mpi_size(), kind=dp)
    n_part_new = mean_n_part + n_add
    if (n_part_new == 0d0) return
    n_part_ideal_local_group = aero_state%n_part_ideal(i_group, i_class) &
         / real(pmc_mpi_size(), kind=dp)
    if ((n_part_new < n_part_ideal_local_group / 2d0) &
         .or. (n_part_new > n_part_ideal_local_group * 2d0)) &
         then
       weight_ratio = n_part_new / n_part_ideal_local_group
       call aero_state_scale_weight(aero_state, aero_data, i_group, &
            i_class, weight_ratio, allow_doubling, allow_halving)
    end if

  end subroutine aero_state_prepare_weight_for_add

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Generates a Poisson sample of an \c aero_dist, adding to \c
  !> aero_state, with the given sample proportion.
  subroutine aero_state_add_aero_dist_sample(aero_state, aero_data, &
       aero_dist, sample_prop, create_time, allow_doubling, allow_halving, &
       n_part_add)

    !> Aero state to add to.
    type(aero_state_t), intent(inout) :: aero_state
    !> Aero data values.
    type(aero_data_t), intent(in) :: aero_data
    !> Distribution to sample.
    type(aero_dist_t), intent(in) :: aero_dist
    !> Fraction to sample (1).
    real(kind=dp), intent(in) :: sample_prop
    !> Creation time for new particles (s).
    real(kind=dp), intent(in) :: create_time
    !> Whether to allow doubling of the population.
    logical, intent(in) :: allow_doubling
    !> Whether to allow halving of the population.
    logical, intent(in) :: allow_halving
    !> Number of particles added.
    integer, intent(out), optional :: n_part_add

    real(kind=dp) :: n_samp_avg, radius, total_vol
    real(kind=dp) :: vols(aero_data_n_spec(aero_data))
    integer :: n_samp, i_mode, i_samp, i_group, i_class, n_group, n_class
    type(aero_particle_t) :: aero_particle

    n_group = size(aero_state%awa%weight, 1)
    n_class = size(aero_state%awa%weight, 2)
    if (present(n_part_add)) then
       n_part_add = 0
    end if
    do i_group = 1,n_group
       do i_mode = 1,aero_dist_n_mode(aero_dist)
          i_class = aero_state_weight_class_for_source(aero_state, &
               aero_dist%mode(i_mode)%source)

          ! adjust weight if necessary
          n_samp_avg = sample_prop * aero_mode_number(aero_dist%mode(i_mode), &
               aero_state%awa%weight(i_group, i_class))
          call aero_state_prepare_weight_for_add(aero_state, aero_data, &
               i_group, i_class, n_samp_avg, allow_doubling, allow_halving)
          if (n_samp_avg == 0d0) cycle

          ! sample particles
          n_samp_avg = sample_prop * aero_mode_number(aero_dist%mode(i_mode), &
               aero_state%awa%weight(i_group, i_class))
          n_samp = rand_poisson(n_samp_avg)
          if (present(n_part_add)) then
             n_part_add = n_part_add + n_samp
          end if
          do i_samp = 1,n_samp
             call aero_particle_zero(aero_particle, aero_data)
             call aero_mode_sample_radius(aero_dist%mode(i_mode), aero_data, &
                  aero_state%awa%weight(i_group, i_class), radius)
             total_vol = aero_data_rad2vol(aero_data, radius)
             call aero_mode_sample_vols(aero_dist%mode(i_mode), total_vol, &
                  vols)
             call aero_particle_set_vols(aero_particle, vols)
             call aero_particle_new_id(aero_particle)
             call aero_particle_set_weight(aero_particle, i_group, i_class)
             call aero_particle_set_create_time(aero_particle, create_time)
             call aero_particle_set_source(aero_particle, &
                  aero_dist%mode(i_mode)%source)
             call aero_state_add_particle(aero_state, aero_particle, aero_data)
          end do
       end do
    end do

  end subroutine aero_state_add_aero_dist_sample

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Choose a random particle from the aero_state.
  subroutine aero_state_rand_particle(aero_state, i_part)

    !> Original state.
    type(aero_state_t), intent(in) :: aero_state
    !> Chosen random particle number.
    integer, intent(out) :: i_part

    call assert(950725003, aero_state_n_part(aero_state) > 0)
    i_part = pmc_rand_int(aero_state_n_part(aero_state))

  end subroutine aero_state_rand_particle

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Generates a random sample by removing particles from
  !> aero_state_from and adding them to aero_state_to, which must be
  !> already allocated (and should have its weight set).
  !!
  !! None of the weights are altered by this sampling, making this the
  !! equivalent of aero_state_add_particles().
  subroutine aero_state_sample_particles(aero_state_from, aero_state_to, &
       aero_data, sample_prob, removal_action)

    !> Original state.
    type(aero_state_t), intent(inout) :: aero_state_from
    !> Destination state.
    type(aero_state_t), intent(inout) :: aero_state_to
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
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
    call aero_state_zero(aero_state_to)
    call aero_state_copy_weight(aero_state_from, aero_state_to)
    n_transfer = rand_binomial(aero_state_total_particles(aero_state_from), &
         sample_prob)
    i_transfer = 0
    do while (i_transfer < n_transfer)
       if (aero_state_total_particles(aero_state_from) <= 0) exit
       call aero_state_rand_particle(aero_state_from, i_part)
       num_conc_from = aero_weight_array_num_conc(aero_state_from%awa, &
            aero_state_from%apa%particle(i_part), aero_data)
       num_conc_to = aero_weight_array_num_conc(aero_state_to%awa, &
            aero_state_from%apa%particle(i_part), aero_data)

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
               aero_state_from%apa%particle(i_part), aero_data)
       end if
       if (do_remove) then
          if (removal_action /= AERO_INFO_NONE) then
             aero_info%id = aero_state_from%apa%particle(i_part)%id
             aero_info%action = removal_action
             call aero_state_remove_particle_with_info(aero_state_from, &
                  i_part, aero_info)
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
  !> weight as well. This is the equivalent of aero_state_add().
  subroutine aero_state_sample(aero_state_from, aero_state_to, &
       aero_data, sample_prob, removal_action)

    !> Original state.
    type(aero_state_t), intent(inout) :: aero_state_from
    !> Destination state (previous contents will be lost).
    type(aero_state_t), intent(inout) :: aero_state_to
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Probability of sampling each particle.
    real(kind=dp), intent(in) :: sample_prob
    !> Action for removal (see pmc_aero_info module for action
    !> parameters). Set to AERO_INFO_NONE to not log removal.
    integer, intent(in) :: removal_action

    integer :: n_transfer, i_transfer, i_part
    logical :: do_add, do_remove, overwrite_to
    real(kind=dp) :: num_conc_from, num_conc_to
    type(aero_info_t) :: aero_info

    call assert(393205561, (sample_prob >= 0d0) .and. (sample_prob <= 1d0))
    call aero_state_zero(aero_state_to)
    call aero_state_copy_weight(aero_state_from, aero_state_to)
    call aero_weight_array_normalize(aero_state_to%awa)
    n_transfer = rand_binomial(aero_state_total_particles(aero_state_from), &
         sample_prob)
    do i_transfer = 1,n_transfer
       if (aero_state_total_particles(aero_state_from) <= 0) exit
       call aero_state_rand_particle(aero_state_from, i_part)

       call aero_state_add_particle(aero_state_to, &
            aero_state_from%apa%particle(i_part), aero_data)
       if (removal_action /= AERO_INFO_NONE) then
          aero_info%id = aero_state_from%apa%particle(i_part)%id
          aero_info%action = removal_action
          call aero_state_remove_particle_with_info(aero_state_from, &
               i_part, aero_info)
       else
          call aero_state_remove_particle_no_info(aero_state_from, &
               i_part)
       end if
    end do
    overwrite_to = .true.
    call aero_weight_array_shift(aero_state_from%awa, aero_state_to%awa, &
         sample_prob, overwrite_to)

  end subroutine aero_state_sample

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Create binned number and mass arrays.
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
    real(kind=dp) :: factor

    aero_binned%num_conc = 0d0
    aero_binned%vol_conc = 0d0
    do i_part = 1,aero_state_n_part(aero_state)
       i_bin = bin_grid_find(bin_grid, &
            aero_particle_radius(aero_state%apa%particle(i_part), aero_data))
       if ((i_bin < 1) .or. (i_bin > bin_grid_size(bin_grid))) then
          call warn_msg(980232449, "particle ID " &
               // trim(integer_to_string(aero_state%apa%particle(i_part)%id)) &
               // " outside of bin_grid, discarding")
       else
          factor = aero_weight_array_num_conc(aero_state%awa, &
               aero_state%apa%particle(i_part), aero_data) &
               / bin_grid%widths(i_bin)
          aero_binned%vol_conc(i_bin,:) = aero_binned%vol_conc(i_bin,:) &
               + aero_state%apa%particle(i_part)%vol * factor
          aero_binned%num_conc(i_bin) = aero_binned%num_conc(i_bin) &
               + factor
       end if
    end do

  end subroutine aero_state_to_binned

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the IDs of all particles.
  function aero_state_ids(aero_state)

    !> Aerosol state.
    type(aero_state_t), intent(in) :: aero_state

    !> Return value.
    integer :: aero_state_ids(aero_state_n_part(aero_state))

    integer :: i_part

    do i_part = 1,aero_state_n_part(aero_state)
       aero_state_ids(i_part) = aero_state%apa%particle(i_part)%id
    end do

  end function aero_state_ids

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the diameters of all particles.
  function aero_state_diameters(aero_state, aero_data, include, exclude)

    !> Aerosol state.
    type(aero_state_t), intent(in) :: aero_state
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Species names to include in the diameter.
    character(len=*), optional, intent(in) :: include(:)
    !> Species names to exclude in the diameter.
    character(len=*), optional, intent(in) :: exclude(:)

    !> Return diameters array (m).
    real(kind=dp) :: aero_state_diameters(aero_state_n_part(aero_state))

    !> Per-particle volume of included components
    real(kind=dp) :: volumes(aero_state_n_part(aero_state))

    volumes = aero_state_volumes(aero_state, aero_data, include, exclude)
    aero_state_diameters = rad2diam(sphere_vol2rad(volumes))

  end function aero_state_diameters

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the dry diameters of all particles.
  function aero_state_dry_diameters(aero_state, aero_data)

    !> Aerosol state.
    type(aero_state_t), intent(in) :: aero_state
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data

    !> Return value (m).
    real(kind=dp) :: aero_state_dry_diameters(aero_state_n_part(aero_state))

    aero_state_dry_diameters = aero_particle_dry_diameter( &
         aero_state%apa%particle(1:aero_state_n_part(aero_state)), aero_data)

  end function aero_state_dry_diameters

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the mobility diameters of all particles.
  function aero_state_mobility_diameters(aero_state, aero_data, env_state)

    !> Aerosol state.
    type(aero_state_t), intent(in) :: aero_state
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Environment state.
    type(env_state_t), intent(in) :: env_state

    !> Return value (m).
    real(kind=dp) :: aero_state_mobility_diameters( &
         aero_state_n_part(aero_state))

    integer :: i_part

    do i_part = 1,aero_state_n_part(aero_state)
       aero_state_mobility_diameters(i_part) &
            = aero_particle_mobility_diameter( &
            aero_state%apa%particle(i_part), aero_data, env_state)
    end do

  end function aero_state_mobility_diameters

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the volumes of all particles.
  !!
  !! If \c include is specified then only those species are included
  !! in computing the volumes. If \c exclude is specified then all
  !! species except those species are included. If both \c include and
  !! \c exclude arguments are specified then only those species in \c
  !! include but not in \c exclude are included.
  function aero_state_volumes(aero_state, aero_data, include, exclude)

    !> Aerosol state.
    type(aero_state_t), intent(in) :: aero_state
    !> Aerosol data.
    type(aero_data_t), optional, intent(in) :: aero_data
    !> Species names to include in the mass.
    character(len=*), optional, intent(in) :: include(:)
    !> Species names to exclude in the mass.
    character(len=*), optional, intent(in) :: exclude(:)

    !> Return volumes array (m^3).
    real(kind=dp) :: aero_state_volumes(aero_state_n_part(aero_state))

    logical, allocatable :: use_species(:)
    integer :: i_name, i_spec

    if ((.not. present(include)) .and. (.not. present(exclude))) then
       aero_state_volumes = aero_particle_volume( &
            aero_state%apa%particle(1:aero_state_n_part(aero_state)))
    else
       call assert_msg(599558703, present(aero_data), &
            "must provide 'aero_data' if using 'include' or 'exclude'")
       allocate(use_species(aero_data_n_spec(aero_data)))
       if (present(include)) then
          use_species = .false.
          do i_name = 1, size(include)
             i_spec = aero_data_spec_by_name(aero_data, include(i_name))
             call assert_msg(111852070, i_spec > 0, &
                  "unknown species: " // trim(include(i_name)))
             use_species(i_spec) = .true.
          end do
       else
          use_species = .true.
       end if
       if (present(exclude)) then
          do i_name = 1, size(exclude)
             i_spec = aero_data_spec_by_name(aero_data, exclude(i_name))
             call assert_msg(182075590, i_spec > 0, &
                  "unknown species: " // trim(exclude(i_name)))
             use_species(i_spec) = .false.
          end do
       end if
       aero_state_volumes = 0d0
       do i_spec = 1,aero_data_n_spec(aero_data)
          if (use_species(i_spec)) then
             aero_state_volumes = aero_state_volumes &
                  + aero_particle_species_volume( &
                  aero_state%apa%particle(1:aero_state_n_part(aero_state)), &
                  i_spec)
          end if
       end do
    end if

  end function aero_state_volumes

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the masses of all particles.
  !!
  !! If \c include is specified then only those species are included
  !! in computing the masses. If \c exclude is specified then all
  !! species except those species are included. If both \c include and
  !! \c exclude arguments are specified then only those species in \c
  !! include but not in \c exclude are included.
  function aero_state_masses(aero_state, aero_data, include, exclude)

    !> Aerosol state.
    type(aero_state_t), intent(in) :: aero_state
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Species names to include in the mass.
    character(len=*), optional, intent(in) :: include(:)
    !> Species names to exclude in the mass.
    character(len=*), optional, intent(in) :: exclude(:)

    !> Return masses array (kg).
    real(kind=dp) :: aero_state_masses(aero_state_n_part(aero_state))

    logical :: use_species(aero_data_n_spec(aero_data))
    integer :: i_name, i_spec

    if ((.not. present(include)) .and. (.not. present(exclude))) then
       aero_state_masses = aero_particle_mass( &
            aero_state%apa%particle(1:aero_state_n_part(aero_state)), &
            aero_data)
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
       aero_state_masses = 0d0
       do i_spec = 1,aero_data_n_spec(aero_data)
          if (use_species(i_spec)) then
             aero_state_masses = aero_state_masses &
                  + aero_particle_species_mass( &
                  aero_state%apa%particle(1:aero_state_n_part(aero_state)), &
                  i_spec, aero_data)
          end if
       end do
    end if

  end function aero_state_masses

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the number concentrations of all particles.
  function aero_state_num_concs(aero_state, aero_data)

    !> Aerosol state.
    type(aero_state_t), intent(in) :: aero_state
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data

    !> Return number concentrations array (m^{-3}).
    real(kind=dp) :: aero_state_num_concs(aero_state_n_part(aero_state))

    integer :: i_part

    do i_part = 1,aero_state_n_part(aero_state)
       aero_state_num_concs(i_part) &
            = aero_state_particle_num_conc(aero_state, &
            aero_state%apa%particle(i_part), aero_data)
    end do

  end function aero_state_num_concs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the total number concentration.
  real(kind=dp) function aero_state_total_num_conc(aero_state, aero_data)

    !> Aerosol state.
    type(aero_state_t), intent(in) :: aero_state
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data

    integer :: i_part

    aero_state_total_num_conc = 0d0
    do i_part = 1,aero_state_n_part(aero_state)
       aero_state_total_num_conc = aero_state_total_num_conc &
            + aero_state_particle_num_conc(aero_state, &
            aero_state%apa%particle(i_part), aero_data)
    end do

  end function aero_state_total_num_conc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the mass-entropies of all particles.
  !!
  !! If \c include is specified then only those species are included
  !! in computing the entropy. If \c exclude is specified then all
  !! species except those species are included. If both \c include and
  !! \c exclude arguments are specified then only those species in \c
  !! include but not in \c exclude are included. If \c group is
  !! specified then the species are divided into two sets, given by
  !! those in the group and those not in the group. The entropy is
  !! then computed using the total mass of each set. Alternatively \c
  !! groups can be specified, which lists several groups of species.
  function aero_state_mass_entropies(aero_state, aero_data, include, exclude, &
       group, groups)

    !> Aerosol state.
    type(aero_state_t), intent(in) :: aero_state
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Species names to include in the mass.
    character(len=*), optional :: include(:)
    !> Species names to exclude in the mass.
    character(len=*), optional :: exclude(:)
    !> Species names to group together.
    character(len=*), optional :: group(:)
    !> Sets of species names to group together.
    character(len=*), optional :: groups(:,:)

    !> Return value.
    real(kind=dp) :: aero_state_mass_entropies(aero_state_n_part(aero_state))

    logical :: use_species(aero_data_n_spec(aero_data))
    logical :: group_species(aero_data_n_spec(aero_data))
    integer :: i_name, i_spec, i_part, i_group, n_group
    integer :: species_group_numbers(aero_data_n_spec(aero_data))
    real(kind=dp) :: group_mass, non_group_mass, mass
    real(kind=dp), allocatable :: group_masses(:)

    if (present(include)) then
       use_species = .false.
       do i_name = 1, size(include)
          i_spec = aero_data_spec_by_name(aero_data, include(i_name))
          call assert_msg(890212002, i_spec > 0, &
               "unknown species: " // trim(include(i_name)))
          use_species(i_spec) = .true.
       end do
    else
       use_species = .true.
    end if
    if (present(exclude)) then
       do i_name = 1, size(exclude)
          i_spec = aero_data_spec_by_name(aero_data, exclude(i_name))
          call assert_msg(859945006, i_spec > 0, &
               "unknown species: " // trim(exclude(i_name)))
          use_species(i_spec) = .false.
       end do
    end if
    if (present(group)) then
       group_species = .false.
       do i_name = 1, size(group)
          i_spec = aero_data_spec_by_name(aero_data, group(i_name))
          call assert_msg(376359046, i_spec > 0, &
               "unknown species: " // trim(group(i_name)))
          group_species(i_spec) = .true.
       end do
       do i_part = 1,aero_state_n_part(aero_state)
          group_mass = 0d0
          non_group_mass = 0d0
          do i_spec = 1,aero_data_n_spec(aero_data)
             if (use_species(i_spec)) then
                mass = aero_particle_species_mass( &
                     aero_state%apa%particle(i_part), i_spec, aero_data)
                if (group_species(i_spec)) then
                   group_mass = group_mass + mass
                else
                   non_group_mass = non_group_mass + mass
                end if
             end if
          end do
          aero_state_mass_entropies(i_part) &
               = entropy([group_mass, non_group_mass])
       end do
    else if (present(groups)) then
       call assert_msg(161633285, .not. present(include), &
            "cannot specify both 'include' and 'groups' arguments")
       call assert_msg(273540426, .not. present(exclude), &
            "cannot specify both 'exclude' and 'groups' arguments")
       call assert_msg(499993914, .not. present(group), &
            "cannot specify both 'group' and 'groups' arguments")

       n_group = size(groups, 1)
       ! species_group_numbers(i_spec) will give the group number for
       ! each species, or zero for non-included
       species_group_numbers = 0 ! extra species go into zero (unuesd) group
       do i_group = 1, n_group
          do i_name = 1, size(groups, 2)
             if (len_trim(groups(i_group, i_name)) > 0) then
                i_spec = aero_data_spec_by_name(aero_data, &
                     groups(i_group, i_name))
                call assert_msg(926066826, i_spec > 0, &
                     "unknown species: " // trim(groups(i_group, i_name)))
                if (use_species(i_spec)) then
                   species_group_numbers(i_spec) = i_group
                end if
             end if
          end do
       end do

       allocate(group_masses(n_group))
       do i_part = 1,aero_state_n_part(aero_state)
          group_masses = 0d0
          do i_spec = 1,aero_data_n_spec(aero_data)
             mass = aero_particle_species_mass( &
                  aero_state%apa%particle(i_part), i_spec, aero_data)
             i_group = species_group_numbers(i_spec)
             if (i_group > 0) then
                group_masses(i_group) = group_masses(i_group) + mass
             end if
          end do
          aero_state_mass_entropies(i_part) = entropy(group_masses)
       end do
    else
       do i_part = 1,aero_state_n_part(aero_state)
          aero_state_mass_entropies(i_part) = entropy(pack( &
               aero_particle_species_masses(aero_state%apa%particle(i_part), &
               aero_data), use_species))
       end do
    end if

  end function aero_state_mass_entropies

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the mixing state metrics of the population.
  !!
  !! If \c include is specified then only those species are included
  !! in computing the entropies. If \c exclude is specified then all
  !! species except those species are included. If both \c include and
  !! \c exclude arguments are specified then only those species in \c
  !! include but not in \c exclude are included. If \c group is
  !! specified then the species are divided into two sets, given by
  !! those in the group and those not in the group. The entropies are
  !! then computed using the total mass of each set. Alternatively \c
  !! groups can be specified, which lists several groups of
  !! species. If \c groups is provided, only species explicitly listed
  !! will be included.
  subroutine aero_state_mixing_state_metrics(aero_state, aero_data, d_alpha, &
       d_gamma, chi, include, exclude, group, groups)

    !> Aerosol state.
    type(aero_state_t), intent(in) :: aero_state
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Average particle diversity.
    real(kind=dp), intent(out) :: d_alpha
    !> Bulk diversity.
    real(kind=dp), intent(out) :: d_gamma
    !> Mixing state index.
    real(kind=dp), intent(out) :: chi
    !> Species names to include in the mass.
    character(len=*), optional :: include(:)
    !> Species names to exclude in the mass.
    character(len=*), optional :: exclude(:)
    !> Species names to group together.
    character(len=*), optional :: group(:)
    !> Sets of species names to group together.
    character(len=*), optional :: groups(:,:)

    real(kind=dp), allocatable :: entropies(:), entropies_of_avg_part(:)
    real(kind=dp), allocatable :: masses(:), num_concs(:), &
         num_concs_of_avg_part(:), masses_of_avg_part(:)
    type(aero_state_t) :: aero_state_averaged
    type(bin_grid_t) :: avg_bin_grid

    ! per-particle masses need to take groups into account

    if (present(groups)) then
       call assert_msg(726652236, .not. present(include), &
            "cannot specify both 'include' and 'groups' arguments")
       call assert_msg(891097454, .not. present(exclude), &
            "cannot specify both 'exclude' and 'groups' arguments")
       call assert_msg(938789093, .not. present(group), &
            "cannot specify both 'group' and 'groups' arguments")
       masses = aero_state_masses(aero_state, aero_data, &
            include=pack(groups, len_trim(groups) > 0))
    else
       masses = aero_state_masses(aero_state, aero_data, include, exclude)
    end if

    ! other per-particle properties
    num_concs = aero_state_num_concs(aero_state, aero_data)
    entropies = aero_state_mass_entropies(aero_state, aero_data, &
         include, exclude, group, groups)

    d_alpha = exp(sum(entropies * masses * num_concs) &
             / sum(masses * num_concs))

    ! per-particle properties of averaged particles
    call bin_grid_make(avg_bin_grid, BIN_GRID_TYPE_LOG, 1, 1d-30, 1d10)
    aero_state_averaged = aero_state
    call aero_state_bin_average_comp(aero_state_averaged, avg_bin_grid, &
         aero_data)
    num_concs_of_avg_part = aero_state_num_concs(aero_state_averaged, &
         aero_data)
    if (present(groups)) then
       masses_of_avg_part = aero_state_masses(aero_state_averaged, aero_data, &
            include=pack(groups, len_trim(groups) > 0))
    else
       masses_of_avg_part = aero_state_masses(aero_state_averaged, aero_data, &
            include, exclude)
    end if
    entropies_of_avg_part = aero_state_mass_entropies(aero_state_averaged, &
         aero_data, include, exclude, group, groups)

    d_gamma = exp(sum(entropies_of_avg_part * masses_of_avg_part &
         * num_concs_of_avg_part) &
         / sum(masses_of_avg_part * num_concs_of_avg_part))

    chi = (d_alpha - 1) / (d_gamma - 1)

  end subroutine aero_state_mixing_state_metrics

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the approximate critical relative humidity for all particles (1).
  function aero_state_approx_crit_rel_humids(aero_state, aero_data, env_state)

    !> Aerosol state.
    type(aero_state_t), intent(in) :: aero_state
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Environment state.
    type(env_state_t), intent(in) :: env_state

    !> Return value.
    real(kind=dp) &
         :: aero_state_approx_crit_rel_humids(aero_state_n_part(aero_state))

    integer :: i_part

    do i_part = 1,aero_state_n_part(aero_state)
       aero_state_approx_crit_rel_humids(i_part) = &
            aero_particle_approx_crit_rel_humid( &
            aero_state%apa%particle(i_part), aero_data, env_state)
    end do

  end function aero_state_approx_crit_rel_humids

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the critical relative humidity for all particles (1).
  function aero_state_crit_rel_humids(aero_state, aero_data, env_state)

    !> Aerosol state.
    type(aero_state_t), intent(in) :: aero_state
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Environment state.
    type(env_state_t), intent(in) :: env_state

    !> Return value.
    real(kind=dp) :: aero_state_crit_rel_humids(aero_state_n_part(aero_state))

    integer :: i_part

    do i_part = 1,aero_state_n_part(aero_state)
       aero_state_crit_rel_humids(i_part) = aero_particle_crit_rel_humid( &
            aero_state%apa%particle(i_part), aero_data, env_state)
    end do

  end function aero_state_crit_rel_humids

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
    real(kind=dp) :: factor

    aero_binned%num_conc = 0d0
    aero_binned%vol_conc = 0d0
    do i_part = 1,aero_state_n_part(aero_state)
       i_bin = bin_grid_find(bin_grid, &
            aero_particle_solute_radius(aero_state%apa%particle(i_part), &
            aero_data))
       if ((i_bin < 1) .or. (i_bin > bin_grid_size(bin_grid))) then
          call warn_msg(503871022, "particle ID " &
               // trim(integer_to_string(aero_state%apa%particle(i_part)%id)) &
               // " outside of bin_grid, discarding")
       else
          factor = aero_weight_array_num_conc(aero_state%awa, &
               aero_state%apa%particle(i_part), aero_data) &
               / bin_grid%widths(i_bin)
          aero_binned%vol_conc(i_bin,:) = aero_binned%vol_conc(i_bin,:) &
               + aero_state%apa%particle(i_part)%vol * factor
          aero_binned%num_conc(i_bin) = aero_binned%num_conc(i_bin) &
               + factor
       end if
    end do

  end subroutine aero_state_to_binned_dry

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Doubles number of particles in the given weight group.
  subroutine aero_state_double(aero_state, aero_data, i_group, i_class)

    !> Aerosol state.
    type(aero_state_t), intent(inout) :: aero_state
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Weight group to double.
    integer, intent(in) :: i_group
    !> Weight class to double.
    integer, intent(in) :: i_class

    integer :: i_part
    type(aero_particle_t) :: aero_particle

    do i_part = 1,aero_state_n_part(aero_state)
       if ((aero_state%apa%particle(i_part)%weight_group == i_group) &
            .and. (aero_state%apa%particle(i_part)%weight_class == i_class)) &
            then
          aero_particle = aero_state%apa%particle(i_part)
          call aero_particle_new_id(aero_particle)
          call aero_state_add_particle(aero_state, aero_particle, aero_data)
       end if
    end do
    aero_state%valid_sort = .false.
    call aero_weight_scale(aero_state%awa%weight(i_group, i_class), 0.5d0)

  end subroutine aero_state_double

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Remove approximately half of the particles in the given weight group.
  subroutine aero_state_halve(aero_state, i_group, i_class)

    !> Aerosol state.
    type(aero_state_t), intent(inout) :: aero_state
    !> Weight group to halve.
    integer, intent(in) :: i_group
    !> Weight class to halve.
    integer, intent(in) :: i_class

    integer :: i_part
    type(aero_info_t) :: aero_info

    do i_part = aero_state_n_part(aero_state),1,-1
       if ((aero_state%apa%particle(i_part)%weight_group == i_group) &
            .and. (aero_state%apa%particle(i_part)%weight_class == i_class)) &
            then
          if (pmc_random() < 0.5d0) then
             aero_info%id = aero_state%apa%particle(i_part)%id
             aero_info%action = AERO_INFO_HALVED
             call aero_state_remove_particle_with_info(aero_state, i_part, &
                  aero_info)
          end if
       end if
    end do
    call aero_weight_scale(aero_state%awa%weight(i_group, i_class), 2d0)

  end subroutine aero_state_halve

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Double or halve the particle population in each weight group to
  !> maintain close to \c n_part_ideal particles per process,
  !> allocated equally amongst the weight groups.
  subroutine aero_state_rebalance(aero_state, aero_data, allow_doubling, &
       allow_halving, initial_state_warning)

    !> Aerosol state.
    type(aero_state_t), intent(inout) :: aero_state
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Whether to allow doubling of the population.
    logical, intent(in) :: allow_doubling
    !> Whether to allow halving of the population.
    logical, intent(in) :: allow_halving
    !> Whether to warn due to initial state doubling/halving.
    logical, intent(in) :: initial_state_warning

    integer :: i_group, i_class, n_group, n_class, global_n_part

    n_group = size(aero_state%awa%weight, 1)
    n_class = size(aero_state%awa%weight, 2)

    ! if we have less than half the maximum number of particles then
    ! double until we fill up the array
    if (allow_doubling) then
       do i_group = 1,n_group
          do i_class = 1,n_class
             global_n_part &
                  = aero_state_total_particles_all_procs(aero_state, i_group, &
                  i_class)
             do while ((real(global_n_part, kind=dp) &
                  < aero_state%n_part_ideal(i_group, i_class) / 2d0) &
                  .and. (global_n_part > 0))
                if (initial_state_warning) then
                   call warn_msg(716882783, "doubling particles in initial " &
                        // "condition")
                end if
                call aero_state_double(aero_state, aero_data, i_group, i_class)
                global_n_part &
                     = aero_state_total_particles_all_procs(aero_state, &
                     i_group, i_class)
             end do
          end do
       end do
    end if

    ! same for halving if we have too many particles
    if (allow_halving) then
       do i_group = 1,n_group
          do i_class = 1,n_class
             do while (real(aero_state_total_particles_all_procs(aero_state, &
                  i_group, i_class), kind=dp) &
                  > aero_state%n_part_ideal(i_group, i_class) * 2d0)
                if (initial_state_warning) then
                   call warn_msg(661936373, &
                        "halving particles in initial condition")
                end if
                call aero_state_halve(aero_state, i_group, i_class)
             end do
          end do
       end do
    end if

  end subroutine aero_state_rebalance

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Scale the weighting of the given group/class by the given ratio,
  !> altering particle number as necessary to preserve the number
  !> concentration.
  subroutine aero_state_scale_weight(aero_state, aero_data, i_group, &
       i_class, weight_ratio, allow_doubling, allow_halving)

    !> Aerosol state.
    type(aero_state_t), intent(inout) :: aero_state
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Weight group number.
    integer, intent(in) :: i_group
    !> Weight class number.
    integer, intent(in) :: i_class
    !> Ratio of <tt>new_weight / old_weight</tt>.
    real(kind=dp), intent(in) :: weight_ratio
    !> Whether to allow doubling of the population.
    logical, intent(in) :: allow_doubling
    !> Whether to allow halving of the population.
    logical, intent(in) :: allow_halving

    real(kind=dp) :: ratio
    integer :: i_part, i_remove, n_remove, i_entry, n_part
    type(aero_info_t) :: aero_info

    ! We could use the ratio < 1 case unconditionally, but that would
    ! have higher variance for the ratio > 1 case than the current
    ! scheme.

    call aero_state_sort(aero_state, aero_data)
    n_part = integer_varray_n_entry( &
         aero_state%aero_sorted%group_class%inverse(i_group, i_class))

    if ((weight_ratio > 1d0) .and. (allow_halving .or. (n_part == 0))) then
       call aero_weight_scale(aero_state%awa%weight(i_group, i_class), &
            weight_ratio)
       n_remove = prob_round(real(n_part, kind=dp) &
            * (1d0 - 1d0 / weight_ratio))
       do i_remove = 1,n_remove
          i_entry = pmc_rand_int(integer_varray_n_entry( &
               aero_state%aero_sorted%group_class%inverse(i_group, i_class)))
          i_part = aero_state%aero_sorted%group_class%inverse(i_group, &
               i_class)%entry(i_entry)
          aero_info%id = aero_state%apa%particle(i_part)%id
          aero_info%action = AERO_INFO_HALVED
          call aero_state_remove_particle(aero_state, i_part, .true., &
               aero_info)
       end do
    elseif ((weight_ratio < 1d0) &
         .and. (allow_doubling .or. (n_part == 0))) then
       call aero_weight_scale(aero_state%awa%weight(i_group, i_class), &
            weight_ratio)
       do i_entry = n_part,1,-1
          i_part = aero_state%aero_sorted%group_class%inverse(i_group, &
               i_class)%entry(i_entry)
          call aero_state_dup_particle(aero_state, aero_data, i_part, &
               1d0 / weight_ratio)
       end do
    end if

  end subroutine aero_state_scale_weight

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
               aero_state_sends(i_proc + 1), aero_data, &
               prob_transfer_given_not_transferred, AERO_INFO_NONE)
          prob_not_transferred = prob_not_transferred - prob_transfer
       end if
    end do

    ! exchange the particles
    call aero_state_mpi_alltoall(aero_state_sends, aero_state_recvs)

    ! process the received particles
    do i_proc = 0,(n_proc - 1)
       if (i_proc /= rank) then
          call aero_state_add(aero_state, aero_state_recvs(i_proc + 1), &
               aero_data)
       end if
    end do

    ! cleanup
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
  subroutine aero_state_bin_average_comp(aero_state, bin_grid, aero_data)

    !> Aerosol state to average.
    type(aero_state_t), intent(inout) :: aero_state
    !> Bin grid to average within.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data

    real(kind=dp) :: species_volume_conc(aero_data_n_spec(aero_data))
    real(kind=dp) :: total_volume_conc, particle_volume, num_conc
    integer :: i_bin, i_class, i_entry, i_part, i_spec

    call aero_state_sort(aero_state, aero_data, bin_grid)

    do i_bin = 1,bin_grid_size(bin_grid)
       species_volume_conc = 0d0
       total_volume_conc = 0d0
       do i_class = 1,size(aero_state%awa%weight, 2)
          do i_entry = 1,integer_varray_n_entry( &
               aero_state%aero_sorted%size_class%inverse(i_bin, i_class))
             i_part = aero_state%aero_sorted%size_class%inverse(i_bin, &
                  i_class)%entry(i_entry)
             num_conc = aero_weight_array_num_conc(aero_state%awa, &
                  aero_state%apa%particle(i_part), aero_data)
             particle_volume = aero_particle_volume( &
                  aero_state%apa%particle(i_part))
             species_volume_conc = species_volume_conc &
                  + num_conc * aero_state%apa%particle(i_part)%vol
             total_volume_conc = total_volume_conc + num_conc * particle_volume
          end do
       end do
       do i_class = 1,size(aero_state%awa%weight, 2)
          do i_entry = 1,integer_varray_n_entry( &
               aero_state%aero_sorted%size_class%inverse(i_bin, i_class))
             i_part = aero_state%aero_sorted%size_class%inverse(i_bin, &
                  i_class)%entry(i_entry)
             particle_volume = aero_particle_volume( &
                  aero_state%apa%particle(i_part))
             aero_state%apa%particle(i_part)%vol &
                  = particle_volume * species_volume_conc / total_volume_conc
          end do
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
       bin_center, preserve_number)

    !> Aerosol state to average.
    type(aero_state_t), intent(inout) :: aero_state
    !> Bin grid to average within.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
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
    integer :: i_bin, i_class, i_entry, i_part, i_bisect, n_part
    logical :: monotone_increasing, monotone_decreasing

    call aero_state_sort(aero_state, aero_data, bin_grid)

    do i_bin = 1,bin_grid_size(bin_grid)
       do i_class = 1,size(aero_state%awa%weight, 2)
          if (integer_varray_n_entry( &
               aero_state%aero_sorted%size_class%inverse(i_bin, i_class)) &
               == 0) then
             cycle
          end if

          n_part = integer_varray_n_entry( &
               aero_state%aero_sorted%size_class%inverse(i_bin, i_class))
          total_num_conc = 0d0
          total_volume_conc = 0d0
          do i_entry = 1,integer_varray_n_entry( &
               aero_state%aero_sorted%size_class%inverse(i_bin, i_class))
             i_part = aero_state%aero_sorted%size_class%inverse(i_bin, &
                  i_class)%entry(i_entry)
             num_conc = aero_weight_array_num_conc(aero_state%awa, &
                  aero_state%apa%particle(i_part), aero_data)
             total_num_conc = total_num_conc + num_conc
             particle_volume = aero_particle_volume( &
                  aero_state%apa%particle(i_part))
             total_volume_conc = total_volume_conc &
                  + num_conc * particle_volume
          end do

          ! determine the new_particle_volume for all particles in this bin
          if (bin_center) then
             new_particle_volume = aero_data_rad2vol(aero_data, &
                  bin_grid%centers(i_bin))
          elseif (aero_weight_array_check_flat(aero_state%awa)) then
             num_conc & ! any radius will have the same num_conc
                  = aero_weight_array_num_conc_at_radius(aero_state%awa, &
                  i_class, 1d0)
             new_particle_volume = total_volume_conc / num_conc &
                  / real(integer_varray_n_entry( &
                  aero_state%aero_sorted%size_class%inverse(i_bin, i_class)), &
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
                i_part = aero_state%aero_sorted%size_class%inverse(i_bin, &
                     i_class)%entry(i_entry)
                particle_volume = aero_particle_volume( &
                     aero_state%apa%particle(i_part))
                if (i_part == 1) then
                   lower_volume = particle_volume
                   upper_volume = particle_volume
                end if
                lower_volume = min(lower_volume, particle_volume)
                upper_volume = max(upper_volume, particle_volume)
             end do

             lower_function = real(n_part, kind=dp) &
                  * aero_weight_array_num_conc_at_radius(aero_state%awa, &
                  i_class, aero_data_vol2rad(aero_data, lower_volume)) &
                  - total_num_conc
             upper_function = real(n_part, kind=dp) &
                  * aero_weight_array_num_conc_at_radius(aero_state%awa, &
                  i_class, aero_data_vol2rad(aero_data, upper_volume)) &
                  - total_num_conc

             ! do 50 rounds of bisection (2^50 = 10^15)
             do i_bisect = 1,50
                center_volume = (lower_volume + upper_volume) / 2d0
                center_function = real(n_part, kind=dp) &
                     * aero_weight_array_num_conc_at_radius(aero_state%awa, &
                     i_class, aero_data_vol2rad(aero_data, center_volume)) &
                     - total_num_conc
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
                i_part = aero_state%aero_sorted%size_class%inverse(i_bin, &
                     i_class)%entry(i_entry)
                particle_volume = aero_particle_volume( &
                     aero_state%apa%particle(i_part))
                if (i_part == 1) then
                   lower_volume = particle_volume
                   upper_volume = particle_volume
                end if
                lower_volume = min(lower_volume, particle_volume)
                upper_volume = max(upper_volume, particle_volume)
             end do

             lower_function = real(n_part, kind=dp) &
                  * aero_weight_array_num_conc_at_radius( &
                  aero_state%awa, i_class, aero_data_vol2rad(aero_data, &
                  lower_volume)) * lower_volume - total_volume_conc
             upper_function = real(n_part, kind=dp) &
                  * aero_weight_array_num_conc_at_radius( &
                  aero_state%awa, i_class, aero_data_vol2rad(aero_data, &
                  upper_volume)) * upper_volume - total_volume_conc

             ! do 50 rounds of bisection (2^50 = 10^15)
             do i_bisect = 1,50
                center_volume = (lower_volume + upper_volume) / 2d0
                center_function = real(n_part, kind=dp) &
                     * aero_weight_array_num_conc_at_radius( &
                     aero_state%awa, i_class, aero_data_vol2rad(aero_data, &
                     center_volume)) * center_volume - total_volume_conc
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
             i_part = aero_state%aero_sorted%size_class%inverse(i_bin, &
                  i_class)%entry(i_entry)
             particle_volume = aero_particle_volume( &
                  aero_state%apa%particle(i_part))
             aero_state%apa%particle(i_part)%vol &
                  = aero_state%apa%particle(i_part)%vol &
                  / particle_volume * new_particle_volume
          end do
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
    real(kind=dp) :: reweight_num_conc(aero_state%apa%n_part)

    ! We're modifying particle diameters, so bin sorting is now invalid
    aero_state%valid_sort = .false.

    call aero_state_num_conc_for_reweight(aero_state, aero_data, &
         reweight_num_conc)
    if (aero_data%i_water > 0) then
       do i_part = 1,aero_state_n_part(aero_state)
          aero_state%apa%particle(i_part)%vol(aero_data%i_water) = 0d0
       end do
       aero_state%valid_sort = .false.
    end if
    ! adjust particles to account for weight changes
    call aero_state_reweight(aero_state, aero_data, reweight_num_conc)

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
    total_size = total_size + pmc_mpi_pack_size_real_array_2d(val%n_part_ideal)
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
    call pmc_mpi_pack_real_array_2d(buffer, position, val%n_part_ideal)
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
    call pmc_mpi_unpack_real_array_2d(buffer, position, val%n_part_ideal)
    call pmc_mpi_unpack_aero_info_array(buffer, position, val%aero_info_array)
    call assert(132104747, &
         position - prev_position <= pmc_mpi_pack_size_aero_state(val))
#endif

  end subroutine pmc_mpi_unpack_aero_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Gathers data from all processes into one aero_state on process 0.
  subroutine aero_state_mpi_gather(aero_state, aero_state_total, aero_data)

    !> Local aero_state.
    type(aero_state_t), intent(in) :: aero_state
    !> Centralized aero_state (only on process 0).
    type(aero_state_t), intent(inout) :: aero_state_total
    !> Aero data values.
    type(aero_data_t), intent(in) :: aero_data

#ifdef PMC_USE_MPI
    type(aero_state_t) :: aero_state_transfer
    integer :: n_proc, ierr, status(MPI_STATUS_SIZE)
    integer :: buffer_size, max_buffer_size, i_proc, position
    character, allocatable :: buffer(:)
#endif

    if (pmc_mpi_rank() == 0) then
       aero_state_total = aero_state
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
          call aero_state_zero(aero_state_transfer)
          position = 0
          call pmc_mpi_unpack_aero_state(buffer, position, &
               aero_state_transfer)
          call assert(518174881, position == buffer_size)
          deallocate(buffer)

          call aero_state_add(aero_state_total, aero_state_transfer, &
               aero_data)
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
    integer :: aero_particle_centers(aero_state_n_part(aero_state))

    ! try to get the dimension ID
    status = nf90_inq_dimid(ncid, "aero_particle", dimid_aero_particle)
    if (status == NF90_NOERR) return
    if (status /= NF90_EBADDIM) call pmc_nc_check(status)

    ! dimension not defined, so define now define it
    call pmc_nc_check(nf90_redef(ncid))

    call pmc_nc_check(nf90_def_dim(ncid, "aero_particle", &
         aero_state_n_part(aero_state), dimid_aero_particle))

    call pmc_nc_check(nf90_enddef(ncid))

    do i_part = 1,aero_state_n_part(aero_state)
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
    integer :: aero_removed_centers( &
         max(1, aero_info_array_n_item(aero_state%aero_info_array)))

    ! try to get the dimension ID
    status = nf90_inq_dimid(ncid, "aero_removed", dimid_aero_removed)
    if (status == NF90_NOERR) return
    if (status /= NF90_EBADDIM) call pmc_nc_check(status)

    ! dimension not defined, so define now define it
    call pmc_nc_check(nf90_redef(ncid))

    dim_size = max(1, aero_info_array_n_item(aero_state%aero_info_array))
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
    real(kind=dp) :: aero_particle_mass(aero_state_n_part(aero_state), &
         aero_data_n_spec(aero_data))
    integer :: aero_n_orig_part(aero_state_n_part(aero_state), &
         aero_data_n_source(aero_data))
    integer :: aero_particle_weight_group(aero_state_n_part(aero_state))
    integer :: aero_particle_weight_class(aero_state_n_part(aero_state))
    real(kind=dp) :: aero_absorb_cross_sect(aero_state_n_part(aero_state))
    real(kind=dp) :: aero_scatter_cross_sect(aero_state_n_part(aero_state))
    real(kind=dp) :: aero_asymmetry(aero_state_n_part(aero_state))
    real(kind=dp) :: aero_refract_shell_real(aero_state_n_part(aero_state))
    real(kind=dp) :: aero_refract_shell_imag(aero_state_n_part(aero_state))
    real(kind=dp) :: aero_refract_core_real(aero_state_n_part(aero_state))
    real(kind=dp) :: aero_refract_core_imag(aero_state_n_part(aero_state))
    real(kind=dp) :: aero_core_vol(aero_state_n_part(aero_state))
    integer :: aero_water_hyst_leg(aero_state_n_part(aero_state))
    real(kind=dp) :: aero_num_conc(aero_state_n_part(aero_state))
    integer :: aero_id(aero_state_n_part(aero_state))
    real(kind=dp) :: aero_least_create_time(aero_state_n_part(aero_state))
    real(kind=dp) :: aero_greatest_create_time(aero_state_n_part(aero_state))
    integer :: aero_removed_id( &
         max(1, aero_info_array_n_item(aero_state%aero_info_array)))
    integer :: aero_removed_action( &
         max(1, aero_info_array_n_item(aero_state%aero_info_array)))
    integer :: aero_removed_other_id( &
         max(1, aero_info_array_n_item(aero_state%aero_info_array)))

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
    !! Each aerosol particle \c i represents a number concentration
    !! given by <tt>aero_num_conc(i)</tt>. Multiplying a per-particle
    !! quantity by the respective number concentration gives the
    !! concentration of that quantity contributed by the particle. For
    !! example, summing <tt>aero_particle_mass(i,s) *
    !! aero_num_conc(i)</tt> over all \c i gives the total mass
    !! concentration of species \c s in kg/m^3. Similarly, summing
    !! <tt>aero_absorb_cross_sect(i) * aero_num_conc(i)</tt> over all
    !! \c i will give the concentration of scattering cross section in
    !! m^2/m^3.
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
    !!   - \b aero_particle_weight_group (dim <tt>aero_particle</tt>):
    !!     weight group number (1 to <tt>aero_weight_group</tt>) of
    !!     each aerosol particle
    !!   - \b aero_particle_weight_class (dim <tt>aero_particle</tt>):
    !!     weight class number (1 to <tt>aero_weight_class</tt>) of each
    !!     aerosol particle
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
    !!   - \b aero_num_conc (unit m^{-3}, dim \c aero_particle): number
    !!     concentration associated with each particle
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

    if (aero_state_n_part(aero_state) > 0) then
       call aero_state_netcdf_dim_aero_particle(aero_state, ncid, &
            dimid_aero_particle)

       do i_part = 1,aero_state_n_part(aero_state)
          aero_particle_mass(i_part, :) &
               = aero_state%apa%particle(i_part)%vol * aero_data%density
          aero_n_orig_part(i_part, :) &
               = aero_state%apa%particle(i_part)%n_orig_part
          aero_particle_weight_group(i_part) &
               = aero_state%apa%particle(i_part)%weight_group
          aero_particle_weight_class(i_part) &
               = aero_state%apa%particle(i_part)%weight_class
          aero_water_hyst_leg(i_part) &
               = aero_state%apa%particle(i_part)%water_hyst_leg
          aero_num_conc(i_part) &
               = aero_state_particle_num_conc(aero_state, &
               aero_state%apa%particle(i_part), aero_data)
          aero_id(i_part) = aero_state%apa%particle(i_part)%id
          aero_least_create_time(i_part) &
               = aero_state%apa%particle(i_part)%least_create_time
          aero_greatest_create_time(i_part) &
               = aero_state%apa%particle(i_part)%greatest_create_time
          if (record_optical) then
             aero_absorb_cross_sect(i_part) &
                  = aero_state%apa%particle(i_part)%absorb_cross_sect
             aero_scatter_cross_sect(i_part) &
                  = aero_state%apa%particle(i_part)%scatter_cross_sect
             aero_asymmetry(i_part) = aero_state%apa%particle(i_part)%asymmetry
             aero_refract_shell_real(i_part) &
                  = real(aero_state%apa%particle(i_part)%refract_shell)
             aero_refract_shell_imag(i_part) &
                  = aimag(aero_state%apa%particle(i_part)%refract_shell)
             aero_refract_core_real(i_part) &
                  = real(aero_state%apa%particle(i_part)%refract_core)
             aero_refract_core_imag(i_part) &
                  = aimag(aero_state%apa%particle(i_part)%refract_core)
             aero_core_vol(i_part) = aero_state%apa%particle(i_part)%core_vol
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
       call pmc_nc_write_integer_1d(ncid, aero_particle_weight_group, &
            "aero_particle_weight_group", (/ dimid_aero_particle /), &
            long_name="weight group number of each aerosol particle")
       call pmc_nc_write_integer_1d(ncid, aero_particle_weight_class, &
            "aero_particle_weight_class", (/ dimid_aero_particle /), &
            long_name="weight class number of each aerosol particle")
       call pmc_nc_write_integer_1d(ncid, aero_water_hyst_leg, &
            "aero_water_hyst_leg", (/ dimid_aero_particle /), &
            long_name="leg of the water hysteresis curve leg of each "&
            // "aerosol particle")
       call pmc_nc_write_real_1d(ncid, aero_num_conc, &
            "aero_num_conc", (/ dimid_aero_particle /), unit="m^{-3}", &
            long_name="number concentration for each particle")
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

    ! FIXME: move this to aero_info_array.F90, together with
    ! aero_state_netcdf_dim_aero_removed() ?
    if (record_removals) then
       call aero_state_netcdf_dim_aero_removed(aero_state, ncid, &
            dimid_aero_removed)
       if (aero_info_array_n_item(aero_state%aero_info_array) >= 1) then
          do i_remove = 1,aero_info_array_n_item(aero_state%aero_info_array)
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

  ! this belongs in the subroutine above, but is outside because
  ! Doxygen 1.8.7 doesn't resolve references when multiple \page
  ! blocks are in one subroutine

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
    integer, allocatable :: aero_particle_weight_group(:)
    integer, allocatable :: aero_particle_weight_class(:)
    real(kind=dp), allocatable :: aero_absorb_cross_sect(:)
    real(kind=dp), allocatable :: aero_scatter_cross_sect(:)
    real(kind=dp), allocatable :: aero_asymmetry(:)
    real(kind=dp), allocatable :: aero_refract_shell_real(:)
    real(kind=dp), allocatable :: aero_refract_shell_imag(:)
    real(kind=dp), allocatable :: aero_refract_core_real(:)
    real(kind=dp), allocatable :: aero_refract_core_imag(:)
    real(kind=dp), allocatable :: aero_core_vol(:)
    integer, allocatable :: aero_water_hyst_leg(:)
    real(kind=dp), allocatable :: aero_num_conc(:)
    integer, allocatable :: aero_id(:)
    real(kind=dp), allocatable :: aero_least_create_time(:)
    real(kind=dp), allocatable :: aero_greatest_create_time(:)
    integer, allocatable :: aero_removed_id(:)
    integer, allocatable :: aero_removed_action(:)
    integer, allocatable :: aero_removed_other_id(:)

    call aero_state_zero(aero_state)

    status = nf90_inq_dimid(ncid, "aero_particle", dimid_aero_particle)
    if (status == NF90_EBADDIM) then
       ! no aero_particle dimension means no particles present
       call aero_weight_array_input_netcdf(aero_state%awa, ncid)
       return
    end if
    call pmc_nc_check(status)
    call pmc_nc_check(nf90_Inquire_Dimension(ncid, dimid_aero_particle, &
         name, n_part))

    call pmc_nc_read_real_2d(ncid, aero_particle_mass, &
         "aero_particle_mass")
    call pmc_nc_read_integer_2d(ncid, aero_n_orig_part, &
         "aero_n_orig_part")
    call pmc_nc_read_integer_1d(ncid, aero_particle_weight_group, &
         "aero_particle_weight_group")
    call pmc_nc_read_integer_1d(ncid, aero_particle_weight_class, &
         "aero_particle_weight_class")
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
    call pmc_nc_read_real_1d(ncid, aero_num_conc, &
         "aero_num_conc")
    call pmc_nc_read_integer_1d(ncid, aero_id, &
         "aero_id")
    call pmc_nc_read_real_1d(ncid, aero_least_create_time, &
         "aero_least_create_time")
    call pmc_nc_read_real_1d(ncid, aero_greatest_create_time, &
         "aero_greatest_create_time")

    call aero_weight_array_input_netcdf(aero_state%awa, ncid)
    call aero_state_set_n_part_ideal(aero_state, 0d0)

    do i_part = 1,n_part
       aero_particle%vol = aero_particle_mass(i_part, :) / aero_data%density
       aero_particle%n_orig_part = aero_n_orig_part(i_part, :)
       aero_particle%weight_group = aero_particle_weight_group(i_part)
       aero_particle%weight_class = aero_particle_weight_class(i_part)
       if (size(aero_absorb_cross_sect) == n_part) then
          aero_particle%absorb_cross_sect = aero_absorb_cross_sect(i_part)
       end if
       if (size(aero_scatter_cross_sect) == n_part) then
            aero_particle%scatter_cross_sect = aero_scatter_cross_sect(i_part)
         end if
       if (size(aero_asymmetry) == n_part) then
            aero_particle%asymmetry = aero_asymmetry(i_part)
         end if
       if ((size(aero_refract_shell_real) == n_part) &
            .and. (size(aero_refract_shell_imag) == n_part)) then
          aero_particle%refract_shell = &
               cmplx(aero_refract_shell_real(i_part), &
               aero_refract_shell_imag(i_part), kind=dc)
       end if
       if ((size(aero_refract_core_real) == n_part) &
            .and. (size(aero_refract_core_imag) == n_part)) then
          aero_particle%refract_core = cmplx(aero_refract_core_real(i_part), &
               aero_refract_core_imag(i_part), kind=dc)
       end if
       if (size(aero_core_vol) == n_part) then
          aero_particle%core_vol = aero_core_vol(i_part)
       end if
       aero_particle%water_hyst_leg = aero_water_hyst_leg(i_part)
       aero_particle%id = aero_id(i_part)
       aero_particle%least_create_time = aero_least_create_time(i_part)
       aero_particle%greatest_create_time = aero_greatest_create_time(i_part)

       call assert(314368871, almost_equal(aero_num_conc(i_part), &
            aero_weight_array_num_conc(aero_state%awa, aero_particle, &
            aero_data)))

       call aero_state_add_particle(aero_state, aero_particle, aero_data)
    end do

    next_id = maxval(aero_id) + 1

    call pmc_nc_read_integer_1d(ncid, aero_removed_id, &
         "aero_removed_id", must_be_present=.false.)
    call pmc_nc_read_integer_1d(ncid, aero_removed_action, &
         "aero_removed_action", must_be_present=.false.)
    call pmc_nc_read_integer_1d(ncid, aero_removed_other_id, &
         "aero_removed_other_id", must_be_present=.false.)

    n_info_item = size(aero_removed_id)
    if (n_info_item >= 1) then
       if ((n_info_item > 1) &
            .or. ((n_info_item == 1) .and. (aero_removed_id(1) /= 0))) then
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
    end if

  end subroutine aero_state_input_netcdf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Sorts the particles if necessary.
  subroutine aero_state_sort(aero_state, aero_data, bin_grid, all_procs_same)

    !> Aerosol state to sort.
    type(aero_state_t), intent(inout) :: aero_state
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Bin grid.
    type(bin_grid_t), optional, intent(in) :: bin_grid
    !> Whether all processors should use the same bin grid.
    logical, optional, intent(in) :: all_procs_same

    call aero_sorted_remake_if_needed(aero_state%aero_sorted, aero_state%apa, &
         aero_data, aero_state%valid_sort, size(aero_state%awa%weight, 1), &
         size(aero_state%awa%weight, 2), bin_grid, all_procs_same)
    aero_state%valid_sort = .true.

  end subroutine aero_state_sort

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Check that aerosol state data is consistent.
  subroutine aero_state_check(aero_state, aero_data)

    !> Aerosol state to check.
    type(aero_state_t), intent(in) :: aero_state
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data

    logical, parameter :: continue_on_error = .false.

    call aero_particle_array_check(aero_state%apa, aero_data, &
         continue_on_error)
    if (aero_state%valid_sort) then
       call aero_sorted_check(aero_state%aero_sorted, aero_state%apa, &
            aero_data, size(aero_state%awa%weight, 1), &
            size(aero_state%awa%weight, 2), continue_on_error)
    end if

  end subroutine aero_state_check

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef PMC_USE_CAMP
  !> Initialize the aero_state_t variable with camp chem data
  subroutine aero_state_initialize(aero_state, aero_data, camp_core)

    !> Aerosol state.
    type(aero_state_t), intent(inout) :: aero_state
    !> Aerosol data.
    class(aero_data_t), intent(in) :: aero_data
    !> CAMP core.
    type(camp_core_t), intent(in) :: camp_core

    select type( aero_rep => aero_data%aero_rep_ptr)
       type is(aero_rep_single_particle_t)
          ! Set up the update data objects for number
          call camp_core%initialize_update_object(aero_rep, &
                                                 aero_state%update_number)
       class default
          call die_msg(927605681, "Wrong aerosol representation type")
    end select

  end subroutine aero_state_initialize
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_aero_state
