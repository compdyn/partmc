! Copyright (C) 2005-2010 Nicole Riemer and Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_aero_state module.

!> The aero_state_t structure and assocated subroutines.
module pmc_aero_state

  use pmc_aero_particle_array
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
#ifdef PMC_USE_MPI
  use mpi
#endif

  !> MPI tag for mixing particles between processors.
  integer, parameter :: AERO_STATE_TAG_MIX = 4989

  !> The current collection of aerosol particles.
  !!
  !! The particles in aero_state_t are stored sorted per-bin, to
  !! improve efficiency of access and sampling. If a particle has
  !! total radius \c r then calling <tt> i_bin =
  !! bin_grid_particle_in_bin(bin_grid, r)</tt> finds the bin number
  !! i_bin where that particle should go. That particle is then stored
  !! as \c aero_state%%bin(i_bin)%%particle(i_part), where \c i_part
  !! is the index within the bin. \c
  !! aero_state%%v(i_bin)%%p(i_part)%%vol(i_spec) is thus the volume
  !! of the \c i_spec-th species in the \c i_part-th particle in the
  !! \c i_bin-th bin.
  !!
  !! Typically most of the bins have only a few particles, while a
  !! small number of bins have many particles. To avoid having too
  !! much storage allocated for the bins with only a few particles, we
  !! do dynamic allocation and deallocation of the storage
  !! per-bin. With Fortran 90 we can't have arrays of arrays, so we
  !! have to use an array of pointers, and then allocate each pointer.
  !!
  !! To avoid doing allocation and deallocation every time we add or
  !! remove a particle to a bin, we always double or halve the bin
  !! storage as necessary. The actual number of particles stored in a
  !! bin will generally be less than the actual memory allocated for
  !! that bin, so we store the current number of particles in a bin in
  !! \c aero_state%%bin(i_bin)%%n_part. The allocated size of bin
  !! storage in \c aero_state%%bin(i_bin)%%particle is not stored
  !! explicitly, but can be obtained with the Fortran 90 SIZE()
  !! intrinsic function.
  !!
  !! Every time we remove particles we keep track of the particle ID
  !! and the action performed in the aero_info_array_t structure. This
  !! is typically cleared each time we output data to disk.
  type aero_state_t
     !> Bin arrays.
     type(aero_particle_array_t), pointer :: bin(:)
     !> Computational volume (m^3).
     real(kind=dp) :: comp_vol
     !> Total number of particles.
     integer :: n_part
     !> Information on removed particles.
     type(aero_info_array_t) :: aero_info_array
  end type aero_state_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Allocates aerosol arrays.
  subroutine aero_state_allocate(aero_state)

    !> Aerosol to initialize.
    type(aero_state_t), intent(out) :: aero_state
    
    allocate(aero_state%bin(0))
    call aero_info_array_allocate(aero_state%aero_info_array)

  end subroutine aero_state_allocate
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Allocates aerosol arrays with the given sizes.
  subroutine aero_state_allocate_size(aero_state, n_bin, n_spec)

    !> Aerosol to initialize.
    type(aero_state_t), intent(out) :: aero_state
    !> Number of bins.
    integer, intent(in) :: n_bin
    !> Number of species.
    integer, intent(in) :: n_spec
    
    integer i

    allocate(aero_state%bin(n_bin))
    do i = 1,n_bin
       call aero_particle_array_allocate_size(aero_state%bin(i), 0, n_spec)
    end do
    aero_state%comp_vol = 0d0
    aero_state%n_part = 0
    call aero_info_array_allocate(aero_state%aero_info_array)

  end subroutine aero_state_allocate_size
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Deallocates a previously allocated aerosol.
  subroutine aero_state_deallocate(aero_state)

    !> Aerosol to initialize.
    type(aero_state_t), intent(inout) :: aero_state
    
    integer :: n_bin, i

    n_bin = size(aero_state%bin)
    do i = 1,n_bin
       call aero_particle_array_deallocate(aero_state%bin(i))
    end do
    deallocate(aero_state%bin)
    call aero_info_array_deallocate(aero_state%aero_info_array)

  end subroutine aero_state_deallocate
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Copies aerosol to a destination that has already had
  !> aero_state_allocate() called on it.
  subroutine aero_state_copy(aero_state_from, aero_state_to)

    !> Reference aerosol.
    type(aero_state_t), intent(in) :: aero_state_from
    !> Already allocated.
    type(aero_state_t), intent(inout) :: aero_state_to
    
    integer :: n_bin, i

    n_bin = size(aero_state_from%bin)

    call aero_state_deallocate(aero_state_to)
    call aero_state_allocate_size(aero_state_to, n_bin, 0)

    do i = 1,n_bin
       call aero_particle_array_copy(aero_state_from%bin(i), &
            aero_state_to%bin(i))
    end do

    aero_state_to%comp_vol = aero_state_from%comp_vol
    aero_state_to%n_part = aero_state_from%n_part

    call aero_info_array_copy(aero_state_from%aero_info_array, &
         aero_state_to%aero_info_array)

  end subroutine aero_state_copy
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the total number of particles in an aerosol distribution.
  integer function aero_state_total_particles(aero_state)

    !> Aerosol state.
    type(aero_state_t), intent(in) :: aero_state

    aero_state_total_particles = aero_state%n_part

  end function aero_state_total_particles

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Resets an aero_state to have zero particles per bin. This must
  !> already have had aero_state_allocate() called on it. This
  !> function can be called more than once on the same state.
  subroutine aero_state_zero(aero_state)

    !> State to zero.
    type(aero_state_t), intent(inout) :: aero_state
    
    integer :: i, n_bin

    n_bin = size(aero_state%bin)
    do i = 1,n_bin
       call aero_particle_array_zero(aero_state%bin(i))
    end do
    aero_state%n_part = 0
    call aero_info_array_zero(aero_state%aero_info_array)

  end subroutine aero_state_zero
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Add the given particle.
  subroutine aero_state_add_particle(aero_state, i_bin, aero_particle)

    !> Aerosol state.
    type(aero_state_t), intent(inout) :: aero_state
    !> Bin number of particle to add.
    integer, intent(in) :: i_bin
    !> Particle to add.
    type(aero_particle_t), intent(in) :: aero_particle

    call aero_particle_array_add_particle(aero_state%bin(i_bin), &
         aero_particle)
    aero_state%n_part = aero_state%n_part + 1

  end subroutine aero_state_add_particle

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Remove the given particle without recording it.
  subroutine aero_state_remove_particle_no_info(aero_state, i_bin, index)

    !> Aerosol state.
    type(aero_state_t), intent(inout) :: aero_state
    !> Bin number of particle to remove.
    integer, intent(in) :: i_bin
    !> Index in bin of particle to remove.
    integer, intent(in) :: index

    call aero_particle_array_remove_particle(aero_state%bin(i_bin), index)
    aero_state%n_part = aero_state%n_part - 1

  end subroutine aero_state_remove_particle_no_info

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Remove the given particle and record the removal.
  subroutine aero_state_remove_particle_with_info(aero_state, i_bin, &
       index, aero_info)

    !> Aerosol state.
    type(aero_state_t), intent(inout) :: aero_state
    !> Bin number of particle to remove.
    integer, intent(in) :: i_bin
    !> Index in bin of particle to remove.
    integer, intent(in) :: index
    !> Removal info.
    type(aero_info_t), intent(in) :: aero_info

    call aero_state_remove_particle_no_info(aero_state, i_bin, index)
    call aero_info_array_add_aero_info(aero_state%aero_info_array, &
         aero_info)

  end subroutine aero_state_remove_particle_with_info

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Remove the given particle and possibly record the removal.
  subroutine aero_state_remove_particle(aero_state, i_bin, index, &
       record_removal, aero_info)

    !> Aerosol state.
    type(aero_state_t), intent(inout) :: aero_state
    !> Bin number of particle to remove.
    integer, intent(in) :: i_bin
    !> Index in bin of particle to remove.
    integer, intent(in) :: index
    !> Whether to record the removal in the aero_info_array.
    logical, intent(in) :: record_removal
    !> Removal info.
    type(aero_info_t), intent(in) :: aero_info

    if (record_removal) then
       call aero_state_remove_particle_with_info(aero_state, i_bin, &
            index, aero_info)
    else
       call aero_state_remove_particle_no_info(aero_state, i_bin, index)
    end if

  end subroutine aero_state_remove_particle

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine aero_state_remove_rand_particle_from_bin(aero_state, &
       i_bin, aero_particle)

    !> Aerosol state.
    type(aero_state_t), intent(inout) :: aero_state
    !> Bin number to remove particle from.
    integer, intent(in) :: i_bin
    !> Removed particle.
    type(aero_particle_t), intent(inout) :: aero_particle

    integer :: i_part

    call assert(392182617, aero_state%bin(i_bin)%n_part > 0)
    i_part = pmc_rand_int(aero_state%bin(i_bin)%n_part)
    call aero_particle_copy(aero_state%bin(i_bin)%particle(i_part), &
         aero_particle)
    call aero_state_remove_particle_no_info(aero_state, i_bin, i_part)

  end subroutine aero_state_remove_rand_particle_from_bin

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
    aero_state%comp_vol = aero_state%comp_vol + aero_state_delta%comp_vol

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

    integer :: i_bin, i_part

    call assert(265083067, &
         size(aero_state%bin) == size(aero_state_delta%bin))
    do i_bin = 1,size(aero_state_delta%bin)
       do i_part = 1,aero_state_delta%bin(i_bin)%n_part
          call aero_state_add_particle(aero_state, i_bin, &
               aero_state_delta%bin(i_bin)%particle(i_part))
       end do
    end do
    call aero_info_array_add(aero_state%aero_info_array, &
         aero_state_delta%aero_info_array)
    
  end subroutine aero_state_add_particles

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Generates a Poisson sample of an aero_dist, adding to
  !> aero_state. The sampled amount is sample_prop *
  !> aero_state%comp_vol.
  subroutine aero_state_add_aero_dist_sample(aero_state, bin_grid, &
       aero_data, aero_weight, aero_dist, sample_prop, create_time)

    !> Aero state to add to.
    type(aero_state_t), intent(inout) :: aero_state
    !> Bin grid.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Aero data values.
    type(aero_data_t), intent(in) :: aero_data
    !> Aerosol weight.
    type(aero_weight_t), intent(in) :: aero_weight
    !> Distribution to sample.
    type(aero_dist_t), intent(in) :: aero_dist
    !> Volume fraction to sample (1).
    real(kind=dp), intent(in) :: sample_prop
    !> Creation time for new particles (s).
    real(kind=dp), intent(in) :: create_time

    real(kind=dp) :: n_samp_avg, sample_vol, radius, vol
    integer :: n_samp, i_mode, i_samp, i_bin
    integer :: num_per_bin(bin_grid%n_bin)
    type(aero_mode_t), pointer :: aero_mode
    type(aero_particle_t) :: aero_particle

    call aero_particle_allocate_size(aero_particle, aero_data%n_spec)
    sample_vol = sample_prop * aero_state%comp_vol
    do i_mode = 1,aero_dist%n_mode
       aero_mode => aero_dist%mode(i_mode)
       n_samp_avg = sample_vol &
            * aero_mode_weighted_num_conc(aero_mode, aero_weight)
       n_samp = rand_poisson(n_samp_avg)
       do i_samp = 1,n_samp
          call aero_mode_sample_radius(aero_mode, aero_weight, radius)
          vol = rad2vol(radius)
          call aero_particle_set_vols(aero_particle, aero_mode%vol_frac * vol)
          call aero_particle_new_id(aero_particle)
          call aero_particle_set_create_time(aero_particle, create_time)
          i_bin = aero_particle_in_bin(aero_particle, bin_grid)
          call aero_state_add_particle(aero_state, i_bin, aero_particle)
       end do
    end do
    call aero_particle_deallocate(aero_particle)

  end subroutine aero_state_add_aero_dist_sample
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Choose a random particle from the aero_state.
  subroutine aero_state_rand_particle(aero_state, i_bin, i_part)

    !> Original state.
    type(aero_state_t), intent(in) :: aero_state
    !> Bin number of particle.
    integer, intent(out) :: i_bin
    !> Particle number within bin.
    integer, intent(out) :: i_part

    integer :: n_bin, disc_pdf(size(aero_state%bin))

    call assert(950725003, aero_state%n_part > 0)
    n_bin = size(aero_state%bin)
    disc_pdf = (/(aero_state%bin(i_bin)%n_part, i_bin = 1,n_bin)/)
    i_bin = sample_disc_pdf(disc_pdf)
    i_part = pmc_rand_int(aero_state%bin(i_bin)%n_part)

  end subroutine aero_state_rand_particle

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Generates a Poisson sample by removing particles from
  !> aero_state_from and adding them to aero_state_to, which must
  !> be already allocated (and should have its comp_vol set).
  subroutine aero_state_sample(aero_state_from, aero_state_to, &
       sample_prop, removal_action)

    !> Original state.
    type(aero_state_t), intent(inout) :: aero_state_from
    !> Destination state.
    type(aero_state_t), intent(inout) :: aero_state_to
    !> Proportion to sample.
    real(kind=dp), intent(in) :: sample_prop
    !> Action for removal (see pmc_aero_info module for action
    !> parameters). Set to AERO_INFO_NONE to not log removal.
    integer, intent(in) :: removal_action
    
    integer :: n_transfer, i_transfer, n_bin, i_bin, i_part
    logical :: do_add, do_remove
    real(kind=dp) :: vol_ratio
    type(aero_info_t) :: aero_info

    call assert(721006962, (sample_prop >= 0d0) .and. (sample_prop <= 1d0))
    n_transfer = rand_poisson(sample_prop &
         * real(aero_state_total_particles(aero_state_from), kind=dp))
    n_bin = size(aero_state_from%bin)
    vol_ratio = aero_state_to%comp_vol / aero_state_from%comp_vol
    i_transfer = 0
    do while (i_transfer < n_transfer)
       if (aero_state_total_particles(aero_state_from) <= 0) exit
       call aero_state_rand_particle(aero_state_from, i_bin, i_part)
       if (aero_state_to%comp_vol == aero_state_from%comp_vol) then
          ! to_comp_vol == from_comp_vol so just move the particle
          do_add = .true.
          do_remove = .true.
       elseif (aero_state_to%comp_vol > aero_state_from%comp_vol) then
          ! to_comp_vol is bigger than from_comp_vol, so only maybe
          ! remove the particle but always add it (vol_ratio > 1)
          do_add = .true.
          do_remove = .false.
          if (pmc_random() < 1d0 / vol_ratio) then
             do_remove = .true.
          end if
       else ! aero_state_to%comp_vol < aero_state_from%comp_vol
          ! to_comp_vol is smaller than from_comp_vol, so always
          ! remove the particle but only maybe add it (vol_ratio < 1)
          do_add = .false.
          if (pmc_random() < vol_ratio) then
             do_add = .true.
          end if
          do_remove = .true.
       end if
       if (do_add) then
          call aero_state_add_particle(aero_state_to, i_bin, &
               aero_state_from%bin(i_bin)%particle(i_part))
       end if
       if (do_remove) then
          if (removal_action /= AERO_INFO_NONE) then
             call aero_info_allocate(aero_info)
             aero_info%id = &
                  aero_state_from%bin(i_bin)%particle(i_part)%id
             aero_info%action = removal_action
             call aero_state_remove_particle(aero_state_from, i_bin, &
                  i_part, .true., aero_info)
             call aero_info_deallocate(aero_info)
          else
             call aero_state_remove_particle(aero_state_from, i_bin, &
                  i_part, .false., aero_info)
          end if
          i_transfer = i_transfer + 1
       end if
    end do
    
  end subroutine aero_state_sample
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Generates a rough bin-wise sample by removing particles from
  !> aero_state_from and adding them to aero_state_to, which must
  !> be already allocated (and should have its comp_vol set).
  subroutine aero_state_sample_rough(aero_state_from, aero_state_to, &
       sample_prop, removal_action)

    !> Original state.
    type(aero_state_t), intent(inout) :: aero_state_from
    !> Destination state.
    type(aero_state_t), intent(inout) :: aero_state_to
    !> Proportion to sample.
    real(kind=dp), intent(in) :: sample_prop
    !> Action for removal (see pmc_aero_info module for action
    !> parameters). Set to AERO_INFO_NONE to not log removal.
    integer, intent(in) :: removal_action
    
    integer :: n_transfer, i_transfer, n_bin, i_bin, i_part
    logical :: do_add, do_remove
    real(kind=dp) :: vol_ratio
    type(aero_info_t) :: aero_info

    n_bin = size(aero_state_from%bin)
    vol_ratio = aero_state_to%comp_vol / aero_state_from%comp_vol
    do i_bin = 1,n_bin
       n_transfer = prob_round(sample_prop &
            * real(aero_state_from%bin(i_bin)%n_part, kind=dp))
       i_transfer = 0
       do while (i_transfer < n_transfer)
          if (aero_state_from%bin(i_bin)%n_part <= 0) exit
          i_part = pmc_rand_int(aero_state_from%bin(i_bin)%n_part)
          if (aero_state_to%comp_vol == aero_state_from%comp_vol) then
             ! to_comp_vol == from_comp_vol so just move the particle
             do_add = .true.
             do_remove = .true.
          elseif (aero_state_to%comp_vol > aero_state_from%comp_vol) then
             ! to_comp_vol is bigger than from_comp_vol, so only maybe
             ! remove the particle but always add it (vol_ratio > 1)
             do_add = .true.
             do_remove = .false.
             if (pmc_random() < 1d0 / vol_ratio) then
                do_remove = .true.
             end if
          else ! aero_state_to%comp_vol < aero_state_from%comp_vol
             ! to_comp_vol is smaller than from_comp_vol, so always
             ! remove the particle but only maybe add it (vol_ratio < 1)
             do_add = .false.
             if (pmc_random() < vol_ratio) then
                do_add = .true.
             end if
             do_remove = .true.
          end if
          if (do_add) then
             call aero_state_add_particle(aero_state_to, i_bin, &
                  aero_state_from%bin(i_bin)%particle(i_part))
          end if
          if (do_remove) then
             if (removal_action /= AERO_INFO_NONE) then
                call aero_info_allocate(aero_info)
                aero_info%id = &
                     aero_state_from%bin(i_bin)%particle(i_part)%id
                aero_info%action = removal_action
                call aero_state_remove_particle(aero_state_from, i_bin, &
                     i_part, .true., aero_info)
                call aero_info_deallocate(aero_info)
             else
                call aero_state_remove_particle(aero_state_from, i_bin, &
                     i_part, .false., aero_info)
             end if
             i_transfer = i_transfer + 1
          end if
       end do
    end do
    
  end subroutine aero_state_sample_rough
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Create the bin number and mass arrays from aero_state%v.
  subroutine aero_state_to_binned(bin_grid, aero_data, aero_weight, &
       aero_state, aero_binned)
    
    !> Bin grid.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Aerosol weight.
    type(aero_weight_t), intent(in) :: aero_weight
    !> Aerosol state.
    type(aero_state_t), intent(in) :: aero_state
    !> Binned distributions.
    type(aero_binned_t), intent(inout) :: aero_binned
    
    integer :: b, j, s
    type(aero_particle_t), pointer :: aero_particle
    
    aero_binned%num_conc = 0d0
    aero_binned%vol_conc = 0d0
    do b = 1,bin_grid%n_bin
       do j = 1,aero_state%bin(b)%n_part
          aero_particle => aero_state%bin(b)%particle(j)
          aero_binned%vol_conc(b,:) = aero_binned%vol_conc(b,:) &
               + aero_particle%vol / aero_state%comp_vol &
               * aero_weight_value(aero_weight, &
               aero_particle_radius(aero_particle)) / bin_grid%log_width
          aero_binned%num_conc(b) = aero_binned%num_conc(b) &
               + 1d0 / aero_state%comp_vol &
               * aero_weight_value(aero_weight, &
               aero_particle_radius(aero_particle)) / bin_grid%log_width
       end do
    end do
    
  end subroutine aero_state_to_binned
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Does the same thing as aero_state_to_bin() but based on dry radius.
  subroutine aero_state_to_binned_dry(bin_grid, aero_data, aero_weight, &
       aero_state, aero_binned)
    
    !> Bin grid.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Aerosol weight.
    type(aero_weight_t), intent(in) :: aero_weight
    !> Aerosol state.
    type(aero_state_t), intent(in) :: aero_state
    !> Binned distributions.
    type(aero_binned_t), intent(inout) :: aero_binned
    
    integer :: b, j, s, b_dry
    type(aero_particle_t), pointer :: aero_particle
    
    aero_binned%num_conc = 0d0
    aero_binned%vol_conc = 0d0
    do b = 1,bin_grid%n_bin
       do j = 1,aero_state%bin(b)%n_part
          aero_particle => aero_state%bin(b)%particle(j)
          b_dry = bin_grid_particle_in_bin(bin_grid, &
               aero_particle_solute_radius(aero_particle, aero_data))
          aero_binned%vol_conc(b_dry,:) = aero_binned%vol_conc(b_dry,:) &
               + aero_particle%vol / aero_state%comp_vol &
               * aero_weight_value(aero_weight, &
               aero_particle_radius(aero_particle)) / bin_grid%log_width
          aero_binned%num_conc(b_dry) = aero_binned%num_conc(b_dry) &
               + 1d0 / aero_state%comp_vol &
               * aero_weight_value(aero_weight, &
               aero_particle_radius(aero_particle)) / bin_grid%log_width
       end do
    end do
    
  end subroutine aero_state_to_binned_dry
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Doubles number of particles.
  subroutine aero_state_double(aero_state)
    
    !> Aerosol state.
    type(aero_state_t), intent(inout) :: aero_state
    
    integer :: i, n_bin
    
    n_bin = size(aero_state%bin)
    do i = 1,n_bin
       call aero_particle_array_double(aero_state%bin(i))
    end do
    aero_state%comp_vol = 2d0 * aero_state%comp_vol
    aero_state%n_part = 2 * aero_state%n_part

  end subroutine aero_state_double
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Remove approximately half of the particles in each bin.
  subroutine aero_state_halve(aero_state, bin_grid)
    
    !> Aerosol state.
    type(aero_state_t), intent(inout) :: aero_state
    !> Bin grid.
    type(bin_grid_t), intent(in) :: bin_grid
    
    integer :: i_bin, i_part, n_part_orig, i_remove, n_remove
    type(aero_info_t) :: aero_info

    n_part_orig = aero_state%n_part
    do i_bin = 1,bin_grid%n_bin
       n_remove = prob_round(real(aero_state%bin(i_bin)%n_part, kind=dp) &
            / 2d0)
       do i_remove = 1,n_remove
          i_part = pmc_rand_int(aero_state%bin(i_bin)%n_part)
          call aero_info_allocate(aero_info)
          aero_info%id = &
               aero_state%bin(i_bin)%particle(i_part)%id
          aero_info%action = AERO_INFO_HALVED
          call aero_state_remove_particle(aero_state, i_bin, &
               i_part, .true., aero_info)
          call aero_info_deallocate(aero_info)
       end do
    end do
    aero_state%comp_vol = aero_state%comp_vol &
         * real(aero_state%n_part, kind=dp) / real(n_part_orig, kind=dp)

  end subroutine aero_state_halve
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Takes an aero_state_t where the particles might no longer be in
  !> the correct bins and resorts it so that every particle is in the
  !> correct bin.
  subroutine aero_state_resort(bin_grid, aero_state)
    
    !> Bin_grid.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Aerosol state.
    type(aero_state_t), intent(inout) :: aero_state
    
    integer :: i_bin, j, i_new_bin, k
    type(aero_info_t) :: aero_info ! dummy variable, never used
    
    ! The approach here is inefficient because we might reprocess
    ! particles. For example, if we are doing bin 1 and we shift a
    ! particle up to bin 2, when we do bin 2 we will reprocess it. It
    ! seems to be more trouble than it's worth to worry about this,
    ! however.
    
    do i_bin = 1,bin_grid%n_bin
       j = 1
       do while (j .le. aero_state%bin(i_bin)%n_part)
          ! find the new bin
          i_new_bin = aero_particle_in_bin( &
               aero_state%bin(i_bin)%particle(j), bin_grid)
          
          ! if the bin number has changed, move the particle
          if (i_bin .ne. i_new_bin) then
             call aero_state_add_particle(aero_state, i_new_bin, &
                  aero_state%bin(i_bin)%particle(j))
             call aero_state_remove_particle(aero_state, i_bin, j, &
                  .false., aero_info)

             ! in this case, don't advance j, so that we will still
             ! process the particle we just moved into the hole
          else
             ! if we didn't move the particle, advance j to process
             ! the next particle
             j = j + 1
          end if
       end do
    end do
    
    ! now shrink the bin storage if necessary
    do i_bin = 1,bin_grid%n_bin
       call aero_particle_array_shrink(aero_state%bin(i_bin))
    end do
    
  end subroutine aero_state_resort
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Mix the aero_states between all processes. Currently uses a
  !> simple all-to-all diffusion.
  subroutine aero_state_mix(aero_state, del_t, mix_timescale, &
       aero_data, bin_grid)

    !> Aerosol state.
    type(aero_state_t), intent(inout) :: aero_state
    !> Timestep (s).
    real(kind=dp), intent(in) :: del_t
    !> Mixing timescale (s).
    real(kind=dp), intent(in) :: mix_timescale
    !> Aero data values.
    type(aero_data_t), intent(in) :: aero_data
    !> Bin grid.
    type(bin_grid_t), intent(in) :: bin_grid

#ifdef PMC_USE_MPI
    integer :: rank, n_proc, dest_proc, ierr
    integer :: buffer_size, buffer_size_check
    character, allocatable :: buffer(:)
    type(aero_state_t), allocatable :: aero_state_sends(:)
    real(kind=dp) :: prob_transfer, prob_not_transferred
    real(kind=dp) :: prob_transfer_given_not_transferred
    real(kind=dp), allocatable :: comp_vols(:)

    ! processor information
    rank = pmc_mpi_rank()
    n_proc = pmc_mpi_size()
    if (n_proc == 1) then
       ! buffer allocation below fails if n_proc == 1
       ! so bail out early (nothing to mix anyway)
       return
    end if
    allocate(comp_vols(n_proc))
    call mpi_allgather(aero_state%comp_vol, 1, MPI_REAL8, &
         comp_vols, 1, MPI_REAL8, MPI_COMM_WORLD, ierr)

    ! extract particles to send
    allocate(aero_state_sends(n_proc))
    prob_not_transferred = 1d0
    do dest_proc = 0,(n_proc - 1)
       if (dest_proc /= rank) then
          call aero_state_allocate_size(aero_state_sends(dest_proc + 1), &
               bin_grid%n_bin, aero_data%n_spec)
          aero_state_sends(dest_proc + 1)%comp_vol = aero_state%comp_vol
          prob_transfer = (1d0 - exp(- del_t / mix_timescale)) &
               / real(n_proc, kind=dp) &
               * min(1d0, comp_vols(dest_proc + 1) / comp_vols(rank + 1))
          ! because we are doing sequential sampling from the aero_state
          ! we need to scale up the later transfer probabilities, because
          ! the later particles are being transferred conditioned on the
          ! fact that they did not transfer already
          prob_transfer_given_not_transferred = prob_transfer &
               / prob_not_transferred
          call aero_state_sample(aero_state, &
               aero_state_sends(dest_proc + 1), &
               prob_transfer_given_not_transferred, AERO_INFO_NONE)
          prob_not_transferred = prob_not_transferred - prob_transfer
       end if
    end do

    ! send the particles
    buffer_size = 0
    do dest_proc = 0,(n_proc - 1)
       if (dest_proc /= rank) then
          buffer_size = buffer_size &
               + pmc_mpi_pack_size_aero_state(aero_state_sends(dest_proc + 1))
          buffer_size = buffer_size + MPI_BSEND_OVERHEAD + 1024
       end if
    end do
    allocate(buffer(buffer_size))
    call mpi_buffer_attach(buffer, buffer_size, ierr)
    call pmc_mpi_check_ierr(ierr)
    do dest_proc = 0,(n_proc - 1)
       if (dest_proc /= rank) then
          call send_aero_state_mix(dest_proc, aero_state_sends(dest_proc + 1))
          call aero_state_deallocate(aero_state_sends(dest_proc + 1))
       end if
    end do
    deallocate(aero_state_sends)

    ! receive the particles
    do dest_proc = 0,(n_proc - 1)
       if (dest_proc /= rank) then
          call recv_aero_state_mix(aero_state)
       end if
    end do

    ! clean up
    call pmc_mpi_barrier() ! syncronize to make sure all buffers are done
    call mpi_buffer_detach(buffer, buffer_size_check, ierr)
    call pmc_mpi_check_ierr(ierr)
    call assert(995066222, buffer_size == buffer_size_check)
    deallocate(comp_vols)
#endif

  end subroutine aero_state_mix

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  !> Send the given \c aero_state to the destination processor for
  !> mixing.
  subroutine send_aero_state_mix(dest_proc, aero_state)

    !> Destination processor number.
    integer, intent(in) :: dest_proc
    !> Aero_state to send.
    type(aero_state_t), intent(in) :: aero_state

#ifdef PMC_USE_MPI
    character, allocatable :: buffer(:)
    integer :: buffer_size, position, ierr

    buffer_size = 0
    buffer_size = buffer_size + pmc_mpi_pack_size_aero_state(aero_state)
    allocate(buffer(buffer_size))
    position = 0
    call pmc_mpi_pack_aero_state(buffer, position, aero_state)
    call assert(311031745, buffer_size == position)
    call mpi_bsend(buffer, buffer_size, MPI_CHARACTER, dest_proc, &
         AERO_STATE_TAG_MIX, MPI_COMM_WORLD, ierr)
    call pmc_mpi_check_ierr(ierr)
    deallocate(buffer)
#endif

  end subroutine send_aero_state_mix

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Receive exactly one \c aero_state for mixing and add it on to the
  !> given \c aero_state.
  subroutine recv_aero_state_mix(aero_state)

    !> Base \c aero_state to add new particles to.
    type(aero_state_t), intent(inout) :: aero_state

#ifdef PMC_USE_MPI
    integer :: buffer_size, check_buffer_size, position, sent_proc, ierr
    integer :: status(MPI_STATUS_SIZE)
    character, allocatable :: buffer(:)
    type(aero_state_t) :: aero_state_mix

    ! get the message size
    call mpi_probe(MPI_ANY_SOURCE, AERO_STATE_TAG_MIX, MPI_COMM_WORLD, &
         status, ierr)
    call pmc_mpi_check_ierr(ierr)
    sent_proc = status(MPI_SOURCE)
    call mpi_get_count(status, MPI_CHARACTER, buffer_size, ierr)
    call pmc_mpi_check_ierr(ierr)
    allocate(buffer(buffer_size))

    ! get the message
    call mpi_recv(buffer, buffer_size, MPI_CHARACTER, &
         sent_proc, AERO_STATE_TAG_MIX, MPI_COMM_WORLD, status, ierr)
    call pmc_mpi_check_ierr(ierr)
    call mpi_get_count(status, MPI_CHARACTER, check_buffer_size, ierr)
    call pmc_mpi_check_ierr(ierr)
    call assert(967116476, buffer_size == check_buffer_size)
    call assert(377784810, sent_proc == status(MPI_SOURCE))

    ! unpack it
    position = 0
    call aero_state_allocate(aero_state_mix)
    call pmc_mpi_unpack_aero_state(buffer, position, aero_state_mix)
    call assert(908434410, position == buffer_size)
    deallocate(buffer)

    ! add the particles
    call aero_state_add_particles(aero_state, aero_state_mix)
    call aero_state_deallocate(aero_state_mix)
#endif

  end subroutine recv_aero_state_mix

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Check that all particles are in the correct bins and that the
  !> bin numbers and masses are correct. This is for debugging only.
  subroutine aero_state_check(bin_grid, aero_data, aero_state)
    
    !> Bin_grid.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Aerosol state.
    type(aero_state_t), intent(inout) :: aero_state
    
    real(kind=dp) :: check_bin_v, check_vol_conc, vol_tol
    real(kind=dp) :: num_tol, state_num_conc
    integer :: i, k, k_check, s, n_part_check, id, max_id
    logical :: error
    logical, allocatable :: id_present(:)
    
    error = .false.

    ! check that the total number of particles is correct
    i = 0 ! HACK to shut up gfortran warning
    n_part_check = sum((/(aero_state%bin(i)%n_part, i = 1,bin_grid%n_bin)/))
    if (aero_state%n_part /= n_part_check) then
       write(0,'(a20,a20)') 'aero_state%n_part', 'n_part_check'
       write(0,'(i20,i20)') aero_state%n_part, n_part_check
       error = .true.
    end if
    
    ! check that all particles are in the correct bins
    do k = 1,bin_grid%n_bin
       do i = 1,aero_state%bin(k)%n_part
          k_check = aero_particle_in_bin(aero_state%bin(k)%particle(i), &
               bin_grid)
          if (k .ne. k_check) then
             write(0,'(a10,a10,a10)') 'i', 'k', 'k_check'
             write(0,'(i10,i10,i10)') i, k, k_check
             error = .true.
          end if
       end do
    end do
    
    ! check we don't have duplicate IDs
    max_id = 0
    do k = 1,bin_grid%n_bin
       do i = 1,aero_state%bin(k)%n_part
          id = aero_state%bin(k)%particle(i)%id
          if (id > max_id) max_id = id
       end do
    end do
    allocate(id_present(max_id))
    do i = 1,max_id
       id_present(i) = .false.
    end do
    do k = 1,bin_grid%n_bin
       do i = 1,aero_state%bin(k)%n_part
          id = aero_state%bin(k)%particle(i)%id
          if (id_present(id)) then
             write(0,'(a15,a10,a10)') 'duplicate id', 'bin', 'index'
             write(0,'(i15,i10,i10)') id, k, i
             error = .true.
          end if
          id_present(id) = .true.
       end do
    end do
    deallocate(id_present)
    
    if (error) then
       call die_msg(371923719, 'aero_state_check() failed')
    end if
    
  end subroutine aero_state_check
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Set each aerosol particle to have its original total volume, but
  !> species volume ratios given by the total species volume ratio
  !> within each bin. This preserves the (weighted) total species
  !> volume per bin as well as per-particle total volumes.
  subroutine aero_state_bin_average_comp(aero_state, bin_grid, aero_data, &
       aero_weight, dry_volume)

    !> Aerosol state to average.
    type(aero_state_t), intent(inout) :: aero_state
    !> Bin grid to average within.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Aerosol weight.
    type(aero_weight_t), intent(in) :: aero_weight
    !> Whether to use dry volume (rather than wet).
    logical, intent(in) :: dry_volume

    real(kind=dp) :: weighted_species_volumes(aero_data%n_spec)
    real(kind=dp) :: weighted_total_volume, particle_volume, weight
    integer :: i_bin, i_part, i_spec
    type(aero_particle_t), pointer :: aero_particle

    do i_bin = 1,bin_grid%n_bin
       weighted_species_volumes = 0d0
       weighted_total_volume = 0d0
       do i_part = 1,aero_state%bin(i_bin)%n_part
          aero_particle => aero_state%bin(i_bin)%particle(i_part)
          weight = aero_weight_value(aero_weight, &
               aero_particle_radius(aero_particle))
          particle_volume = aero_particle_volume_maybe_dry(aero_particle, &
               aero_data, dry_volume)
          weighted_species_volumes = weighted_species_volumes &
               + weight * aero_particle%vol
          weighted_total_volume = weighted_total_volume &
               + weight * particle_volume
       end do
       do i_part = 1,aero_state%bin(i_bin)%n_part
          aero_particle => aero_state%bin(i_bin)%particle(i_part)
          particle_volume = aero_particle_volume_maybe_dry(aero_particle, &
               aero_data, dry_volume)
          aero_particle%vol = particle_volume * weighted_species_volumes &
               / weighted_total_volume
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
       aero_weight, dry_volume, bin_center, preserve_number)

    !> Aerosol state to average.
    type(aero_state_t), intent(inout) :: aero_state
    !> Bin grid to average within.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Aerosol weight.
    type(aero_weight_t), intent(in) :: aero_weight
    !> Whether to use dry volume (rather than wet).
    logical, intent(in) :: dry_volume
    !> Whether to assign the bin center volume (rather than the average
    !> volume).
    logical, intent(in) :: bin_center
    !> Whether to use the number-preserving scheme (otherwise will use
    !> the volume-preserving scheme). This parameter has no effect if
    !> \c bin_center is \c .true.
    logical, intent(in) :: preserve_number

    real(kind=dp) :: weighted_total_volume, particle_volume
    real(kind=dp) :: new_particle_volume, weight, total_weight
    real(kind=dp) :: lower_volume, upper_volume, center_volume
    real(kind=dp) :: lower_function, upper_function, center_function
    integer :: i_bin, i_part, i_bisect
    type(aero_particle_t), pointer :: aero_particle

    do i_bin = 1,bin_grid%n_bin
       if (aero_state%bin(i_bin)%n_part == 0) then
          cycle
       end if

       total_weight = 0d0
       weighted_total_volume = 0d0
       do i_part = 1,aero_state%bin(i_bin)%n_part
          aero_particle => aero_state%bin(i_bin)%particle(i_part)
          weight = aero_weight_value(aero_weight, &
               aero_particle_radius(aero_particle))
          total_weight = total_weight + weight
          particle_volume = aero_particle_volume_maybe_dry(aero_particle, &
               aero_data, dry_volume)
          weighted_total_volume = weighted_total_volume &
               + weight * particle_volume
       end do

       ! determine the new_particle_volume for all particles in this bin
       if (bin_center) then
          new_particle_volume = rad2vol(bin_grid%center_radius(i_bin))
       elseif (aero_weight%type == AERO_WEIGHT_TYPE_NONE) then
          new_particle_volume = weighted_total_volume &
               / real(aero_state%bin(i_bin)%n_part, kind=dp)
       elseif (preserve_number) then
          ! number-preserving scheme: Solve the implicit equation:
          ! n_part * W(new_vol) = total_weight
          !
          ! We assume that the weighting function is strictly monotone
          ! so this equation has a unique solution and the solution
          ! lies between the min and max particle volumes. We use
          ! bisection as this doesn't really need to be fast, just
          ! robust.

          ! initialize to min and max particle volumes
          do i_part = 1,aero_state%bin(i_bin)%n_part
             aero_particle => aero_state%bin(i_bin)%particle(i_part)
             particle_volume = aero_particle_volume_maybe_dry(aero_particle, &
                  aero_data, dry_volume)
             if (i_part == 1) then
                lower_volume = particle_volume
                upper_volume = particle_volume
             end if
             lower_volume = min(lower_volume, particle_volume)
             upper_volume = max(upper_volume, particle_volume)
          end do

          lower_function = real(aero_state%bin(i_bin)%n_part, kind=dp) &
               * aero_weight_value(aero_weight, vol2rad(lower_volume)) &
               - total_weight
          upper_function = real(aero_state%bin(i_bin)%n_part, kind=dp) &
               * aero_weight_value(aero_weight, vol2rad(upper_volume)) &
               - total_weight

          ! do 50 rounds of bisection (2^50 = 10^15)
          do i_bisect = 1,50
             center_volume = (lower_volume + upper_volume) / 2d0
             center_function = real(aero_state%bin(i_bin)%n_part, kind=dp) &
                  * aero_weight_value(aero_weight, vol2rad(center_volume)) &
                  - total_weight
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
          ! n_part * W(new_vol) * new_vol = weighted_total_volume
          !
          ! We assume that the weighting function is strictly monotone
          ! so this equation has a unique solution and the solution
          ! lies between the min and max particle volumes. We use
          ! bisection as this doesn't really need to be fast, just
          ! robust.

          ! initialize to min and max particle volumes
          do i_part = 1,aero_state%bin(i_bin)%n_part
             aero_particle => aero_state%bin(i_bin)%particle(i_part)
             particle_volume = aero_particle_volume_maybe_dry(aero_particle, &
                  aero_data, dry_volume)
             if (i_part == 1) then
                lower_volume = particle_volume
                upper_volume = particle_volume
             end if
             lower_volume = min(lower_volume, particle_volume)
             upper_volume = max(upper_volume, particle_volume)
          end do

          lower_function = real(aero_state%bin(i_bin)%n_part, kind=dp) &
               * aero_weight_value(aero_weight, vol2rad(lower_volume)) &
               * lower_volume - weighted_total_volume
          upper_function = real(aero_state%bin(i_bin)%n_part, kind=dp) &
               * aero_weight_value(aero_weight, vol2rad(upper_volume)) &
               * upper_volume - weighted_total_volume

          ! do 50 rounds of bisection (2^50 = 10^15)
          do i_bisect = 1,50
             center_volume = (lower_volume + upper_volume) / 2d0
             center_function = real(aero_state%bin(i_bin)%n_part, kind=dp) &
                  * aero_weight_value(aero_weight, vol2rad(center_volume)) &
                  * center_volume - weighted_total_volume
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

       do i_part = 1,aero_state%bin(i_bin)%n_part
          aero_particle => aero_state%bin(i_bin)%particle(i_part)
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
  subroutine aero_state_make_dry(aero_state, bin_grid, aero_data)

    !> Aerosol state to make dry.
    type(aero_state_t), intent(inout) :: aero_state
    !> Bin grid.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data

    integer :: i_bin, i_part

    if (aero_data%i_water > 0) then
       do i_bin = 1,bin_grid%n_bin
          do i_part = 1,aero_state%bin(i_bin)%n_part
             aero_state%bin(i_bin)%particle(i_part)%vol(aero_data%i_water) &
                  = 0d0
          end do
       end do
       call aero_state_resort(bin_grid, aero_state)
    end if

   end subroutine aero_state_make_dry

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Determines the number of bytes required to pack the given value.
  integer function pmc_mpi_pack_size_aero_state(val)

    !> Value to pack.
    type(aero_state_t), intent(in) :: val

    integer :: i, total_size

    total_size = 0
    total_size = total_size + pmc_mpi_pack_size_real(val%comp_vol)
    total_size = total_size + pmc_mpi_pack_size_integer(val%n_part)
    total_size = total_size + pmc_mpi_pack_size_integer(size(val%bin))
    do i = 1,size(val%bin)
       total_size = total_size + pmc_mpi_pack_size_apa(val%bin(i))
    end do
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
    integer :: prev_position, i

    prev_position = position
    call pmc_mpi_pack_real(buffer, position, val%comp_vol)
    call pmc_mpi_pack_integer(buffer, position, val%n_part)
    call pmc_mpi_pack_integer(buffer, position, size(val%bin))
    do i = 1,size(val%bin)
       call pmc_mpi_pack_aero_particle_array(buffer, position, val%bin(i))
    end do
    call pmc_mpi_pack_aero_info_array(buffer, position, val%aero_info_array)
    call assert(850997402, &
         position - prev_position == pmc_mpi_pack_size_aero_state(val))
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
    integer :: prev_position, i, n

    call aero_state_deallocate(val)
    prev_position = position
    call pmc_mpi_unpack_real(buffer, position, val%comp_vol)
    call pmc_mpi_unpack_integer(buffer, position, val%n_part)
    call pmc_mpi_unpack_integer(buffer, position, n)
    allocate(val%bin(n))
    do i = 1,size(val%bin)
       call aero_particle_array_allocate(val%bin(i))
       call pmc_mpi_unpack_aero_particle_array(buffer, position, val%bin(i))
    end do
    call aero_info_array_allocate(val%aero_info_array)
    call pmc_mpi_unpack_aero_info_array(buffer, position, val%aero_info_array)
    call assert(132104747, &
         position - prev_position == pmc_mpi_pack_size_aero_state(val))
#endif

  end subroutine pmc_mpi_unpack_aero_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Gathers data from all processors into one aero_state on node 0.
  subroutine aero_state_mpi_gather(aero_state, aero_state_total)

    !> Local aero_state.
    type(aero_state_t), intent(in) :: aero_state
    !> Centralized aero_state (only on node 0).
    type(aero_state_t), intent(inout) :: aero_state_total

    integer, parameter :: tag1 = 1053
    integer, parameter :: tag2 = 5658
    integer :: rank, n_proc
#ifdef PMC_USE_MPI
    type(aero_state_t) :: aero_state_transfer
    integer :: ierr, status(MPI_STATUS_SIZE), buffer_size, i_proc, position
    character, allocatable :: buffer(:)
#endif

    rank = pmc_mpi_rank()
    n_proc = pmc_mpi_size()
    if (pmc_mpi_rank() == 0) then
       call aero_state_copy(aero_state, aero_state_total)
    end if

#ifdef PMC_USE_MPI
    
    if (rank /= 0) then
       ! compute and send buffer_size from remote node
       buffer_size = pmc_mpi_pack_size_aero_state(aero_state)
       call mpi_send(buffer_size, 1, MPI_INTEGER, 0, tag1, &
            MPI_COMM_WORLD, ierr)
       call pmc_mpi_check_ierr(ierr)

       ! send buffer from remote node
       allocate(buffer(buffer_size))
       position = 0
       call pmc_mpi_pack_aero_state(buffer, position, aero_state)
       call assert(542772170, position == buffer_size)
       call mpi_send(buffer, buffer_size, MPI_CHARACTER, 0, tag2, &
            MPI_COMM_WORLD, ierr)
       call pmc_mpi_check_ierr(ierr)
       deallocate(buffer)
    else
       ! root node receives data
       
       ! transfer from other nodes
       do i_proc = 1,(n_proc - 1)
          ! get buffer_size at root node
          call mpi_recv(buffer_size, 1, MPI_INTEGER, i_proc, tag1, &
               MPI_COMM_WORLD, status, ierr)
          call pmc_mpi_check_ierr(ierr)

          ! get buffer at root node and unpack it
          allocate(buffer(buffer_size))
          call mpi_recv(buffer, buffer_size, MPI_CHARACTER, i_proc, tag2, &
               MPI_COMM_WORLD, status, ierr)
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

  !> Scatters data from node 0 to all processors by assigning each
  !> particle to a random node.
  subroutine aero_state_mpi_scatter(aero_state_total, aero_state, &
       aero_data)

    !> Centralized aero_state (only on node 0.
    type(aero_state_t), intent(in) :: aero_state_total
    !> Local aero_state.
    type(aero_state_t), intent(inout) :: aero_state
    !> Aero_data.
    type(aero_data_t), intent(in) :: aero_data

#ifndef PMC_USE_MPI
    call aero_state_copy(aero_state_total, aero_state)
#else

    integer, parameter :: tag1 = 137433114
    integer, parameter :: tag2 = 384004297
    integer :: rank, n_proc, i_bin, i_part
    type(aero_state_t), allocatable :: aero_state_transfers(:)
    integer :: ierr, status(MPI_STATUS_SIZE), buffer_size, i_proc, position
    character, allocatable :: buffer(:)
    type(aero_particle_t), pointer :: aero_particle

    rank = pmc_mpi_rank()
    n_proc = pmc_mpi_size()

    if (rank == 0) then
       ! assign particles at random to other processors
       allocate(aero_state_transfers(n_proc))
       do i_proc = 1,n_proc
          call aero_state_allocate_size(aero_state_transfers(i_proc), &
               size(aero_state_total%bin), aero_data%n_spec)
          aero_state_transfers(i_proc)%comp_vol = &
               aero_state_total%comp_vol / real(n_proc, kind=dp)
       end do
       do i_bin = 1,size(aero_state_total%bin)
          do i_part = 1,aero_state_total%bin(i_bin)%n_part
             aero_particle => aero_state_total%bin(i_bin)%particle(i_part)
             i_proc = pmc_rand_int(n_proc)
             call aero_state_add_particle(aero_state_transfers(i_proc), &
                  i_bin, aero_particle)
          end do
       end do

       ! just copy our own transfer aero_state
       call aero_state_copy(aero_state_transfers(1), aero_state)

       ! send the transfer aero_states
       do i_proc = 1,(n_proc - 1)
          ! compute and send buffer_size to remote node
          buffer_size = pmc_mpi_pack_size_aero_state( &
               aero_state_transfers(i_proc + 1))
          call mpi_send(buffer_size, 1, MPI_INTEGER, i_proc, tag1, &
               MPI_COMM_WORLD, ierr)
          call pmc_mpi_check_ierr(ierr)

          ! send buffer to remote node
          allocate(buffer(buffer_size))
          position = 0
          call pmc_mpi_pack_aero_state(buffer, position, &
               aero_state_transfers(i_proc + 1))
          call assert(781876859, position == buffer_size)
          call mpi_send(buffer, buffer_size, MPI_CHARACTER, i_proc, tag2, &
               MPI_COMM_WORLD, ierr)
          call pmc_mpi_check_ierr(ierr)
          deallocate(buffer)
       end do

       ! all done, free data
       do i_proc = 1,n_proc
          call aero_state_deallocate(aero_state_transfers(i_proc))
       end do
       deallocate(aero_state_transfers)
    else
       ! remote nodes just receive data
       
       ! get buffer_size at root node
       call mpi_recv(buffer_size, 1, MPI_INTEGER, 0, tag1, &
            MPI_COMM_WORLD, status, ierr)
       call pmc_mpi_check_ierr(ierr)
       
       ! get buffer at root node and unpack it
       allocate(buffer(buffer_size))
       call mpi_recv(buffer, buffer_size, MPI_CHARACTER, 0, tag2, &
            MPI_COMM_WORLD, status, ierr)
       position = 0
       
       call pmc_mpi_unpack_aero_state(buffer, position, aero_state)
       call assert(492983285, position == buffer_size)
       deallocate(buffer)
    end if

#endif

  end subroutine aero_state_mpi_scatter

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
    integer :: aero_particle_centers(aero_state%n_part)

    ! try to get the dimension ID
    status = nf90_inq_dimid(ncid, "aero_particle", dimid_aero_particle)
    if (status == NF90_NOERR) return
    if (status /= NF90_EBADDIM) call pmc_nc_check(status)

    ! dimension not defined, so define now define it
    call pmc_nc_check(nf90_redef(ncid))

    call pmc_nc_check(nf90_def_dim(ncid, "aero_particle", &
         aero_state%n_part, dimid_aero_particle))

    call pmc_nc_check(nf90_enddef(ncid))

    do i_part = 1,aero_state%n_part
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
  subroutine aero_state_output_netcdf(aero_state, ncid, bin_grid, &
       aero_data, aero_weight, record_removals, record_optical)
    
    !> aero_state to write.
    type(aero_state_t), intent(in) :: aero_state
    !> NetCDF file ID, in data mode.
    integer, intent(in) :: ncid
    !> bin_grid structure.
    type(bin_grid_t), intent(in) :: bin_grid
    !> aero_data structure.
    type(aero_data_t), intent(in) :: aero_data
    !> aero_weight structure.
    type(aero_weight_t), intent(in) :: aero_weight
    !> Whether to output particle removal info.
    logical, intent(in) :: record_removals
    !> Whether to output aerosol optical properties.
    logical, intent(in) :: record_optical

    integer :: dimid_aero_particle, dimid_aero_species
    integer :: dimid_aero_removed
    integer :: i_bin, i_part_in_bin, i_part, i_remove
    type(aero_particle_t), pointer :: particle
    real(kind=dp) :: aero_particle_mass(aero_state%n_part, aero_data%n_spec)
    integer :: aero_n_orig_part(aero_state%n_part)
    real(kind=dp) :: aero_absorb_cross_sect(aero_state%n_part)
    real(kind=dp) :: aero_scatter_cross_sect(aero_state%n_part)
    real(kind=dp) :: aero_asymmetry(aero_state%n_part)
    real(kind=dp) :: aero_refract_shell_real(aero_state%n_part)
    real(kind=dp) :: aero_refract_shell_imag(aero_state%n_part)
    real(kind=dp) :: aero_refract_core_real(aero_state%n_part)
    real(kind=dp) :: aero_refract_core_imag(aero_state%n_part)
    real(kind=dp) :: aero_core_vol(aero_state%n_part)
    integer :: aero_water_hyst_leg(aero_state%n_part)
    real(kind=dp) :: aero_comp_vol(aero_state%n_part)
    integer :: aero_id(aero_state%n_part)
    real(kind=dp) :: aero_least_create_time(aero_state%n_part)
    real(kind=dp) :: aero_greatest_create_time(aero_state%n_part)
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
    !!   - \b aero_n_orig_part (dim \c aero_particle): number of original
    !!     particles that coagulated to form each aerosol particle - when
    !!     particles first enter the simulation (by emissions, dilution,
    !!     etc.) they have <tt>aero_n_orig_part = 1</tt> and when two
    !!     particles coagulate, their values of \c aero_n_orig_part are
    !!     added (the number of coagulation events that formed each particle
    !!     is thus <tt>aero_n_orig_part - 1</tt>)
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

    call aero_data_netcdf_dim_aero_species(aero_data, ncid, &
         dimid_aero_species)
    
    if (aero_state%n_part > 0) then
       call aero_state_netcdf_dim_aero_particle(aero_state, ncid, &
            dimid_aero_particle)
       
       i_part = 0
       do i_bin = 1,bin_grid%n_bin
          do i_part_in_bin = 1,aero_state%bin(i_bin)%n_part
             i_part = i_part + 1
             particle => aero_state%bin(i_bin)%particle(i_part_in_bin)
             aero_particle_mass(i_part, :) = particle%vol * aero_data%density
             aero_n_orig_part(i_part) = particle%n_orig_part
             aero_water_hyst_leg(i_part) = particle%water_hyst_leg
             aero_comp_vol(i_part) = aero_state%comp_vol &
                  / aero_weight_value(aero_weight, &
                  aero_particle_radius(particle))
             aero_id(i_part) = particle%id
             aero_least_create_time(i_part) = particle%least_create_time
             aero_greatest_create_time(i_part) = particle%greatest_create_time
             if (record_optical) then
                aero_absorb_cross_sect(i_part) = particle%absorb_cross_sect
                aero_scatter_cross_sect(i_part) = particle%scatter_cross_sect
                aero_asymmetry(i_part) = particle%asymmetry
                aero_refract_shell_real(i_part) = real(particle%refract_shell)
                aero_refract_shell_imag(i_part) = &
                     aimag(particle%refract_shell)
                aero_refract_core_real(i_part) = real(particle%refract_core)
                aero_refract_core_imag(i_part) = aimag(particle%refract_core)
                aero_core_vol(i_part) = particle%core_vol
             end if
          end do
       end do
       call pmc_nc_write_real_2d(ncid, aero_particle_mass, &
            "aero_particle_mass", (/ dimid_aero_particle, &
            dimid_aero_species /), unit="kg", &
            long_name="constituent masses of each aerosol particle")
       call pmc_nc_write_integer_1d(ncid, aero_n_orig_part, &
            "aero_n_orig_part", (/ dimid_aero_particle /), &
            long_name="number of original constituent particles that " &
            // "coagulated to form each aerosol particle")
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
  subroutine aero_state_input_netcdf(aero_state, ncid, bin_grid, &
       aero_data, aero_weight)
    
    !> aero_state to read.
    type(aero_state_t), intent(inout) :: aero_state
    !> NetCDF file ID, in data mode.
    integer, intent(in) :: ncid
    !> bin_grid structure.
    type(bin_grid_t), intent(in) :: bin_grid
    !> aero_data structure.
    type(aero_data_t), intent(in) :: aero_data
    !> aero_weight structure.
    type(aero_weight_t), intent(in) :: aero_weight

    integer :: dimid_aero_particle, dimid_aero_removed, n_info_item, n_part
    integer :: i_bin, i_part_in_bin, i_part, i_remove, status
    type(aero_particle_t) :: aero_particle
    character(len=1000) :: name

    real(kind=dp), allocatable :: aero_particle_mass(:,:)
    integer, allocatable :: aero_n_orig_part(:)
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
       call aero_state_allocate_size(aero_state, bin_grid%n_bin, &
            aero_data%n_spec)
       return
    end if
    call pmc_nc_check(status)
    call pmc_nc_check(nf90_Inquire_Dimension(ncid, dimid_aero_particle, &
         name, n_part))

    allocate(aero_particle_mass(n_part, aero_data%n_spec))
    allocate(aero_n_orig_part(n_part))
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
    call pmc_nc_read_integer_1d(ncid, aero_n_orig_part, &
         "aero_n_orig_part")
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
    call aero_state_allocate_size(aero_state, bin_grid%n_bin, &
         aero_data%n_spec)

    call aero_particle_allocate_size(aero_particle, aero_data%n_spec)
    do i_part = 1,n_part
       aero_particle%vol = aero_particle_mass(i_part, :) / aero_data%density
       aero_particle%n_orig_part = aero_n_orig_part(i_part)
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
       aero_state%comp_vol = aero_comp_vol(i_part) &
            * aero_weight_value(aero_weight, &
            aero_particle_radius(aero_particle))
       aero_particle%id = aero_id(i_part)
       aero_particle%least_create_time = aero_least_create_time(i_part)
       aero_particle%greatest_create_time = aero_greatest_create_time(i_part)

       i_bin = aero_particle_in_bin(aero_particle, bin_grid)
       call aero_state_add_particle(aero_state, i_bin, aero_particle)
    end do
    call aero_particle_deallocate(aero_particle)

    deallocate(aero_particle_mass)
    deallocate(aero_n_orig_part)
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
  
end module pmc_aero_state
