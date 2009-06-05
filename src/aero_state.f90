! Copyright (C) 2005-2009 Nicole Riemer and Matthew West
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
  use pmc_inout
  use pmc_aero_info
  use pmc_aero_info_array
#ifdef PMC_USE_MPI
  use mpi
#endif

  !> The current collection of aerosol particles.
  !!
  !! The particles in aero_state_t are stored sorted per-bin, to
  !! improve efficiency of access and sampling. If a particle has
  !! total volume \c v then calling <tt> i_bin =
  !! bin_grid_particle_in_bin(bin_grid, v)</tt> finds the bin number
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
     real*8 :: comp_vol
     !> Total number of particles.
     integer :: n_part
     !> Information on removed particles.
     type(aero_info_array_t) :: aero_info_array
  end type aero_state_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Initializes aerosol arrays to have zero particles in each
  !> bin. Do not call this more than once on a given aerosol, use
  !> aero_state_zero() instead to reset to zero.
  subroutine aero_state_alloc(n_bin, n_spec, aero_state)

    !> Number of bins.
    integer, intent(in) :: n_bin
    !> Number of species.
    integer, intent(in) :: n_spec
    !> Aerosol to initialize.
    type(aero_state_t), intent(inout) :: aero_state
    
    integer i

    allocate(aero_state%bin(n_bin))
    do i = 1,n_bin
       call aero_particle_array_alloc(aero_state%bin(i), 0, n_spec)
    end do
    aero_state%comp_vol = 0d0
    aero_state%n_part = 0
    call aero_info_array_alloc(aero_state%aero_info_array, 0)

  end subroutine aero_state_alloc
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Deallocates a previously allocated aerosol.
  subroutine aero_state_free(aero_state)

    !> Aerosol to initialize.
    type(aero_state_t), intent(inout) :: aero_state
    
    integer :: n_bin, i

    n_bin = size(aero_state%bin)
    do i = 1,n_bin
       call aero_particle_array_free(aero_state%bin(i))
    end do
    deallocate(aero_state%bin)
    call aero_info_array_free(aero_state%aero_info_array)

  end subroutine aero_state_free
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Copies aerosol to a destination that has already had
  !> aero_state_alloc() called on it.
  subroutine aero_state_copy(aero_state_from, aero_state_to)

    !> Reference aerosol.
    type(aero_state_t), intent(in) :: aero_state_from
    !> Already allocated.
    type(aero_state_t), intent(inout) :: aero_state_to
    
    integer :: n_bin, i

    n_bin = size(aero_state_from%bin)

    call aero_state_free(aero_state_to)
    call aero_state_alloc(n_bin, 0, aero_state_to)

    do i = 1,n_bin
       call aero_particle_array_copy(aero_state_from%bin(i), &
            aero_state_to%bin(i))
    end do

    aero_state_to%comp_vol = aero_state_from%comp_vol
    aero_state_to%n_part = aero_state_from%n_part

    call aero_info_array_copy(aero_state_from%aero_info_array, &
         aero_state_to%aero_info_array)

  end subroutine aero_state_copy
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the total number of particles in an aerosol distribution.
  integer function aero_state_total_particles(aero_state)

    !> Aerosol state.
    type(aero_state_t), intent(in) :: aero_state

    aero_state_total_particles = aero_state%n_part

  end function aero_state_total_particles

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Resets an aero_state to have zero particles per bin. This must
  !> already have had aero_state_alloc() called on it. This
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
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Remove the given particle.
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

    call aero_particle_array_remove_particle(aero_state%bin(i_bin), index)
    aero_state%n_part = aero_state%n_part - 1
    if (record_removal) then
       call aero_info_array_add_aero_info(aero_state%aero_info_array, &
            aero_info)
    end if

  end subroutine aero_state_remove_particle

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> aero_state += aero_state_delta
  subroutine aero_state_add(aero_state, aero_state_delta)

    !> Aerosol state.
    type(aero_state_t), intent(inout) :: aero_state
    !> Increment.
    type(aero_state_t), intent(in) :: aero_state_delta

    integer :: i_bin, i_part

    do i_bin = 1,size(aero_state_delta%bin)
       do i_part = 1,aero_state_delta%bin(i_bin)%n_part
          call aero_state_add_particle(aero_state, i_bin, &
               aero_state_delta%bin(i_bin)%particle(i_part))
       end do
    end do
    aero_state%comp_vol = aero_state%comp_vol + aero_state_delta%comp_vol

  end subroutine aero_state_add

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Generates a Poisson sample of an aero_dist, adding to
  !> aero_state. The sampled amount is sample_prop *
  !> aero_state%comp_vol.
  subroutine aero_state_add_aero_dist_sample(aero_state, bin_grid, &
       aero_data, aero_dist, sample_prop, create_time)

    !> Aero state to add to.
    type(aero_state_t), intent(inout) :: aero_state
    !> Bin grid.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Aero data values.
    type(aero_data_t), intent(in) :: aero_data
    !> Distribution to sample.
    type(aero_dist_t), intent(in) :: aero_dist
    !> Volume fraction to sample (1).
    real*8, intent(in) :: sample_prop
    !> Creation time for new particles (s).
    real*8, intent(in) :: create_time

    real*8 :: n_samp_avg, sample_vol, radius, vol
    integer :: n_samp, i_mode, i_samp, i_bin
    integer :: num_per_bin(bin_grid%n_bin)
    type(aero_mode_t), pointer :: aero_mode
    type(aero_particle_t) :: aero_particle

    call aero_particle_alloc(aero_particle, aero_data%n_spec)
    sample_vol = sample_prop * aero_state%comp_vol
    do i_mode = 1,aero_dist%n_mode
       aero_mode => aero_dist%mode(i_mode)
       n_samp_avg = sample_vol * aero_mode%num_den
       n_samp = rand_poisson(n_samp_avg)
       do i_samp = 1,n_samp
          call aero_mode_sample_radius(aero_mode, radius)
          vol = rad2vol(radius)
          call aero_particle_set_vols(aero_particle, aero_mode%vol_frac * vol)
          call aero_particle_new_id(aero_particle)
          call aero_particle_set_create_time(aero_particle, create_time)
          i_bin = aero_particle_in_bin(aero_particle, bin_grid)
          call aero_state_add_particle(aero_state, i_bin, aero_particle)
       end do
    end do
    call aero_particle_free(aero_particle)

  end subroutine aero_state_add_aero_dist_sample
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
    i_bin = sample_disc_pdf(n_bin, disc_pdf)
    i_part = pmc_rand_int(aero_state%bin(i_bin)%n_part)

  end subroutine aero_state_rand_particle

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
    real*8, intent(in) :: sample_prop
    !> Action for removal (see pmc_aero_info module for action
    !> parameters). Set to AERO_INFO_NONE to not log removal.
    integer, intent(in) :: removal_action
    
    integer :: n_transfer, i_transfer, n_bin, i_bin, i_part
    logical :: do_add, do_remove
    real*8 :: vol_ratio
    type(aero_info_t) :: aero_info

    call assert(721006962, (sample_prop >= 0d0) .and. (sample_prop <= 1d0))
    n_transfer = rand_poisson(sample_prop &
         * dble(aero_state_total_particles(aero_state_from)))
    n_bin = size(aero_state_from%bin)
    vol_ratio = aero_state_to%comp_vol / aero_state_from%comp_vol
    i_transfer = 0
    do while (i_transfer < n_transfer)
       if (aero_state_total_particles(aero_state_from) <= 0) exit
       call aero_state_rand_particle(aero_state_from, i_bin, i_part)
       if (vol_ratio == 1d0) then
          ! to_comp_vol == from_comp_vol so just move the particle
          do_add = .true.
          do_remove = .true.
       elseif (vol_ratio > 1d0) then
          ! to_comp_vol is bigger than from_comp_vol, so only maybe
          ! remove the particle but always add it
          do_add = .true.
          do_remove = .false.
          if (pmc_random() < 1d0 / vol_ratio) then
             do_remove = .true.
          end if
       else ! vol_ratio < 1d0
          ! to_comp_vol is smaller than from_comp_vol, so always
          ! remove the particle but only maybe add it
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
             call aero_info_alloc(aero_info)
             aero_info%id = &
                  aero_state_from%bin(i_bin)%particle(i_part)%id
             aero_info%action = removal_action
             call aero_state_remove_particle(aero_state_from, i_bin, &
                  i_part, .true., aero_info)
             call aero_info_free(aero_info)
          else
             call aero_state_remove_particle(aero_state_from, i_bin, &
                  i_part, .false., aero_info)
          end if
          i_transfer = i_transfer + 1
       end if
    end do
    
  end subroutine aero_state_sample
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
    real*8, intent(in) :: sample_prop
    !> Action for removal (see pmc_aero_info module for action
    !> parameters). Set to AERO_INFO_NONE to not log removal.
    integer, intent(in) :: removal_action
    
    integer :: n_transfer, i_transfer, n_bin, i_bin, i_part
    logical :: do_add, do_remove
    real*8 :: vol_ratio
    type(aero_info_t) :: aero_info

    n_bin = size(aero_state_from%bin)
    vol_ratio = aero_state_to%comp_vol / aero_state_from%comp_vol
    do i_bin = 1,n_bin
       n_transfer = prob_round(sample_prop &
            * dble(aero_state_from%bin(i_bin)%n_part))
       i_transfer = 0
       do while (i_transfer < n_transfer)
          if (aero_state_from%bin(i_bin)%n_part <= 0) exit
          i_part = pmc_rand_int(aero_state_from%bin(i_bin)%n_part)
          if (vol_ratio == 1d0) then
             ! to_comp_vol == from_comp_vol so just move the particle
             do_add = .true.
             do_remove = .true.
          elseif (vol_ratio > 1d0) then
             ! to_comp_vol is bigger than from_comp_vol, so only maybe
             ! remove the particle but always add it
             do_add = .true.
             do_remove = .false.
             if (pmc_random() < 1d0 / vol_ratio) then
                do_remove = .true.
             end if
          else ! vol_ratio < 1d0
             ! to_comp_vol is smaller than from_comp_vol, so always
             ! remove the particle but only maybe add it
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
                call aero_info_alloc(aero_info)
                aero_info%id = &
                     aero_state_from%bin(i_bin)%particle(i_part)%id
                aero_info%action = removal_action
                call aero_state_remove_particle(aero_state_from, i_bin, &
                     i_part, .true., aero_info)
                call aero_info_free(aero_info)
             else
                call aero_state_remove_particle(aero_state_from, i_bin, &
                     i_part, .false., aero_info)
             end if
             i_transfer = i_transfer + 1
          end if
       end do
    end do
    
  end subroutine aero_state_sample_rough
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Adds aero_state_delta particles to aero_state. The number of
  !> particles added depends on the computational volume ratio, so
  !> either more or less particles might be added than are actually
  !> in aero_state_delta.
  subroutine aero_state_add_particles(aero_state, aero_state_delta)

    !> Base state.
    type(aero_state_t), intent(inout) :: aero_state
    !> State to add.
    type(aero_state_t), intent(in) :: aero_state_delta
    
    integer :: n_bin, i_bin, i_part, n_new_part, i_new_part
    real*8 :: vol_ratio
    
    n_bin = size(aero_state%bin)
    vol_ratio = aero_state%comp_vol / aero_state_delta%comp_vol
    do i_bin = 1,n_bin
       do i_part = 1,aero_state_delta%bin(i_bin)%n_part
          n_new_part = prob_round(vol_ratio)
          do i_new_part = 1,n_new_part
             call aero_state_add_particle(aero_state, i_bin, &
                  aero_state_delta%bin(i_bin)%particle(i_part))
          end do
       end do
    end do
    
  end subroutine aero_state_add_particles
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
    type(aero_binned_t), intent(out) :: aero_binned
    
    integer :: b, j, s
    type(aero_particle_t), pointer :: aero_particle
    
    aero_binned%num_den = 0d0
    aero_binned%vol_den = 0d0
    do b = 1,bin_grid%n_bin
       do j = 1,aero_state%bin(b)%n_part
          aero_particle => aero_state%bin(b)%particle(j)
          aero_binned%vol_den(b,:) = aero_binned%vol_den(b,:) &
               + aero_particle%vol / aero_state%comp_vol / bin_grid%dlnr
          aero_binned%num_den(b) = aero_binned%num_den(b) &
               + 1d0 / aero_state%comp_vol / bin_grid%dlnr
       end do
    end do
    
  end subroutine aero_state_to_binned
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
    type(aero_binned_t), intent(out) :: aero_binned
    
    integer :: b, j, s, b_dry
    type(aero_particle_t), pointer :: aero_particle
    
    aero_binned%num_den = 0d0
    aero_binned%vol_den = 0d0
    do b = 1,bin_grid%n_bin
       do j = 1,aero_state%bin(b)%n_part
          aero_particle => aero_state%bin(b)%particle(j)
          b_dry = bin_grid_particle_in_bin(bin_grid, &
               aero_particle_solute_volume(aero_particle, aero_data))
          aero_binned%vol_den(b_dry,:) = aero_binned%vol_den(b_dry,:) &
               + aero_particle%vol / aero_state%comp_vol / bin_grid%dlnr
          aero_binned%num_den(b_dry) = aero_binned%num_den(b_dry) &
               + 1d0 / aero_state%comp_vol / bin_grid%dlnr
       end do
    end do
    
  end subroutine aero_state_to_binned_dry
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Remove approximately half of the particles in each bin.
  subroutine aero_state_halve(aero_state, aero_binned, bin_grid)
    
    !> Aerosol state.
    type(aero_state_t), intent(inout) :: aero_state
    !> Aero binned.
    type(aero_binned_t), intent(inout) :: aero_binned
    !> Bin grid.
    type(bin_grid_t), intent(in) :: bin_grid
    
    integer :: i_bin, i_part, n_part_orig, i_remove, n_remove
    type(aero_info_t) :: aero_info

    n_part_orig = aero_state%n_part
    do i_bin = 1,bin_grid%n_bin
       n_remove = prob_round(dble(aero_state%bin(i_bin)%n_part) / 2d0)
       do i_remove = 1,n_remove
          i_part = pmc_rand_int(aero_state%bin(i_bin)%n_part)
          call aero_binned_remove_particle_in_bin(aero_binned, bin_grid, &
               i_bin, aero_state%comp_vol, &
               aero_state%bin(i_bin)%particle(i_part))
          call aero_info_alloc(aero_info)
          aero_info%id = &
               aero_state%bin(i_bin)%particle(i_part)%id
          aero_info%action = AERO_INFO_HALVED
          call aero_state_remove_particle(aero_state, i_bin, &
               i_part, .true., aero_info)
          call aero_info_free(aero_info)
       end do
    end do
    aero_state%comp_vol = aero_state%comp_vol &
         * dble(aero_state%n_part) / dble(n_part_orig)

  end subroutine aero_state_halve
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Send a sample to the given process, and receive exactly one
  !> sample from an unspecified source.
  subroutine aero_state_mix_to(aero_state, mix_rate, dest, &
       aero_binned, aero_data, bin_grid)

    !> Aerosol state.
    type(aero_state_t), intent(inout) :: aero_state
    !> Mixing rate (0 to 1).
    real*8, intent(in) :: mix_rate
    !> Process to send to.
    integer, intent(in) :: dest
    !> Aero binned to update.
    type(aero_binned_t), intent(inout) :: aero_binned
    !> Aero data values.
    type(aero_data_t), intent(in) :: aero_data
    !> Bin grid.
    type(bin_grid_t), intent(in) :: bin_grid

#ifdef PMC_USE_MPI
    integer :: ierr, status(MPI_STATUS_SIZE)
    type(aero_state_t) :: aero_state_send, aero_state_recv
    character, allocatable :: buffer_send(:), buffer_recv(:)
    integer :: buffer_size_send, buffer_size_recv, position
    type(aero_binned_t) :: aero_binned_delta

    if (pmc_mpi_rank() == dest) then
       return
    end if

    ! allocate memory
    call aero_binned_alloc(aero_binned_delta, bin_grid%n_bin, &
         aero_data%n_spec)

    ! extract particles to send
    call aero_state_alloc(bin_grid%n_bin, aero_data%n_spec, &
         aero_state_send)
    aero_state_send%comp_vol = aero_state%comp_vol
    ! FIXME: would probably be slightly better for sampling purposes
    ! to use the destination comp_vol here
    call aero_state_sample_rough(aero_state, aero_state_send, &
         mix_rate, AERO_INFO_NONE)
    call aero_state_to_binned(bin_grid, aero_data, aero_state_send, &
         aero_binned_delta)
    call aero_binned_sub(aero_binned, aero_binned_delta)

    ! pack up data to send
    buffer_size_send = pmc_mpi_pack_size_aero_state(aero_state_send)
    allocate(buffer_send(buffer_size_send))
    position = 0
    call pmc_mpi_pack_aero_state(buffer_send, position, &
         aero_state_send)
    call assert(420102502, position == buffer_size_send)

    ! transfer buffer size info and allocate receive buffers
    call mpi_sendrecv(buffer_size_send, 1, MPI_INTEGER, dest, 0, &
         buffer_size_recv, 1, MPI_INTEGER, MPI_ANY_SOURCE, &
         MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
    call pmc_mpi_check_ierr(ierr)
    allocate(buffer_recv(buffer_size_recv))

    ! do data transfer
    call mpi_sendrecv(buffer_send, buffer_size_send, MPI_PACKED, &
         dest, 0, buffer_recv, buffer_size_recv, MPI_PACKED, &
         MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
    call pmc_mpi_check_ierr(ierr)

    ! unpack received data and process it
    position = 0
    call pmc_mpi_unpack_aero_state(buffer_recv, position, &
         aero_state_recv)
    call assert(593694264, position == buffer_size_recv)
    call aero_state_add_particles(aero_state, aero_state_recv)
    call aero_state_to_binned(bin_grid, aero_data, aero_state_recv, &
         aero_binned_delta)
    call aero_binned_add(aero_binned, aero_binned_delta)

    ! cleanup
    call aero_binned_free(aero_binned_delta)
    call aero_state_free(aero_state_send)
    call aero_state_free(aero_state_recv)
    deallocate(buffer_send)
    deallocate(buffer_recv)
#endif

  end subroutine aero_state_mix_to

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Mix the aero_states between all processes. Currently uses a
  !> simple periodic-1D diffusion.
  subroutine aero_state_mix(aero_state, mix_rate, &
       aero_binned, aero_data, bin_grid)

    !> Aerosol state.
    type(aero_state_t), intent(inout) :: aero_state
    !> Mixing rate (0 to 1).
    real*8, intent(in) :: mix_rate
    !> Aero binned to update.
    type(aero_binned_t), intent(inout) :: aero_binned
    !> Aero data values.
    type(aero_data_t), intent(in) :: aero_data
    !> Bin grid.
    type(bin_grid_t), intent(in) :: bin_grid

    integer :: rank, n_proc, dest

    rank = pmc_mpi_rank()
    n_proc = pmc_mpi_size()
    
    ! mix up
    dest = rank + 1
    if (dest == n_proc) then
       dest = 0
    end if
    call aero_state_mix_to(aero_state, mix_rate, dest, &
       aero_binned, aero_data, bin_grid)

    ! synchronize
    call pmc_mpi_barrier()

    ! mix down
    dest = rank - 1
    if (dest == -1) then
       dest = n_proc - 1
    end if
    call aero_state_mix_to(aero_state, mix_rate, dest, &
       aero_binned, aero_data, bin_grid)

  end subroutine aero_state_mix

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Check that all particles are in the correct bins and that the
  !> bin numbers and masses are correct. This is for debugging only.
  subroutine aero_state_check(bin_grid, aero_binned, aero_data, aero_state)
    
    !> Bin_grid.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Binned distributions.
    type(aero_binned_t), intent(out) :: aero_binned
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Aerosol state.
    type(aero_state_t), intent(inout) :: aero_state
    
    real*8 :: check_bin_v, check_vol_den, vol_tol
    real*8 :: num_tol, state_num_den
    integer :: i, k, k_check, s, n_part_check, id, max_id
    logical :: error
    logical, allocatable :: id_present(:)
    
    error = .false.

    ! check that the total number of particles is correct
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
    
    ! check the aero_binned%num_den array
    do k = 1,bin_grid%n_bin
       num_tol = 0.01d0 / aero_state%comp_vol / bin_grid%dlnr
       state_num_den = dble(aero_state%bin(k)%n_part) / aero_state%comp_vol &
            / bin_grid%dlnr
       if (.not. almost_equal_abs(state_num_den, &
            aero_binned%num_den(k), num_tol)) then
          write(0,'(a10,a20,a20,a20,a20)') 'k', 'bins(k)%n_part', &
               'state_num_den', 'num_den(k)', 'comp_vol'
          write(0,'(i10,i20,e20.10,e20.10,e20.10)') k, &
               aero_state%bin(k)%n_part, state_num_den, &
               aero_binned%num_den(k), aero_state%comp_vol
          error = .true.
       end if
    end do
    
    ! check the aero_binned%vol_den array
    do k = 1,bin_grid%n_bin
       vol_tol = bin_grid%v(k) / 1d3 / bin_grid%dlnr
       do s = 1,aero_data%n_spec
          check_vol_den = sum((/(aero_state%bin(k)%particle(i)%vol(s), &
               i = 1,aero_state%bin(k)%n_part)/)) &
               / aero_state%comp_vol / bin_grid%dlnr
          if (.not. almost_equal_abs(check_vol_den, &
               aero_binned%vol_den(k,s), vol_tol)) then
             write(0,'(a10,a10,a25,a25)') 'k', 's', 'check_vol_den', &
                  'vol_den(k,s)'
             write(0,'(i10,i10,e25.10,e25.10)') k, s, check_vol_den, &
                  aero_binned%vol_den(k,s)
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
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Write full state.
  subroutine inout_write_aero_state(file, aero_state)
    
    !> File to write to.
    type(inout_file_t), intent(inout) :: file
    !> Aero_state to write.
    type(aero_state_t), intent(in) :: aero_state

    integer :: n_bin, i
    
    n_bin = size(aero_state%bin)
    call inout_write_comment(file, "begin aero_state")
    call inout_write_real(file, "comp_vol(m^3)", aero_state%comp_vol)
    call inout_write_integer(file, "n_part", aero_state%n_part)
    call inout_write_integer(file, "n_bin", n_bin)
    do i = 1,n_bin
       call inout_write_integer(file, "bin_number", i)
       call inout_write_aero_particle_array(file, aero_state%bin(i))
    end do
    call inout_write_comment(file, "end aero_state")
    
  end subroutine inout_write_aero_state
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read full state.
  subroutine inout_read_aero_state(file, aero_state)
    
    !> File to write to.
    type(inout_file_t), intent(inout) :: file
    !> Aero_state to read.
    type(aero_state_t), intent(out) :: aero_state
    
    integer :: n_bin, i, check_i
    
    call inout_check_comment(file, "begin aero_state")
    call inout_read_real(file, "comp_vol(m^3)", aero_state%comp_vol)
    call inout_read_integer(file, "n_part", aero_state%n_part)
    call inout_read_integer(file, "n_bin", n_bin)
    allocate(aero_state%bin(n_bin))
    do i = 1,n_bin
       call inout_read_integer(file, "bin_number", check_i)
       call inout_check_index(file, i, check_i)
       call inout_read_aero_particle_array(file, aero_state%bin(i))
    end do
    call inout_check_comment(file, "end aero_state")
    
  end subroutine inout_read_aero_state
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
    pmc_mpi_pack_size_aero_state = total_size

  end function pmc_mpi_pack_size_aero_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
    call assert(850997402, &
         position - prev_position == pmc_mpi_pack_size_aero_state(val))
#endif

  end subroutine pmc_mpi_pack_aero_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Unpacks the given value from the buffer, advancing position.
  subroutine pmc_mpi_unpack_aero_state(buffer, position, val)

    !> Memory buffer.
    character, intent(inout) :: buffer(:)
    !> Current buffer position.
    integer, intent(inout) :: position
    !> Value to pack.
    type(aero_state_t), intent(out) :: val

#ifdef PMC_USE_MPI
    integer :: prev_position, i, n

    prev_position = position
    call pmc_mpi_unpack_real(buffer, position, val%comp_vol)
    call pmc_mpi_unpack_integer(buffer, position, val%n_part)
    call pmc_mpi_unpack_integer(buffer, position, n)
    allocate(val%bin(n))
    do i = 1,size(val%bin)
       call pmc_mpi_unpack_aero_particle_array(buffer, position, val%bin(i))
    end do
    call assert(132104747, &
         position - prev_position == pmc_mpi_pack_size_aero_state(val))
#endif

  end subroutine pmc_mpi_unpack_aero_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
    call pmc_nc_check(nf90_def_var(ncid, "aero_particle", NF90_INT, &
         dimid_aero_particle, varid_aero_particle))
    call pmc_nc_check(nf90_put_att(ncid, varid_aero_particle, "unit", "1"))

    call pmc_nc_check(nf90_enddef(ncid))

    do i_part = 1,aero_state%n_part
       aero_particle_centers(i_part) = i_part
    end do
    call pmc_nc_check(nf90_put_var(ncid, varid_aero_particle, &
         aero_particle_centers))

  end subroutine aero_state_netcdf_dim_aero_particle

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
    call pmc_nc_check(nf90_def_var(ncid, "aero_removed", NF90_INT, &
         dimid_aero_removed, varid_aero_removed))
    call pmc_nc_check(nf90_put_att(ncid, varid_aero_removed, "unit", "1"))

    call pmc_nc_check(nf90_enddef(ncid))

    do i_remove = 1,dim_size
       aero_removed_centers(i_remove) = i_remove
    end do
    call pmc_nc_check(nf90_put_var(ncid, varid_aero_removed, &
         aero_removed_centers))

  end subroutine aero_state_netcdf_dim_aero_removed

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Write full state.
  subroutine aero_state_output_netcdf(aero_state, ncid, bin_grid, &
       aero_data, record_removals)
    
    !> aero_state to write.
    type(aero_state_t), intent(in) :: aero_state
    !> NetCDF file ID, in data mode.
    integer, intent(in) :: ncid
    !> bin_grid structure.
    type(bin_grid_t), intent(in) :: bin_grid
    !> aero_data structure.
    type(aero_data_t), intent(in) :: aero_data
    !> Whether to output particle removal info.
    logical, intent(in) :: record_removals

    integer :: dimid_aero_particle, dimid_aero_species
    integer :: dimid_aero_removed
    integer :: i_bin, i_part_in_bin, i_part, i_remove
    type(aero_particle_t), pointer :: particle
    real*8 :: aero_comp_mass(aero_state%n_part, aero_data%n_spec)
    integer :: aero_n_orig_part(aero_state%n_part)
    real*8 :: aero_absorb_cross_sect(aero_state%n_part)
    real*8 :: aero_scatter_cross_sect(aero_state%n_part)
    real*8 :: aero_asymmetry(aero_state%n_part)
    real*8 :: aero_refract_shell_real(aero_state%n_part)
    real*8 :: aero_refract_shell_imag(aero_state%n_part)
    real*8 :: aero_refract_core_real(aero_state%n_part)
    real*8 :: aero_refract_core_imag(aero_state%n_part)
    real*8 :: aero_core_vol(aero_state%n_part)
    integer :: aero_water_hyst_leg(aero_state%n_part)
    real*8 :: aero_comp_vol(aero_state%n_part)
    integer :: aero_id(aero_state%n_part)
    real*8 :: aero_least_create_time(aero_state%n_part)
    real*8 :: aero_greatest_create_time(aero_state%n_part)
    integer :: aero_removed_id(max(aero_state%aero_info_array%n_item,1))
    integer :: aero_removed_action(max(aero_state%aero_info_array%n_item,1))
    integer :: aero_removed_other_id(max(aero_state%aero_info_array%n_item,1))

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
             aero_comp_mass(i_part, :) = particle%vol * aero_data%density
             aero_n_orig_part(i_part) = particle%n_orig_part
             aero_absorb_cross_sect(i_part) = particle%absorb_cross_sect
             aero_scatter_cross_sect(i_part) = particle%scatter_cross_sect
             aero_asymmetry(i_part) = particle%asymmetry
             aero_refract_shell_real(i_part) = real(particle%refract_shell)
             aero_refract_shell_imag(i_part) = aimag(particle%refract_shell)
             aero_refract_core_real(i_part) = real(particle%refract_core)
             aero_refract_core_imag(i_part) = aimag(particle%refract_core)
             aero_core_vol(i_part) = particle%core_vol
             aero_water_hyst_leg(i_part) = particle%water_hyst_leg
             aero_comp_vol(i_part) = aero_state%comp_vol
             aero_id(i_part) = particle%id
             aero_least_create_time(i_part) = particle%least_create_time
             aero_greatest_create_time(i_part) = particle%greatest_create_time
          end do
       end do
       call pmc_nc_write_real_2d(ncid, aero_comp_mass, &
            "aero_comp_mass", "kg", (/ dimid_aero_particle, dimid_aero_species /))
       call pmc_nc_write_integer_1d(ncid, aero_n_orig_part, &
            "aero_n_orig_part", "1", (/ dimid_aero_particle /))
       call pmc_nc_write_real_1d(ncid, aero_absorb_cross_sect, &
            "aero_absorb_cross_sect", "m^2", (/ dimid_aero_particle /))
       call pmc_nc_write_real_1d(ncid, aero_scatter_cross_sect, &
            "aero_scatter_cross_sect", "m^2", (/ dimid_aero_particle /))
       call pmc_nc_write_real_1d(ncid, aero_asymmetry, &
            "aero_asymmetry", "1", (/ dimid_aero_particle /))
       call pmc_nc_write_real_1d(ncid, aero_refract_shell_real, &
            "aero_refract_shell_real", "1", (/ dimid_aero_particle /))
       call pmc_nc_write_real_1d(ncid, aero_refract_shell_imag, &
            "aero_refract_shell_imag", "1", (/ dimid_aero_particle /))
       call pmc_nc_write_real_1d(ncid, aero_refract_core_real, &
            "aero_refract_core_real", "1", (/ dimid_aero_particle /))
       call pmc_nc_write_real_1d(ncid, aero_refract_core_imag, &
            "aero_refract_core_imag", "1", (/ dimid_aero_particle /))
       call pmc_nc_write_real_1d(ncid, aero_core_vol, &
            "aero_core_vol", "m^3", (/ dimid_aero_particle /))
       call pmc_nc_write_integer_1d(ncid, aero_water_hyst_leg, &
            "aero_water_hyst_leg", "1", (/ dimid_aero_particle /))
       call pmc_nc_write_real_1d(ncid, aero_comp_vol, &
            "aero_comp_vol", "m^3", (/ dimid_aero_particle /))
       call pmc_nc_write_integer_1d(ncid, aero_id, &
            "aero_id", "1", (/ dimid_aero_particle /))
       call pmc_nc_write_real_1d(ncid, aero_least_create_time, &
            "aero_least_create_time", "s", (/ dimid_aero_particle /))
       call pmc_nc_write_real_1d(ncid, aero_greatest_create_time, &
            "aero_greatest_create_time", "s", (/ dimid_aero_particle /))
    end if
    
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
            "aero_removed_id", "1", (/ dimid_aero_removed /))
       call pmc_nc_write_integer_1d(ncid, aero_removed_action, &
            "aero_removed_action", "1", (/ dimid_aero_removed /))
       call pmc_nc_write_integer_1d(ncid, aero_removed_other_id, &
            "aero_removed_other_id", "1", (/ dimid_aero_removed /))
    end if

  end subroutine aero_state_output_netcdf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read full state.
  subroutine aero_state_input_netcdf(aero_state, ncid, bin_grid, &
       aero_data)
    
    !> aero_state to read.
    type(aero_state_t), intent(inout) :: aero_state
    !> NetCDF file ID, in data mode.
    integer, intent(in) :: ncid
    !> bin_grid structure.
    type(bin_grid_t), intent(in) :: bin_grid
    !> aero_data structure.
    type(aero_data_t), intent(in) :: aero_data

    integer :: dimid_aero_particle, dimid_aero_removed, n_info_item, n_part
    integer :: i_bin, i_part_in_bin, i_part, i_remove, status
    type(aero_particle_t) :: aero_particle
    character(len=1000) :: unit, name

    real*8, allocatable :: aero_comp_mass(:,:)
    integer, allocatable :: aero_n_orig_part(:)
    real*8, allocatable :: aero_absorb_cross_sect(:)
    real*8, allocatable :: aero_scatter_cross_sect(:)
    real*8, allocatable :: aero_asymmetry(:)
    real*8, allocatable :: aero_refract_shell_real(:)
    real*8, allocatable :: aero_refract_shell_imag(:)
    real*8, allocatable :: aero_refract_core_real(:)
    real*8, allocatable :: aero_refract_core_imag(:)
    real*8, allocatable :: aero_core_vol(:)
    integer, allocatable :: aero_water_hyst_leg(:)
    real*8, allocatable :: aero_comp_vol(:)
    integer, allocatable :: aero_id(:)
    real*8, allocatable :: aero_least_create_time(:)
    real*8, allocatable :: aero_greatest_create_time(:)
    integer, allocatable :: aero_removed_id(:)
    integer, allocatable :: aero_removed_action(:)
    integer, allocatable :: aero_removed_other_id(:)

    call pmc_nc_check(nf90_inq_dimid(ncid, "aero_particle", dimid_aero_particle))
    call pmc_nc_check(nf90_Inquire_Dimension(ncid, dimid_aero_particle, name, n_part))

    allocate(aero_comp_mass(n_part, aero_data%n_spec))
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

    call pmc_nc_read_real_2d(ncid, aero_comp_mass, &
         "aero_comp_mass", unit)
    call pmc_nc_read_integer_1d(ncid, aero_n_orig_part, &
         "aero_n_orig_part", unit)
    call pmc_nc_read_real_1d(ncid, aero_absorb_cross_sect, &
         "aero_absorb_cross_sect", unit)
    call pmc_nc_read_real_1d(ncid, aero_scatter_cross_sect, &
         "aero_scatter_cross_sect", unit)
    call pmc_nc_read_real_1d(ncid, aero_asymmetry, &
         "aero_asymmetry", unit)
    call pmc_nc_read_real_1d(ncid, aero_refract_shell_real, &
         "aero_refract_shell_real", unit)
    call pmc_nc_read_real_1d(ncid, aero_refract_shell_imag, &
         "aero_refract_shell_imag", unit)
    call pmc_nc_read_real_1d(ncid, aero_refract_core_real, &
         "aero_refract_core_real", unit)
    call pmc_nc_read_real_1d(ncid, aero_refract_core_imag, &
         "aero_refract_core_imag", unit)
    call pmc_nc_read_real_1d(ncid, aero_core_vol, &
         "aero_core_vol", unit)
    call pmc_nc_read_integer_1d(ncid, aero_water_hyst_leg, &
         "aero_water_hyst_leg", unit)
    call pmc_nc_read_real_1d(ncid, aero_comp_vol, &
         "aero_comp_vol", unit)
    call pmc_nc_read_integer_1d(ncid, aero_id, &
         "aero_id", unit)
    call pmc_nc_read_real_1d(ncid, aero_least_create_time, &
         "aero_least_create_time", unit)
    call pmc_nc_read_real_1d(ncid, aero_greatest_create_time, &
         "aero_greatest_create_time", unit)

    call aero_state_free(aero_state)
    call aero_state_alloc(bin_grid%n_bin, aero_data%n_spec, aero_state)

    call aero_particle_alloc(aero_particle, aero_data%n_spec)
    do i_part = 1,n_part
       aero_particle%vol = aero_comp_mass(i_part, :) / aero_data%density
       aero_particle%n_orig_part = aero_n_orig_part(i_part)
       aero_particle%absorb_cross_sect = aero_absorb_cross_sect(i_part)
       aero_particle%scatter_cross_sect = aero_scatter_cross_sect(i_part)
       aero_particle%asymmetry = aero_asymmetry(i_part)
       aero_particle%refract_shell = complex(aero_refract_shell_real(i_part), &
            aero_refract_shell_imag(i_part))
       aero_particle%refract_core = complex(aero_refract_core_real(i_part), &
            aero_refract_core_imag(i_part))
       aero_particle%core_vol = aero_core_vol(i_part)
       aero_particle%water_hyst_leg = aero_water_hyst_leg(i_part)
       aero_state%comp_vol = aero_comp_vol(i_part)
       aero_particle%id = aero_id(i_part)
       aero_particle%least_create_time = aero_least_create_time(i_part)
       aero_particle%greatest_create_time = aero_greatest_create_time(i_part)

       i_bin = aero_particle_in_bin(aero_particle, bin_grid)
       call aero_state_add_particle(aero_state, i_bin, aero_particle)
    end do
    call aero_particle_free(aero_particle)

    deallocate(aero_comp_mass)
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
       call pmc_nc_check(nf90_Inquire_Dimension(ncid, dimid_aero_removed, name, n_info_item))

       allocate(aero_removed_id(max(n_info_item,1)))
       allocate(aero_removed_action(max(n_info_item,1)))
       allocate(aero_removed_other_id(max(n_info_item,1)))

       call pmc_nc_read_integer_1d(ncid, aero_removed_id, &
            "aero_removed_id", unit)
       call pmc_nc_read_integer_1d(ncid, aero_removed_action, &
            "aero_removed_action", unit)
       call pmc_nc_read_integer_1d(ncid, aero_removed_other_id, &
            "aero_removed_other_id", unit)

       if ((n_info_item > 1) .or. (aero_removed_id(1) /= 0)) then
          call aero_info_array_enlarge_to(aero_state%aero_info_array, n_info_item)
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
end module pmc_aero_state
