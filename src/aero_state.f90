! Copyright (C) 2005-2007 Nicole Riemer and Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.
!
! The array aero_state%v is the main storage of the particle sizes and
! compositions, together with its sizing array aero_state%n. The
! particles in aero_state%v are stored sorted per-bin, to improve
! efficiency of sampling. If a particle has total volume pv then
! calling particle_in_bin(pv, n_bin, v_bin, i_bin) finds the bin
! number i_bin that that particle should go in. That particle is then
! stored as aero_state%v(i_bin)%p(i_part,:), where i_part is the index
! within the bin. aero_state%v(i_bin)%p(i_part,i_spec) is the volume
! of the i_spec-th species in the i_part-th particle in the i_bin-th
! bin.
!
! Typically most of the bins have only a few particles, while a small
! number of bins have many particles. To avoid having too much storage
! allocated for the bins with only a few particles, we do dynamic
! allocation/deallocation of the storage per-bin.
!
! With Fortran 90 we can't have arrays of arrays, so we have to use an
! array of pointers, and then allocate each pointer. We really want a
! 3D structure, with indices (i_bin, i_part, i_spec) specifiying
! species i_spec in particle number i_part stored in bin i_bin. This
! is stored as an array of pointers, one per bin, pointing to 2D
! arrays for which each row is a single particle (with the columns
! giving the volumes of the individual species).
!
! To avoid doing allocation and deallocation every time we add or
! remove a particle to a bin, we always double or halve the bin
! storage as necessary. The actual number of particles stored in a bin
! will generally be less than the actual memory allocated for that
! bin, so we store the current number of particles in a bin in the
! array aero_state%n. The allocated size of bin storage in
! aero_state%v(i_bin) is not stored explicitly, but can be obtained
! with the Fortran 90 SIZE() intrinsic function.

module pmc_aero_state

  use pmc_aero_particle_array
  use pmc_bin_grid
  use pmc_aero_data
  use pmc_aero_particle
  use pmc_aero_dist
  use pmc_util
  use pmc_rand_poisson
  use pmc_aero_binned
  use pmc_mpi
  use pmc_inout
#ifdef PMC_USE_MPI
#ifndef PMC_EVEREST
  use mpi
#endif
#endif

  type aero_state_t
     type(aero_particle_array_t), pointer :: bins(:) ! bin arrays
     real*8 :: comp_vol                 ! computational volume (m^3)
     integer :: n_part                  ! total number of particles
  end type aero_state_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine aero_state_alloc(n_bin, n_spec, aero_state)

    ! Initializes aerosol arrays to have zero particles in each
    ! bin. Do not call this more than once on a given aerosol, use
    ! aero_state_zero() instead to reset to zero.

    integer, intent(in) :: n_bin        ! number of bins
    integer, intent(in) :: n_spec       ! number of species
    type(aero_state_t), intent(inout) :: aero_state ! aerosol to initialize
    
    integer i

    allocate(aero_state%bins(n_bin))
    do i = 1,n_bin
       call aero_particle_array_alloc(aero_state%bins(i), 0, n_spec)
    end do
    aero_state%comp_vol = 0d0
    aero_state%n_part = 0

  end subroutine aero_state_alloc
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine aero_state_free(aero_state)

    ! Deallocates a previously allocated aerosol.

    type(aero_state_t), intent(inout) :: aero_state ! aerosol to initialize
    
    integer :: n_bin, i

    n_bin = size(aero_state%bins)
    do i = 1,n_bin
       call aero_particle_array_free(aero_state%bins(i))
    end do
    deallocate(aero_state%bins)

  end subroutine aero_state_free
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine aero_state_copy(aero_state_from, aero_state_to)

    ! Copies aerosol to a destination that has already had
    ! aero_state_alloc() called on it.

    type(aero_state_t), intent(in) :: aero_state_from ! reference aerosol
    type(aero_state_t), intent(inout) :: aero_state_to ! already allocated
    
    integer :: n_bin, i

    n_bin = size(aero_state_from%bins)

    call aero_state_free(aero_state_to)
    call aero_state_alloc(n_bin, 0, aero_state_to)

    do i = 1,n_bin
       call aero_particle_array_copy(aero_state_from%bins(i), &
            aero_state_to%bins(i))
    end do

    aero_state_to%comp_vol = aero_state_from%comp_vol
    aero_state_to%n_part = aero_state_from%n_part

  end subroutine aero_state_copy
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function aero_state_total_particles(aero_state)

    ! Returns the total number of particles in an aerosol distribution.

    type(aero_state_t), intent(in) :: aero_state ! aerosol state

    aero_state_total_particles = aero_state%n_part

  end function aero_state_total_particles

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine aero_state_zero(aero_state)

    ! Resets an aero_state to have zero particles per bin. This must
    ! already have had aero_state_alloc() called on it. This
    ! function can be called more than once on the same state.

    type(aero_state_t), intent(inout) :: aero_state ! state to zero
    
    integer :: i, n_bin

    n_bin = size(aero_state%bins)
    do i = 1,n_bin
       call aero_particle_array_zero(aero_state%bins(i))
    end do
    aero_state%n_part = 0

  end subroutine aero_state_zero
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine aero_state_add_particle(aero_state, bin, aero_particle)

    ! Add the given particle.

    type(aero_state_t), intent(inout) :: aero_state ! aerosol state
    integer, intent(in) :: bin          ! bin number of particle to remove
    type(aero_particle_t), intent(in) :: aero_particle ! particle to add

    call aero_particle_array_add_particle(aero_state%bins(bin), &
         aero_particle)
    aero_state%n_part = aero_state%n_part + 1

  end subroutine aero_state_add_particle

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine aero_state_remove_particle(aero_state, bin, index)

    ! Remove the given particle.

    type(aero_state_t), intent(inout) :: aero_state ! aerosol state
    integer, intent(in) :: bin          ! bin number of particle to remove
    integer, intent(in) :: index        ! index in bin of particle to remove

    call aero_particle_array_remove_particle(aero_state%bins(bin), index)
    aero_state%n_part = aero_state%n_part - 1

  end subroutine aero_state_remove_particle

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine aero_state_add(aero_state, aero_state_delta)

    ! aero_state += aero_state_delta

    type(aero_state_t), intent(inout) :: aero_state ! aerosol state
    type(aero_state_t), intent(in) :: aero_state_delta ! increment

    integer :: i_bin, i_part

    do i_bin = 1,size(aero_state_delta%bins)
       do i_part = 1,aero_state_delta%bins(i_bin)%n_part
          call aero_state_add_particle(aero_state, i_bin, &
               aero_state_delta%bins(i_bin)%particle(i_part))
       end do
    end do
    aero_state%comp_vol = aero_state%comp_vol + aero_state_delta%comp_vol

  end subroutine aero_state_add

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine aero_state_add_disc_mode(bin_grid, aero_data, vol_frac, &
       bin_n, aero_state)

    ! Makes particles from the given number distribution and appends
    ! them to the aero_state%v array.
    
    type(bin_grid_t), intent(in) :: bin_grid ! bin grid
    type(aero_data_t), intent(in) :: aero_data ! aero_data data
    real*8, intent(in) :: vol_frac(aero_data%n_spec) ! composition of particles
    integer, intent(in) :: bin_n(bin_grid%n_bin) ! number in bins
    type(aero_state_t), intent(inout) :: aero_state ! aerosol, must be
                                              ! allocated already
    
    real*8 v_low, v_high, pv
    integer k, i
    type(aero_particle_t) :: aero_particle

    call aero_particle_alloc(aero_particle, aero_data%n_spec)
    do k = 1,bin_grid%n_bin
       v_low = bin_edge(bin_grid, k)
       v_high = bin_edge(bin_grid, k + 1)
       do i = 1,bin_n(k)
          ! we used to do:
          ! pv = dble(i) / dble(bin_n(k) + 1) * (v_high - v_low) + v_low
          ! but this doesn't actually work as well as:
          pv = bin_grid%v(k)
          call aero_particle_set_vols(aero_particle, vol_frac * pv)
          call aero_state_add_particle(aero_state, k, aero_particle)
       end do
    end do
    call aero_particle_free(aero_particle)

  end subroutine aero_state_add_disc_mode
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine aero_dist_to_state(bin_grid, aero_data, aero_dist, &
       n_part, aero_state)

    ! Convert a continuous distribution into particles.
    
    type(bin_grid_t), intent(in) :: bin_grid ! bin grid
    type(aero_data_t), intent(in) :: aero_data ! aero_data data
    type(aero_dist_t), intent(in) :: aero_dist ! aerosol distribution
    integer, intent(in) :: n_part       ! total number of particles
    type(aero_state_t), intent(out) :: aero_state ! aerosol distribution,
                                                  ! will be alloced

    integer :: i
    real*8 :: total_n_den
    real*8 :: mode_n_dens(aero_dist%n_mode)
    integer :: mode_n_parts(aero_dist%n_mode)
    integer :: num_per_bin(bin_grid%n_bin)

    ! find the total number density of each mode
    total_n_den = 0d0
    do i = 1,aero_dist%n_mode
       mode_n_dens(i) = sum(aero_dist%mode(i)%num_den)
    end do
    total_n_den = sum(mode_n_dens)

    ! allocate particles to modes proportional to their number densities
    call vec_cts_to_disc(aero_dist%n_mode, mode_n_dens, n_part, mode_n_parts)

    ! allocate particles within each mode in proportion to mode shape
    call aero_state_alloc(bin_grid%n_bin, aero_data%n_spec, aero_state)
    do i = 1,aero_dist%n_mode
       call vec_cts_to_disc(bin_grid%n_bin, aero_dist%mode(i)%num_den, &
            mode_n_parts(i), num_per_bin)
       call aero_state_add_disc_mode(bin_grid, aero_data, &
            aero_dist%mode(i)%vol_frac, num_per_bin, aero_state)
    end do

    aero_state%comp_vol = dble(n_part) &
         / aero_dist_total_num_den(bin_grid, aero_dist)

  end subroutine aero_dist_to_state
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine aero_dist_sample(bin_grid, aero_data, aero_dist, &
       sample_vol, aero_state)

    ! Generates a Poisson sample of an aero_dist, adding to aero_state.

    ! FIXME: we should not take sample_vol, but rather sample_amount
    ! that we then multiply by aero_state%comp_vol

    type(bin_grid_t), intent(in) :: bin_grid ! bin grid
    type(aero_data_t), intent(in) :: aero_data ! aero data values
    type(aero_dist_t), intent(in) :: aero_dist ! distribution to sample
    real*8, intent(in) :: sample_vol    ! volume to sample
    type(aero_state_t), intent(inout) :: aero_state ! aero state to add to

    real*8 :: n_samp_avg
    integer :: n_samp, i_mode
    integer :: num_per_bin(bin_grid%n_bin)

    do i_mode = 1,aero_dist%n_mode
       n_samp_avg = sample_vol * sum(aero_dist%mode(i_mode)%num_den) &
            * bin_grid%dlnr
       n_samp = rand_poisson(n_samp_avg)
       call sample_vec_cts_to_disc(bin_grid%n_bin, &
            aero_dist%mode(i_mode)%num_den, n_samp, num_per_bin)
       call aero_state_add_disc_mode(bin_grid, aero_data, &
            aero_dist%mode(i_mode)%vol_frac, num_per_bin, aero_state)
    end do

  end subroutine aero_dist_sample
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine aero_state_rand_particle(aero_state, i_bin, i_part)

    ! Choose a random particle from the aero_state.

    type(aero_state_t), intent(in) :: aero_state ! original state
    integer, intent(out) :: i_bin       ! bin number of particle
    integer, intent(out) :: i_part      ! particle number within bin

    integer :: n_bin, disc_pdf(size(aero_state%bins))

    call assert(950725003, aero_state%n_part > 0)
    n_bin = size(aero_state%bins)
    disc_pdf = (/(aero_state%bins(i_bin)%n_part, i_bin = 1,n_bin)/)
    i_bin = sample_disc_pdf(n_bin, disc_pdf)
    i_part = util_rand_int(aero_state%bins(i_bin)%n_part)

  end subroutine aero_state_rand_particle

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine aero_state_sample(aero_state_from, aero_state_to, &
       sample_prop)
    
    ! Generates a Poisson sample by removing particles from
    ! aero_state_from and adding them to aero_state_to, which must
    ! be already allocated (and should have its comp_vol set).

    type(aero_state_t), intent(inout) :: aero_state_from ! original state
    type(aero_state_t), intent(inout) :: aero_state_to ! destination state
    real*8, intent(in) :: sample_prop   ! proportion to sample
    
    integer :: n_transfer, i_transfer, n_bin, i_bin, i_part
    logical :: do_add, do_remove
    real*8 :: vol_ratio

    call assert(721006962, (sample_prop >= 0d0) .and. (sample_prop <= 1d0))
    n_transfer = rand_poisson(sample_prop &
         * dble(aero_state_total_particles(aero_state_from)))
    n_bin = size(aero_state_from%bins)
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
          if (util_rand() < 1d0 / vol_ratio) then
             do_remove = .true.
          end if
       else ! vol_ratio < 1d0
          ! to_comp_vol is smaller than from_comp_vol, so always
          ! remove the particle but only maybe add it
          do_add = .false.
          if (util_rand() < vol_ratio) then
             do_add = .true.
          end if
          do_remove = .true.
       end if
       if (do_add) then
          call aero_state_add_particle(aero_state_to, i_bin, &
               aero_state_from%bins(i_bin)%particle(i_part))
       end if
       if (do_remove) then
          call aero_state_remove_particle(aero_state_from, i_bin, i_part)
          i_transfer = i_transfer + 1
       end if
    end do
    
  end subroutine aero_state_sample
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine aero_state_sample_rough(aero_state_from, aero_state_to, &
       sample_prop)
    
    ! Generates a rough bin-wise sample by removing particles from
    ! aero_state_from and adding them to aero_state_to, which must
    ! be already allocated (and should have its comp_vol set).

    type(aero_state_t), intent(inout) :: aero_state_from ! original state
    type(aero_state_t), intent(inout) :: aero_state_to ! destination state
    real*8, intent(in) :: sample_prop   ! proportion to sample
    
    integer :: n_transfer, i_transfer, n_bin, i_bin, i_part
    logical :: do_add, do_remove
    real*8 :: vol_ratio

    n_bin = size(aero_state_from%bins)
    vol_ratio = aero_state_to%comp_vol / aero_state_from%comp_vol
    do i_bin = 1,n_bin
       n_transfer = prob_round(sample_prop &
            * dble(aero_state_from%bins(i_bin)%n_part))
       i_transfer = 0
       do while (i_transfer < n_transfer)
          if (aero_state_from%bins(i_bin)%n_part <= 0) exit
          i_part = util_rand_int(aero_state_from%bins(i_bin)%n_part)
          if (vol_ratio == 1d0) then
             ! to_comp_vol == from_comp_vol so just move the particle
             do_add = .true.
             do_remove = .true.
          elseif (vol_ratio > 1d0) then
             ! to_comp_vol is bigger than from_comp_vol, so only maybe
             ! remove the particle but always add it
             do_add = .true.
             do_remove = .false.
             if (util_rand() < 1d0 / vol_ratio) then
                do_remove = .true.
             end if
          else ! vol_ratio < 1d0
             ! to_comp_vol is smaller than from_comp_vol, so always
             ! remove the particle but only maybe add it
             do_add = .false.
             if (util_rand() < vol_ratio) then
                do_add = .true.
             end if
             do_remove = .true.
          end if
          if (do_add) then
             call aero_state_add_particle(aero_state_to, i_bin, &
                  aero_state_from%bins(i_bin)%particle(i_part))
          end if
          if (do_remove) then
             call aero_state_remove_particle(aero_state_from, i_bin, i_part)
             i_transfer = i_transfer + 1
          end if
       end do
    end do
    
  end subroutine aero_state_sample_rough
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine aero_state_add_particles(aero_state, aero_state_delta)
    
    ! Adds aero_state_delta particles to aero_state. The number of
    ! particles added depends on the computational volume ratio, so
    ! either more or less particles might be added than are actually
    ! in aero_state_delta.

    type(aero_state_t), intent(inout) :: aero_state ! base state
    type(aero_state_t), intent(in) :: aero_state_delta ! state to add
    
    integer :: n_bin, i_bin, i_part, n_new_part, i_new_part
    real*8 :: vol_ratio
    
    n_bin = size(aero_state%bins)
    vol_ratio = aero_state%comp_vol / aero_state_delta%comp_vol
    do i_bin = 1,n_bin
       do i_part = 1,aero_state_delta%bins(i_bin)%n_part
          n_new_part = prob_round(vol_ratio)
          do i_new_part = 1,n_new_part
             call aero_state_add_particle(aero_state, i_bin, &
                  aero_state_delta%bins(i_bin)%particle(i_part))
          end do
       end do
    end do
    
  end subroutine aero_state_add_particles
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine aero_state_to_binned(bin_grid, aero_data, aero_state, &
       aero_binned)
    
    ! Create the bin number and mass arrays from aero_state%v.
    
    type(bin_grid_t), intent(in) :: bin_grid ! bin grid
    type(aero_data_t), intent(in) :: aero_data ! aerosol data
    type(aero_state_t), intent(in) :: aero_state ! aerosol state
    type(aero_binned_t), intent(out) :: aero_binned ! binned distributions
    
    integer b, j, s
    
    aero_binned%num_den = 0d0
    aero_binned%vol_den = 0d0
    do b = 1,bin_grid%n_bin
       do j = 1,aero_state%bins(b)%n_part
          aero_binned%vol_den(b,:) = aero_binned%vol_den(b,:) &
               + aero_state%bins(b)%particle(j)%vol / aero_state%comp_vol &
               / bin_grid%dlnr
       end do
       aero_binned%num_den(b) = dble(aero_state%bins(b)%n_part) &
            / aero_state%comp_vol / bin_grid%dlnr
    end do
    
  end subroutine aero_state_to_binned
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine aero_state_double(aero_state)
    
    ! Doubles number of particles.
    
    type(aero_state_t), intent(inout) :: aero_state ! aerosol state
    
    integer :: i, n_bin
    
    n_bin = size(aero_state%bins)
    do i = 1,n_bin
       call aero_particle_array_double(aero_state%bins(i))
    end do
    aero_state%comp_vol = 2d0 * aero_state%comp_vol
    aero_state%n_part = 2 * aero_state%n_part

  end subroutine aero_state_double
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine aero_state_halve(aero_state, aero_binned, bin_grid)
    
    ! Remove approximately half of the particles in each bin.
    
    type(aero_state_t), intent(inout) :: aero_state ! aerosol state
    type(aero_binned_t), intent(inout) :: aero_binned ! aero binned
    type(bin_grid_t), intent(in) :: bin_grid ! bin grid
    
    integer :: i_bin, i_part, n_part_orig, i_remove, n_remove

    n_part_orig = aero_state%n_part
    do i_bin = 1,bin_grid%n_bin
       n_remove = prob_round(dble(aero_state%bins(i_bin)%n_part) / 2d0)
       do i_remove = 1,n_remove
          i_part = util_rand_int(aero_state%bins(i_bin)%n_part)
          call aero_binned_remove_particle_in_bin(aero_binned, bin_grid, &
               i_bin, aero_state%comp_vol, &
               aero_state%bins(i_bin)%particle(i_part))
          call aero_state_remove_particle(aero_state, i_bin, i_part)
       end do
    end do
    aero_state%comp_vol = aero_state%comp_vol &
         * dble(aero_state%n_part) / dble(n_part_orig)

  end subroutine aero_state_halve
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine aero_state_resort(bin_grid, aero_state)
    
    ! Takes a VH array where the particle volumes might no longer be
    ! correct for the bins they are in and resorts it so that every
    ! particle is in the correct bin.
    
    type(bin_grid_t), intent(in) :: bin_grid ! bin_grid
    type(aero_state_t), intent(inout) :: aero_state ! aerosol state
    
    integer :: bin, j, new_bin, k
    
    ! The approach here is inefficient because we might reprocess
    ! particles. For example, if we are doing bin 1 and we shift a
    ! particle up to bin 2, when we do bin 2 we will reprocess it. It
    ! seems to be more trouble than it's worth to worry about this
    ! yet, however.
    
    do bin = 1,bin_grid%n_bin
       j = 1
       do while (j .le. aero_state%bins(bin)%n_part)
          ! find the new bin
          new_bin = aero_particle_in_bin(aero_state%bins(bin)%particle(j), &
               bin_grid)
          
          ! if the bin number has changed, move the particle
          if (bin .ne. new_bin) then
             call aero_state_add_particle(aero_state, new_bin, &
                  aero_state%bins(bin)%particle(j))
             call aero_state_remove_particle(aero_state, bin, j)

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
    do bin = 1,bin_grid%n_bin
       call aero_particle_array_shrink(aero_state%bins(bin))
    end do
    
  end subroutine aero_state_resort
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine aero_state_mix_to(aero_state, mix_rate, dest, &
       aero_binned, aero_data, bin_grid)

    ! Send a sample to the given process, and receive exactly one
    ! sample from an unspecified source.

    type(aero_state_t), intent(inout) :: aero_state ! aerosol state
    real*8, intent(in) :: mix_rate      ! mixing rate (0 to 1)
    integer, intent(in) :: dest         ! process to send to
    type(aero_binned_t), intent(inout) :: aero_binned ! aero binned to update
    type(aero_data_t), intent(in) :: aero_data ! aero data values
    type(bin_grid_t), intent(in) :: bin_grid ! bin grid

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
    call aero_binned_alloc(aero_binned_delta, bin_grid%n_bin, aero_data%n_spec)

    ! extract particles to send
    call aero_state_alloc(bin_grid%n_bin, aero_data%n_spec, aero_state_send)
    aero_state_send%comp_vol = aero_state%comp_vol
    ! FIXME: would probably be slightly better for sampling purposes
    ! to use the destination comp_vol here
    call aero_state_sample_rough(aero_state, aero_state_send, mix_rate)
    call aero_state_to_binned(bin_grid, aero_data, aero_state_send, &
         aero_binned_delta)
    call aero_binned_sub(aero_binned, aero_binned_delta)

    ! pack up data to send
    buffer_size_send = pmc_mpi_pack_size_aero_state(aero_state_send)
    allocate(buffer_send(buffer_size_send))
    position = 0
    call pmc_mpi_pack_aero_state(buffer_send, position, aero_state_send)
    call assert(420102502, position == buffer_size_send)

    ! transfer buffer size info and allocate receive buffers
    call mpi_sendrecv(buffer_size_send, 1, MPI_INTEGER, dest, 0, &
         buffer_size_recv, 1, MPI_INTEGER, MPI_ANY_SOURCE, &
         MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
    call pmc_mpi_check_ierr(ierr)
    allocate(buffer_recv(buffer_size_recv))

    ! do data transfer
    call mpi_sendrecv(buffer_send, buffer_size_send, MPI_PACKED, dest, 0, &
         buffer_recv, buffer_size_recv, MPI_PACKED, MPI_ANY_SOURCE, &
         MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
    call pmc_mpi_check_ierr(ierr)

    ! unpack received data and process it
    position = 0
    call pmc_mpi_unpack_aero_state(buffer_recv, position, aero_state_recv)
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

  subroutine aero_state_mix(aero_state, mix_rate, &
       aero_binned, aero_data, bin_grid)

    ! Mix the aero_states between all processes. Currently uses a
    ! simple periodic-1D diffusion.

    type(aero_state_t), intent(inout) :: aero_state ! aerosol state
    real*8, intent(in) :: mix_rate      ! mixing rate (0 to 1)
    type(aero_binned_t), intent(inout) :: aero_binned ! aero binned to update
    type(aero_data_t), intent(in) :: aero_data ! aero data values
    type(bin_grid_t), intent(in) :: bin_grid ! bin grid

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
  
  subroutine aero_state_check(bin_grid, aero_binned, aero_data, aero_state)
    
    ! Check that all particles are in the correct bins and that the
    ! bin numbers and masses are correct. This is for debugging only.
    
    type(bin_grid_t), intent(in) :: bin_grid ! bin_grid
    type(aero_binned_t), intent(out) :: aero_binned ! binned distributions
    type(aero_data_t), intent(in) :: aero_data ! aerosol data
    type(aero_state_t), intent(inout) :: aero_state ! aerosol state
    
    real*8 :: check_bin_v, check_vol_den, vol_tol
    real*8 :: num_tol, state_num_den
    integer :: i, k, k_check, s, n_part_check
    logical :: error
    
    error = .false.

    ! check that the total number of particles is correct
    n_part_check = sum((/(aero_state%bins(i)%n_part, i = 1,bin_grid%n_bin)/))
    if (aero_state%n_part /= n_part_check) then
       write(0,'(a20,a20)') 'aero_state%n_part', 'n_part_check'
       write(0,'(i20,i20)') aero_state%n_part, n_part_check
       error = .true.
    end if
    
    ! check that all particles are in the correct bins
    do k = 1,bin_grid%n_bin
       do i = 1,aero_state%bins(k)%n_part
          k_check = aero_particle_in_bin(aero_state%bins(k)%particle(i), &
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
       state_num_den = dble(aero_state%bins(k)%n_part) / aero_state%comp_vol &
            / bin_grid%dlnr
       if (.not. almost_equal_abs(state_num_den, &
            aero_binned%num_den(k), num_tol)) then
          write(0,'(a10,a20,a20,a20,a20)') 'k', 'bins(k)%n_part', &
               'state_num_den', 'num_den(k)', 'comp_vol'
          write(0,'(i10,i20,e20.10,e20.10,e20.10)') k, &
               aero_state%bins(k)%n_part, state_num_den, &
               aero_binned%num_den(k), aero_state%comp_vol
          error = .true.
       end if
    end do
    
    ! check the aero_binned%vol_den array
    do k = 1,bin_grid%n_bin
       vol_tol = bin_grid%v(k) / 1d3 / bin_grid%dlnr
       do s = 1,aero_data%n_spec
          check_vol_den = sum((/(aero_state%bins(k)%particle(i)%vol(s), &
               i = 1,aero_state%bins(k)%n_part)/)) &
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
    
    if (error) then
       write(0,*) 'ERROR: aero_state_check() failed'
       call exit(2)
    end if
    
  end subroutine aero_state_check
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine inout_write_aero_state(file, aero_state)
    
    ! Write full state.
    
    type(inout_file_t), intent(inout) :: file ! file to write to
    type(aero_state_t), intent(in) :: aero_state ! aero_state to write

    integer :: n_bin, i
    
    n_bin = size(aero_state%bins)
    call inout_write_comment(file, "begin aero_state")
    call inout_write_real(file, "comp_vol(m^3)", aero_state%comp_vol)
    call inout_write_integer(file, "n_part", aero_state%n_part)
    call inout_write_integer(file, "n_bin", n_bin)
    do i = 1,n_bin
       call inout_write_integer(file, "bin_number", i)
       call inout_write_aero_particle_array(file, aero_state%bins(i))
    end do
    call inout_write_comment(file, "end aero_state")
    
  end subroutine inout_write_aero_state
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine inout_read_aero_state(file, aero_state)
    
    ! Read full state.
    
    type(inout_file_t), intent(inout) :: file ! file to write to
    type(aero_state_t), intent(out) :: aero_state ! aero_state to read
    
    integer :: n_bin, i, check_i
    
    call inout_check_comment(file, "begin aero_state")
    call inout_read_real(file, "comp_vol(m^3)", aero_state%comp_vol)
    call inout_read_integer(file, "n_part", aero_state%n_part)
    call inout_read_integer(file, "n_bin", n_bin)
    allocate(aero_state%bins(n_bin))
    do i = 1,n_bin
       call inout_read_integer(file, "bin_number", check_i)
       call inout_check_index(file, i, check_i)
       call inout_read_aero_particle_array(file, aero_state%bins(i))
    end do
    call inout_check_comment(file, "end aero_state")
    
  end subroutine inout_read_aero_state
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function pmc_mpi_pack_size_aero_state(val)

    ! Determines the number of bytes required to pack the given value.

    type(aero_state_t), intent(in) :: val ! value to pack

    integer :: i, total_size

    total_size = 0
    total_size = total_size + pmc_mpi_pack_size_real(val%comp_vol)
    total_size = total_size + pmc_mpi_pack_size_integer(val%n_part)
    total_size = total_size + pmc_mpi_pack_size_integer(size(val%bins))
    do i = 1,size(val%bins)
       total_size = total_size + pmc_mpi_pack_size_apa(val%bins(i))
    end do
    pmc_mpi_pack_size_aero_state = total_size

  end function pmc_mpi_pack_size_aero_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine pmc_mpi_pack_aero_state(buffer, position, val)

    ! Packs the given value into the buffer, advancing position.

    character, intent(inout) :: buffer(:) ! memory buffer
    integer, intent(inout) :: position  ! current buffer position
    type(aero_state_t), intent(in) :: val ! value to pack

#ifdef PMC_USE_MPI
    integer :: prev_position, i

    prev_position = position
    call pmc_mpi_pack_real(buffer, position, val%comp_vol)
    call pmc_mpi_pack_integer(buffer, position, val%n_part)
    call pmc_mpi_pack_integer(buffer, position, size(val%bins))
    do i = 1,size(val%bins)
       call pmc_mpi_pack_aero_particle_array(buffer, position, val%bins(i))
    end do
    call assert(850997402, position - prev_position == pmc_mpi_pack_size_aero_state(val))
#endif

  end subroutine pmc_mpi_pack_aero_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine pmc_mpi_unpack_aero_state(buffer, position, val)

    ! Unpacks the given value from the buffer, advancing position.

    character, intent(inout) :: buffer(:) ! memory buffer
    integer, intent(inout) :: position  ! current buffer position
    type(aero_state_t), intent(out) :: val ! value to pack

#ifdef PMC_USE_MPI
    integer :: prev_position, i, n

    prev_position = position
    call pmc_mpi_unpack_real(buffer, position, val%comp_vol)
    call pmc_mpi_unpack_integer(buffer, position, val%n_part)
    call pmc_mpi_unpack_integer(buffer, position, n)
    allocate(val%bins(n))
    do i = 1,size(val%bins)
       call pmc_mpi_unpack_aero_particle_array(buffer, position, val%bins(i))
    end do
    call assert(132104747, position - prev_position == pmc_mpi_pack_size_aero_state(val))
#endif

  end subroutine pmc_mpi_unpack_aero_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
end module pmc_aero_state
