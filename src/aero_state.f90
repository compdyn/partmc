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

module mod_aero_state

  type bin_p_t
     real*8, dimension(:,:), pointer :: p ! particle volumes (m^3)
     ! dimension of p is (# particles in bin) x n_spec
  end type bin_p_t

  type aero_state_t
     integer, dimension(:), pointer :: n ! number of particles in each bin
     type(bin_p_t), dimension(:), pointer :: v ! particle volumes (m^3)
     real*8 :: comp_vol                 ! computational volume (m^3)
  end type aero_state_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine alloc_aero_state(n_bin, n_spec, aero)

    ! Initializes aerosol arrays to have zero particles in each
    ! bin. Do not call this more than once on a given aerosol, use
    ! zero_aero_state() instead to reset to zero.

    integer, intent(in) :: n_bin        ! number of bins
    integer, intent(in) :: n_spec       ! number of species
    type(aero_state_t), intent(inout) :: aero ! aerosol to initialize
    
    integer i

    allocate(aero%n(n_bin))
    aero%n = 0

    allocate(aero%v(n_bin))
    do i = 1,n_bin
       allocate(aero%v(i)%p(0, n_spec))
    end do

  end subroutine alloc_aero_state
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine free_aero_state(aero)

    ! Deallocates a previously allocated aerosol.

    type(aero_state_t), intent(inout) :: aero ! aerosol to initialize
    
    integer :: n_bin, i

    n_bin = size(aero%n)
    deallocate(aero%n)
    do i = 1,n_bin
       deallocate(aero%v(i)%p)
    end do
    deallocate(aero%v)

  end subroutine free_aero_state
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine copy_aero_state(from_aero, to_aero)

    ! Copies aerosol to a destination that has already had
    ! alloc_aero_state() called on it.

    type(aero_state_t), intent(in) :: from_aero ! reference aerosol
    type(aero_state_t), intent(inout) :: to_aero ! must already be allocated
    
    integer :: n_bin, n_spec, n_part, i
    integer :: arr_shape(2)

    n_bin = size(from_aero%n)
    arr_shape = shape(from_aero%v(1)%p)
    n_spec = arr_shape(2)

    call free_aero_state(to_aero)
    call alloc_aero_state(n_bin, n_spec, to_aero)

    to_aero%n = from_aero%n
    do i = 1,n_bin
       arr_shape = shape(from_aero%v(i)%p)
       n_part = arr_shape(1)
       call enlarge_bin_to(to_aero%v(i), n_part)
       to_aero%v(i)%p = from_aero%v(i)%p(1:n_part,:)
    end do

    to_aero%comp_vol = from_aero%comp_vol

  end subroutine copy_aero_state
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function total_particles(aero_state)

    ! Returns the total number of particles in an aerosol distribution.

    type(aero_state_t), intent(in) :: aero_state ! aerosol state

    total_particles = sum(aero_state%n)

  end function total_particles

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine zero_aero_state(aero_state)

    ! Resets an aero_state to have zero particles per bin. This must
    ! already have had alloc_aero_state() called on it. This
    ! function can be called more than once on the same state.

    type(aero_state_t), intent(inout) :: aero_state ! state to zero
    
    integer :: i, n_bin, n_spec, p_shape(2)

    n_bin = size(aero_state%v)
    p_shape = shape(aero_state%v(1)%p)
    n_spec = p_shape(2)
    do i = 1,n_bin
       deallocate(aero_state%v(i)%p)
       allocate(aero_state%v(i)%p(0, n_spec))
    end do
    aero_state%n = 0

  end subroutine zero_aero_state
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine enlarge_bin(bin)

    ! Enlarges the given bin (which must be allocated) by at least one
    ! element (currently doubles the length).

    type(bin_p_t), intent(inout) :: bin   ! bin data

    integer :: n_part, n_spec, new_n_part
    real*8, dimension(:,:), pointer :: new_p

    ! FIXME: should use SHAPE instead of SIZE here?
    n_part = size(bin%p, 1)
    n_spec = size(bin%p, 2)
    new_n_part = max(n_part * 2, n_part + 1)
    allocate(new_p(new_n_part, n_spec))
    new_p(1:n_part,:) = bin%p
    deallocate(bin%p)
    bin%p => new_p
    
  end subroutine enlarge_bin

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine enlarge_bin_to(bin, n)

    ! Enlarges the given bin so that it is at least of size n.

    type(bin_p_t), intent(inout) :: bin   ! bin data
    integer, intent(in) :: n            ! minimum new size of bin

    do while (size(bin%p,1) < n)
       call enlarge_bin(bin)
    end do

  end subroutine enlarge_bin_to

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine shrink_bin(n_used, bin)

    ! Possibly shrinks the storage of the given bin, ensuring that it
    ! is at least of length n_used.

    integer, intent(in) :: n_used       ! number of used entries in bin
    type(bin_p_t), intent(inout) :: bin   ! bin data

    integer :: n_part, n_spec, new_n_part
    real*8, dimension(:,:), pointer :: new_p

    ! FIXME: should use SHAPE instead of SIZE here?
    n_part = size(bin%p, 1)
    n_spec = size(bin%p, 2)
    new_n_part = n_part / 2
    do while (n_used <= new_n_part)
       allocate(new_p(new_n_part, n_spec))
       new_p(:,:) = bin%p(1:new_n_part,:)
       deallocate(bin%p)
       bin%p => new_p
       n_part = new_n_part
       new_n_part = n_part / 2
       ! FIXME: gfortran 4.1.1 requires the "then" in the following
       ! statement, rather than using a single-line "if" statement.
       if (new_n_part == 0) then
          exit
       end if
    end do

  end subroutine shrink_bin

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine add_particles(bin_grid, aero_data, vol_frac, bin_n, aero)

    ! Makes particles from the given number distribution and appends
    ! them to the aero_state%v array.
    
    use mod_bin_grid
    use mod_aero_data

    type(bin_grid_t), intent(in) :: bin_grid ! bin grid
    type(aero_data_t), intent(in) :: aero_data   ! aero_data data
    real*8, intent(in) :: vol_frac(aero_data%n_spec) ! composition of particles
    integer, intent(in) :: bin_n(bin_grid%n_bin) ! number in bins
    type(aero_state_t), intent(inout) :: aero ! aerosol, must be
                                              ! allocated already
    
    real*8 total_vol_frac, v_low, v_high, pv
    integer k, i

    total_vol_frac = sum(vol_frac)
    do k = 1,bin_grid%n_bin
       call bin_edge(bin_grid, k, v_low)
       call bin_edge(bin_grid, k + 1, v_high)
       do i = 1,bin_n(k)
          ! we used to do:
          ! pv = dble(i) / dble(bin_n(k) + 1) * (v_high - v_low) + v_low
          ! but this doesn't actually work as well as:
          pv = bin_grid%v(k)
          aero%n(k) = aero%n(k) + 1
          call enlarge_bin_to(aero%v(k), aero%n(k))
          aero%v(k)%p(aero%n(k),:) = vol_frac / total_vol_frac * pv
       end do
    end do

  end subroutine add_particles
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine aero_dist_to_part(bin_grid, aero_data, aero_dist, &
       n_part, aero_state)

    ! Convert a continuous distribution into particles.
    
    use mod_aero_data
    use mod_bin_grid
    use mod_aero_dist
    use mod_util

    type(bin_grid_t), intent(in) :: bin_grid ! bin grid
    type(aero_data_t), intent(in) :: aero_data ! aero_data data
    type(aero_dist_t), intent(in) :: aero_dist ! aerosol distribution
    integer, intent(in) :: n_part       ! total number of particles
    type(aero_state_t), intent(out) :: aero_state ! aerosol distribution,
                                                  ! will be alloced

    integer :: i
    real*8 :: total_n_den
    real*8 :: mode_n_dens(aero_dist%n_modes)
    integer :: mode_n_parts(aero_dist%n_modes)
    integer :: num_per_bin(bin_grid%n_bin)

    ! find the total number density of each mode
    total_n_den = 0d0
    do i = 1,aero_dist%n_modes
       mode_n_dens(i) = sum(aero_dist%modes(i)%n_den)
    end do
    total_n_den = sum(mode_n_dens)

    ! allocate particles to modes proportional to their number densities
    call vec_cts_to_disc(aero_dist%n_modes, mode_n_dens, n_part, mode_n_parts)

    ! allocate particles within each mode in proportion to mode shape
    call alloc_aero_state(bin_grid%n_bin, aero_data%n_spec, aero_state)
    do i = 1,aero_dist%n_modes
       call vec_cts_to_disc(bin_grid%n_bin, aero_dist%modes(i)%n_den, &
            mode_n_parts(i), num_per_bin)
       call add_particles(bin_grid, aero_data, aero_dist%modes(i)%vol_frac, &
            num_per_bin, aero_state)
    end do

    aero_state%comp_vol = dble(n_part) / dist_num_conc(bin_grid, aero_dist)

  end subroutine aero_dist_to_part
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine aero_state_to_binned(bin_grid, aero_binned, aero_data, aero_state)
    
    ! Create the bin number and mass arrays from aero_state%v.

    use mod_bin_grid
    use mod_aero_data
    use mod_aero_binned
    
    type(bin_grid_t), intent(in) :: bin_grid ! bin grid
    type(aero_binned_t), intent(out) :: aero_binned ! binned distributions
    type(aero_data_t), intent(in) :: aero_data ! aerosol data
    type(aero_state_t), intent(inout) :: aero_state ! aerosol state
    
    integer b, j, s
    
    aero_binned%num_den = 0d0
    aero_binned%vol_den = 0d0
    do b = 1,bin_grid%n_bin
       do j = 1,aero_state%n(b)
          aero_binned%vol_den(b,:) = aero_binned%vol_den(b,:) &
               + aero_state%v(b)%p(j,:) / aero_state%comp_vol / bin_grid%dlnr
       end do
    end do
    aero_binned%num_den = dble(aero_state%n) / aero_state%comp_vol &
         / bin_grid%dlnr
   
  end subroutine aero_state_to_binned
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine double(aero_state)
    
    ! Doubles number of particles.

    type(aero_state_t), intent(inout) :: aero_state ! aerosol state
    
    integer :: n_bin, k, n
    
    n_bin = size(aero_state%n)
    do k = 1,n_bin
       n = aero_state%n(k)
       call enlarge_bin_to(aero_state%v(k), 2 * n)
       aero_state%v(k)%p((n+1):(2*n), :) = aero_state%v(k)%p(1:n, :)
       aero_state%n(k) = 2 * aero_state%n(k)
    end do
    aero_state%comp_vol = 2d0 * aero_state%comp_vol
    
  end subroutine double
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine resort_aero_state(bin_grid, aero_state)
    
    ! Takes a VH array where the particle volumes might no longer be
    ! correct for the bins they are in and resorts it so that every
    ! particle is in the correct bin.

    use mod_aero_data
    use mod_bin_grid
    
    type(bin_grid_t), intent(in) :: bin_grid ! bin_grid
    type(aero_state_t), intent(inout) :: aero_state ! aerosol state
    
    integer :: bin, j, new_bin, k
    real*8 :: pv
    
    ! The approach here is inefficient because we might reprocess
    ! particles. For example, if we are doing bin 1 and we shift a
    ! particle up to bin 2, when we do bin 2 we will reprocess it. It
    ! seems to be more trouble than it's worth to worry about this
    ! yet, however.
    
    do bin = 1,bin_grid%n_bin
       j = 1
       do while (j .le. aero_state%n(bin))
          ! find the new volume and new bin
          pv = particle_volume(aero_state%v(bin)%p(j,:))
          call particle_in_bin(pv, bin_grid, new_bin)
          
          ! if the bin number has changed, move the particle
          if (bin .ne. new_bin) then
             ! move the particle to the new bin, leaving a hole
             aero_state%n(new_bin) = aero_state%n(new_bin) + 1
             call enlarge_bin_to(aero_state%v(new_bin), aero_state%n(new_bin))
             aero_state%v(new_bin)%p(aero_state%n(new_bin),:) &
                  = aero_state%v(bin)%p(j,:)
             
             ! copy the last particle in the current bin into the hole
             ! if the hole isn't in fact the last particle
             if (j .lt. aero_state%n(bin)) then
                aero_state%v(bin)%p(j,:) &
                     = aero_state%v(bin)%p(aero_state%n(bin),:)
             end if
             aero_state%n(bin) = aero_state%n(bin) - 1
             if (aero_state%n(bin) .lt. 0) then
                write(0,*) 'ERROR: invalid aero_state%n in bin ', bin
                call exit(2)
             end if
             
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
       call shrink_bin(aero_state%n(bin), aero_state%v(bin))
    end do
    
  end subroutine resort_aero_state
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine check_aero_state(bin_grid, aero_binned, aero_data, aero_state)
    
    ! Check that all particles are in the correct bins and that the
    ! bin numbers and masses are correct. This is for debugging only.

    use mod_util
    use mod_bin_grid
    use mod_aero_data
    use mod_aero_binned
    
    type(bin_grid_t), intent(in) :: bin_grid ! bin_grid
    type(aero_binned_t), intent(out) :: aero_binned ! binned distributions
    type(aero_data_t), intent(in) :: aero_data ! aerosol data
    type(aero_state_t), intent(inout) :: aero_state ! aerosol state
    
    real*8 pv, check_bin_v, check_vol_den(aero_data%n_spec), vol_tol
    real*8 num_tol, state_num_den
    integer i, k, k_check, s
    logical error
    
    error = .false.
    
    ! check that all particles are in the correct bins
    do k = 1,bin_grid%n_bin
       do i = 1,aero_state%n(k)
          pv = particle_volume(aero_state%v(k)%p(i,:))
          call particle_in_bin(pv, bin_grid, k_check)
          if (k .ne. k_check) then
             write(0,'(a10,a10,a20,a10)') 'k', 'i', 'aero_state%v(k, i)', &
                  'k_check'
             write(0,'(i10,i10,e20.8,i10)') k, i, pv, k_check
             error = .true.
          end if
       end do
    end do
    
    ! check the aero_binned%num_den array
    do k = 1,bin_grid%n_bin
       num_tol = 0.01d0 / aero_state%comp_vol / bin_grid%dlnr
       state_num_den = dble(aero_state%n(k)) / aero_state%comp_vol &
            / bin_grid%dlnr
       if (.not. almost_equal_abs(state_num_den, &
            aero_binned%num_den(k), num_tol)) then
          write(0,'(a10,a20,a20,a20,a20)') 'k', 'aero_state%n(k)', &
               'state_num_den', 'a_binned%num_den(k)', 'comp_vol'
          write(0,'(i10,i20,e20.10,e20.10,e20.10)') k, aero_state%n(k), &
               state_num_den, aero_binned%num_den(k), aero_state%comp_vol
          error = .true.
       end if
    end do
    
    ! check the aero_binned%vol_den array
    do k = 1,bin_grid%n_bin
       check_vol_den = sum(aero_state%v(k)%p(1:aero_state%n(k),:), 1) &
            / aero_state%comp_vol / bin_grid%dlnr
       vol_tol = bin_grid%v(k) / 1d3 / bin_grid%dlnr
       do s = 1,aero_data%n_spec
          if (.not. almost_equal_abs(check_vol_den(s), &
               aero_binned%vol_den(k,s), vol_tol)) then
             write(0,'(a10,a10,a25,a25)') 'k', 's', 'check_vol_den(s)', &
                  'vol_den(k,s)'
             write(0,'(i10,i10,e25.10,e25.10)') k, s, check_vol_den(s), &
                  aero_binned%vol_den(k,s)
             error = .true.
          end if
       end do
    end do
    
    if (error) then
       write(0,*) 'ERROR: check_aero_state() failed'
       call exit(2)
    end if
    
  end subroutine check_aero_state
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine inout_write_bin_p(file, bin_p)
    
    ! Write full state.
    
    use mod_inout
    
    type(inout_file_t), intent(inout) :: file ! file to write to
    type(bin_p_t), intent(in) :: bin_p ! bin_p to write

    call inout_write_real_array_2d(file, 'bin_p', bin_p%p)
    
  end subroutine inout_write_bin_p

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine inout_write_bin_p_array(file, bin_ps)
    
    ! Write full state.
    
    use mod_inout
    
    type(inout_file_t), intent(inout) :: file ! file to write to
    type(bin_p_t), intent(in) :: bin_ps(:) ! bin_p array to write

    integer :: length, i

    length = size(bin_ps)
    call inout_write_integer(file, 'bin_p_array_len', length)
    do i = 1,length
       call inout_write_integer(file, 'bin_p_array_entry', i)
       call inout_write_bin_p(file, bin_ps(i))
    end do
    
  end subroutine inout_write_bin_p_array

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine inout_write_aero_state(file, aero_state)
    
    ! Write full state.
    
    use mod_inout
    
    type(inout_file_t), intent(inout) :: file ! file to write to
    type(aero_state_t), intent(in) :: aero_state ! aero_state to write

    call inout_write_real(file, "comp_vol(m^3)", aero_state%comp_vol)
    call inout_write_integer_array(file, "number_per_bin", aero_state%n)
    call inout_write_bin_p_array(file, aero_state%v)
    
  end subroutine inout_write_aero_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine inout_read_bin_p(file, bin_p)
    
    ! Read full state.
    
    use mod_inout
    
    type(inout_file_t), intent(inout) :: file ! file to write to
    type(bin_p_t), intent(out) :: bin_p ! bin_p to read

    call inout_read_real_array_2d(file, 'bin_p', bin_p%p)
    
  end subroutine inout_read_bin_p

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine inout_read_bin_p_array(file, bin_ps)
    
    ! Read full state.
    
    use mod_inout
    
    type(inout_file_t), intent(inout) :: file ! file to write to
    type(bin_p_t), pointer :: bin_ps(:) ! bin_p array to read

    integer :: length, i, check_i

    call inout_read_integer(file, 'bin_p_array_len', length)
    allocate(bin_ps(length))
    do i = 1,length
       call inout_read_integer(file, 'bin_p_array_entry', check_i)
       call inout_check_index(file, i, check_i)
       call inout_read_bin_p(file, bin_ps(i))
    end do
    
  end subroutine inout_read_bin_p_array

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine inout_read_aero_state(file, aero_state)
    
    ! Read full state.
    
    use mod_inout
    
    type(inout_file_t), intent(inout) :: file ! file to write to
    type(aero_state_t), intent(out) :: aero_state ! aero_state to read

    call inout_read_real(file, "comp_vol(m^3)", aero_state%comp_vol)
    call inout_read_integer_array(file, "number_per_bin", aero_state%n)
    call inout_read_bin_p_array(file, aero_state%v)
    
  end subroutine inout_read_aero_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
end module mod_aero_state
