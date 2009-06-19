! Copyright (C) 2005-2009 Nicole Riemer and Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_coagulation module.

!> Aerosol particle coagulation.
module pmc_coagulation

  use pmc_bin_grid
  use pmc_aero_data
  use pmc_util
  use pmc_env_state
  use pmc_aero_state
  use pmc_mpi
#ifdef PMC_USE_MPI
  use mpi
#endif

  integer, parameter :: COAG_MAX_BUFFER_SIZE      = 1000000 ! 1 MB
  integer, parameter :: COAG_TAG_REQUEST_PARTICLE = 885784964
  integer, parameter :: COAG_TAG_REMOVE_PARTICLE  = 223971199
  integer, parameter :: COAG_TAG_SEND_PARTICLE    = 302095212
  integer, parameter :: COAG_TAG_DONE             = 116195367
  
contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Do coagulation for time del_t.
  subroutine mc_coag_mpi_controlled(kernel, bin_grid, env_state, aero_data, &
       aero_state, del_t, k_max, tot_n_samp, tot_n_coag)

    !> Bin grid.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Environment state.
    type(env_state_t), intent(in) :: env_state
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Aerosol state.
    type(aero_state_t), intent(inout) :: aero_state
    !> Timestep.
    real(kind=dp), intent(in) :: del_t
    !> Maximum kernel.
    real(kind=dp), intent(in) :: k_max(bin_grid%n_bin,bin_grid%n_bin)
    !> Total number of samples tested.
    integer, intent(out) :: tot_n_samp
    !> Number of coagulation events.
    integer, intent(out) :: tot_n_coag

#ifndef DOXYGEN_SKIP_DOC
    interface
       subroutine kernel(aero_particle_1, aero_particle_2, aero_data, &
            env_state, k)
         use pmc_aero_particle
         use pmc_aero_data
         use pmc_env_state
         type(aero_particle_t), intent(in) :: aero_particle_1
         type(aero_particle_t), intent(in) :: aero_particle_2
         type(aero_data_t), intent(in) :: aero_data
         type(env_state_t), intent(in) :: env_state  
         real(kind=dp), intent(out) :: k
       end subroutine kernel
    end interface
#endif

    logical :: did_coag
    integer :: i, j, n_samp, i_samp, rank, i_proc, n_proc, i_bin
    real(kind=dp) :: n_samp_real
    integer, allocatable :: n_parts(:,:)
    real(kind=dp), allocatable :: comp_vols(:)
    
    !>DEBUG
    !call pmc_mpi_barrier()
    !write(*,*) 'mc_coag_mpi_controlled: entry ', pmc_mpi_rank()
    !<DEBUG
    rank = pmc_mpi_rank()
    n_proc = pmc_mpi_size()
    allocate(n_parts(bin_grid%n_bin, n_proc))
    allocate(comp_vols(n_proc))
    do i_proc = 0,(n_proc - 1)
       !>DEBUG
       !write(*,*) 'mc_coag_mpi_controlled: transfer comp_vol ', i_proc, 0, pmc_mpi_rank()
       !<DEBUG
       call pmc_mpi_transfer_real(aero_state%comp_vol, &
            comp_vols(i_proc + 1), i_proc, 0)
       do i_bin = 1,bin_grid%n_bin
          !>DEBUG
          !write(*,*) 'mc_coag_mpi_controlled: transfer n_part ', i_bin, i_proc, 0, pmc_mpi_rank()
          !<DEBUG
          call pmc_mpi_transfer_integer(aero_state%bin(i_bin)%n_part, &
               n_parts(i_bin, i_proc + 1), i_proc, 0)
       end do
    end do
    !>DEBUG
    !call pmc_mpi_barrier()
    !write(*,*) 'mc_coag_mpi_controlled: done with transfer ', pmc_mpi_rank()
    !<DEBUG

    tot_n_samp = 0
    tot_n_coag = 0
    if (rank == 0) then
       !>DEBUG
       !call sleep(2)
       !write(*,*) 'mc_coag_mpi_controlled: comp_vols = ', comp_vols
       !do i_bin = 1,bin_grid%n_bin
       !   write(*,'(i20)',advance='no') i_bin
       !   do i_proc = 0,(n_proc - 1)
       !      write(*,'(i20)',advance='no') n_parts(i_bin, i_proc + 1)
       !   end do
       !   write(*,*) ''
       !end do
       !write(*,*) 'total particle counts = ', sum(n_parts, 1)
       !<DEBUG
       ! root node actually does the coagulation
       do i = 1,bin_grid%n_bin
          do j = 1,bin_grid%n_bin
             !>DEBUG
             !write(*,*) 'mc_coag_mpi_controlled: bins = ', i, j
             !write(*,*) 'mc_coag_mpi_controlled: n_parts = ', sum(n_parts(i,:)), sum(n_parts(j,:))
             !write(*,*) 'mc_coag_mpi_controlled: comp_vol = ', sum(comp_vols)
             !<DEBUG
             call compute_n_samp(sum(n_parts(i,:)), &
                  sum(n_parts(j,:)), i == j, k_max(i,j), &
                  sum(comp_vols), del_t, n_samp_real)
             ! probabalistically determine n_samp to cope with < 1 case
             n_samp = prob_round(n_samp_real)
             tot_n_samp = tot_n_samp + n_samp
             !>DEBUG
             !if (n_samp > 0) then
             !   write(*,*) 'mc_coag_mpi_controlled: i,j, n_samp = ', i, j, n_samp
             !end if
             !<DEBUG
             do i_samp = 1,n_samp
                ! check we still have enough particles to coagulate
                if ((sum(n_parts(i,:)) < 1) &
                     .or. (sum(n_parts(j,:)) < 1) &
                     .or. ((i == j) .and. (sum(n_parts(i,:)) < 2))) then
                   exit
                end if
                call maybe_coag_pair_mpi_controlled(bin_grid, env_state, &
                     aero_data, aero_state, i, j, del_t, k_max(i,j), &
                     kernel, did_coag, n_parts, comp_vols)
                if (did_coag) tot_n_coag = tot_n_coag + 1
             enddo
          enddo
       enddo
       ! terminate remote helper node loops
       do i_proc = 0,(n_proc - 1)
          call coag_remote_done(i_proc)
       end do
    else
       ! remote nodes just implement protocol
       call coag_remote_agent(bin_grid, aero_data, aero_state)
    end if
    deallocate(n_parts)
    deallocate(comp_vols)

  end subroutine mc_coag_mpi_controlled

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Compute the number of samples required for the pair of bins.
  subroutine compute_n_samp(ni, nj, same_bin, k_max, comp_vol, &
       del_t, n_samp_real)

    !> Number particles in first bin .
    integer, intent(in) :: ni
    !> Number particles in second bin.
    integer, intent(in) :: nj
    !> Whether first bin is second bin.
    logical, intent(in) :: same_bin
    !> Maximum kernel value.
    real(kind=dp), intent(in) :: k_max
    !> Computational volume (m^3).
    real(kind=dp), intent(in) :: comp_vol
    !> Timestep (s).
    real(kind=dp), intent(in) :: del_t
    !> Number of samples per timestep.
    real(kind=dp), intent(out) :: n_samp_real
    
    real(kind=dp) r_samp
    real(kind=dp) n_possible ! use real(kind=dp) to avoid integer overflow
    ! FIXME: should use integer*8 or integer(kind = 8)
    
    if (same_bin) then
       n_possible = real(ni, kind=dp) * (real(nj, kind=dp) - 1d0) / 2d0
    else
       n_possible = real(ni, kind=dp) * real(nj, kind=dp) / 2d0
    endif
    
    r_samp = k_max / comp_vol * del_t
    n_samp_real = r_samp * n_possible
    
  end subroutine compute_n_samp
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Implement remote protocol.
  subroutine coag_remote_agent(bin_grid, aero_data, aero_state)

    !> Bin grid.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Aerosol state.
    type(aero_state_t), intent(inout) :: aero_state

#ifdef PMC_USE_MPI
    logical :: done, record_removal
    integer :: status(MPI_STATUS_SIZE), ierr, i_bin, i_part
    integer :: buffer_size, position
    character :: buffer(COAG_MAX_BUFFER_SIZE)
    type(aero_info_t) :: aero_info
    type(aero_particle_t) :: aero_particle

    done = .false.
    do while (.not. done)
       !>DEBUG
       !write(*,*) 'coag_remote_agent: waiting... ', pmc_mpi_rank()
       !<DEBUG
       call mpi_recv(buffer, COAG_MAX_BUFFER_SIZE, MPI_CHARACTER, 0, &
            MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
       call pmc_mpi_check_ierr(ierr)
       call mpi_get_count(status, MPI_CHARACTER, buffer_size, ierr)
       call pmc_mpi_check_ierr(ierr)
       if (status(MPI_TAG) == COAG_TAG_REQUEST_PARTICLE) then
          !>DEBUG
          !write(*,*) 'coag_remote_agent: COAG_TAG_REQUEST_PARTICLE ', pmc_mpi_rank()
          !<DEBUG
          position = 0
          call pmc_mpi_unpack_integer(buffer, position, i_bin)
          call pmc_mpi_unpack_integer(buffer, position, i_part)
          call assert(743978983, position == buffer_size)
          !>DEBUG
          !write(*,*) 'coag_remote_agent: i_bin = ', i_bin, pmc_mpi_rank()
          !write(*,*) 'coag_remote_agent: i_part = ', i_part, pmc_mpi_rank()
          !<DEBUG

          call assert(256332719, i_bin >= 1)
          call assert(887849380, i_bin <= bin_grid%n_bin)
          call assert(998347627, i_part >= 1)
          call assert(774017621, i_part <= aero_state%bin(i_bin)%n_part)

          position = 0
          call pmc_mpi_pack_aero_particle(buffer, position, &
               aero_state%bin(i_bin)%particle(i_part))
          buffer_size = position
          !>DEBUG
          !write(*,*) 'coag_remote_agent: sending particle back', pmc_mpi_rank()
          !<DEBUG
          call mpi_send(buffer, buffer_size, MPI_CHARACTER, 0, &
               COAG_TAG_SEND_PARTICLE, MPI_COMM_WORLD, ierr)
          call pmc_mpi_check_ierr(ierr)
       elseif (status(MPI_TAG) == COAG_TAG_REMOVE_PARTICLE) then
          !>DEBUG
          !write(*,*) 'coag_remote_agent: COAG_TAG_REMOVE_PARTICLE ', pmc_mpi_rank()
          !<DEBUG
          call aero_info_allocate(aero_info)
          position = 0
          call pmc_mpi_unpack_integer(buffer, position, i_bin)
          call pmc_mpi_unpack_integer(buffer, position, i_part)
          call pmc_mpi_unpack_logical(buffer, position, record_removal)
          call pmc_mpi_unpack_aero_info(buffer, position, aero_info)
          call assert(822092586, position == buffer_size)
          !>DEBUG
          !write(*,*) 'coag_remote_agent: i_bin = ', i_bin, pmc_mpi_rank()
          !write(*,*) 'coag_remote_agent: i_part = ', i_part, pmc_mpi_rank()
          !write(*,*) 'coag_remote_agent: record_removal = ', record_removal, pmc_mpi_rank()
          !<DEBUG
          
          call assert(703875248, i_bin >= 1)
          call assert(254613726, i_bin <= bin_grid%n_bin)
          call assert(652695715, i_part >= 1)
          call assert(987564730, i_part <= aero_state%bin(i_bin)%n_part)

          call aero_state_remove_particle(aero_state, i_bin, i_part, &
               record_removal, aero_info)
          call aero_info_deallocate(aero_info)
       elseif (status(MPI_TAG) == COAG_TAG_SEND_PARTICLE) then
          !>DEBUG
          !write(*,*) 'coag_remote_agent: COAG_TAG_SEND_PARTICLE ', pmc_mpi_rank()
          !<DEBUG
          call aero_particle_allocate(aero_particle)
          position = 0
          call pmc_mpi_unpack_aero_particle(buffer, position, &
               aero_particle)
          call assert(517349299, position == buffer_size)

          i_bin = aero_particle_in_bin(aero_particle, bin_grid)
          !>DEBUG
          !write(*,*) 'coag_remote_agent: i_bin = ', i_bin, pmc_mpi_rank()
          !<DEBUG
          call aero_state_add_particle(aero_state, i_bin, aero_particle)
          call aero_particle_deallocate(aero_particle)
       elseif (status(MPI_TAG) == COAG_TAG_DONE) then
          !>DEBUG
          !write(*,*) 'coag_remote_agent: COAG_TAG_DONE ', pmc_mpi_rank()
          !<DEBUG
          done = .true.
       else
          call die_msg(734832087, "unknown tag")
       end if
    end do
#endif

  end subroutine coag_remote_agent

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine coag_remote_fetch_particle(aero_state, i_proc, i_bin, i_part, &
       aero_particle)

    !> Aerosol state on node 0.
    type(aero_state_t), intent(in) :: aero_state
    !> Processor to fetch from.
    integer, intent(in) :: i_proc
    !> Bin to fetch from.
    integer, intent(in) :: i_bin
    !> Particle number to fetch.
    integer, intent(in) :: i_part
    !> Particle to fetch into.
    type(aero_particle_t), intent(inout) :: aero_particle
    
#ifdef PMC_USE_MPI
    integer :: buffer_size, position, ierr, status(MPI_STATUS_SIZE)
    character :: buffer(COAG_MAX_BUFFER_SIZE)

    !>DEBUG
    !write(*,*) 'coag_remote_fetch_particle: i_proc = ', i_proc, pmc_mpi_rank()
    !write(*,*) 'coag_remote_fetch_particle: i_bin = ', i_bin, pmc_mpi_rank()
    !write(*,*) 'coag_remote_fetch_particle: i_part = ', i_part, pmc_mpi_rank()
    !<DEBUG
    if (i_proc == 0) then
       ! just read it out of aero_state locally on node 0
       !>DEBUG
       !write(*,*) 'coag_remote_fetch_particle: local ', pmc_mpi_rank()
       !<DEBUG
       call assert(977735278, i_bin >= 1)
       call assert(974530513, i_bin <= size(aero_state%bin))
       call assert(156583013, i_part >= 1)
       call assert(578709336, i_part <= aero_state%bin(i_bin)%n_part)
       call aero_particle_copy(aero_state%bin(i_bin)%particle(i_part), &
            aero_particle)
    else
       ! request particle
       position = 0
       call pmc_mpi_pack_integer(buffer, position, i_bin)
       call pmc_mpi_pack_integer(buffer, position, i_part)
       buffer_size = position
       !>DEBUG
       !write(*,*) 'coag_remote_fetch_particle: sending COAG_TAG_REQUEST_PARTICLE ', pmc_mpi_rank()
       !<DEBUG
       call mpi_send(buffer, buffer_size, MPI_CHARACTER, i_proc, &
            COAG_TAG_REQUEST_PARTICLE, MPI_COMM_WORLD, ierr)
       call pmc_mpi_check_ierr(ierr)

       ! get particle
       !>DEBUG
       !write(*,*) 'coag_remote_fetch_particle: waiting for COAG_TAG_SEND_PARTICLE ', pmc_mpi_rank()
       !<DEBUG
       call mpi_recv(buffer, COAG_MAX_BUFFER_SIZE, MPI_CHARACTER, i_proc, &
            COAG_TAG_SEND_PARTICLE, MPI_COMM_WORLD, status, ierr)
       call pmc_mpi_check_ierr(ierr)
       call mpi_get_count(status, MPI_CHARACTER, buffer_size, ierr)
       call pmc_mpi_check_ierr(ierr)

       position = 0
       call pmc_mpi_unpack_aero_particle(buffer, position, aero_particle)
       call assert(739976989, position == buffer_size)
    end if
#endif

  end subroutine coag_remote_fetch_particle

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine coag_remote_remove_particle(aero_state, i_proc, i_bin, i_part, &
       record_removal, aero_info)

    !> Aerosol state on node 0.
    type(aero_state_t), intent(inout) :: aero_state
    !> Processor to remove from.
    integer, intent(in) :: i_proc
    !> Bin to remove from.
    integer, intent(in) :: i_bin
    !> Particle number to remove.
    integer, intent(in) :: i_part
    !> Whether to record the removal in the aero_info_array.
    logical, intent(in) :: record_removal
    !> Removal info.
    type(aero_info_t), intent(in) :: aero_info
    
#ifdef PMC_USE_MPI
    integer :: buffer_size, position, ierr
    character :: buffer(COAG_MAX_BUFFER_SIZE)

    !>DEBUG
    !write(*,*) 'coag_remote_remove_particle: i_proc = ', i_proc, pmc_mpi_rank()
    !write(*,*) 'coag_remote_remove_particle: i_bin = ', i_bin, pmc_mpi_rank()
    !write(*,*) 'coag_remote_remove_particle: i_part = ', i_part, pmc_mpi_rank()
    !write(*,*) 'coag_remote_remove_particle: record_removal = ', record_removal, pmc_mpi_rank()
    !<DEBUG
    if (i_proc == 0) then
       ! just do it directly on the local aero_state
       !>DEBUG
       !write(*,*) 'coag_remote_remove_particle: local ', pmc_mpi_rank()
       !<DEBUG
       call aero_state_remove_particle(aero_state, i_bin, i_part, &
            record_removal, aero_info)
    else
       position = 0
       call pmc_mpi_pack_integer(buffer, position, i_bin)
       call pmc_mpi_pack_integer(buffer, position, i_part)
       call pmc_mpi_pack_logical(buffer, position, record_removal)
       call pmc_mpi_pack_aero_info(buffer, position, aero_info)
       buffer_size = position
       
       !>DEBUG
       !write(*,*) 'coag_remote_remove_particle: sending COAG_TAG_REMOVE_PARTICLE ', pmc_mpi_rank()
       !<DEBUG
       call mpi_send(buffer, buffer_size, MPI_CHARACTER, i_proc, &
            COAG_TAG_REMOVE_PARTICLE, MPI_COMM_WORLD, ierr)
       call pmc_mpi_check_ierr(ierr)
    end if
#endif

  end subroutine coag_remote_remove_particle

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine coag_remote_add_particle(bin_grid, aero_state, i_proc, &
       aero_particle)

    !> Bin grid.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Aerosol state on node 0.
    type(aero_state_t), intent(inout) :: aero_state
    !> Processor to add to.
    integer, intent(in) :: i_proc
    !> Particle to add.
    type(aero_particle_t), intent(in) :: aero_particle
    
#ifdef PMC_USE_MPI
    integer :: buffer_size, position, i_bin, ierr
    character :: buffer(COAG_MAX_BUFFER_SIZE)

    !>DEBUG
    !write(*,*) 'coag_remote_add_particle: i_proc = ', i_proc, pmc_mpi_rank()
    !<DEBUG
    if (i_proc == 0) then
       ! just do it directly on the local aero_state
       !>DEBUG
       !write(*,*) 'coag_remote_add_particle: local ', pmc_mpi_rank()
       !<DEBUG
       i_bin = aero_particle_in_bin(aero_particle, bin_grid)
       call aero_state_add_particle(aero_state, i_bin, aero_particle)
    else
       position = 0
       call pmc_mpi_pack_aero_particle(buffer, position, aero_particle)
       buffer_size = position
       
       !>DEBUG
       !write(*,*) 'coag_remote_add_particle: sending COAG_TAG_SEND_PARTICLE ', pmc_mpi_rank()
       !<DEBUG
       call mpi_send(buffer, buffer_size, MPI_CHARACTER, i_proc, &
            COAG_TAG_SEND_PARTICLE, MPI_COMM_WORLD, ierr)
       call pmc_mpi_check_ierr(ierr)
    end if
#endif

  end subroutine coag_remote_add_particle

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine coag_remote_done(i_proc)

    !> Processor to send "done" to.
    integer, intent(in) :: i_proc
    
#ifdef PMC_USE_MPI
    integer :: buffer_size, ierr
    character :: buffer(COAG_MAX_BUFFER_SIZE)

    !>DEBUG
    !write(*,*) 'coag_remote_done: i_proc = ', i_proc, pmc_mpi_rank()
    !<DEBUG
    if (i_proc > 0) then
       buffer_size = 0
       !>DEBUG
       !write(*,*) 'coag_remote_done: sending COAG_TAG_DONE ', pmc_mpi_rank()
       !<DEBUG
       call mpi_send(buffer, buffer_size, MPI_CHARACTER, i_proc, &
            COAG_TAG_DONE, MPI_COMM_WORLD, ierr)
       call pmc_mpi_check_ierr(ierr)
    end if
#endif

  end subroutine coag_remote_done

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Choose a random pair for potential coagulation and test its
  !> probability of coagulation. If it happens, do the coagulation and
  !> update all structures.
  !!
  !! The probability of a coagulation will be taken as <tt>(kernel /
  !! k_max)</tt>.
  subroutine maybe_coag_pair_mpi_controlled(bin_grid, env_state, &
       aero_data, aero_state, b1, b2, del_t, k_max, kernel, did_coag, &
       n_parts, comp_vols)

    !> Bin grid.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Environment state.
    type(env_state_t), intent(in) :: env_state
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Aerosol state.
    type(aero_state_t), intent(inout) :: aero_state
    !> Bin of first particle.
    integer, intent(in) :: b1
    !> Bin of second particle.
    integer, intent(in) :: b2
    !> Timestep.
    real(kind=dp), intent(in) :: del_t
    !> K_max scale factor.
    real(kind=dp), intent(in) :: k_max
    !> Whether a coagulation occured.
    logical, intent(out) :: did_coag
    !> Number of particles per-bin and per-processor.
    integer :: n_parts(:, :)
    !> Computational volumes for each processor.
    real(kind=dp) :: comp_vols(:)
    
#ifndef DOXYGEN_SKIP_DOC
    interface
       subroutine kernel(aero_particle_1, aero_particle_2, aero_data, &
            env_state, k)
         use pmc_aero_particle
         use pmc_aero_data
         use pmc_env_state
         type(aero_particle_t), intent(in) :: aero_particle_1
         type(aero_particle_t), intent(in) :: aero_particle_2
         type(aero_data_t), intent(in) :: aero_data
         type(env_state_t), intent(in) :: env_state  
         real(kind=dp), intent(out) :: k
       end subroutine kernel
    end interface
#endif
    
    integer :: p1, p2, s1, s2
    real(kind=dp) :: p, k
    type(aero_particle_t) :: particle_1, particle_2
    
    !>DEBUG
    !write(*,*) '************************************************************************************'
    !write(*,*) 'maybe_coag_pair_mpi_controlled: entry ', pmc_mpi_rank()
    !<DEBUG
    call aero_particle_allocate(particle_1)
    call aero_particle_allocate(particle_2)

    call assert(717627403, sum(n_parts(b1,:)) >= 1)
    call assert(994541886, sum(n_parts(b2,:)) >= 1)
    if (b1 == b2) then
       call assert(968872192, sum(n_parts(b1,:)) >= 2)
    end if
    
    did_coag = .false.
    
    call find_rand_pair_mpi_controlled(n_parts, b1, b2, p1, p2, s1, s2)
    !>DEBUG
    !write(*,*) 'maybe_coag_pair_mpi_controlled: find_rand_pair: ', p1, b1, s1, p2, b2, s2
    !write(*,*) 'maybe_coag_pair_mpi_controlled: fetching particle 1 ', pmc_mpi_rank()
    !<DEBUG
    call coag_remote_fetch_particle(aero_state, p1, b1, s1, particle_1)
    !>DEBUG
    !write(*,*) 'maybe_coag_pair_mpi_controlled: fetching particle 2 ', pmc_mpi_rank()
    !<DEBUG
    call coag_remote_fetch_particle(aero_state, p2, b2, s2, particle_2)
    call kernel(particle_1, particle_2, aero_data, env_state, k)
    p = k / k_max

    !>DEBUG
    !write(*,*) 'maybe_coag_pair_mpi_controlled: p1,p2 = ', p1, p2
    !<DEBUG
    if (pmc_random() .lt. p) then
       !>DEBUG
       !write(*,*) 'maybe_coag_pair_mpi_controlled: coag happened ', pmc_mpi_rank()
       !<DEBUG
       call coagulate_mpi_controlled(bin_grid, aero_data, aero_state, &
            p1, b1, s1, p2, b2, s2, particle_1, particle_2, n_parts, comp_vols)
       did_coag = .true.
    end if

    call aero_particle_deallocate(particle_1)
    call aero_particle_deallocate(particle_2)
    
  end subroutine maybe_coag_pair_mpi_controlled
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Choose a random pair for potential coagulation and test its
  !> probability of coagulation. If it happens, do the coagulation and
  !> update all structures.
  !!
  !! The probability of a coagulation will be taken as <tt>(kernel /
  !! k_max)</tt>.
  subroutine maybe_coag_pair(bin_grid, env_state, aero_data, &
       aero_state, b1, b2, del_t, k_max, kernel, did_coag)

    !> Bin grid.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Environment state.
    type(env_state_t), intent(in) :: env_state
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Aerosol state.
    type(aero_state_t), intent(inout) :: aero_state
    !> Bin of first particle.
    integer, intent(in) :: b1
    !> Bin of second particle.
    integer, intent(in) :: b2
    !> Timestep.
    real(kind=dp), intent(in) :: del_t
    !> K_max scale factor.
    real(kind=dp), intent(in) :: k_max
    !> Whether a coagulation occured.
    logical, intent(out) :: did_coag
    
#ifndef DOXYGEN_SKIP_DOC
    interface
       subroutine kernel(aero_particle_1, aero_particle_2, aero_data, &
            env_state, k)
         use pmc_aero_particle
         use pmc_aero_data
         use pmc_env_state
         type(aero_particle_t), intent(in) :: aero_particle_1
         type(aero_particle_t), intent(in) :: aero_particle_2
         type(aero_data_t), intent(in) :: aero_data
         type(env_state_t), intent(in) :: env_state  
         real(kind=dp), intent(out) :: k
       end subroutine kernel
    end interface
#endif
    
    integer :: s1, s2
    real(kind=dp) :: p, k
    
    did_coag = .false.
    
    call assert(210827476, aero_state%bin(b1)%n_part >= 1)
    call assert(368973460, aero_state%bin(b2)%n_part >= 1)
    if (b1 == b2) then
       call assert(528541565, aero_state%bin(b1)%n_part >= 2)
    end if
    
    call find_rand_pair(aero_state, b1, b2, s1, s2)
    call kernel(aero_state%bin(b1)%particle(s1), &
         aero_state%bin(b2)%particle(s2), aero_data, env_state, k)
    p = k / k_max
    
    if (pmc_random() .lt. p) then
       call coagulate(bin_grid, aero_data, aero_state, &
            b1, s1, b2, s2)
       did_coag = .true.
    end if
    
  end subroutine maybe_coag_pair
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Given bins b1 and b2, find a random pair of particles (b1, s1)
  !> and (b2, s2) that are not the same particle particle as each
  !> other.
  subroutine find_rand_pair(aero_state, b1, b2, s1, s2)
    
    !> Aerosol state.
    type(aero_state_t), intent(in) :: aero_state
    !> Bin number of first particle.
    integer, intent(in) :: b1
    !> Bin number of second particle.
    integer, intent(in) :: b2
    !> First rand particle.
    integer, intent(out) :: s1
    !> Second rand particle.
    integer, intent(out) :: s2

    ! check we have enough particles to avoid being stuck in an
    ! infinite loop below
    call assert(362349482, aero_state%bin(b1)%n_part >= 1)
    call assert(479121681, aero_state%bin(b2)%n_part >= 1)
    if (b1 == b2) then
       call assert(161928491, aero_state%bin(b1)%n_part >= 2)
    end if
    
    ! FIXME: rand() only returns a REAL*4, so we might not be able to
    ! generate all integers between 1 and M if M is too big.

100 s1 = int(pmc_random() * real(aero_state%bin(b1)%n_part, kind=dp)) + 1
    if ((s1 .lt. 1) .or. (s1 .gt. aero_state%bin(b1)%n_part)) goto 100
101 s2 = int(pmc_random() * real(aero_state%bin(b2)%n_part, kind=dp)) + 1
    if ((s2 .lt. 1) .or. (s2 .gt. aero_state%bin(b2)%n_part)) goto 101
    if ((b1 .eq. b2) .and. (s1 .eq. s2)) goto 101

! FIXME: enable this and delete the above junk
!    do
!       s1 = pmc_rand_int(aero_state%bin(b1)%n_part)
!       s2 = pmc_rand_int(aero_state%bin(b2)%n_part)
!       if (.not. ((b1 .eq. b2) .and. (s1 .eq. s2))) then
!          exit
!       end if
!    end do
    
  end subroutine find_rand_pair
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Given bins b1 and b2, find a random pair of particles (b1, s1)
  !> and (b2, s2) that are not the same particle particle as each
  !> other.
  subroutine find_rand_pair_mpi_controlled(n_parts, b1, b2, p1, p2, s1, s2)
    
    !> Number of particles per-bin and per-processor.
    integer :: n_parts(:,:)
    !> Bin number of first particle.
    integer, intent(in) :: b1
    !> Bin number of second particle.
    integer, intent(in) :: b2
    !> Processor of first particle.
    integer, intent(out) :: p1
    !> Processor of second particle.
    integer, intent(out) :: p2
    !> First rand particle.
    integer, intent(out) :: s1
    !> Second rand particle.
    integer, intent(out) :: s2

    ! check we have enough particles to avoid being stuck in an
    ! infinite loop below
    call assert(329209888, sum(n_parts(b1,:)) >= 1)
    call assert(745799706, sum(n_parts(b2,:)) >= 1)
    if (b1 == b2) then
       call assert(755331221, sum(n_parts(b1,:)) >= 2)
    end if
    
    ! FIXME: rand() only returns a REAL*4, so we might not be able to
    ! generate all integers between 1 and M if M is too big.

100 s1 = int(pmc_random() * dble(sum(n_parts(b1,:)))) + 1
    if ((s1 .lt. 1) .or. (s1 .gt. sum(n_parts(b1,:)))) goto 100
101 s2 = int(pmc_random() * dble(sum(n_parts(b2,:)))) + 1
    if ((s2 .lt. 1) .or. (s2 .gt. sum(n_parts(b2,:)))) goto 101
    if ((b1 .eq. b2) .and. (s1 .eq. s2)) goto 101

    !>DEBUG
    !write(*,*) 'find_rand_pair_mpi_controlled: s1/np1 = ', real(s1, kind=dp) / dble(sum(n_parts(b1,:)))
    !write(*,*) 'find_rand_pair_mpi_controlled: s2/np2 = ', real(s2, kind=dp) / dble(sum(n_parts(b2,:)))
    !write(*,*) 'find_rand_pair_mpi_controlled: b1,np1,s1,nps1 = ', b1, sum(n_parts(b1,:)), s1, n_parts(b1,:)
    !write(*,*) 'find_rand_pair_mpi_controlled: b2,np2,s2,nps2 = ', b2, sum(n_parts(b2,:)), s2, n_parts(b2,:)
    !<DEBUG
    do p1 = 0,(pmc_mpi_size() - 1)
       if (s1 <= n_parts(b1, p1 + 1)) then
          exit
       end if
       s1 = s1 - n_parts(b1, p1 + 1)
    end do
    do p2 = 0,(pmc_mpi_size() + 1)
       if (s2 <= n_parts(b2, p2 + 1)) then
          exit
       end if
       s2 = s2 - n_parts(b2, p2 + 1)
    end do
    !>DEBUG
    !write(*,*) 'find_rand_pair_mpi_controlled: b1,p1,s1 = ', b1, p1, s1
    !write(*,*) 'find_rand_pair_mpi_controlled: b2,p2,s2 = ', b2, p2, s2
    !<DEBUG
    
  end subroutine find_rand_pair_mpi_controlled
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Join together particles (b1, s1) and (b2, s2), updating all
  !> particle and bin structures to reflect the change.
  subroutine coagulate(bin_grid, aero_data, aero_state, b1, s1, b2, s2)
 
    !> Bin grid.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Aerosol state.
    type(aero_state_t), intent(inout) :: aero_state
    !> First particle (bin number).
    integer, intent(in) :: b1
    !> First particle (number in bin).
    integer, intent(in) :: s1
    !> Second particle (bin number).
    integer, intent(in) :: b2
    !> Second particle (number in bin).
    integer, intent(in) :: s2
    
    type(aero_particle_t), pointer :: particle_1, particle_2
    type(aero_particle_t) :: new_particle
    integer :: bn
    type(aero_info_t) :: aero_info
    logical :: p1_removed, p2_removed

    call aero_particle_allocate_size(new_particle, aero_data%n_spec)
    particle_1 => aero_state%bin(b1)%particle(s1)
    particle_2 => aero_state%bin(b2)%particle(s2)
    call assert(371947172, particle_1%id /= particle_2%id)

    ! coagulate particles
    call aero_particle_coagulate(particle_1, particle_2, new_particle)
    bn = aero_particle_in_bin(new_particle, bin_grid)

    ! remove old particles
    call aero_info_allocate(aero_info)
    if (new_particle%id /= particle_1%id) then
       ! particle_1 is the removed particle
       call assert(361912382, new_particle%id == particle_2%id)
       aero_info%id = particle_1%id
       aero_info%action = AERO_INFO_COAG
       aero_info%other_id = particle_2%id
       p1_removed = .true.
       p2_removed = .false.
    else
       ! particle_2 is the removed particle
       call assert(742917292, new_particle%id /= particle_2%id)
       aero_info%id = particle_2%id
       aero_info%action = AERO_INFO_COAG
       aero_info%other_id = particle_1%id
       p1_removed = .false.
       p2_removed = .true.
    end if
    if ((b1 == b2) .and. (s2 > s1)) then
       ! handle a tricky corner case where we have to watch for s2 or
       ! s1 being the last entry in the array and being repacked when
       ! the other one is removed
       call aero_state_remove_particle(aero_state, b2, s2, &
            p2_removed, aero_info)
       call aero_state_remove_particle(aero_state, b1, s1, &
            p1_removed, aero_info)
    else
       call aero_state_remove_particle(aero_state, b1, s1, &
            p1_removed, aero_info)
       call aero_state_remove_particle(aero_state, b2, s2, &
            p2_removed, aero_info)
    end if
    call aero_info_deallocate(aero_info)

    ! add new particle
    call aero_state_add_particle(aero_state, bn, new_particle)
    call aero_particle_deallocate(new_particle)
    
  end subroutine coagulate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Join together particles (b1, s1) and (b2, s2), updating all
  !> particle and bin structures to reflect the change.
  subroutine coagulate_mpi_controlled(bin_grid, aero_data, aero_state, &
       p1, b1, s1, p2, b2, s2, particle_1, particle_2, n_parts, comp_vols)

    !> Bin grid.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Aerosol state.
    type(aero_state_t), intent(inout) :: aero_state
    !> Processor of first particle.
    integer, intent(out) :: p1
    !> First particle (bin number).
    integer, intent(in) :: b1
    !> First particle (number in bin).
    integer, intent(in) :: s1
    !> Processor of second particle.
    integer, intent(out) :: p2
    !> Second particle (bin number).
    integer, intent(in) :: b2
    !> Second particle (number in bin).
    integer, intent(in) :: s2
    !> Copy of particle_1
    type(aero_particle_t) :: particle_1
    !> Copy of particle_2
    type(aero_particle_t) :: particle_2
    !> Number of particles per-bin and per-processor.
    integer :: n_parts(:, :)
    !> Computational volumes for each processor.
    real(kind=dp) :: comp_vols(:)
    
    type(aero_particle_t) :: new_particle
    type(aero_info_t) :: aero_info
    logical :: p1_removed, p2_removed
    integer :: pn, bn

    !>DEBUG
    !write(*,*) 'coagulate_mpi_controlled: entry ', pmc_mpi_rank()
    !<DEBUG
    call assert(907178147, particle_1%id /= particle_2%id)

    ! coagulate particles
    call aero_particle_allocate_size(new_particle, aero_data%n_spec)
    call aero_particle_coagulate(particle_1, particle_2, new_particle)

    ! remove old particles
    call aero_info_allocate(aero_info)
    if (new_particle%id /= particle_1%id) then
       ! particle_1 is the removed particle
       call assert(681297276, new_particle%id == particle_2%id)
       aero_info%id = particle_1%id
       aero_info%action = AERO_INFO_COAG
       aero_info%other_id = particle_2%id
       p1_removed = .true.
       p2_removed = .false.
    else
       ! particle_2 is the removed particle
       call assert(605428059, new_particle%id /= particle_2%id)
       aero_info%id = particle_2%id
       aero_info%action = AERO_INFO_COAG
       aero_info%other_id = particle_1%id
       p1_removed = .false.
       p2_removed = .true.
    end if
    if ((p1 == p2) .and. (b1 == b2) .and. (s2 > s1)) then
       ! handle a tricky corner case where we have to watch for s2 or
       ! s1 being the last entry in the array and being repacked when
       ! the other one is removed
       !>DEBUG
       !write(*,*) 'coagulate_mpi_controlled: removing ', p2, b2, s2, p2_removed, pmc_mpi_rank()
       !<DEBUG
       call coag_remote_remove_particle(aero_state, p2, b2, s2, &
            p2_removed, aero_info)
       !>DEBUG
       !write(*,*) 'coagulate_mpi_controlled: removing ', p1, b1, s1, p1_removed, pmc_mpi_rank()
       !<DEBUG
       call coag_remote_remove_particle(aero_state, p1, b1, s1, &
            p1_removed, aero_info)
    else
       !>DEBUG
       !write(*,*) 'coagulate_mpi_controlled: removing ', p1, b1, s1, p1_removed, pmc_mpi_rank()
       !<DEBUG
       call coag_remote_remove_particle(aero_state, p1, b1, s1, &
            p1_removed, aero_info)
       !>DEBUG
       !write(*,*) 'coagulate_mpi_controlled: removing ', p2, b2, s2, p2_removed, pmc_mpi_rank()
       !<DEBUG
       call coag_remote_remove_particle(aero_state, p2, b2, s2, &
            p2_removed, aero_info)
    end if
    call aero_info_deallocate(aero_info)

    ! where does the new particle go?
    bn = aero_particle_in_bin(new_particle, bin_grid)
    if (pmc_random() &
         < comp_vols(p1 + 1) / (comp_vols(p1 + 1) + comp_vols(p2 + 1))) then
       pn = p1
    else
       pn = p2
    end if
    
    ! add new particle
    !>DEBUG
    !write(*,*) 'coagulate_mpi_controlled: adding to ', pn, pmc_mpi_rank()
    !<DEBUG
    call coag_remote_add_particle(bin_grid, aero_state, pn, &
         new_particle)
    call aero_particle_deallocate(new_particle)

    ! fix up n_parts
    !>DEBUG
    !write(*,*) 'coagulate_mpi_controlled: b1,p1 = ', b1, p1, pmc_mpi_rank()
    !write(*,*) 'coagulate_mpi_controlled: b2,p2 = ', b2, p2, pmc_mpi_rank()
    !write(*,*) 'coagulate_mpi_controlled: bn,pn = ', bn, pn, pmc_mpi_rank()
    !write(*,*) 'coagulate_mpi_controlled: n_parts = ', size(n_parts,1), size(n_parts,2), pmc_mpi_rank()
    !<DEBUG
    n_parts(b1, p1 + 1) = n_parts(b1, p1 + 1) - 1
    n_parts(b2, p2 + 1) = n_parts(b2, p2 + 1) - 1
    n_parts(bn, pn + 1) = n_parts(bn, pn + 1) + 1
    !>DEBUG
    !write(*,*) 'coagulate_mpi_controlled: p1,p2,pn = ', p1, p2, pn
    !<DEBUG
    
  end subroutine coagulate_mpi_controlled

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
end module pmc_coagulation
