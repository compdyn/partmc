! Copyright (C) 2005-2009 Nicole Riemer and Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_coagulation_mpi_equal module.

!> Aerosol particle coagulation with MPI where each node has its own
!> aero_state and all nodes perform coagulation equally.
module pmc_coagulation_mpi_equal

  use pmc_bin_grid
  use pmc_aero_data
  use pmc_util
  use pmc_env_state
  use pmc_aero_state
  use pmc_coagulation
  use pmc_mpi
#ifdef PMC_USE_MPI
  use mpi
#endif

  !> Size of the outgoing buffer for \c bsend (bytes).
  !!
  !! FIXME: check that this size is big enough. It must be large
  !! enough to handle the required number of messages of the given
  !! sizes, plus MPI_BSEND_OVERHEAD per message, plus some extra room
  !! because it's only kind of a circular buffer --- the messages
  !! themselves aren't allowed to wrap around then end, so we might
  !! need extra space up to the size of the largest message type.
  integer, parameter :: COAG_EQUAL_OUTGOING_BUFFER_SIZE      = 1000000
  !> Size of send and receive buffer for each message (bytes).
  !!
  !! FIXME: check that this size is big enough. The biggest message
  !! type will be one of the particle-sending types, for which we need
  !! pmc_mpi_pack_size_aero_particle(aero_particle), plus a couple of
  !! integers or something. At the moment this means something like (10
  !! + n_spec) reals, (3 + 2) integers, which for n_spec = 20 gives a
  !! size of 260 bytes.
  integer, parameter :: COAG_EQUAL_MAX_BUFFER_SIZE           = 10000
  integer, parameter :: COAG_EQUAL_MAX_REQUESTS              = 1
  integer, parameter :: COAG_EQUAL_TAG_REQUEST_PARTICLE      = 5321
  integer, parameter :: COAG_EQUAL_TAG_RETURN_REQ_PARTICLE   = 5322
  integer, parameter :: COAG_EQUAL_TAG_RETURN_UNREQ_PARTICLE = 5323
  integer, parameter :: COAG_EQUAL_TAG_RETURN_NO_PARTICLE    = 5324
  integer, parameter :: COAG_EQUAL_TAG_DONE                  = 5325

  type request_t
     !> Local \c aero_particle to maybe coagulate with the received
     !> particle.
     type(aero_particle_t) :: local_aero_particle
     !> Remote processor number that we sent the request to
     !> (-1 means this request is currently not used).
     integer :: remote_proc
     !> Local bin number from which we took \c local_aero_particle.
     integer :: local_bin
     !> Remote bin number from which we requested an \c aero_particle.
     integer :: remote_bin
     !> Whether this request is currently active
     logical :: active
  end type request_t
  
contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Allocate a request object and set its state to 
  subroutine request_allocate(request)

    !> Request object to allocate.
    type(request_t), intent(out) :: request

    call aero_particle_allocate(request%local_aero_particle)
    request%active = .false.

  end subroutine request_allocate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Deallocate a request object and set it to be invalid.
  subroutine request_deallocate(request)
   
    !> Request object to deallocate
    type(request_t), intent(inout) :: request

    call aero_particle_deallocate(request%local_aero_particle)
    request%active = .false.

  end subroutine request_deallocate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Whether the given reqest object is currectly active.
  logical function request_is_active(request)

    !> Request object to test for activeness.
    type(request_t), intent(in) :: request

    request_is_active = request%active

  end function request_is_active

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Do coagulation for time del_t.
  subroutine mc_coag_mpi_equal(kernel, bin_grid, env_state, aero_data, &
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

#ifdef PMC_USE_MPI
    logical :: samps_remaining, sent_dones
    integer :: i, j, n_samp, i_samp, rank, i_proc, n_proc, i_bin
    integer :: ierr, status(MPI_STATUS_SIZE), current_i, current_j, i_req
    real(kind=dp) :: n_samp_real
    integer, allocatable :: n_parts(:,:)
    real(kind=dp), allocatable :: comp_vols(:)
    type(request_t) :: requests(COAG_EQUAL_MAX_REQUESTS)
    integer :: n_samps(bin_grid%n_bin, bin_grid%n_bin)
    logical, allocatable :: procs_done(:)
    integer :: outgoing_buffer(COAG_EQUAL_OUTGOING_BUFFER_SIZE)
    integer :: outgoing_buffer_size_check
    
    rank = pmc_mpi_rank()
    n_proc = pmc_mpi_size()

    call pmc_mpi_barrier()

    allocate(n_parts(bin_grid%n_bin, n_proc))
    allocate(comp_vols(n_proc))
    call sync_info(aero_state%bin(:)%n_part, aero_state%comp_vol, &
         n_parts, comp_vols)
    !>DEBUG
    !call sleep(rank * 2)
    !<DEBUG

    call generate_n_samps(bin_grid, n_parts, comp_vols, del_t, k_max, &
         n_samps)
    tot_n_samp = sum(n_samps)

    ! main loop
    do i_req = 1,COAG_EQUAL_MAX_REQUESTS
       call request_allocate(requests(i_req))
    end do
    samps_remaining = .true.
    current_i = 1
    current_j = 1
    allocate(procs_done(n_proc))
    procs_done = .false.
    sent_dones = .false.
    call mpi_buffer_attach(outgoing_buffer, &
         COAG_EQUAL_OUTGOING_BUFFER_SIZE, ierr)
    call pmc_mpi_check_ierr(ierr)
    do while (.not. all(procs_done))
       ! add requests if we have any slots available call
       call add_coagulation_requests(bin_grid, env_state, aero_data, &
            aero_state, requests, n_parts, current_i, current_j, n_samps, &
            samps_remaining, del_t, k_max, comp_vols, procs_done)

       ! if we have no outstanding requests, send done messages
       if (.not. sent_dones) then
          if (.not. any_requests_active(requests)) then
             sent_dones = .true.
             do i_proc = 0, (n_proc - 1)
                call send_done(i_proc)
             end do
          end if
       end if

       !>DEBUG
       !call sleep(1 + rank * 2)
       !<DEBUG
       ! receive exactly one message
       call coag_equal_recv(requests, bin_grid, &
            env_state, aero_data, aero_state, del_t, k_max, kernel, &
            tot_n_coag, comp_vols, procs_done)
    end do

    do i_req = 1,COAG_EQUAL_MAX_REQUESTS
       call assert(502009333, .not. request_is_active(requests(i_req)))
       call request_deallocate(requests(i_req))
    end do
    deallocate(procs_done)
    deallocate(n_parts)
    deallocate(comp_vols)
    call mpi_buffer_detach(outgoing_buffer, &
         outgoing_buffer_size_check, ierr)
    call pmc_mpi_check_ierr(ierr)
    call assert(577822730, &
         COAG_EQUAL_OUTGOING_BUFFER_SIZE == outgoing_buffer_size_check)
#endif

  end subroutine mc_coag_mpi_equal

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine coag_equal_recv(requests, bin_grid, &
       env_state, aero_data, aero_state, del_t, k_max, kernel, &
       tot_n_coag, comp_vols, procs_done)

    !> Array of outstanding requests.
    type(request_t), intent(inout) :: requests(COAG_EQUAL_MAX_REQUESTS)
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
    !> Number of coagulation events.
    integer, intent(inout) :: tot_n_coag
    !> Computational volumes on all processors.
    real(kind=dp), intent(in) :: comp_vols(:)
    !> Which processors are finished with coagulation.
    logical, intent(inout) :: procs_done(:)

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

#ifdef PMC_USE_MPI
    character(len=100) :: error_msg
    integer :: status(MPI_STATUS_SIZE), ierr
    !>DEBUG
    integer :: i_req
    !<DEBUG

    !>DEBUG
    !write(*,*) pmc_mpi_rank(), 'coag_equal_recv: entry'
    !do i_req = 1,COAG_EQUAL_MAX_REQUESTS
    !   if (request_is_active(requests(i_req))) then
    !      write(*,*) pmc_mpi_rank(), ' request/active/remote = ', i_req, &
    !           request_is_active(requests(i_req)), requests(i_req)%remote_proc, &
    !           requests(i_req)%remote_bin
    !   else
    !      write(*,*) pmc_mpi_rank(), ' request/active = ', i_req, &
    !           request_is_active(requests(i_req))
    !   end if
    !end do
    !<DEBUG
    call mpi_probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &
         status, ierr)
    call pmc_mpi_check_ierr(ierr)
    !>DEBUG
    !write(*,*) pmc_mpi_rank(), 'coag_equal_recv: tag = ', status(MPI_TAG)
    !write(*,*) pmc_mpi_rank(), 'coag_equal_recv: from = ', status(MPI_SOURCE)
    !<DEBUG
    if (status(MPI_TAG) == COAG_EQUAL_TAG_REQUEST_PARTICLE) then
       call recv_request_particle(aero_state)
    elseif (status(MPI_TAG) == COAG_EQUAL_TAG_RETURN_REQ_PARTICLE) then
       call recv_return_req_particle(requests, bin_grid, &
            env_state, aero_data, aero_state, del_t, k_max, kernel, &
            tot_n_coag, comp_vols)
    elseif (status(MPI_TAG) == COAG_EQUAL_TAG_RETURN_UNREQ_PARTICLE) then
       call recv_return_unreq_particle(aero_state, bin_grid)
    elseif (status(MPI_TAG) == COAG_EQUAL_TAG_RETURN_NO_PARTICLE) then
       call recv_return_no_particle(requests, bin_grid, &
            aero_data, aero_state)
    elseif (status(MPI_TAG) == COAG_EQUAL_TAG_DONE) then
       call recv_done(procs_done)
    else
       write(error_msg, '(a,i20)') 'unknown tag', status(MPI_TAG)
       call die_msg(856123972, error_msg)
    end if
#endif
    
  end subroutine coag_equal_recv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine add_coagulation_requests(bin_grid, env_state, aero_data, &
       aero_state, requests, n_parts, local_bin, remote_bin, &
       n_samps, samps_remaining, del_t, k_max, comp_vols, procs_done)

    !> Bin grid.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Environment state.
    type(env_state_t), intent(in) :: env_state
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Aerosol state.
    type(aero_state_t), intent(inout) :: aero_state
    !> Array of outstanding requests.
    type(request_t), intent(inout) :: requests(COAG_EQUAL_MAX_REQUESTS)
    !> Number of particles per bin per processor.
    integer, intent(in) :: n_parts(:,:)
    !> Bin index of first particle we need to coagulate.
    integer, intent(inout) :: local_bin
    !> Bin index of second particle we need to coagulate.
    integer, intent(inout) :: remote_bin
    !> Whether there are still coagulation samples that need to be done.
    logical, intent(inout) :: samps_remaining
    !> Number of samples remaining per bin pair
    integer, intent(inout) :: n_samps(bin_grid%n_bin,bin_grid%n_bin)
    !> Timestep.
    real(kind=dp), intent(in) :: del_t
    !> Maximum kernel.
    real(kind=dp), intent(in) :: k_max(bin_grid%n_bin,bin_grid%n_bin)
    !> Computational volumes on all processors.
    real(kind=dp), intent(in) :: comp_vols(:)
    !> Which processors are finished with coagulation.
    logical, intent(inout) :: procs_done(:)

    integer :: i_req

    if (.not. samps_remaining) return

    outer: do i_req = 1,COAG_EQUAL_MAX_REQUESTS
       if (.not. request_is_active(requests(i_req))) then
          inner: do
             call update_n_samps(bin_grid, n_samps, local_bin, &
                  remote_bin, samps_remaining)
             if (.not. samps_remaining) exit outer
             if (aero_state%bin(local_bin)%n_part > 0) then
                call find_rand_remote_proc(bin_grid, n_parts, &
                     remote_bin, requests(i_req)%remote_proc)
                requests(i_req)%active = .true.
                requests(i_req)%local_bin = local_bin
                requests(i_req)%remote_bin = remote_bin
                !>DEBUG
                !write(*,*) pmc_mpi_rank(), 'add request, local,remote = ', pmc_mpi_rank(), local_bin, &
                !     requests(i_req)%remote_proc, remote_bin
                !<DEBUG
                call aero_state_remove_rand_particle_from_bin(aero_state, &
                     local_bin, requests(i_req)%local_aero_particle)
                call send_request_particle(requests(i_req)%remote_proc, &
                     requests(i_req)%remote_bin)
                exit inner
             end if
          end do inner
       end if
    end do outer

  end subroutine add_coagulation_requests

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns \c .true. if any of the requests are active, otherwise
  !> returns \c .false.
  logical function any_requests_active(requests)

    !> Array of outstanding requests.
    type(request_t), intent(inout) :: requests(COAG_EQUAL_MAX_REQUESTS)
    
    integer :: i_req

    do i_req = 1,COAG_EQUAL_MAX_REQUESTS
       if (request_is_active(requests(i_req))) then
          any_requests_active = .true.
          return
       end if
    end do
    any_requests_active = .false.
    
  end function any_requests_active

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine find_rand_remote_proc(bin_grid, n_parts, remote_bin, remote_proc)

    !> Bin grid.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Number of particles per bin per processor.
    integer, intent(in) :: n_parts(:,:)
    !> Remote bin number.
    integer, intent(in) :: remote_bin
    !> Remote processor number chosen at random.
    integer, intent(out) :: remote_proc

#ifdef PMC_USE_MPI
    integer :: rank, n_proc

    rank = pmc_mpi_rank()
    n_proc = pmc_mpi_size()
    call assert(542705260, size(n_parts, 1) == bin_grid%n_bin)
    call assert(770964285, size(n_parts, 2) == n_proc)
    remote_proc = sample_disc_pdf(n_proc, n_parts(remote_bin,:)) - 1
#endif

  end subroutine find_rand_remote_proc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine update_n_samps(bin_grid, n_samps, local_bin, remote_bin, &
       samps_remaining)

    !> Bin grid.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Number of samples remaining per bin pair
    integer, intent(inout) :: n_samps(bin_grid%n_bin,bin_grid%n_bin)
    !> Bin index of first particle we need to coagulate.
    integer, intent(inout) :: local_bin
    !> Bin index of second particle we need to coagulate.
    integer, intent(inout) :: remote_bin
    !> Whether there are still coagulation samples that need to be done.
    logical, intent(inout) :: samps_remaining

    if (.not. samps_remaining) return

    do
       if (n_samps(local_bin, remote_bin) > 0) exit

       remote_bin = remote_bin + 1
       if (remote_bin > bin_grid%n_bin) then
          remote_bin = 1
          local_bin = local_bin + 1
       end if
       if (local_bin > bin_grid%n_bin) exit
    end do
    
    if (local_bin > bin_grid%n_bin) then
       samps_remaining = .false.
    else
       n_samps(local_bin, remote_bin) = n_samps(local_bin, remote_bin) - 1
    end if

  end subroutine update_n_samps

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine send_request_particle(remote_proc, remote_bin)

    !> Remote processor number.
    integer, intent(in) :: remote_proc
    !> Remote bin number.
    integer, intent(in) :: remote_bin

#ifdef PMC_USE_MPI
    character :: buffer(COAG_EQUAL_MAX_BUFFER_SIZE)
    integer :: buffer_size, position, ierr

    !>DEBUG
    !write(*,*) pmc_mpi_rank(), 'send_request_particle: entry'
    !write(*,*) pmc_mpi_rank(), 'send_request_particle: dest proc/bin = ', remote_proc, remote_bin
    !<DEBUG
    position = 0
    call pmc_mpi_pack_integer(buffer, position, remote_bin)
    buffer_size = position
    call assert(490250818, buffer_size < COAG_EQUAL_MAX_BUFFER_SIZE)
    call mpi_bsend(buffer, buffer_size, MPI_CHARACTER, remote_proc, &
         COAG_EQUAL_TAG_REQUEST_PARTICLE, MPI_COMM_WORLD, ierr)
    call pmc_mpi_check_ierr(ierr)
#endif

  end subroutine send_request_particle

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine recv_request_particle(aero_state)

    !> Aero state.
    type(aero_state_t), intent(inout) :: aero_state

#ifdef PMC_USE_MPI
    integer :: buffer_size, position, request_bin, sent_proc, ierr
    integer :: remote_proc, status(MPI_STATUS_SIZE)
    character :: buffer(COAG_EQUAL_MAX_BUFFER_SIZE)
    type(aero_particle_t) :: aero_particle

    !>DEBUG
    !write(*,*) pmc_mpi_rank(), 'recv_request_particle: entry'
    !<DEBUG
    ! get the message
    call mpi_recv(buffer, COAG_EQUAL_MAX_BUFFER_SIZE, MPI_CHARACTER, &
         MPI_ANY_SOURCE, COAG_EQUAL_TAG_REQUEST_PARTICLE, MPI_COMM_WORLD, &
         status, ierr)
    call pmc_mpi_check_ierr(ierr)
    call assert(920139874, status(MPI_TAG) &
         == COAG_EQUAL_TAG_REQUEST_PARTICLE)
    call mpi_get_count(status, MPI_CHARACTER, buffer_size, ierr)
    call pmc_mpi_check_ierr(ierr)
    call assert(190658659, buffer_size < COAG_EQUAL_MAX_BUFFER_SIZE)
    remote_proc = status(MPI_SOURCE)

    ! unpack it
    position = 0
    call pmc_mpi_unpack_integer(buffer, position, request_bin)
    call assert(895128380, position == buffer_size)
    !>DEBUG
    !write(*,*) pmc_mpi_rank(), 'recv_request_particle: from proc/bin = ', remote_proc, request_bin
    !<DEBUG

    ! send the particle back if we have one
    if (aero_state%bin(request_bin)%n_part == 0) then
       call send_return_no_particle(remote_proc, request_bin)
    else
       call aero_particle_allocate(aero_particle)
       call aero_state_remove_rand_particle_from_bin(aero_state, &
            request_bin, aero_particle)
       call send_return_req_particle(aero_particle, request_bin, &
            remote_proc)
       call aero_particle_deallocate(aero_particle)
    end if
#endif

  end subroutine recv_request_particle

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine send_return_no_particle(dest_proc, i_bin)

    !> Processor number to send message to.
    integer, intent(in) :: dest_proc
    !> Bin number where there was no particle.
    integer, intent(in) :: i_bin

#ifdef PMC_USE_MPI
    character :: buffer(COAG_EQUAL_MAX_BUFFER_SIZE)
    integer :: buffer_size, position, ierr

    !>DEBUG
    !write(*,*) pmc_mpi_rank(), 'send_return_no_particle: entry'
    !write(*,*) pmc_mpi_rank(), 'send_return_no_particle: dest proc/bin = ', dest_proc, i_bin
    !<DEBUG
    position = 0
    call pmc_mpi_pack_integer(buffer, position, i_bin)
    buffer_size = position
    call assert(473131470, buffer_size < COAG_EQUAL_MAX_BUFFER_SIZE)
    call mpi_bsend(buffer, buffer_size, MPI_CHARACTER, dest_proc, &
         COAG_EQUAL_TAG_RETURN_NO_PARTICLE, MPI_COMM_WORLD, ierr)
    call pmc_mpi_check_ierr(ierr)
#endif

  end subroutine send_return_no_particle

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine recv_return_no_particle(requests, bin_grid, &
       aero_data, aero_state)

    !> Array of outstanding requests.
    type(request_t), intent(inout) :: requests(COAG_EQUAL_MAX_REQUESTS)
    !> Bin grid.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Aerosol state.
    type(aero_state_t), intent(inout) :: aero_state

#ifdef PMC_USE_MPI
    logical :: found_request
    integer :: buffer_size, position, sent_bin, sent_proc, i_req, ierr
    integer :: status(MPI_STATUS_SIZE)
    character :: buffer(COAG_EQUAL_MAX_BUFFER_SIZE)

    !>DEBUG
    !write(*,*) pmc_mpi_rank(), 'recv_return_no_particle: entry'
    !<DEBUG
    ! get the message
    call mpi_recv(buffer, COAG_EQUAL_MAX_BUFFER_SIZE, MPI_CHARACTER, &
         MPI_ANY_SOURCE, COAG_EQUAL_TAG_RETURN_NO_PARTICLE, &
         MPI_COMM_WORLD, status, ierr)
    call pmc_mpi_check_ierr(ierr)
    call assert(918153221, status(MPI_TAG) &
         == COAG_EQUAL_TAG_RETURN_NO_PARTICLE)
    call mpi_get_count(status, MPI_CHARACTER, buffer_size, ierr)
    call pmc_mpi_check_ierr(ierr)
    call assert(461111487, buffer_size < COAG_EQUAL_MAX_BUFFER_SIZE)
    sent_proc = status(MPI_SOURCE)

    ! unpack it
    position = 0
    call pmc_mpi_unpack_integer(buffer, position, sent_bin)
    call assert(518172999, position == buffer_size)
    !>DEBUG
    !write(*,*) pmc_mpi_rank(), 'recv_return_no_particle: from proc/bin = ', sent_proc, sent_bin
    !<DEBUG

    ! find the matching request
    found_request = .false.
    do i_req = 1,COAG_EQUAL_MAX_REQUESTS
       if ((requests(i_req)%remote_proc == sent_proc) &
            .and. (requests(i_req)%remote_bin == sent_bin)) then
          found_request = .true.
          exit
       end if
    end do
    call assert(215612776, found_request)

    ! we can't do coagulation with the local particle, so store it back
    call aero_state_add_particle(aero_state, requests(i_req)%local_bin, &
         requests(i_req)%local_aero_particle)
    call request_deallocate(requests(i_req))
    call request_allocate(requests(i_req))
#endif

  end subroutine recv_return_no_particle

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine send_return_req_particle(aero_particle, i_bin, dest_proc)

    !> Aero particle to send.
    type(aero_particle_t), intent(in) :: aero_particle
    !> Bin that the particle is in.
    integer, intent(in) :: i_bin
    !> Processor number to send particle to.
    integer, intent(in) :: dest_proc

#ifdef PMC_USE_MPI
    character :: buffer(COAG_EQUAL_MAX_BUFFER_SIZE)
    integer :: buffer_size, position, ierr

    !>DEBUG
    !write(*,*) pmc_mpi_rank(), 'send_return_req_particle: entry'
    !write(*,*) pmc_mpi_rank(), 'send_return_req_particle: to proc/bin/id = ', dest_proc, i_bin, aero_particle%id
    !<DEBUG
    position = 0
    call pmc_mpi_pack_integer(buffer, position, i_bin)
    call pmc_mpi_pack_aero_particle(buffer, position, aero_particle)
    buffer_size = position
    call assert(622476860, buffer_size < COAG_EQUAL_MAX_BUFFER_SIZE)
    call mpi_bsend(buffer, buffer_size, MPI_CHARACTER, dest_proc, &
         COAG_EQUAL_TAG_RETURN_REQ_PARTICLE, MPI_COMM_WORLD, ierr)
    call pmc_mpi_check_ierr(ierr)
#endif

  end subroutine send_return_req_particle

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine recv_return_req_particle(requests, bin_grid, &
       env_state, aero_data, aero_state, del_t, k_max, kernel, &
       tot_n_coag, comp_vols)

    !> Array of outstanding requests.
    type(request_t), intent(inout) :: requests(COAG_EQUAL_MAX_REQUESTS)
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
    !> Number of coagulation events.
    integer, intent(inout) :: tot_n_coag
    !> Computational volumes on all processors.
    real(kind=dp), intent(in) :: comp_vols(:)

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

#ifdef PMC_USE_MPI
    logical :: found_request
    integer :: buffer_size, position, sent_bin, sent_proc, i_req, ierr
    integer :: status(MPI_STATUS_SIZE)
    character :: buffer(COAG_EQUAL_MAX_BUFFER_SIZE)
    type(aero_particle_t) :: sent_aero_particle
    real(kind=dp) :: k, p

    !>DEBUG
    !write(*,*) pmc_mpi_rank(), 'recv_return_req_particle: entry'
    !<DEBUG
    ! get the message
    call mpi_recv(buffer, COAG_EQUAL_MAX_BUFFER_SIZE, MPI_CHARACTER, &
         MPI_ANY_SOURCE, COAG_EQUAL_TAG_RETURN_REQ_PARTICLE, &
         MPI_COMM_WORLD, status, ierr)
    call pmc_mpi_check_ierr(ierr)
    !>DEBUG
    !write(*,*) pmc_mpi_rank(), 'recv_return_req_particle: tag/from = ', &
    !     status(MPI_TAG), status(MPI_SOURCE)
    !<DEBUG
    call assert(133285061, status(MPI_TAG) &
         == COAG_EQUAL_TAG_RETURN_REQ_PARTICLE)
    call mpi_get_count(status, MPI_CHARACTER, buffer_size, ierr)
    call pmc_mpi_check_ierr(ierr)
    call assert(461111487, buffer_size < COAG_EQUAL_MAX_BUFFER_SIZE)
    sent_proc = status(MPI_SOURCE)

    ! unpack it
    position = 0
    call pmc_mpi_unpack_integer(buffer, position, sent_bin)
    call aero_particle_allocate(sent_aero_particle)
    call pmc_mpi_unpack_aero_particle(buffer, position, sent_aero_particle)
    call assert(753356021, position == buffer_size)
    !>DEBUG
    !write(*,*) pmc_mpi_rank(), 'recv_return_req_particle: from proc/bin/id = ', &
    !     sent_proc, sent_bin, sent_aero_particle%id
    !<DEBUG

    ! find the matching request
    found_request = .false.
    do i_req = 1,COAG_EQUAL_MAX_REQUESTS
       if ((requests(i_req)%remote_proc == sent_proc) &
            .and. (requests(i_req)%remote_bin == sent_bin)) then
          found_request = .true.
          exit
       end if
    end do
    call assert(579308475, found_request)
    !>DEBUG
    !write(*,*) pmc_mpi_rank(), 'recv_return_req_particle: matches request # ', i_req
    !<DEBUG
    
    ! maybe do coagulation
    call kernel(requests(i_req)%local_aero_particle, sent_aero_particle, &
         aero_data, env_state, k)
    p = k / k_max(requests(i_req)%local_bin, sent_bin)

    if (pmc_random() .lt. p) then
       !>DEBUG
       !write(*,*) pmc_mpi_rank(), 'recv_return_req_particle: coagulation occured'
       !<DEBUG
       ! coagulation happened, do it
       tot_n_coag = tot_n_coag + 1
       call coagulate_mpi_equal(bin_grid, aero_data, aero_state, &
            requests(i_req)%local_aero_particle, sent_aero_particle, &
            sent_proc, comp_vols)
    else
       !>DEBUG
       !write(*,*) pmc_mpi_rank(), 'recv_return_req_particle: coagulation did not occur'
       !<DEBUG
       ! coagulation didn't happen, send the particles back
       call aero_state_add_particle(aero_state, requests(i_req)%local_bin, &
            requests(i_req)%local_aero_particle)
       call send_return_unreq_particle(sent_aero_particle, sent_proc)
    end if

    call request_deallocate(requests(i_req))
    call request_allocate(requests(i_req))
    call aero_particle_deallocate(sent_aero_particle)
#endif

  end subroutine recv_return_req_particle

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine send_return_unreq_particle(aero_particle, dest_proc)

    !> Aero particle to send.
    type(aero_particle_t), intent(in) :: aero_particle
    !> Process to send the particle to.
    integer, intent(in) :: dest_proc

#ifdef PMC_USE_MPI
    character :: buffer(COAG_EQUAL_MAX_BUFFER_SIZE)
    integer :: buffer_size, position, ierr

    !>DEBUG
    !write(*,*) pmc_mpi_rank(), 'send_return_unreq_particle: entry'
    !write(*,*) pmc_mpi_rank(), 'send_return_unreq_particle: dest proc/id = ', dest_proc, aero_particle%id
    !<DEBUG
    position = 0
    call pmc_mpi_pack_aero_particle(buffer, position, aero_particle)
    buffer_size = position
    call assert(622476860, buffer_size < COAG_EQUAL_MAX_BUFFER_SIZE)
    call mpi_bsend(buffer, buffer_size, MPI_CHARACTER, dest_proc, &
         COAG_EQUAL_TAG_RETURN_UNREQ_PARTICLE, MPI_COMM_WORLD, ierr)
    call pmc_mpi_check_ierr(ierr)
#endif

  end subroutine send_return_unreq_particle

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine recv_return_unreq_particle(aero_state, bin_grid)

    !> Aerosol state.
    type(aero_state_t), intent(inout) :: aero_state
    !> Bin grid.
    type(bin_grid_t), intent(in) :: bin_grid

#ifdef PMC_USE_MPI
    logical :: found_request
    integer :: buffer_size, position, sent_proc, ierr
    character :: buffer(COAG_EQUAL_MAX_BUFFER_SIZE)
    type(aero_particle_t) :: aero_particle
    integer :: i_bin, status(MPI_STATUS_SIZE), send_proc

    !>DEBUG
    !write(*,*) pmc_mpi_rank(), 'recv_return_unreq_particle: entry'
    !<DEBUG
    ! get the message
    call mpi_recv(buffer, COAG_EQUAL_MAX_BUFFER_SIZE, MPI_CHARACTER, &
         MPI_ANY_SOURCE, COAG_EQUAL_TAG_RETURN_UNREQ_PARTICLE, &
         MPI_COMM_WORLD, status, ierr)
    call pmc_mpi_check_ierr(ierr)
    call assert(496247788, status(MPI_TAG) &
         == COAG_EQUAL_TAG_RETURN_UNREQ_PARTICLE)
    call mpi_get_count(status, MPI_CHARACTER, buffer_size, ierr)
    call pmc_mpi_check_ierr(ierr)
    call assert(590644042, buffer_size < COAG_EQUAL_MAX_BUFFER_SIZE)
    sent_proc = status(MPI_SOURCE)

    ! unpack it
    position = 0
    call aero_particle_allocate(aero_particle)
    call pmc_mpi_unpack_aero_particle(buffer, position, aero_particle)
    call assert(833588594, position == buffer_size)
    !>DEBUG
    !write(*,*) pmc_mpi_rank(), 'recv_return_unreq_particle: from proc/id = ', sent_proc, aero_particle%id
    !<DEBUG

    ! put it back
    i_bin = aero_particle_in_bin(aero_particle, bin_grid)
    call aero_state_add_particle(aero_state, i_bin, aero_particle)
    call aero_particle_deallocate(aero_particle)
#endif

  end subroutine recv_return_unreq_particle

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Send a message saying that this processor is finished with its
  !> coagulation.
  subroutine send_done(dest_proc)

    !> Process to send the message to.
    integer, intent(in) :: dest_proc

#ifdef PMC_USE_MPI
    character :: buffer(COAG_EQUAL_MAX_BUFFER_SIZE)
    integer :: buffer_size, ierr

    !>DEBUG
    !write(*,*) pmc_mpi_rank(), 'sent_done: entry'
    !write(*,*) pmc_mpi_rank(), 'sent_done: dest proc = ', dest_proc
    !<DEBUG
    buffer_size = 0
    call mpi_bsend(buffer, buffer_size, MPI_CHARACTER, dest_proc, &
         COAG_EQUAL_TAG_DONE, MPI_COMM_WORLD, ierr)
    call pmc_mpi_check_ierr(ierr)
#endif

  end subroutine send_done

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Receive a done message.
  subroutine recv_done(procs_done)

    !> Which processors are finished with coagulation.
    logical, intent(inout) :: procs_done(:)
    
#ifdef PMC_USE_MPI
    integer :: buffer_size, sent_proc, ierr
    character :: buffer(COAG_EQUAL_MAX_BUFFER_SIZE)
    integer :: status(MPI_STATUS_SIZE)

    !>DEBUG
    !write(*,*) pmc_mpi_rank(), 'recv_done: entry'
    !<DEBUG
    ! get the message
    call mpi_recv(buffer, COAG_EQUAL_MAX_BUFFER_SIZE, MPI_CHARACTER, &
         MPI_ANY_SOURCE, COAG_EQUAL_TAG_DONE, MPI_COMM_WORLD, &
         status, ierr)
    call pmc_mpi_check_ierr(ierr)
    call assert(348737947, status(MPI_TAG) &
         == COAG_EQUAL_TAG_DONE)
    call mpi_get_count(status, MPI_CHARACTER, buffer_size, ierr)
    call pmc_mpi_check_ierr(ierr)
    call assert(214904056, buffer_size == 0)
    sent_proc = status(MPI_SOURCE)
    !>DEBUG
    !write(*,*) pmc_mpi_rank(), 'recv_done: from proc = ', sent_proc
    !<DEBUG

    ! process it
    procs_done(sent_proc + 1) = .true.
#endif

  end subroutine recv_done

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Do an allgather to exchange number of particles and computational
  !> volume information between all processors.
  subroutine sync_info(local_n_parts, local_comp_vol, &
       global_n_parts, global_comp_vols)

    !> Number of particles per bin on the local processor.
    integer, intent(in) :: local_n_parts(:)
    !> Computational volume on the local processor.
    real(kind=dp), intent(in) :: local_comp_vol
    !> Number of particles per bin on all processors.
    integer, intent(out) :: global_n_parts(:,:)
    !> Computational volumes on all processors (m^3).
    real(kind=dp), intent(out) :: global_comp_vols(:)

#ifdef PMC_USE_MPI
    integer :: n_bin, n_proc, ierr
    integer, allocatable :: send_buf(:), recv_buf(:)
    !> DEBUG
    integer :: i_bin, i_proc
    !< DEBUG

    n_bin = size(local_n_parts)
    n_proc = pmc_mpi_size()
    call assert(816230609, all(shape(global_n_parts) == (/n_bin, n_proc/)))
    call assert(883861456, all(shape(global_comp_vols) == (/n_proc/)))

    ! use a new send_buf to make sure the memory is contiguous
    allocate(send_buf(n_bin))
    allocate(recv_buf(n_bin * n_proc))
    send_buf = local_n_parts
    call mpi_allgather(send_buf, n_bin, MPI_INTEGER, &
         recv_buf, n_bin, MPI_INTEGER, MPI_COMM_WORLD, ierr)
    call pmc_mpi_check_ierr(ierr)
    global_n_parts = reshape(recv_buf, (/n_bin, n_proc/))
    deallocate(send_buf)
    deallocate(recv_buf)

    call mpi_allgather(local_comp_vol, 1, MPI_REAL8, &
         global_comp_vols, 1, MPI_REAL8, MPI_COMM_WORLD, ierr)
    call pmc_mpi_check_ierr(ierr)

    !> DEBUG
    if (pmc_mpi_rank() == 0) then
       !write(*,*) pmc_mpi_rank(), 'sync_info: global_comp_vols = ', global_comp_vols
       !write(*,*) pmc_mpi_rank(), 'sync_info: global_n_parts = '
       do i_bin = 1,n_bin
          !write(*,'(i20,i20)',advance='no') pmc_mpi_rank(), i_bin
          do i_proc = 0,(n_proc - 1)
             !write(*,'(i20)',advance='no') global_n_parts(i_bin, i_proc + 1)
          end do
          !write(*,*) ''
       end do
       !write(*,*) pmc_mpi_rank(), 'sync_info: total n_parts = ', sum(global_n_parts, 1)
    end if
    !< DEBUG
#endif

  end subroutine sync_info

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> generate the number of samples to do per bin pair.
  subroutine generate_n_samps(bin_grid, n_parts, comp_vols, del_t, &
       k_max, n_samps)
    
    !> Bin grid.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Number of particles per bin on all processors.
    integer, intent(in) :: n_parts(:,:)
    !> Computational volumes on all processors..
    real(kind=dp), intent(in) :: comp_vols(:)
    !> Timestep.
    real(kind=dp), intent(in) :: del_t
    !> Maximum kernel.
    real(kind=dp), intent(in) :: k_max(bin_grid%n_bin,bin_grid%n_bin)
    !> Number of samples to do per bin pair.
    integer, intent(out) :: n_samps(bin_grid%n_bin,bin_grid%n_bin)

    integer :: i, j, rank
    real(kind=dp) :: n_samp_real

    rank = pmc_mpi_rank()
    do i = 1,bin_grid%n_bin
       do j = 1,bin_grid%n_bin
          call coag_equal_compute_n_samp(n_parts(i, rank + 1), &
               sum(n_parts(j,:)), i == j, k_max(i,j), &
               sum(comp_vols), del_t, n_samp_real)
          ! probabalistically determine n_samp to cope with < 1 case
          n_samps(i,j) = prob_round(n_samp_real)
       end do
    end do

  end subroutine generate_n_samps

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine coag_equal_compute_n_samp(n_parts_local, &
       n_parts_total, same_bin, k_max, &
       comp_vol_total, del_t, n_samp_real)

    !> Number of particle on the local processor.
    integer, intent(in) :: n_parts_local
    !> Number of particles on all processors.
    integer, intent(in) :: n_parts_total
    !> Whether the bins are the same.
    logical, intent(in) :: same_bin
    !> Maximum kernel value (s^{-1}).
    real(kind=dp), intent(in) :: k_max
    !> Total computational volume
    real(kind=dp), intent(in) :: comp_vol_total
    !> Timestep (s).
    real(kind=dp), intent(in) :: del_t
    !> Number of samples to do per timestep.
    real(kind=dp), intent(out) :: n_samp_real

    real(kind=dp) :: r_samp
    real(kind=dp) :: n_possible ! use real(kind=dp) to avoid integer overflow
    ! FIXME: should use integer*8 or integer(kind = 8)
    ! or even better, di = selected_int_kind(18), integer(kind=di)
    ! to represent 10^{-18} to 10^{18}
    
    if (same_bin) then
       n_possible = real(n_parts_local, kind=dp) &
            * (real(n_parts_total, kind=dp) - 1d0) / 2d0
    else
       n_possible = real(n_parts_local, kind=dp) &
            * real(n_parts_total, kind=dp) / 2d0
    endif
    
    r_samp = k_max / comp_vol_total * del_t
    n_samp_real = r_samp * n_possible
  
  end subroutine coag_equal_compute_n_samp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine coagulate_mpi_equal(bin_grid, aero_data, aero_state, &
       aero_particle_1, aero_particle_2, remote_proc, comp_vols)

    !> Bin grid.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Aerosol state.
    type(aero_state_t), intent(inout) :: aero_state
    !> First particle to coagulate.
    type(aero_particle_t), intent(in) :: aero_particle_1
    !> Second particle to coagulate.
    type(aero_particle_t), intent(in) :: aero_particle_2
    !> Remote processor that the particle came from.
    integer, intent(in) :: remote_proc
    !> Computational volumes on all processors (m^3).
    real(kind=dp), intent(in) :: comp_vols(:)
    
    type(aero_particle_t) :: aero_particle_new
    integer :: bin_new, rank, proc_new
    type(aero_info_t) :: aero_info

    !>DEBUG
    !write(*,*) pmc_mpi_rank(), 'coagulate_mpi_equal: entry'
    !<DEBUG
    rank = pmc_mpi_rank()

    call aero_particle_allocate_size(aero_particle_new, aero_data%n_spec)
    call assert(400463852, aero_particle_1%id /= aero_particle_2%id)

    ! coagulate particles
    call aero_particle_coagulate(aero_particle_1, aero_particle_2, &
         aero_particle_new)
    bin_new = aero_particle_in_bin(aero_particle_new, bin_grid)

    ! removal information
    call aero_info_allocate(aero_info)
    if (aero_particle_new%id /= aero_particle_1%id) then
       ! aero_particle_1 is the removed particle
       call assert(209246658, aero_particle_new%id == aero_particle_2%id)
       aero_info%id = aero_particle_1%id
       aero_info%action = AERO_INFO_COAG
       aero_info%other_id = aero_particle_2%id
    else
       ! aero_particle_2 is the removed particle
       call assert(967187262, aero_particle_new%id /= aero_particle_2%id)
       aero_info%id = aero_particle_2%id
       aero_info%action = AERO_INFO_COAG
       aero_info%other_id = aero_particle_1%id
    end if
    call aero_info_array_add_aero_info(aero_state%aero_info_array, &
         aero_info)
    call aero_info_deallocate(aero_info)

    ! add new particle
    proc_new = sample_cts_pdf(size(comp_vols), comp_vols) - 1
    !>DEBUG
    !write(*,*) pmc_mpi_rank(), 'coagulate_mpi_equal: new proc = ', proc_new
    !<DEBUG
    call send_return_unreq_particle(aero_particle_new, proc_new)

    call aero_particle_deallocate(aero_particle_new)

  end subroutine coagulate_mpi_equal

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
end module pmc_coagulation_mpi_equal
