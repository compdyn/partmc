! Copyright (C) 2005-2012 Nicole Riemer and Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_coagulation_dist module.

!> Parallel aerosol particle coagulation with MPI
module pmc_coagulation_dist

  use pmc_bin_grid
  use pmc_aero_data
  use pmc_util
  use pmc_env_state
  use pmc_aero_state
  use pmc_aero_weight_array
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
  integer, parameter :: COAG_DIST_OUTGOING_BUFFER_SIZE      = 1000000
  !> Size of send and receive buffer for each message (bytes).
  !!
  !! The biggest message type will be one of the particle-sending
  !! types, for which we need pmc_mpi_pack_size_aero_particle(), plus
  !! a couple of integers or something. At the moment this means
  !! something like (10 + n_spec) reals, (3 + 2) integers, which for
  !! n_spec = 20 gives a size of 260 bytes.
  integer, parameter :: COAG_DIST_MAX_BUFFER_SIZE           = 10000
  integer, parameter :: COAG_DIST_MAX_REQUESTS              = 1
  integer, parameter :: COAG_DIST_TAG_REQUEST_PARTICLE      = 5321
  integer, parameter :: COAG_DIST_TAG_RETURN_REQ_PARTICLE   = 5322
  integer, parameter :: COAG_DIST_TAG_RETURN_UNREQ_PARTICLE = 5323
  integer, parameter :: COAG_DIST_TAG_RETURN_NO_PARTICLE    = 5324
  integer, parameter :: COAG_DIST_TAG_DONE                  = 5325

  !> A single outstanding request for a remote particle.
  type request_t
     !> Local \c aero_particle to maybe coagulate with the received
     !> particle.
     type(aero_particle_t) :: local_aero_particle
     !> Remote process number that we sent the request to
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

  ! Allocate a request object and set its state to invalid.
  subroutine request_allocate(request)

    !> Request object to allocate.
    type(request_t), intent(out) :: request

    request%active = .false.

  end subroutine request_allocate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Deallocate a request object and set it to be invalid.
  subroutine request_deallocate(request)

    !> Request object to deallocate
    type(request_t), intent(inout) :: request

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
  subroutine mc_coag_dist(coag_kernel_type, env_state, aero_data, &
       aero_state, del_t, tot_n_samp, tot_n_coag)

    !> Coagulation kernel type.
    integer, intent(in) :: coag_kernel_type
    !> Environment state.
    type(env_state_t), intent(in) :: env_state
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Aerosol state.
    type(aero_state_t), intent(inout) :: aero_state
    !> Timestep.
    real(kind=dp), intent(in) :: del_t
    !> Total number of samples tested.
    integer, intent(out) :: tot_n_samp
    !> Number of coagulation events.
    integer, intent(out) :: tot_n_coag

    integer, parameter :: s1 = 1
    integer, parameter :: s2 = 1
    integer, parameter :: sc = 1

#ifdef PMC_USE_MPI
    logical :: samps_remaining, sent_dones
    integer :: i_bin, j_bin, n_samp, i_samp, i_proc, n_proc
    integer :: ierr, status(MPI_STATUS_SIZE), current_i, current_j, i_req
    real(kind=dp) :: n_samp_real, f_max
    integer, allocatable :: n_parts(:,:)
    real(kind=dp), allocatable :: magnitudes(:,:)
    type(request_t) :: requests(COAG_DIST_MAX_REQUESTS)
    integer, allocatable :: n_samps(:,:)
    real(kind=dp), allocatable :: accept_factors(:,:), k_max(:,:)
    logical, allocatable :: procs_done(:)
    integer, allocatable :: outgoing_buffer(:)
    integer :: outgoing_buffer_size_check
    type(aero_weight_array_t) :: aero_weight_total

    call assert_msg(667898741, &
         aero_sorted_n_class(aero_state%aero_sorted) == 1, &
         "FIXME: mc_coag_dist() can only handle one weight class")

    allocate(outgoing_buffer(COAG_DIST_OUTGOING_BUFFER_SIZE))

    n_proc = pmc_mpi_size()

    call pmc_mpi_barrier()

    call aero_state_sort(aero_state, aero_data, all_procs_same=.true.)
    if (.not. aero_state%aero_sorted%coag_kernel_bounds_valid) then
       call est_k_minmax_binned_unweighted(aero_state%aero_sorted%bin_grid, &
            coag_kernel_type, aero_data, env_state, &
            aero_state%aero_sorted%coag_kernel_min, &
            aero_state%aero_sorted%coag_kernel_max)
       aero_state%aero_sorted%coag_kernel_bounds_valid = .true.
    end if

    allocate(n_samps(bin_grid_size(aero_state%aero_sorted%bin_grid), &
         bin_grid_size(aero_state%aero_sorted%bin_grid)))
    allocate(accept_factors(bin_grid_size(aero_state%aero_sorted%bin_grid), &
         bin_grid_size(aero_state%aero_sorted%bin_grid)))

    allocate(n_parts(bin_grid_size(aero_state%aero_sorted%bin_grid), n_proc))
    call pmc_mpi_allgather_integer_array(integer_varray_n_entry( &
         aero_state%aero_sorted%size_class%inverse(:, s1)), n_parts)

    allocate(magnitudes(size(aero_state%awa%weight), n_proc))
    call pmc_mpi_allgather_real_array(aero_state%awa%weight(:, s1)%magnitude, &
         magnitudes)

    aero_weight_total = aero_state%awa
    aero_weight_total%weight(:, s1)%magnitude = 1d0 / sum(1d0 / magnitudes, 2)

    allocate(k_max(bin_grid_size(aero_state%aero_sorted%bin_grid), &
         bin_grid_size(aero_state%aero_sorted%bin_grid)))
    do i_bin = 1,bin_grid_size(aero_state%aero_sorted%bin_grid)
       do j_bin = 1,bin_grid_size(aero_state%aero_sorted%bin_grid)
          call max_coag_num_conc_factor(aero_weight_total, &
               aero_data, aero_state%aero_sorted%bin_grid, &
               i_bin, j_bin, s1, s2, sc, f_max)
          k_max(i_bin, j_bin) &
               = aero_state%aero_sorted%coag_kernel_max(i_bin, j_bin) * f_max
       end do
    end do

    call generate_n_samps(n_parts, del_t, aero_state%aero_sorted%bin_grid, &
         aero_weight_total, k_max, n_samps, accept_factors)
    tot_n_samp = sum(n_samps)
    tot_n_coag = 0

    ! main loop
    do i_req = 1,COAG_DIST_MAX_REQUESTS
       call request_allocate(requests(i_req))
    end do
    samps_remaining = .true.
    current_i = 1
    current_j = 1
    allocate(procs_done(n_proc))
    procs_done = .false.
    sent_dones = .false.
    call mpi_buffer_attach(outgoing_buffer, &
         COAG_DIST_OUTGOING_BUFFER_SIZE, ierr)
    call pmc_mpi_check_ierr(ierr)
    do while (.not. all(procs_done))
       ! add requests if we have any slots available call
       call add_coagulation_requests(aero_state, requests, n_parts, &
            current_i, current_j, n_samps, samps_remaining)

       ! if we have no outstanding requests, send done messages
       if (.not. sent_dones) then
          if (.not. any_requests_active(requests)) then
             sent_dones = .true.
             do i_proc = 0, (n_proc - 1)
                call send_done(i_proc)
             end do
          end if
       end if

       ! receive exactly one message
       call coag_dist_recv(requests, env_state, aero_weight_total, aero_data, &
            aero_state, accept_factors, coag_kernel_type, tot_n_coag, &
            magnitudes, procs_done)
    end do

    do i_req = 1,COAG_DIST_MAX_REQUESTS
       call assert(502009333, .not. request_is_active(requests(i_req)))
       call request_deallocate(requests(i_req))
    end do
    deallocate(procs_done)
    deallocate(n_samps)
    deallocate(accept_factors)
    deallocate(n_parts)
    deallocate(magnitudes)
    call mpi_buffer_detach(outgoing_buffer, &
         outgoing_buffer_size_check, ierr)
    call pmc_mpi_check_ierr(ierr)
    call assert(577822730, &
         COAG_DIST_OUTGOING_BUFFER_SIZE == outgoing_buffer_size_check)
    deallocate(outgoing_buffer)
#endif

  end subroutine mc_coag_dist

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine coag_dist_recv(requests, env_state, aero_weight_total, &
       aero_data, aero_state, accept_factors, coag_kernel_type, tot_n_coag, &
       magnitudes, procs_done)

    !> Array of outstanding requests.
    type(request_t), intent(inout) :: requests(COAG_DIST_MAX_REQUESTS)
    !> Environment state.
    type(env_state_t), intent(in) :: env_state
    !> Total weighting functions.
    type(aero_weight_array_t), intent(in) :: aero_weight_total
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Aerosol state.
    type(aero_state_t), intent(inout) :: aero_state
    !> Accept scale factors per bin pair (1).
    real(kind=dp), intent(in) :: accept_factors(:,:)
    !> Coagulation kernel type.
    integer, intent(in) :: coag_kernel_type
    !> Number of coagulation events.
    integer, intent(inout) :: tot_n_coag
    !> Computational volumes on all processes.
    real(kind=dp), intent(in) :: magnitudes(:,:)
    !> Which processes are finished with coagulation.
    logical, intent(inout) :: procs_done(:)

#ifdef PMC_USE_MPI
    integer :: status(MPI_STATUS_SIZE), ierr

    call mpi_probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &
         status, ierr)
    call pmc_mpi_check_ierr(ierr)
    if (status(MPI_TAG) == COAG_DIST_TAG_REQUEST_PARTICLE) then
       call recv_request_particle(aero_state)
    elseif (status(MPI_TAG) == COAG_DIST_TAG_RETURN_REQ_PARTICLE) then
       call recv_return_req_particle(requests, env_state, aero_weight_total, &
            aero_data, aero_state, accept_factors, coag_kernel_type, &
            tot_n_coag, magnitudes)
    elseif (status(MPI_TAG) == COAG_DIST_TAG_RETURN_UNREQ_PARTICLE) then
       call recv_return_unreq_particle(aero_state, aero_data)
    elseif (status(MPI_TAG) == COAG_DIST_TAG_RETURN_NO_PARTICLE) then
       call recv_return_no_particle(requests, aero_data, aero_state)
    elseif (status(MPI_TAG) == COAG_DIST_TAG_DONE) then
       call recv_done(procs_done)
    else
       call die_msg(856123972, &
            'unknown tag: ' // trim(integer_to_string(status(MPI_TAG))))
    end if
#endif

  end subroutine coag_dist_recv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine add_coagulation_requests(aero_state, requests, n_parts, &
       local_bin, remote_bin, n_samps, samps_remaining)

    !> Aerosol state.
    type(aero_state_t), intent(inout) :: aero_state
    !> Array of outstanding requests.
    type(request_t), intent(inout) :: requests(COAG_DIST_MAX_REQUESTS)
    !> Number of particles per bin per process.
    integer, intent(in) :: n_parts(:,:)
    !> Bin index of first particle we need to coagulate.
    integer, intent(inout) :: local_bin
    !> Bin index of second particle we need to coagulate.
    integer, intent(inout) :: remote_bin
    !> Number of samples remaining per bin pair
    integer, intent(inout) :: n_samps(:,:)
    !> Whether there are still coagulation samples that need to be done.
    logical, intent(inout) :: samps_remaining

    integer, parameter :: s1 = 1
    integer, parameter :: s2 = 1
    integer, parameter :: sc = 1

    integer :: i_req

    if (.not. samps_remaining) return

    outer: do i_req = 1,COAG_DIST_MAX_REQUESTS
       if (.not. request_is_active(requests(i_req))) then
          inner: do
             call update_n_samps(n_samps, local_bin, remote_bin, &
                  samps_remaining)
             if (.not. samps_remaining) exit outer
             if (integer_varray_n_entry( &
                  aero_state%aero_sorted%size_class%inverse(local_bin, s2)) &
                  > 0) then
                call find_rand_remote_proc(n_parts, remote_bin, &
                     requests(i_req)%remote_proc)
                requests(i_req)%active = .true.
                requests(i_req)%local_bin = local_bin
                requests(i_req)%remote_bin = remote_bin
                call aero_state_remove_rand_particle_from_bin(aero_state, &
                     local_bin, s2, requests(i_req)%local_aero_particle)
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
    type(request_t), intent(inout) :: requests(COAG_DIST_MAX_REQUESTS)

    integer :: i_req

    do i_req = 1,COAG_DIST_MAX_REQUESTS
       if (request_is_active(requests(i_req))) then
          any_requests_active = .true.
          return
       end if
    end do
    any_requests_active = .false.

  end function any_requests_active

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine find_rand_remote_proc(n_parts, remote_bin, remote_proc)

    !> Number of particles per bin per process.
    integer, intent(in) :: n_parts(:,:)
    !> Remote bin number.
    integer, intent(in) :: remote_bin
    !> Remote process number chosen at random.
    integer, intent(out) :: remote_proc

#ifdef PMC_USE_MPI
    call assert(770964285, size(n_parts, 2) == pmc_mpi_size())
    remote_proc = sample_disc_pdf(n_parts(remote_bin, :)) - 1
#endif

  end subroutine find_rand_remote_proc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine update_n_samps(n_samps, local_bin, remote_bin, samps_remaining)

    !> Number of samples remaining per bin pair
    integer, intent(inout) :: n_samps(:,:)
    !> Bin index of first particle we need to coagulate.
    integer, intent(inout) :: local_bin
    !> Bin index of second particle we need to coagulate.
    integer, intent(inout) :: remote_bin
    !> Whether there are still coagulation samples that need to be done.
    logical, intent(inout) :: samps_remaining

    integer :: n_bin

    if (.not. samps_remaining) return

    n_bin = size(n_samps, 1)
    do
       if (n_samps(local_bin, remote_bin) > 0) exit

       remote_bin = remote_bin + 1
       if (remote_bin > n_bin) then
          remote_bin = 1
          local_bin = local_bin + 1
       end if
       if (local_bin > n_bin) exit
    end do

    if (local_bin > n_bin) then
       samps_remaining = .false.
    else
       n_samps(local_bin, remote_bin) = n_samps(local_bin, remote_bin) - 1
    end if

  end subroutine update_n_samps

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine send_request_particle(remote_proc, remote_bin)

    !> Remote process number.
    integer, intent(in) :: remote_proc
    !> Remote bin number.
    integer, intent(in) :: remote_bin

#ifdef PMC_USE_MPI
    character :: buffer(COAG_DIST_MAX_BUFFER_SIZE)
    integer :: buffer_size, max_buffer_size, position, ierr

    max_buffer_size = 0
    max_buffer_size = max_buffer_size + pmc_mpi_pack_size_integer(remote_bin)
    call assert(893545122, max_buffer_size <= COAG_DIST_MAX_BUFFER_SIZE)
    position = 0
    call pmc_mpi_pack_integer(buffer, position, remote_bin)
    call assert(610314213, position <= max_buffer_size)
    buffer_size = position
    call mpi_bsend(buffer, buffer_size, MPI_CHARACTER, remote_proc, &
         COAG_DIST_TAG_REQUEST_PARTICLE, MPI_COMM_WORLD, ierr)
    call pmc_mpi_check_ierr(ierr)
#endif

  end subroutine send_request_particle

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine recv_request_particle(aero_state)

    !> Aero state.
    type(aero_state_t), intent(inout) :: aero_state

    integer, parameter :: s1 = 1
    integer, parameter :: s2 = 1
    integer, parameter :: sc = 1

#ifdef PMC_USE_MPI
    integer :: buffer_size, position, request_bin, sent_proc
    integer :: ierr, remote_proc, status(MPI_STATUS_SIZE)
    character :: buffer(COAG_DIST_MAX_BUFFER_SIZE)
    type(aero_particle_t) :: aero_particle

    ! get the message
    call mpi_recv(buffer, COAG_DIST_MAX_BUFFER_SIZE, MPI_CHARACTER, &
         MPI_ANY_SOURCE, COAG_DIST_TAG_REQUEST_PARTICLE, MPI_COMM_WORLD, &
         status, ierr)
    call pmc_mpi_check_ierr(ierr)
    call assert(920139874, status(MPI_TAG) &
         == COAG_DIST_TAG_REQUEST_PARTICLE)
    call mpi_get_count(status, MPI_CHARACTER, buffer_size, ierr)
    call pmc_mpi_check_ierr(ierr)
    call assert(190658659, buffer_size <= COAG_DIST_MAX_BUFFER_SIZE)
    remote_proc = status(MPI_SOURCE)

    ! unpack it
    position = 0
    call pmc_mpi_unpack_integer(buffer, position, request_bin)
    call assert(895128380, position == buffer_size)

    ! send the particle back if we have one
    if (integer_varray_n_entry( &
         aero_state%aero_sorted%size_class%inverse(request_bin, s1)) == 0) then
       call send_return_no_particle(remote_proc, request_bin)
    else
       call aero_state_remove_rand_particle_from_bin(aero_state, &
            request_bin, s1, aero_particle)
       call send_return_req_particle(aero_particle, request_bin, &
            remote_proc)
    end if
#endif

  end subroutine recv_request_particle

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine send_return_no_particle(dest_proc, i_bin)

    !> Process number to send message to.
    integer, intent(in) :: dest_proc
    !> Bin number where there was no particle.
    integer, intent(in) :: i_bin

#ifdef PMC_USE_MPI
    character :: buffer(COAG_DIST_MAX_BUFFER_SIZE)
    integer :: buffer_size, max_buffer_size, position, ierr

    max_buffer_size = 0
    max_buffer_size = max_buffer_size + pmc_mpi_pack_size_integer(i_bin)
    call assert(744787119, max_buffer_size <= COAG_DIST_MAX_BUFFER_SIZE)
    position = 0
    call pmc_mpi_pack_integer(buffer, position, i_bin)
    call assert(445960340, position <= max_buffer_size)
    buffer_size = position
    call mpi_bsend(buffer, buffer_size, MPI_CHARACTER, dest_proc, &
         COAG_DIST_TAG_RETURN_NO_PARTICLE, MPI_COMM_WORLD, ierr)
    call pmc_mpi_check_ierr(ierr)
#endif

  end subroutine send_return_no_particle

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine recv_return_no_particle(requests, aero_data, aero_state)

    !> Array of outstanding requests.
    type(request_t), intent(inout) :: requests(COAG_DIST_MAX_REQUESTS)
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Aerosol state.
    type(aero_state_t), intent(inout) :: aero_state

#ifdef PMC_USE_MPI
    logical :: found_request
    integer :: buffer_size, position, sent_bin, sent_proc, i_req
    integer :: ierr, status(MPI_STATUS_SIZE)
    character :: buffer(COAG_DIST_MAX_BUFFER_SIZE)

    ! get the message
    call mpi_recv(buffer, COAG_DIST_MAX_BUFFER_SIZE, MPI_CHARACTER, &
         MPI_ANY_SOURCE, COAG_DIST_TAG_RETURN_NO_PARTICLE, &
         MPI_COMM_WORLD, status, ierr)
    call pmc_mpi_check_ierr(ierr)
    call assert(918153221, status(MPI_TAG) &
         == COAG_DIST_TAG_RETURN_NO_PARTICLE)
    call mpi_get_count(status, MPI_CHARACTER, buffer_size, ierr)
    call pmc_mpi_check_ierr(ierr)
    call assert(461111487, buffer_size <= COAG_DIST_MAX_BUFFER_SIZE)
    sent_proc = status(MPI_SOURCE)

    ! unpack it
    position = 0
    call pmc_mpi_unpack_integer(buffer, position, sent_bin)
    call assert(518172999, position == buffer_size)

    ! find the matching request
    found_request = .false.
    do i_req = 1,COAG_DIST_MAX_REQUESTS
       if ((requests(i_req)%remote_proc == sent_proc) &
            .and. (requests(i_req)%remote_bin == sent_bin)) then
          found_request = .true.
          exit
       end if
    end do
    call assert(215612776, found_request)

    ! We can't do coagulation with the local particle, so store it
    ! back. If we wanted to, we could use the knowledge that it should
    ! go into bin requests(i_req)%local_bin
    call aero_state_add_particle(aero_state, &
         requests(i_req)%local_aero_particle, aero_data, &
         allow_resort=.false.)
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
    !> Process number to send particle to.
    integer, intent(in) :: dest_proc

#ifdef PMC_USE_MPI
    character :: buffer(COAG_DIST_MAX_BUFFER_SIZE)
    integer :: buffer_size, max_buffer_size, position, ierr

    max_buffer_size = 0
    max_buffer_size = max_buffer_size + pmc_mpi_pack_size_integer(i_bin)
    max_buffer_size = max_buffer_size &
         + pmc_mpi_pack_size_aero_particle(aero_particle)
    call assert(496283814, max_buffer_size <= COAG_DIST_MAX_BUFFER_SIZE)
    position = 0
    call pmc_mpi_pack_integer(buffer, position, i_bin)
    call pmc_mpi_pack_aero_particle(buffer, position, aero_particle)
    call assert(263666386, position <= max_buffer_size)
    buffer_size = position
    call mpi_bsend(buffer, buffer_size, MPI_CHARACTER, dest_proc, &
         COAG_DIST_TAG_RETURN_REQ_PARTICLE, MPI_COMM_WORLD, ierr)
    call pmc_mpi_check_ierr(ierr)
#endif

  end subroutine send_return_req_particle

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine recv_return_req_particle(requests, env_state, aero_weight_total, &
       aero_data, aero_state, accept_factors, coag_kernel_type, tot_n_coag, &
       magnitudes)

    !> Array of outstanding requests.
    type(request_t), intent(inout) :: requests(COAG_DIST_MAX_REQUESTS)
    !> Environment state.
    type(env_state_t), intent(in) :: env_state
    !> Total weighting array.
    type(aero_weight_array_t), intent(in) :: aero_weight_total
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Aerosol state.
    type(aero_state_t), intent(inout) :: aero_state
    !> Accept scale factors per bin pair (1).
    real(kind=dp), intent(in) :: accept_factors(:,:)
    !> Coagulation kernel type.
    integer, intent(in) :: coag_kernel_type
    !> Number of coagulation events.
    integer, intent(inout) :: tot_n_coag
    !> Computational volumes on all processes.
    real(kind=dp), intent(in) :: magnitudes(:,:)

    integer, parameter :: s1 = 1
    integer, parameter :: s2 = 1
    integer, parameter :: sc = 1

#ifdef PMC_USE_MPI
    logical :: found_request, remove_1, remove_2
    integer :: buffer_size, position, sent_bin, sent_proc, i_req
    integer :: ierr, status(MPI_STATUS_SIZE)
    character :: buffer(COAG_DIST_MAX_BUFFER_SIZE)
    type(aero_particle_t) :: sent_aero_particle
    real(kind=dp) :: k, p

    ! get the message
    call mpi_recv(buffer, COAG_DIST_MAX_BUFFER_SIZE, MPI_CHARACTER, &
         MPI_ANY_SOURCE, COAG_DIST_TAG_RETURN_REQ_PARTICLE, &
         MPI_COMM_WORLD, status, ierr)
    call pmc_mpi_check_ierr(ierr)
    call assert(133285061, status(MPI_TAG) &
         == COAG_DIST_TAG_RETURN_REQ_PARTICLE)
    call mpi_get_count(status, MPI_CHARACTER, buffer_size, ierr)
    call pmc_mpi_check_ierr(ierr)
    call assert(563012836, buffer_size <= COAG_DIST_MAX_BUFFER_SIZE)
    sent_proc = status(MPI_SOURCE)

    ! unpack it
    position = 0
    call pmc_mpi_unpack_integer(buffer, position, sent_bin)
    call pmc_mpi_unpack_aero_particle(buffer, position, sent_aero_particle)
    call assert(753356021, position == buffer_size)

    ! find the matching request
    found_request = .false.
    do i_req = 1,COAG_DIST_MAX_REQUESTS
       if ((requests(i_req)%remote_proc == sent_proc) &
            .and. (requests(i_req)%remote_bin == sent_bin)) then
          found_request = .true.
          exit
       end if
    end do
    call assert(579308475, found_request)

    ! maybe do coagulation
    call num_conc_weighted_kernel(coag_kernel_type, &
         requests(i_req)%local_aero_particle, sent_aero_particle, &
         s1, s2, sc, aero_data, aero_weight_total, env_state, k)
    p = k * accept_factors(requests(i_req)%local_bin, sent_bin)

    if (pmc_random() .lt. p) then
       ! coagulation happened, do it
       tot_n_coag = tot_n_coag + 1
       call coagulate_dist(aero_data, aero_state, &
            requests(i_req)%local_aero_particle, sent_aero_particle, &
            sent_proc, aero_weight_total, magnitudes, remove_1, remove_2)
    else
       remove_1 = .false.
       remove_2 = .false.
    end if

    ! send the particles back
    if (.not. remove_1) then
       ! If we wanted to, we could use the knowledge that this will go
       ! into bin requests(i_req)%local_bin
       call aero_state_add_particle(aero_state, &
            requests(i_req)%local_aero_particle, aero_data, &
            allow_resort=.false.)
    end if
    if (.not. remove_2) then
       call send_return_unreq_particle(sent_aero_particle, sent_proc)
    end if

    call request_deallocate(requests(i_req))
    call request_allocate(requests(i_req))
#endif

  end subroutine recv_return_req_particle

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine send_return_unreq_particle(aero_particle, dest_proc)

    !> Aero particle to send.
    type(aero_particle_t), intent(in) :: aero_particle
    !> Process to send the particle to.
    integer, intent(in) :: dest_proc

#ifdef PMC_USE_MPI
    character :: buffer(COAG_DIST_MAX_BUFFER_SIZE)
    integer :: buffer_size, max_buffer_size, position, ierr

    max_buffer_size = 0
    max_buffer_size = max_buffer_size &
         + pmc_mpi_pack_size_aero_particle(aero_particle)
    call assert(414990602, max_buffer_size <= COAG_DIST_MAX_BUFFER_SIZE)
    position = 0
    call pmc_mpi_pack_aero_particle(buffer, position, aero_particle)
    call assert(898537822, position <= max_buffer_size)
    buffer_size = position
    call mpi_bsend(buffer, buffer_size, MPI_CHARACTER, dest_proc, &
         COAG_DIST_TAG_RETURN_UNREQ_PARTICLE, MPI_COMM_WORLD, ierr)
    call pmc_mpi_check_ierr(ierr)
#endif

  end subroutine send_return_unreq_particle

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine recv_return_unreq_particle(aero_state, aero_data)

    !> Aerosol state.
    type(aero_state_t), intent(inout) :: aero_state
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data

#ifdef PMC_USE_MPI
    logical :: found_request
    integer :: buffer_size, position, sent_proc, ierr
    character :: buffer(COAG_DIST_MAX_BUFFER_SIZE)
    type(aero_particle_t) :: aero_particle
    integer :: status(MPI_STATUS_SIZE), send_proc

    ! get the message
    call mpi_recv(buffer, COAG_DIST_MAX_BUFFER_SIZE, MPI_CHARACTER, &
         MPI_ANY_SOURCE, COAG_DIST_TAG_RETURN_UNREQ_PARTICLE, &
         MPI_COMM_WORLD, status, ierr)
    call pmc_mpi_check_ierr(ierr)
    call assert(496247788, status(MPI_TAG) &
         == COAG_DIST_TAG_RETURN_UNREQ_PARTICLE)
    call mpi_get_count(status, MPI_CHARACTER, buffer_size, ierr)
    call pmc_mpi_check_ierr(ierr)
    call assert(590644042, buffer_size <= COAG_DIST_MAX_BUFFER_SIZE)
    sent_proc = status(MPI_SOURCE)

    ! unpack it
    position = 0
    call pmc_mpi_unpack_aero_particle(buffer, position, aero_particle)
    call assert(833588594, position == buffer_size)

    ! put it back
    call aero_state_add_particle(aero_state, aero_particle, aero_data, &
         allow_resort=.false.)
#endif

  end subroutine recv_return_unreq_particle

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Send a message saying that this process is finished with its
  !> coagulation.
  subroutine send_done(dest_proc)

    !> Process to send the message to.
    integer, intent(in) :: dest_proc

#ifdef PMC_USE_MPI
    character :: buffer(COAG_DIST_MAX_BUFFER_SIZE)
    integer :: buffer_size, ierr

    buffer_size = 0
    call mpi_bsend(buffer, buffer_size, MPI_CHARACTER, dest_proc, &
         COAG_DIST_TAG_DONE, MPI_COMM_WORLD, ierr)
    call pmc_mpi_check_ierr(ierr)
#endif

  end subroutine send_done

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Receive a done message.
  subroutine recv_done(procs_done)

    !> Which processes are finished with coagulation.
    logical, intent(inout) :: procs_done(:)

#ifdef PMC_USE_MPI
    integer :: buffer_size, sent_proc, ierr
    character :: buffer(COAG_DIST_MAX_BUFFER_SIZE)
    integer :: status(MPI_STATUS_SIZE)

    ! get the message
    call mpi_recv(buffer, COAG_DIST_MAX_BUFFER_SIZE, MPI_CHARACTER, &
         MPI_ANY_SOURCE, COAG_DIST_TAG_DONE, MPI_COMM_WORLD, &
         status, ierr)
    call pmc_mpi_check_ierr(ierr)
    call assert(348737947, status(MPI_TAG) &
         == COAG_DIST_TAG_DONE)
    call mpi_get_count(status, MPI_CHARACTER, buffer_size, ierr)
    call pmc_mpi_check_ierr(ierr)
    call assert(214904056, buffer_size == 0)
    sent_proc = status(MPI_SOURCE)

    ! process it
    procs_done(sent_proc + 1) = .true.
#endif

  end subroutine recv_done

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> generate the number of samples to do per bin pair.
  subroutine generate_n_samps(n_parts, del_t, bin_grid, aero_weight_array, &
       k_max, n_samps, accept_factors)

    !> Number of particles per bin on all processes.
    integer, intent(in) :: n_parts(:,:)
    !> Timestep.
    real(kind=dp), intent(in) :: del_t
    !> Bin grid.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Weighting function array.
    type(aero_weight_array_t), intent(in) :: aero_weight_array
    !> Maximum kernel.
    real(kind=dp), intent(in) :: k_max(:,:)
    !> Number of samples to do per bin pair.
    integer, intent(out) :: n_samps(:,:)
    !> Accept scale factors per bin pair (1).
    real(kind=dp), intent(out) :: accept_factors(:,:)

    integer :: i_bin, j_bin, rank, n_bin
    real(kind=dp) :: n_samp_mean

    n_bin = size(k_max, 1)
    rank = pmc_mpi_rank()
    n_samps = 0
    do i_bin = 1,n_bin
       if (n_parts(i_bin, rank + 1) == 0) &
            cycle
       do j_bin = i_bin,n_bin
          call compute_n_samp(n_parts(i_bin, rank + 1), &
               sum(n_parts(j_bin, :)), (i_bin == j_bin), &
               k_max(i_bin, j_bin), del_t, n_samp_mean, &
               n_samps(i_bin, j_bin), accept_factors(i_bin, j_bin))
       end do
    end do

  end subroutine generate_n_samps

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine coagulate_dist(aero_data, aero_state, aero_particle_1, &
       aero_particle_2, remote_proc, aero_weight_total, magnitudes, &
       remove_1, remove_2)

    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Aerosol state.
    type(aero_state_t), intent(inout) :: aero_state
    !> First particle to coagulate.
    type(aero_particle_t), intent(in) :: aero_particle_1
    !> Second particle to coagulate.
    type(aero_particle_t), intent(in) :: aero_particle_2
    !> Remote process that the particle came from.
    integer, intent(in) :: remote_proc
    !> Total weight across all processes.
    type(aero_weight_array_t), intent(in) :: aero_weight_total
    !> Computational volumes on all processes (m^3).
    real(kind=dp), intent(in) :: magnitudes(:,:)
    !> Whether to remove aero_particle_1 after the coagulation.
    logical, intent(out) :: remove_1
    !> Whether to remove aero_particle_2 after the coagulation.
    logical, intent(out) :: remove_2

    integer, parameter :: s1 = 1
    integer, parameter :: s2 = 1
    integer, parameter :: sc = 1

    type(aero_particle_t) :: aero_particle_new
    integer :: new_proc, new_group
    type(aero_info_t) :: aero_info_1, aero_info_2
    logical :: create_new, id_1_lost, id_2_lost

    call coagulate_weighting(aero_particle_1, aero_particle_2, &
         aero_particle_new, s1, s2, sc, aero_data, aero_state%awa, &
         remove_1, remove_2, create_new, id_1_lost, id_2_lost, &
         aero_info_1, aero_info_2)

    if (id_1_lost) then
       call aero_info_array_add_aero_info(aero_state%aero_info_array, &
            aero_info_1)
    end if
    if (id_2_lost) then
       call aero_info_array_add_aero_info(aero_state%aero_info_array, &
            aero_info_2)
    end if

    ! add new particle
    if (create_new) then
       new_group = aero_weight_array_rand_group(aero_weight_total, sc, &
            aero_particle_radius(aero_particle_new, aero_data))
       aero_particle_new%weight_group = new_group
       new_proc = sample_cts_pdf(1d0 / magnitudes(new_group, :)) - 1
       call send_return_unreq_particle(aero_particle_new, new_proc)
    end if

  end subroutine coagulate_dist

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_coagulation_dist
