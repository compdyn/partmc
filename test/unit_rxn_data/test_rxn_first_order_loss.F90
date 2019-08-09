! Copyright (C) 2017-2018 Matt Dawson
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_test_first_order_loss program

!> Test of first_order_loss reaction module
program pmc_test_first_order_loss

  use pmc_util,                         only: i_kind, dp, assert, &
                                              almost_equal, string_t, &
                                              warn_msg
  use pmc_rxn_data
  use pmc_rxn_first_order_loss
  use pmc_rxn_factory
  use pmc_mechanism_data
  use pmc_chem_spec_data
  use pmc_phlex_core
  use pmc_phlex_state
  use pmc_solver_stats
#ifdef PMC_USE_JSON
  use json_module
#endif
  use pmc_mpi

  use iso_c_binding

  implicit none

  ! Number of timesteps to output in mechanisms
  integer(kind=i_kind) :: NUM_TIME_STEP = 100

  ! initialize mpi
  call pmc_mpi_init()

  if (run_first_order_loss_tests()) then
    if (pmc_mpi_rank().eq.0) write(*,*) &
          "First-Order loss reaction tests - PASS"
  else
    if (pmc_mpi_rank().eq.0) write(*,*) &
          "First-Order loss reaction tests - FAIL"
    stop 3
  end if

  ! finalize mpi
  call pmc_mpi_finalize()

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Run all pmc_chem_mech_solver tests
  logical function run_first_order_loss_tests() result(passed)

    use pmc_phlex_solver_data

    type(phlex_solver_data_t), pointer :: phlex_solver_data

    phlex_solver_data => phlex_solver_data_t()

    if (phlex_solver_data%is_solver_available()) then
      passed = run_first_order_loss_test()
    else
      call warn_msg(366534359, "No solver available")
      passed = .true.
    end if

    deallocate(phlex_solver_data)

  end function run_first_order_loss_tests

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Solve a mechanism of consecutive reactions
  !!
  !! The mechanism is of the form:
  !!
  !!   A -k1->
  !!   B -k2->
  !!
  !! where k1 and k2 are first-order loss reaction rate constants.
  logical function run_first_order_loss_test()

    use pmc_constants

    type(phlex_core_t), pointer :: phlex_core
    type(phlex_state_t), pointer :: phlex_state
    character(len=:), allocatable :: input_file_path, key, str_val
    type(string_t), allocatable, dimension(:) :: output_file_path

    real(kind=dp), dimension(0:NUM_TIME_STEP, 2) :: model_conc, true_conc
    integer(kind=i_kind) :: idx_A, idx_B, i_time, i_spec, i_rxn, i_rxn_A, &
                            i_mech_rxn_A
    real(kind=dp) :: time_step, time, k1, k2, temp, pressure, rate_1
    type(chem_spec_data_t), pointer :: chem_spec_data
    class(rxn_data_t), pointer :: rxn
#ifdef PMC_USE_MPI
    character, allocatable :: buffer(:), buffer_copy(:)
    integer(kind=i_kind) :: pack_size, pos, i_elem, results
#endif

    type(solver_stats_t), target :: solver_stats

    ! For setting rates
    type(mechanism_data_t), pointer :: mechanism
    type(rxn_factory_t) :: rxn_factory
    type(rxn_update_data_first_order_loss_rate_t) :: rate_update

    run_first_order_loss_test = .true.

    ! Set the rate constants (for calculating the true value)
    temp = 272.5d0
    pressure = 101253.3d0
    rate_1 = 0.954d0
    k1 = rate_1
    k2 = 1.0d-02 * 12.3d0

    ! Set output time step (s)
    time_step = 1.0

#ifdef PMC_USE_MPI
    ! Load the model data on the root process and pass it to process 1 for solving
    if (pmc_mpi_rank().eq.0) then
#endif

      ! Get the first_order_loss reaction mechanism json file
      input_file_path = 'test_first_order_loss_config.json'

      ! Construct a phlex_core variable
      phlex_core => phlex_core_t(input_file_path)

      deallocate(input_file_path)

      ! Initialize the model
      call phlex_core%initialize()

      ! Find the mechanism
      key = "first order loss"
      call assert(158978657, phlex_core%get_mechanism(key, mechanism))

      ! Find the A loss reaction
      key = "rxn id"
      i_rxn_A = 342
      i_mech_rxn_A = 0
      do i_rxn = 1, mechanism%size()
        rxn => mechanism%get_rxn(i_rxn)
        if (rxn%property_set%get_string(key, str_val)) then
          if (trim(str_val).eq."rxn A") then
            i_mech_rxn_A = i_rxn
            select type (rxn_loss => rxn)
              class is (rxn_first_order_loss_t)
                call rxn_loss%set_rxn_id(i_rxn_A)
            end select
          end if
        end if
      end do
      call assert(873933421, i_mech_rxn_A.eq.1)

      ! Get the chemical species data
      call assert(533619613, phlex_core%get_chem_spec_data(chem_spec_data))

      ! Get species indices
      key = "A"
      idx_A = chem_spec_data%gas_state_id(key);
      key = "B"
      idx_B = chem_spec_data%gas_state_id(key);

      ! Make sure the expected species are in the model
      call assert(870574648, idx_A.gt.0)
      call assert(982892993, idx_B.gt.0)

#ifdef PMC_USE_MPI
      ! pack the phlex core
      pack_size = phlex_core%pack_size()
      allocate(buffer(pack_size))
      pos = 0
      call phlex_core%bin_pack(buffer, pos)
      call assert(642579185, pos.eq.pack_size)
    end if

    ! broadcast the species ids
    call pmc_mpi_bcast_integer(idx_A)
    call pmc_mpi_bcast_integer(idx_B)
    call pmc_mpi_bcast_integer(i_rxn_A)

    ! broadcast the buffer size
    call pmc_mpi_bcast_integer(pack_size)

    if (pmc_mpi_rank().eq.1) then
      ! allocate the buffer to receive data
      allocate(buffer(pack_size))
    end if

    ! broadcast the data
    call pmc_mpi_bcast_packed(buffer)

    if (pmc_mpi_rank().eq.1) then
      ! unpack the data
      phlex_core => phlex_core_t()
      pos = 0
      call phlex_core%bin_unpack(buffer, pos)
      call assert(320136977, pos.eq.pack_size)
      allocate(buffer_copy(pack_size))
      pos = 0
      call phlex_core%bin_pack(buffer_copy, pos)
      call assert(356745163, pos.eq.pack_size)
      do i_elem = 1, pack_size
        call assert_msg(981439754, buffer(i_elem).eq.buffer_copy(i_elem), &
                "Mismatch in element: "//trim(to_string(i_elem)))
      end do

      ! solve and evaluate results on process 1
#endif

      ! Initialize the solver
      call phlex_core%solver_initialize()

      ! Get a model state variable
      phlex_state => phlex_core%new_state()

      ! Set the environmental conditions
      phlex_state%env_state%temp = temp
      phlex_state%env_state%pressure = pressure
      call phlex_state%update_env_state()

      ! Save the initial concentrations
      true_conc(0,idx_A) = 1.0
      true_conc(0,idx_B) = 1.0
      model_conc(0,:) = true_conc(0,:)

      ! Set the initial concentrations in the model
      phlex_state%state_var(:) = model_conc(0,:)

      ! Set the rxn B rate
      call rxn_factory%initialize_update_data(rate_update)
      call rate_update%set_rate(i_rxn_A, rate_1)
      call phlex_core%update_rxn_data(rate_update)

#ifdef PMC_DEBUG
      ! Evaluate the Jacobian during solving
      solver_stats%eval_Jac = .true.
#endif

      ! Integrate the mechanism
      do i_time = 1, NUM_TIME_STEP

        ! Get the modeled conc
        call phlex_core%solve(phlex_state, time_step, &
                              solver_stats = solver_stats)
        model_conc(i_time,:) = phlex_state%state_var(:)

#ifdef PMC_DEBUG
        ! Check the Jacobian evaluations
        call assert_msg(748874843, solver_stats%Jac_eval_fails.eq.0, &
                        trim( to_string( solver_stats%Jac_eval_fails ) )// &
                        " Jacobian evaluation failures at time step "// &
                        trim( to_string( i_time ) ) )
#endif

        ! Get the analytic conc
        time = i_time * time_step
        true_conc(i_time,idx_A) = true_conc(0,idx_A) * exp(-(k1)*time)
        true_conc(i_time,idx_B) = true_conc(0,idx_B) * exp(-(k2)*time)

      end do

      ! Save the results
      open(unit=7, file="out/first_order_loss_results.txt", status="replace", &
              action="write")
      do i_time = 0, NUM_TIME_STEP
        write(7,*) i_time*time_step, &
              ' ', true_conc(i_time, idx_A),' ', model_conc(i_time, idx_A), &
              ' ', true_conc(i_time, idx_B),' ', model_conc(i_time, idx_B)
      end do
      close(7)

      ! Analyze the results
      do i_time = 1, NUM_TIME_STEP
        do i_spec = 1, size(model_conc, 2)
          call assert_msg(704515935, &
            almost_equal(model_conc(i_time, i_spec), &
            true_conc(i_time, i_spec), real(1.0e-2, kind=dp)).or. &
            (model_conc(i_time, i_spec).lt.1e-5*model_conc(1, i_spec).and. &
            true_conc(i_time, i_spec).lt.1e-5*true_conc(1, i_spec)), &
            "time: "//trim(to_string(i_time))//"; species: "// &
            trim(to_string(i_spec))//"; mod: "// &
            trim(to_string(model_conc(i_time, i_spec)))//"; true: "// &
            trim(to_string(true_conc(i_time, i_spec))))
        end do
      end do

      deallocate(phlex_state)

#ifdef PMC_USE_MPI
      ! convert the results to an integer
      if (run_first_order_loss_test) then
        results = 0
      else
        results = 1
      end if
    end if
    
    ! Send the results back to the primary process
    call pmc_mpi_transfer_integer(results, results, 1, 0)

    ! convert the results back to a logical value
    if (pmc_mpi_rank().eq.0) then
      if (results.eq.0) then
        run_first_order_loss_test = .true.
      else
        run_first_order_loss_test = .false.
      end if
    end if

    deallocate(buffer)
#endif

    deallocate(phlex_core)

  end function run_first_order_loss_test

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program pmc_test_first_order_loss
