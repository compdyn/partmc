! Copyright (C) 2017-2018 Matt Dawson
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_test_emission program

!> Test of emission reaction module
program pmc_test_emission

  use pmc_util,                         only: i_kind, dp, assert, &
                                              almost_equal, string_t, &
                                              warn_msg
  use pmc_rxn_data
  use pmc_rxn_emission
  use pmc_rxn_factory
  use pmc_mechanism_data
  use pmc_chem_spec_data
  use pmc_camp_core
  use pmc_camp_state
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

  if (run_emission_tests()) then
    if (pmc_mpi_rank().eq.0) write(*,*) &
          "Emission reaction tests - PASS"
  else
    if (pmc_mpi_rank().eq.0) write(*,*) &
          "Emission reaction tests - FAIL"
    stop 3
  end if

  ! finalize mpi
  call pmc_mpi_finalize()

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Run all pmc_chem_mech_solver tests
  logical function run_emission_tests() result(passed)

    use pmc_camp_solver_data

    type(camp_solver_data_t), pointer :: camp_solver_data

    camp_solver_data => camp_solver_data_t()

    if (camp_solver_data%is_solver_available()) then
      passed = run_emission_test()
    else
      call warn_msg(599253453, "No solver available")
      passed = .true.
    end if

    deallocate(camp_solver_data)

  end function run_emission_tests

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Solve a mechanism of consecutive reactions
  !!
  !! The mechanism is of the form:
  !!
  !!   -k1-> A
  !!   -k2-> B
  !!
  !! where k1 and k2 are emission reaction rate constants.
  logical function run_emission_test()

    use pmc_constants

    type(camp_core_t), pointer :: camp_core
    type(camp_state_t), pointer :: camp_state
    character(len=:), allocatable :: input_file_path, key, str_val
    type(string_t), allocatable, dimension(:) :: output_file_path

    real(kind=dp), dimension(0:NUM_TIME_STEP, 2) :: model_conc, true_conc
    integer(kind=i_kind) :: idx_A, idx_B, i_time, i_spec, i_rxn
    integer(kind=i_kind) ::  i_mech_rxn_A, i_mech_rxn_B
    real(kind=dp) :: time_step, time, k1, k2, temp, pressure, rate_A, rate_B
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
    type(rxn_update_data_emission_t) :: rate_update_A, rate_update_B

    run_emission_test = .true.

    ! Set the rate constants (for calculating the true value)
    temp = 272.5d0
    pressure = 101253.3d0
    rate_A = 0.954d0
    rate_B = 1.0d-2
    k1 = rate_A
    k2 = rate_B * 12.3d0

    ! Set output time step (s)
    time_step = 1.0

#ifdef PMC_USE_MPI
    ! Load the model data on the root process and pass it to process 1 for solving
    if (pmc_mpi_rank().eq.0) then
#endif

      ! Get the emission reaction mechanism json file
      input_file_path = 'test_emission_config.json'

      ! Construct a camp_core variable
      camp_core => camp_core_t(input_file_path)

      deallocate(input_file_path)

      ! Initialize the model
      call camp_core%initialize()

      ! Find the mechanism
      key = "emissions"
      call assert(260845179, camp_core%get_mechanism(key, mechanism))

      ! Find the A emission reaction
      key = "rxn id"
      i_mech_rxn_A = 0
      i_mech_rxn_B = 0
      do i_rxn = 1, mechanism%size()
        rxn => mechanism%get_rxn(i_rxn)
        if (rxn%property_set%get_string(key, str_val)) then
          if (trim(str_val).eq."rxn A") then
            i_mech_rxn_A = i_rxn
            select type (rxn_emis => rxn)
              class is (rxn_emission_t)
                call rxn_factory%initialize_update_data(rxn_emis, &
                                                        rate_update_A)
            end select
          end if
          if (trim(str_val).eq."rxn B") then
            i_mech_rxn_B = i_rxn
            select type (rxn_emis => rxn)
              class is (rxn_emission_t)
                call rxn_factory%initialize_update_data(rxn_emis, &
                                                        rate_update_B)
            end select
          end if
        end if
      end do
      call assert(262750713, i_mech_rxn_A.eq.1)
      call assert(508177387, i_mech_rxn_B.eq.2)

      ! Get the chemical species data
      call assert(310060658, camp_core%get_chem_spec_data(chem_spec_data))

      ! Get species indices
      key = "A"
      idx_A = chem_spec_data%gas_state_id(key);
      key = "B"
      idx_B = chem_spec_data%gas_state_id(key);

      ! Make sure the expected species are in the model
      call assert(369804751, idx_A.gt.0)
      call assert(199647847, idx_B.gt.0)

#ifdef PMC_USE_MPI
      ! pack the camp core
      pack_size = camp_core%pack_size() &
                + rate_update_A%pack_size() &
                + rate_update_B%pack_size()
      allocate(buffer(pack_size))
      pos = 0
      call camp_core%bin_pack(buffer, pos)
      call rate_update_A%bin_pack(buffer, pos)
      call rate_update_B%bin_pack(buffer, pos)
      call assert(194383540, pos.eq.pack_size)
    end if

    ! broadcast the species ids
    call pmc_mpi_bcast_integer(idx_A)
    call pmc_mpi_bcast_integer(idx_B)

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
      camp_core => camp_core_t()
      pos = 0
      call camp_core%bin_unpack(buffer, pos)
      call rate_update_A%bin_unpack(buffer, pos)
      call rate_update_B%bin_unpack(buffer, pos)
      call assert(785350501, pos.eq.pack_size)
      allocate(buffer_copy(pack_size))
      pos = 0
      call camp_core%bin_pack(buffer_copy, pos)
      call rate_update_A%bin_pack(buffer_copy, pos)
      call rate_update_B%bin_pack(buffer_copy, pos)
      call assert(215135696, pos.eq.pack_size)
      do i_elem = 1, pack_size
        call assert_msg(819078131, buffer(i_elem).eq.buffer_copy(i_elem), &
                "Mismatch in element: "//trim(to_string(i_elem)))
      end do

      ! solve and evaluate results on process 1
#endif

      ! Initialize the solver
      call camp_core%solver_initialize()

      ! Get a model state variable
      camp_state => camp_core%new_state()

      ! Set the environmental conditions
      call camp_state%env_states(1)%set_temperature_K(   temp )
      call camp_state%env_states(1)%set_pressure_Pa( pressure )

      ! Save the initial concentrations
      true_conc(0,idx_A) = 1.0
      true_conc(0,idx_B) = 1.0
      model_conc(0,:) = true_conc(0,:)

      ! Set the initial concentrations in the model
      camp_state%state_var(:) = model_conc(0,:)

      ! Set the rxn rates
      call rate_update_A%set_rate(rate_A)
      call rate_update_B%set_rate(42.11d0)
      call camp_core%update_rxn_data(rate_update_A)
      call camp_core%update_rxn_data(rate_update_B)

      ! Test re-setting of the rxn B rate
      call rate_update_B%set_rate(rate_B)
      call camp_core%update_rxn_data(rate_update_B)

#ifdef PMC_DEBUG
      ! Evaluate the Jacobian during solving
      solver_stats%eval_Jac = .true.
#endif

      ! Integrate the mechanism
      do i_time = 1, NUM_TIME_STEP

        ! Get the modeled conc
        call camp_core%solve(camp_state, time_step, &
                              solver_stats = solver_stats)
        model_conc(i_time,:) = camp_state%state_var(:)

#ifdef PMC_DEBUG
        ! Check the Jacobian evaluations
        call assert_msg(823507267, solver_stats%Jac_eval_fails.eq.0, &
                        trim( to_string( solver_stats%Jac_eval_fails ) )// &
                        " Jacobian evaluation failures at time step "// &
                        trim( to_string( i_time ) ) )
#endif

        ! Get the analytic conc
        time = i_time * time_step
        true_conc(i_time,idx_A) = true_conc(0,idx_A) + (k1)*time
        true_conc(i_time,idx_B) = true_conc(0,idx_B) + (k2)*time

      end do

      ! Save the results
      open(unit=7, file="out/emission_results.txt", status="replace", &
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
          call assert_msg(870199144, &
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

      deallocate(camp_state)

#ifdef PMC_USE_MPI
      ! convert the results to an integer
      if (run_emission_test) then
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
        run_emission_test = .true.
      else
        run_emission_test = .false.
      end if
    end if

    deallocate(buffer)
#endif

    deallocate(camp_core)

  end function run_emission_test

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program pmc_test_emission
