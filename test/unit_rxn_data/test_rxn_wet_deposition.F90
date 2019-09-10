! Copyright (C) 2017-2018 Matt Dawson
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_test_wet_deposition program

!> Test of wet_deposition reaction module
program pmc_test_wet_deposition

  use pmc_aero_rep_data
  use pmc_chem_spec_data
  use pmc_mechanism_data
  use pmc_camp_core
  use pmc_camp_state
  use pmc_rxn_data
  use pmc_rxn_wet_deposition
  use pmc_rxn_factory
  use pmc_solver_stats
  use pmc_util,                         only: i_kind, dp, assert, &
                                              almost_equal, string_t, &
                                              warn_msg
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

  if (run_wet_deposition_tests()) then
    if (pmc_mpi_rank().eq.0) write(*,*) &
          "Wet Deposition reaction tests - PASS"
  else
    if (pmc_mpi_rank().eq.0) write(*,*) &
          "Wet Deposition reaction tests - FAIL"
    stop 3
  end if

  ! finalize mpi
  call pmc_mpi_finalize()

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Run all pmc_rxn_wet_deposition tests
  logical function run_wet_deposition_tests() result(passed)

    use pmc_camp_solver_data

    type(camp_solver_data_t), pointer :: camp_solver_data

    camp_solver_data => camp_solver_data_t()

    if (camp_solver_data%is_solver_available()) then
      passed = run_wet_deposition_test()
    else
      call warn_msg(280044966, "No solver available")
      passed = .true.
    end if

    deallocate(camp_solver_data)

  end function run_wet_deposition_tests

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Solve a mechanism of wet deposition reactions
  !!
  !! The mechanism is of the form:
  !!
  !!   A -k->
  !!   B -k->
  !!
  !! where k is the wet_deposition reaction rate constant.
  logical function run_wet_deposition_test()

    use pmc_constants

    type(camp_core_t), pointer :: camp_core
    type(camp_state_t), pointer :: camp_state
    character(len=:), allocatable :: input_file_path, key, str_val
    type(string_t), allocatable, dimension(:) :: output_file_path

    real(kind=dp), dimension(0:NUM_TIME_STEP, 8) :: model_conc, true_conc
    integer(kind=i_kind) :: idx_1RA, idx_1RB, idx_1CB, idx_1CC
    integer(kind=i_kind) :: idx_2RA, idx_2RB, idx_2CB, idx_2CC
    integer(kind=i_kind) :: i_time, i_spec, i_rxn
    integer(kind=i_kind) :: i_mech_rxn_rain, i_mech_rxn_cloud
    real(kind=dp) :: time_step, time, k_rain, k_cloud, temp, pressure, &
                     rate_rain, rate_cloud
    class(rxn_data_t), pointer :: rxn
    class(aero_rep_data_t), pointer :: aero_rep
    type(string_t), allocatable :: unique_names(:)
    character(len=:), allocatable :: spec_name, phase_name, rep_name
#ifdef PMC_USE_MPI
    character, allocatable :: buffer(:), buffer_copy(:)
    integer(kind=i_kind) :: pack_size, pos, i_elem, results
#endif

    type(solver_stats_t), target :: solver_stats

    ! For setting rates
    type(mechanism_data_t), pointer :: mechanism
    type(rxn_factory_t) :: rxn_factory
    type(rxn_update_data_wet_deposition_t) :: rate_update_rain
    type(rxn_update_data_wet_deposition_t) :: rate_update_cloud

    run_wet_deposition_test = .true.

    ! Set the rate constants (for calculating the true value)
    temp = 272.5d0
    pressure = 101253.3d0
    rate_rain = 0.954d0
    rate_cloud = 1.0d-2
    k_rain = rate_rain
    k_cloud = rate_cloud * 12.3d0

    ! Set output time step (s)
    time_step = 1.0

#ifdef PMC_USE_MPI
    ! Load the model data on the root process and pass it to process 1 for solving
    if (pmc_mpi_rank().eq.0) then
#endif

      ! Get the wet_deposition reaction mechanism json file
      input_file_path = 'test_wet_deposition_config.json'

      ! Construct a camp_core variable
      camp_core => camp_core_t(input_file_path)

      deallocate(input_file_path)

      ! Initialize the model
      call camp_core%initialize()

      ! Find the mechanism
      key = "wet deposition"
      call assert(868882379, camp_core%get_mechanism(key, mechanism))

      ! Find the A wet_deposition reaction
      key = "rxn id"
      i_mech_rxn_rain = 0
      i_mech_rxn_cloud = 0
      do i_rxn = 1, mechanism%size()
        rxn => mechanism%get_rxn(i_rxn)
        if (rxn%property_set%get_string(key, str_val)) then
          if (trim(str_val).eq."rxn rain") then
            i_mech_rxn_rain = i_rxn
            select type (rxn_loss => rxn)
              class is (rxn_wet_deposition_t)
                call rxn_factory%initialize_update_data(rxn_loss, &
                                                        rate_update_rain)
            end select
          end if
          if (trim(str_val).eq."rxn cloud") then
            i_mech_rxn_cloud = i_rxn
            select type (rxn_loss => rxn)
              class is (rxn_wet_deposition_t)
                call rxn_factory%initialize_update_data(rxn_loss, &
                                                        rate_update_cloud)
            end select
          end if
        end if
      end do
      call assert(300573108, i_mech_rxn_rain.eq.1)
      call assert(222343365, i_mech_rxn_cloud.eq.2)

      ! Find species in the first aerosol representation
      rep_name = "my first particle"
      call assert(229751752, camp_core%get_aero_rep(rep_name, aero_rep))

      ! Get species indices
      phase_name = "rain"
      spec_name = "A"
      unique_names = aero_rep%unique_names( phase_name = phase_name, &
                                            spec_name = spec_name )
      call assert(635673803, size(unique_names).eq.1)
      idx_1RA = aero_rep%spec_state_id(unique_names(1)%string);

      spec_name = "B"
      unique_names = aero_rep%unique_names( phase_name = phase_name, &
                                            spec_name = spec_name )
      call assert(565853391, size(unique_names).eq.1)
      idx_1RB = aero_rep%spec_state_id(unique_names(1)%string);

      phase_name = "cloud"
      spec_name = "B"
      unique_names = aero_rep%unique_names( phase_name = phase_name, &
                                            spec_name = spec_name )
      call assert(957288212, size(unique_names).eq.1)
      idx_1CB = aero_rep%spec_state_id(unique_names(1)%string);

      spec_name = "C"
      unique_names = aero_rep%unique_names( phase_name = phase_name, &
                                            spec_name = spec_name )
      call assert(388978941, size(unique_names).eq.1)
      idx_1CC = aero_rep%spec_state_id(unique_names(1)%string);

      ! Find species in the second aerosol representation
      rep_name = "my second particle"
      call assert(440099954, camp_core%get_aero_rep(rep_name, aero_rep))

      phase_name = "rain"
      spec_name = "A"
      unique_names = aero_rep%unique_names( phase_name = phase_name, &
                                            spec_name = spec_name )
      call assert(271848584, size(unique_names).eq.1)
      idx_2RA = aero_rep%spec_state_id(unique_names(1)%string);

      spec_name = "B"
      unique_names = aero_rep%unique_names( phase_name = phase_name, &
                                            spec_name = spec_name )
      call assert(384166929, size(unique_names).eq.1)
      idx_2RB = aero_rep%spec_state_id(unique_names(1)%string);

      phase_name = "cloud"
      spec_name = "B"
      unique_names = aero_rep%unique_names( phase_name = phase_name, &
                                            spec_name = spec_name )
      call assert(778960523, size(unique_names).eq.1)
      idx_2CB = aero_rep%spec_state_id(unique_names(1)%string);

      spec_name = "C"
      unique_names = aero_rep%unique_names( phase_name = phase_name, &
                                            spec_name = spec_name )
      call assert(273754118, size(unique_names).eq.1)
      idx_2CC = aero_rep%spec_state_id(unique_names(1)%string);

      ! Make sure the expected species are in the model
      call assert(100238441, idx_1RA.gt.0)
      call assert(830081536, idx_1RB.gt.0)
      call assert(659924632, idx_1CB.gt.0)
      call assert(489767728, idx_1CC.gt.0)
      call assert(319610824, idx_2RA.gt.0)
      call assert(149453920, idx_2RB.gt.0)
      call assert(879297015, idx_2CB.gt.0)
      call assert(709140111, idx_2CC.gt.0)

#ifdef PMC_USE_MPI
      ! pack the camp core
      pack_size = camp_core%pack_size()
      allocate(buffer(pack_size))
      pos = 0
      call camp_core%bin_pack(buffer, pos)
      call assert(768884204, pos.eq.pack_size)
    end if

    ! broadcast the species ids
    call pmc_mpi_bcast_integer(idx_1RA)
    call pmc_mpi_bcast_integer(idx_1RB)
    call pmc_mpi_bcast_integer(idx_1CB)
    call pmc_mpi_bcast_integer(idx_1CC)
    call pmc_mpi_bcast_integer(idx_2RA)
    call pmc_mpi_bcast_integer(idx_2RB)
    call pmc_mpi_bcast_integer(idx_2CB)
    call pmc_mpi_bcast_integer(idx_2CC)
    call pmc_mpi_bcast_integer(i_rxn_rain)

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
      call assert(413229770, pos.eq.pack_size)
      allocate(buffer_copy(pack_size))
      pos = 0
      call camp_core%bin_pack(buffer_copy, pos)
      call assert(243072866, pos.eq.pack_size)
      do i_elem = 1, pack_size
        call assert_msg(809928898, buffer(i_elem).eq.buffer_copy(i_elem), &
                "Mismatch in element: "//trim(to_string(i_elem)))
      end do

      ! solve and evaluate results on process 1
#endif

      ! Initialize the solver
      call camp_core%solver_initialize()

      ! Get a model state variable
      camp_state => camp_core%new_state()

      ! Set the environmental conditions
      camp_state%env_state%temp = temp
      camp_state%env_state%pressure = pressure
      call camp_state%update_env_state()

      ! Save the initial concentrations
      true_conc(0,idx_1RA) = 1.0
      true_conc(0,idx_1RB) = 2.2
      true_conc(0,idx_1CB) = 1.7
      true_conc(0,idx_1CC) = 2.4
      true_conc(0,idx_2RA) = 2.5
      true_conc(0,idx_2RB) = 0.7
      true_conc(0,idx_2CB) = 1.6
      true_conc(0,idx_2CC) = 1.9
      model_conc(0,:) = true_conc(0,:)

      ! Set the initial concentrations in the model
      camp_state%state_var(:) = model_conc(0,:)

      ! Set the rain rxn rate
      call rate_update_rain%set_rate(rate_rain)
      call rate_update_cloud%set_rate(43912.5d0)
      call camp_core%update_rxn_data(rate_update_rain)
      call camp_core%update_rxn_data(rate_update_cloud)

      ! Test re-setting of the rxn B rate
      call rate_update_cloud%set_rate(rate_cloud)
      call camp_core%update_rxn_data(rate_update_cloud)

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
        call assert_msg(172394787, solver_stats%Jac_eval_fails.eq.0, &
                        trim( to_string( solver_stats%Jac_eval_fails ) )// &
                        " Jacobian evaluation failures at time step "// &
                        trim( to_string( i_time ) ) )
#endif

        ! Get the analytic conc
        time = i_time * time_step
        true_conc(i_time,idx_1RA) = true_conc(0,idx_1RA) * exp(-(k_rain)*time)
        true_conc(i_time,idx_1RB) = true_conc(0,idx_1RB) * exp(-(k_rain)*time)
        true_conc(i_time,idx_1CB) = true_conc(0,idx_1CB) * exp(-(k_cloud)*time)
        true_conc(i_time,idx_1CC) = true_conc(0,idx_1CC) * exp(-(k_cloud)*time)
        true_conc(i_time,idx_2RA) = true_conc(0,idx_2RA) * exp(-(k_rain)*time)
        true_conc(i_time,idx_2RB) = true_conc(0,idx_2RB) * exp(-(k_rain)*time)
        true_conc(i_time,idx_2CB) = true_conc(0,idx_2CB) * exp(-(k_cloud)*time)
        true_conc(i_time,idx_2CC) = true_conc(0,idx_2CC) * exp(-(k_cloud)*time)

      end do

      ! Save the results
      open(unit=7, file="out/wet_deposition_results.txt", status="replace", &
              action="write")
      do i_time = 0, NUM_TIME_STEP
        write(7,*) i_time*time_step, &
              ' ', true_conc(i_time, idx_1RA),' ', model_conc(i_time, idx_1RA), &
              ' ', true_conc(i_time, idx_1RB),' ', model_conc(i_time, idx_1RB), &
              ' ', true_conc(i_time, idx_1CB),' ', model_conc(i_time, idx_1CB), &
              ' ', true_conc(i_time, idx_1CC),' ', model_conc(i_time, idx_1CC), &
              ' ', true_conc(i_time, idx_2RA),' ', model_conc(i_time, idx_2RA), &
              ' ', true_conc(i_time, idx_2RB),' ', model_conc(i_time, idx_2RB), &
              ' ', true_conc(i_time, idx_2CB),' ', model_conc(i_time, idx_2CB), &
              ' ', true_conc(i_time, idx_2CC),' ', model_conc(i_time, idx_2CC)
      end do
      close(7)

      ! Analyze the results
      do i_time = 1, NUM_TIME_STEP
        do i_spec = 1, size(model_conc, 2)
          call assert_msg(281211082, &
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
      if (run_wet_deposition_test) then
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
        run_wet_deposition_test = .true.
      else
        run_wet_deposition_test = .false.
      end if
    end if

    deallocate(buffer)
#endif

    deallocate(camp_core)

  end function run_wet_deposition_test

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program pmc_test_wet_deposition
