! Copyright (C) 2017 Matt Dawson
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_test_CMAQ_OH_HNO3 program

!> Test of CMAQ_OH_HNO3 reaction module
program pmc_test_CMAQ_OH_HNO3

  use pmc_util,                         only: i_kind, dp, assert, &
                                              almost_equal, string_t, &
                                              warn_msg
  use pmc_camp_core
  use pmc_camp_state
  use pmc_chem_spec_data
  use pmc_solver_stats
#ifdef PMC_USE_JSON
  use json_module
#endif
  use pmc_mpi

  implicit none

  ! Number of timesteps to output in mechanisms
  integer(kind=i_kind) :: NUM_TIME_STEP = 100

  ! initialize mpi
  call pmc_mpi_init()

  if (run_CMAQ_OH_HNO3_tests()) then
    if (pmc_mpi_rank().eq.0) write(*,*) "CMAQ_OH_HNO3 reaction tests - PASS"
  else
    if (pmc_mpi_rank().eq.0) write(*,*) "CMAQ_OH_HNO3 reaction tests - FAIL"
    stop 3
  end if

  ! finalize mpi
  call pmc_mpi_finalize()

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Run all pmc_chem_mech_solver tests
  logical function run_CMAQ_OH_HNO3_tests() result(passed)

    use pmc_camp_solver_data

    type(camp_solver_data_t), pointer :: camp_solver_data

    camp_solver_data => camp_solver_data_t()

    if (camp_solver_data%is_solver_available()) then
      passed = run_CMAQ_OH_HNO3_test()
    else
      call warn_msg(196869173, "No solver available")
      passed = .true.
    end if

    deallocate(camp_solver_data)

  end function run_CMAQ_OH_HNO3_tests

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Solve a mechanism of consecutive reactions
  !!
  !! The mechanism is of the form:
  !!
  !!   A -k1-> B -k2-> C
  !!
  !! where k1 and k2 are CMAQ_OH_HNO3 reaction rate constants.
  logical function run_CMAQ_OH_HNO3_test()

    use pmc_constants

    type(camp_core_t), pointer :: camp_core
    type(camp_state_t), pointer :: camp_state
    character(len=:), allocatable :: input_file_path
    type(string_t), allocatable, dimension(:) :: output_file_path

    type(chem_spec_data_t), pointer :: chem_spec_data
    real(kind=dp), dimension(0:NUM_TIME_STEP, 3) :: model_conc, true_conc
    integer(kind=i_kind) :: idx_A, idx_B, idx_C
    character(len=:), allocatable :: key
    integer(kind=i_kind) :: i_time, i_spec
    real(kind=dp) :: time_step, time
#ifdef PMC_USE_MPI
    character, allocatable :: buffer(:), buffer_copy(:)
    integer(kind=i_kind) :: pack_size, pos, i_elem, results
#endif

    type(solver_stats_t), target :: solver_stats

    ! Parameters for calculating true concentrations
    real(kind=dp) :: k1, k2, air_conc, temp, pressure, k_0, k_2, k_3, conv

    run_CMAQ_OH_HNO3_test = .true.

    ! Set the rate constants (for calculating the true value)
    temp = 272.5d0
    pressure = 101253.3d0
    air_conc = 1.0e6
    conv = const%avagadro / const%univ_gas_const * 10.0d0**(-12.0d0) * &
            pressure / temp
    k_0 = 1.0d0
    k_2 = 1.0d0
    k_3 = 4.0e-21
    k1 = (k_0 + (k_3*air_conc*conv) / (1 + k_3*air_conc*conv/k_2))
    k_0 = 1476.0d0 * exp(-398.0d0/temp) * (temp/300.0d0)**(60.0d0)
    k_2 = 1350.0d0 * exp(-450.0d0/temp) * (temp/300.0d0)**(58.0d0)
    k_3 = 4.0e-20 * exp(-2.0d0/temp) * (temp/300.0d0)**(1.5d0)
    k2 = (k_0 + (k_3*air_conc*conv) / (1 + k_3*air_conc*conv/k_2)) / 60.0d0

    ! Set output time step (s)
    time_step = 1.0

#ifdef PMC_USE_MPI
    ! Load the model data on the root process and pass it to process 1 for solving
    if (pmc_mpi_rank().eq.0) then
#endif

      ! Get the CMAQ_OH_HNO3 reaction mechanism json file
      input_file_path = 'test_CMAQ_OH_HNO3_config.json'

      ! Construct a camp_core variable
      camp_core => camp_core_t(input_file_path)

      deallocate(input_file_path)

      ! Initialize the model
      call camp_core%initialize()

      ! Get the chemical species data
      call assert(673122252, camp_core%get_chem_spec_data(chem_spec_data))

      ! Get species indices
      key = "A"
      idx_A = chem_spec_data%gas_state_id(key);
      key = "B"
      idx_B = chem_spec_data%gas_state_id(key);
      key = "C"
      idx_C = chem_spec_data%gas_state_id(key);

      ! Make sure the expected species are in the model
      call assert(698716805, idx_A.gt.0)
      call assert(811035150, idx_B.gt.0)
      call assert(640878246, idx_C.gt.0)

#ifdef PMC_USE_MPI
      ! pack the camp core
      pack_size = camp_core%pack_size()
      allocate(buffer(pack_size))
      pos = 0
      call camp_core%bin_pack(buffer, pos)
      call assert(585267143, pos.eq.pack_size)
    end if

    ! broadcast the species ids
    call pmc_mpi_bcast_integer(idx_A)
    call pmc_mpi_bcast_integer(idx_B)
    call pmc_mpi_bcast_integer(idx_C)

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
      call assert(580002836, pos.eq.pack_size)
      allocate(buffer_copy(pack_size))
      pos = 0
      call camp_core%bin_pack(buffer_copy, pos)
      call assert(357271680, pos.eq.pack_size)
      do i_elem = 1, pack_size
        call assert_msg(752065274, buffer(i_elem).eq.buffer_copy(i_elem), &
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
      true_conc(0,idx_B) = 0.0
      true_conc(0,idx_C) = 0.0
      model_conc(0,:) = true_conc(0,:)

      ! Set the initial concentrations in the model
      camp_state%state_var(:) = model_conc(0,:)

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
        call assert_msg(189902473, solver_stats%Jac_eval_fails.eq.0, &
                        trim( to_string( solver_stats%Jac_eval_fails ) )// &
                        " Jacobian evaluation failures at time step "// &
                        trim( to_string( i_time ) ) )
#endif

        ! Get the analytic conc
        time = i_time * time_step
        true_conc(i_time,idx_A) = true_conc(0,idx_A) * exp(-(k1)*time)
        true_conc(i_time,idx_B) = true_conc(0,idx_A) * (k1/(k2-k1)) * &
                (exp(-k1*time) - exp(-k2*time))
        true_conc(i_time,idx_C) = true_conc(0,idx_A) * &
               (1.0 + (k1*exp(-k2*time) - k2*exp(-k1*time))/(k2-k1))

      end do

      ! Save the results
      open(unit=7, file="out/CMAQ_OH_HNO3_results.txt", status="replace", &
              action="write")
      do i_time = 0, NUM_TIME_STEP
        write(7,*) i_time*time_step, &
              ' ', true_conc(i_time, idx_A),' ', model_conc(i_time, idx_A), &
              ' ', true_conc(i_time, idx_B),' ', model_conc(i_time, idx_B), &
              ' ', true_conc(i_time, idx_C),' ', model_conc(i_time, idx_C)
      end do
      close(7)

      ! Analyze the results
      do i_time = 1, NUM_TIME_STEP
        do i_spec = 1, size(model_conc, 2)
          call assert_msg(740083421, &
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
      if (run_CMAQ_OH_HNO3_test) then
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
        run_CMAQ_OH_HNO3_test = .true.
      else
        run_CMAQ_OH_HNO3_test = .false.
      end if
    end if

    deallocate(buffer)
#endif

    deallocate(camp_core)

    run_CMAQ_OH_HNO3_test = .true.

  end function run_CMAQ_OH_HNO3_test

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program pmc_test_CMAQ_OH_HNO3
