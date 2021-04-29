! Copyright (C) 2017 Matt Dawson
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_test_wennberg_no_ro2 program

!> Test of wennberg_no_ro2 reaction module
program pmc_test_wennberg_no_ro2

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

  if (run_wennberg_no_ro2_tests()) then
    if (pmc_mpi_rank().eq.0) write(*,*) "Wennberg NO + RO2 reaction tests - PASS"
  else
    if (pmc_mpi_rank().eq.0) write(*,*) "Wennberg NO + RO2 reaction tests - FAIL"
    stop 3
  end if

  ! finalize mpi
  call pmc_mpi_finalize()

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Run all pmc_chem_mech_solver tests
  logical function run_wennberg_no_ro2_tests() result(passed)

    use pmc_camp_solver_data

    type(camp_solver_data_t), pointer :: camp_solver_data

    camp_solver_data => camp_solver_data_t()

    if (camp_solver_data%is_solver_available()) then
      passed = run_wennberg_no_ro2_test()
    else
      call warn_msg(383750890, "No solver available")
      passed = .true.
    end if

    deallocate(camp_solver_data)

  end function run_wennberg_no_ro2_tests

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Solve a mechanism of consecutive reactions
  !!
  !! The mechanism is of the form:
  !!
  !!   A -k1_alkoxy->  C
  !!     -k1_nitrate-> B -k2_alkoxy->  C
  !!                     -k2_nitrate-> D
  !!
  !! where k1 and k2 are Wennberg NO + RO2 reaction rate constants.
  logical function run_wennberg_no_ro2_test()

    use pmc_constants

    type(camp_core_t), pointer :: camp_core
    type(camp_state_t), pointer :: camp_state
    character(len=:), allocatable :: input_file_path, key
    type(string_t), allocatable, dimension(:) :: output_file_path

    type(chem_spec_data_t), pointer :: chem_spec_data
    real(kind=dp), dimension(0:NUM_TIME_STEP, 4) :: model_conc, true_conc
    integer(kind=i_kind) :: idx_A, idx_B, idx_C, idx_D, i_time, i_spec
    real(kind=dp) :: time_step, time, k1_alkoxy, k1_nitrate, k2_alkoxy, &
                     k1, k2, Z, A, air_density, k2_nitrate, temp, pressure, &
                     conv
#ifdef PMC_USE_MPI
    character, allocatable :: buffer(:), buffer_copy(:)
    integer(kind=i_kind) :: pack_size, pos, i_elem, results
#endif

    type(solver_stats_t), target :: solver_stats

    run_wennberg_no_ro2_test = .true.

    ! Set the rate constants (for calculating the true value)
    temp = 272.5d0 ! [K]
    pressure = 101253.3d0 ! [Pa]
    air_density = pressure / ( const%boltzmann * temp ) &
                  * 1.0d-6 ! [molec/cm3]
    conv = const%avagadro / const%univ_gas_const * 10.0d0**(-12.0d0) * &
            pressure / temp
    Z = 0d0
    A = get_A( temp, air_density, 0 )
    k1_alkoxy =  4.0d-3 * Z / ( Z + A )
    k1_nitrate = 4.0d-3 * A / ( A + Z )
    k1 = 4.0d-3
    Z = get_A( 293.0d0, 2.45d19, 9 ) * ( 1.0d0 - 0.15d0 ) / 0.15d0
    A = get_A( temp, air_density, 9 )
    k2_alkoxy  = 1.2d-4 * exp( -167d0 / temp ) * Z / ( Z + A ) / 60.0d0
    k2_nitrate = 1.2d-4 * exp( -167d0 / temp ) * A / ( A + Z ) / 60.0d0
    k2 = 1.2d-4 * exp( -167d0 / temp )

    ! Set output time step (s)
    time_step = 1.0

#ifdef PMC_USE_MPI
    ! Load the model data on the root process and pass it to process 1 for solving
    if (pmc_mpi_rank().eq.0) then
#endif

      ! Get the wennberg_no_ro2 reaction mechanism json file
      input_file_path = 'test_wennberg_no_ro2_config.json'

      ! Construct a camp_core variable
      camp_core => camp_core_t(input_file_path)

      ! Initialize the model
      call camp_core%initialize()

      ! Get the chemical species data
      call assert(584648230, camp_core%get_chem_spec_data(chem_spec_data))

      ! Get species indices
      key = "A"
      idx_A = chem_spec_data%gas_state_id(key);
      key = "B"
      idx_B = chem_spec_data%gas_state_id(key);
      key = "C"
      idx_C = chem_spec_data%gas_state_id(key);
      key = "D"
      idx_D = chem_spec_data%gas_state_id(key);

      ! Make sure the expected species are in the model
      call assert(244334422, idx_A.gt.0)
      call assert(691702268, idx_B.gt.0)
      call assert(521545364, idx_C.gt.0)
      call assert(251791386, idx_D.gt.0)

#ifdef PMC_USE_MPI
      ! pack the camp core
      pack_size = camp_core%pack_size()
      allocate(buffer(pack_size))
      pos = 0
      call camp_core%bin_pack(buffer, pos)
      call assert(581289457, pos.eq.pack_size)
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
      call assert(188401397, pos.eq.pack_size)
      allocate(buffer_copy(pack_size))
      pos = 0
      call camp_core%bin_pack(buffer_copy, pos)
      call assert(300719742, pos.eq.pack_size)
      do i_elem = 1, pack_size
        call assert_msg(348029687, buffer(i_elem).eq.buffer_copy(i_elem), &
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
      true_conc(0,idx_D) = 0.0
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
        call assert_msg(572666377, solver_stats%Jac_eval_fails.eq.0, &
                        trim( to_string( solver_stats%Jac_eval_fails ) )// &
                        " Jacobian evaluation failures at time step "// &
                        trim( to_string( i_time ) ) )
#endif

        ! Get the analytic conc
        time = i_time * time_step
        true_conc(i_time,idx_A) = true_conc(0,idx_A) * exp(-(k1)*time)
        true_conc(i_time,idx_B) = true_conc(0,idx_A) * (k1_nitrate/(k2-k1_nitrate)) * &
                (exp(-k1_nitrate*time) - exp(-k2*time))
        true_conc(i_time,idx_C) = true_conc(0,idx_A) * &
               (1.0 + (k1_nitrate*exp(-k2_alkoxy*time) - &
               k2_alkoxy*exp(-k1_nitrate*time))/(k2_alkoxy-k1_nitrate))
        true_conc(i_time,idx_D) = true_conc(0,idx_A) * &
               (1.0 + (k1_nitrate*exp(-k2_nitrate*time) - &
               k2_nitrate*exp(-k1_nitrate*time))/(k2_nitrate-k1_nitrate))

      end do

      ! Save the results
      open(unit=7, file="out/wennberg_no_ro2_results.txt", status="replace", &
            action="write")
      do i_time = 0, NUM_TIME_STEP
        write(7,*) i_time*time_step, &
              ' ', true_conc(i_time, idx_A),' ', model_conc(i_time, idx_A), &
              ' ', true_conc(i_time, idx_B),' ', model_conc(i_time, idx_B), &
              ' ', true_conc(i_time, idx_C),' ', model_conc(i_time, idx_C), &
              ' ', true_conc(i_time, idx_D),' ', model_conc(i_time, idx_D)
      end do
      close(7)

      ! Analyze the results
      do i_time = 1, NUM_TIME_STEP
        do i_spec = 1, size(model_conc, 2)
          call assert_msg(232352569, &
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
      if (run_wennberg_no_ro2_test) then
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
        run_wennberg_no_ro2_test = .true.
      else
        run_wennberg_no_ro2_test = .false.
      end if
    end if

    deallocate(buffer)
#endif

    deallocate(camp_core)

  end function run_wennberg_no_ro2_test

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Get A(T,[M],n)
  real(kind=dp) function get_A( T, M, n )

    real(kind=dp),        intent(in) :: T, M
    integer(kind=i_kind), intent(in) :: n
    real(kind=dp) :: k0M, kinf

    k0M = 2d-22 * exp( real( n, kind=dp ) ) * M
    kinf = 0.43 * ( T / 298d0 )**(-8)
    get_A = k0M / (1.0d0 + k0M/kinf) * &
            0.41**(1.0d0/(1.0d0 + (log10(k0M/kinf))**2))

  end function get_A

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program pmc_test_wennberg_no_ro2
