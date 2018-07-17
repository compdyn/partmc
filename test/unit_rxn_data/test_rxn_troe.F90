! Copyright (C) 2017 Matt Dawson
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_test_troe program

!> Test of troe reaction module
program pmc_test_troe

  use pmc_util,                         only: phlex_real, phlex_int, &
                                              assert, almost_equal, &
                                              string_t, warn_msg
  use pmc_phlex_core
  use pmc_phlex_state
  use pmc_chem_spec_data
#ifdef PMC_USE_JSON
  use json_module
#endif
  use pmc_mpi

  implicit none

  ! Number of timesteps to output in mechanisms
  integer(kind=phlex_int) :: NUM_TIME_STEP = 100

  ! initialize mpi
  call pmc_mpi_init()

  if (run_troe_tests()) then
    if (pmc_mpi_rank().eq.0) write(*,*) "Troe reaction tests - PASS"
  else
    if (pmc_mpi_rank().eq.0) write(*,*) "Troe reaction tests - FAIL"
  end if

  ! finalize mpi
  call pmc_mpi_finalize()

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Run all pmc_chem_mech_solver tests
  logical function run_troe_tests() result(passed)

    use pmc_phlex_solver_data

    type(phlex_solver_data_t), pointer :: phlex_solver_data

    phlex_solver_data => phlex_solver_data_t()

    if (phlex_solver_data%is_solver_available()) then
      passed = run_troe_test()
    else
      call warn_msg(109427097, "No solver available")
      passed = .true.
    end if

    deallocate(phlex_solver_data)

  end function run_troe_tests

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Solve a mechanism of consecutive reactions
  !!
  !! The mechanism is of the form:
  !!
  !!   A -k1-> B -k2-> C
  !!
  !! where k1 and k2 are Troe reaction rate constants.
  logical function run_troe_test()

    use pmc_constants

    type(phlex_core_t), pointer :: phlex_core
    type(phlex_state_t), pointer :: phlex_state
    character(len=:), allocatable :: input_file_path, key
    type(string_t), allocatable, dimension(:) :: output_file_path

    type(chem_spec_data_t), pointer :: chem_spec_data
    real(kind=phlex_real), dimension(0:NUM_TIME_STEP, 3) :: model_conc, true_conc
    integer(kind=phlex_int) :: idx_A, idx_B, idx_C, i_time, i_spec
    real(kind=phlex_real) :: time_step, time, k1, k2, air_conc, temp, pressure, k_0, &
            k_inf, conv
#ifdef PMC_USE_MPI
    character, allocatable :: buffer(:), buffer_copy(:)
    integer(kind=phlex_int) :: pack_size, pos, i_elem, results
#endif

    run_troe_test = .true.

    ! Set the rate constants (for calculating the true value)
    temp = 272.5d0
    pressure = 101253.3d0
    air_conc = 1.0d6
    conv = const%avagadro / const%univ_gas_const * 10.0d0**(-12.0d0) * &
            pressure / temp
    k_0 = air_conc * 4.0e-18 * conv
    k_inf = k_0 / 1.0d0
    k1 = k_0/(1+k_inf) * 0.6d0**(1.0d0/(1.0d0 + (log10(k_inf)/1.0)**(2)))
    
    k_0 = air_conc * 1.2d-12 * exp( 3.0d0/temp) * (temp/300.0d0)**(167.0d0) &
            * conv
    k_inf = 136.0d0 * exp( 24.0d0/temp) * (temp/300.0d0)**(5.0d0)
    k_inf = k_0 / k_inf
    k2 = k_0/(1+k_inf) * 0.9d0**(1.0d0/(1.0d0 + &
            (log10(k_inf)/0.8d0)**(2)))/60.0d0

    ! Set output time step (s)
    time_step = 1.0

#ifdef PMC_USE_MPI
    ! Load the model data on the root process and pass it to process 1 for solving
    if (pmc_mpi_rank().eq.0) then
#endif

      ! Get the troe reaction mechanism json file
      input_file_path = 'test_troe_config.json'

      ! Construct a phlex_core variable
      phlex_core => phlex_core_t(input_file_path)

      ! Initialize the model
      call phlex_core%initialize()

      ! Get the chemical species data
      call assert(145118257, phlex_core%get_chem_spec_data(chem_spec_data))

      ! Get species indices
      key = "A"
      idx_A = chem_spec_data%gas_state_id(key);
      key = "B"
      idx_B = chem_spec_data%gas_state_id(key);
      key = "C"
      idx_C = chem_spec_data%gas_state_id(key);

      ! Make sure the expected species are in the model
      call assert(341233628, idx_A.gt.0)
      call assert(171076724, idx_B.gt.0)
      call assert(900919819, idx_C.gt.0)

#ifdef PMC_USE_MPI
      ! pack the phlex core
      pack_size = phlex_core%pack_size()
      allocate(buffer(pack_size))
      pos = 0
      call phlex_core%bin_pack(buffer, pos)
      call assert(363562528, pos.eq.pack_size)
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
      phlex_core => phlex_core_t()
      pos = 0
      call phlex_core%bin_unpack(buffer, pos)
      call assert(423306621, pos.eq.pack_size)
      allocate(buffer_copy(pack_size))
      pos = 0
      call phlex_core%bin_pack(buffer_copy, pos)
      call assert(818100215, pos.eq.pack_size)
      do i_elem = 1, pack_size
        call assert_msg(247885410, buffer(i_elem).eq.buffer_copy(i_elem), &
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
      true_conc(0,idx_B) = 0.0
      true_conc(0,idx_C) = 0.0
      model_conc(0,:) = true_conc(0,:)

      ! Set the initial concentrations in the model
      phlex_state%state_var(:) = model_conc(0,:)

      ! Integrate the mechanism
      do i_time = 1, NUM_TIME_STEP

        ! Get the modeled conc
        call phlex_core%solve(phlex_state, time_step)
        model_conc(i_time,:) = phlex_state%state_var(:)

        ! Get the analytic conc
        time = i_time * time_step
        true_conc(i_time,idx_A) = true_conc(0,idx_A) * exp(-(k1)*time)
        true_conc(i_time,idx_B) = true_conc(0,idx_A) * (k1/(k2-k1)) * &
                (exp(-k1*time) - exp(-k2*time))
        true_conc(i_time,idx_C) = true_conc(0,idx_A) * &
               (1.0 + (k1*exp(-k2*time) - k2*exp(-k1*time))/(k2-k1))

      end do

      ! Save the results
      open(unit=7, file="out/troe_results.txt", status="replace", &
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
          call assert_msg(911807542, &
            almost_equal(model_conc(i_time, i_spec), &
            true_conc(i_time, i_spec), real(1.0e-2, kind=phlex_real)).or. &
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
      if (run_troe_test) then
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
        run_troe_test = .true.
      else
        run_troe_test = .false.
      end if
    end if

    deallocate(buffer)
#endif

    deallocate(phlex_core)

  end function run_troe_test

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program pmc_test_troe
