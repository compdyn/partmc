! Copyright (C) 2017 Matt Dawson
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_test_arrhenius program

!> Test of arrhenius reaction module
program pmc_test_arrhenius

  use pmc_util,                         only: i_kind, dp, assert, &
                                              almost_equal, string_t, &
                                              warn_msg
  use pmc_phlex_core
  use pmc_phlex_state
#ifdef PMC_USE_JSON
  use json_module
#endif
  use pmc_mpi

  implicit none

  ! Number of timesteps to output in mechanisms
  integer(kind=i_kind) :: NUM_TIME_STEP = 100

  ! initialize mpi
  call pmc_mpi_init()

  if (run_arrhenius_tests()) then
    write(*,*) "Arrhenius reaction tests - PASS"
  else
    write(*,*) "Arrhenius reaction tests - FAIL"
  end if

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Run all pmc_chem_mech_solver tests
  logical function run_arrhenius_tests() result(passed)

    use pmc_phlex_solver_data

    type(phlex_solver_data_t), pointer :: phlex_solver_data

    phlex_solver_data => phlex_solver_data_t()

    if (phlex_solver_data%is_solver_available()) then
      passed = run_arrhenius_test()
    else
      call warn_msg(398972036, "No solver available")
      passed = .true.
    end if

  end function run_arrhenius_tests

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Solve a mechanism of consecutive reactions
  !!
  !! The mechanism is of the form:
  !!
  !!   A -k1-> B -k2-> C
  !!
  !! where k1 and k2 are Arrhenius reaction rate constants.
  logical function run_arrhenius_test()

    use pmc_constants

    type(phlex_core_t), pointer :: phlex_core
    type(phlex_state_t), target :: phlex_state
    character(len=:), allocatable :: input_file_path
    type(string_t), allocatable, dimension(:) :: output_file_path

    real(kind=dp), dimension(0:NUM_TIME_STEP, 3) :: model_conc, true_conc
    integer(kind=i_kind) :: idx_A, idx_B, idx_C
    character(len=:), allocatable :: key
    integer(kind=i_kind) :: i_time, i_spec
    real(kind=dp) :: time_step, time

    ! Parameters for calculating true concentrations
    real(kind=dp) :: k1, k2, temp, pressure

    run_arrhenius_test = .true.

    ! Set the rate constants (for calculating the true value)
    temp = 272.5d0
    pressure = 101253.3d0
    k1 = 1.0d0
    k1 = k1 + 1476.0d0 * exp( -5.5d-21 / (const%boltzmann * temp) ) * &
            (temp/300.0d0)**(150.0d0) * (1.0d0 + 0.15d0 * pressure) / 60.0d0
    k2 = 21.0d0 * exp( -4000.0d0/temp ) * (temp/315.0d0)**(11.0d0) * &
            (1.0d0 + 0.05d0 * pressure)

    ! Set output time step (s)
    time_step = 1.0

    ! Get the arrhenius reaction mechanism json file
    input_file_path = 'test_arrhenius_config.json'

    ! Construct a phlex_core variable
    phlex_core => phlex_core_t(input_file_path)

    ! Initialize the model
    call phlex_core%initialize()

    ! Get a model state variable
    phlex_state = phlex_core%new_state()

    ! Set the environmental conditions
    phlex_state%env_state%temp = temp
    phlex_state%env_state%pressure = pressure
    call phlex_state%update_env_state()

    ! Get species indices
    key = "A"
    idx_A = phlex_core%chem_spec_data%gas_state_id(key);
    key = "B"
    idx_B = phlex_core%chem_spec_data%gas_state_id(key);
    key = "C"
    idx_C = phlex_core%chem_spec_data%gas_state_id(key);

    ! Make sure the expected species are in the model
    call assert(629811894, idx_A.gt.0)
    call assert(226395220, idx_B.gt.0)
    call assert(338713565, idx_C.gt.0)

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
    open(unit=7, file="out/arrhenius_results.txt", status="replace", action="write")
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
        call assert_msg(848069355, &
          almost_equal(model_conc(i_time, i_spec), true_conc(i_time, i_spec), &
          real(1.0e-2, kind=dp)), "time: "//to_string(i_time)//"; species: "// &
          to_string(i_spec)//"; mod: "//to_string(model_conc(i_time, i_spec))// &
          "; true: "//to_string(true_conc(i_time, i_spec)))
      end do
    end do

    run_arrhenius_test = .true.

  end function run_arrhenius_test

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program pmc_test_arrhenius
