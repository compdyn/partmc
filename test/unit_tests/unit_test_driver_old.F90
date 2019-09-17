! Copyright (C) 2019 Christian Guzman and Matt Dawson
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The unit test driver program.
program unit_test_driver

  use pmc_util,                          only : assert_msg, almost_equal, &
                                                to_string
  use pmc_chem_spec_data
  use pmc_camp_core
  use pmc_camp_state
  use pmc_solver_stats
  use pmc_mpi
  use pmc_rand
  use pmc_unit_test_data
  !use UNIT_TEST_MODULE_

  !Unit test modules
  use pmc_unit_test_rxn_arrhenius

  !todo mapping of tests (DEFINE ARRHENIUS_ID 1, CMAQ_ID 2, etc), to preserve the correct order and be clear

  implicit none

  !Test mapping
  integer(kind=i_kind), parameter :: ARRHENIUS = 1

  ! Number of grid cells to solve
  integer(kind=i_kind), parameter :: N_CELLS = 1
  integer(kind=i_kind), parameter :: NUM_TIME_STEP = 20
  integer(kind=i_kind), parameter :: N_TESTS = 2
  real(kind=dp), parameter :: TIME_STEP = 1.0

  character(len=200), dimension(N_TESTS):: input_files_path

  !Paths to input files
  !input_files_path(1) = "arrhenius_config.json"


  ! State pointer for building array of states
  !type :: state_ptr
  !  class(camp_state_t), pointer :: val
  !end type state_ptr

  ! Unit test object
  ! CAMP core
  class(camp_core_t), pointer :: camp_core

  !Unit test objects
  type(chem_spec_data_t), pointer :: chem_spec_data
  type(unit_test_rxn_arrhenius_t), pointer :: test_arrhenius
  type(camp_state_t), pointer :: camp_state
  type(solver_stats_t) :: solver_stats

  !type(state_ptr) :: state(N_CELLS)
  !TODO: Store all the states for time steps?
  real(kind=dp), allocatable :: test_state(:)
  real(kind=dp), allocatable :: test_temp(:)
  real(kind=dp), allocatable :: test_pressure(:)

  character(len=:), allocatable :: input_file_path

  integer(kind=i_kind) :: i, id_test, i_time, i_cell = 1

  real(kind=dp) :: time
  logical :: passed

  !Paths input data for all tests !TODO: Move to a better place
  input_files_path(1) = "input_files/arrhenius_config.json"
  input_files_path(2) = "input_files/arrhenius_config.json"

  !the true workflow:

  !get user_data() !from json
  !get monarch_data
  !init cvode & set data on C files
  !run cvode with her functions
  !print state
  !deallocate all data
  !the end

#ifdef PMC_USE_MPI
  call pmc_mpi_init()
#endif

#ifdef PMC_DEBUG
  ! Evaluate the Jacobian during solving
  solver_stats%eval_Jac = .true.
#endif

  do id_test = 1, N_TESTS

    ! Initialize the model
    camp_core => camp_core_t(input_files_path(id_test))
    call camp_core%initialize()

    !Get input data form json files
    !TODO: improve for load more than chem_spec_data (aerosol, gas...)
    call assert(657781626, camp_core%get_chem_spec_data(chem_spec_data))

  !#ifdef PMC_USE_MPI
    !... core from primary process to run on another process ...
  !#endif

    ! Initialize the solver
    call camp_core%solver_initialize()

    !Create state object
    camp_state => camp_core%new_state()

    !Initialize test object and init variables(state,temp,press...)
    call init_test()

    !If N_cells!=1 -> Set_values for the rest of cells (state, temps...)

    !for n_cells

    do i_cell=1, N_CELLS

      ! Set the environmental conditions
      camp_state%env_state%temp = temp(i_cell)
      camp_state%env_state%pressure = pressure(i_cell)
      call camp_state%update_env_state(i_cell)

    end do

    !Solve
    ! Loop over time
    do i_time=1, NUM_TIME_STEP

      call camp_core%solve(this%camp_state, &
              real(time_step, kind=dp), solver_stats = solver_stats)

      !todo check all time steps or only 1?
      !passed = unit_test%analyze_state(i_cell, camp_core, &
      ! grid_cell_state(i_cell)%val, grid_cell_state_id(i_cell), time)

      ! output state to file

    end do


    call deallocate_data()

    end do


  ! deallocate pointers and arrays
  !deallocate(unit_test)
  !deallocate(camp_core)
  !deallocate(state)
  !do i_cell = 1, N_CELLS
  !  deallocate(grid_cell_state(i_cell)%val)
  !end do

  passed = .true.

  ! Finalize the test
  if (passed) then
    if (pmc_mpi_rank().eq.0) write(*,*) "Unit test - PASS"
    call pmc_mpi_finalize()
  else
    if (pmc_mpi_rank().eq.0) write(*,*) "Unit test - FAIL"
    call pmc_mpi_finalize()
    stop 3
  end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

contains

  !init test/model data (state)
  subroutine init_test()

    integer(kind=i_kind) :: size_state_per_cell
    integer(kind=i_kind) :: size_state
    integer(kind=i_kind) :: state_id
    integer(kind=i_kind) :: i_key
    
    print*, "Initializing test", id_test

    select case (id_test)
      case (1)

        !Create object
        test_arrhenius => unit_test_rxn_arrhenius_t()

        !Get sizes
        size_state_per_cell = test_arrhenius%size_state
        size_state = size_state_per_cell * N_CELLS

        ! Allocate variables
        allocate(test_state(size_state))
        allocate(test_temp(N_CELLS))
        allocate(test_pressure(N_CELLS))

        !Link chem_spec_ids with test_ids
        do i_key=1, size_state_per_cell
          !Get state_id for the key
          state_id = chem_spec_data%gas_state_id(test_arrhenius%keys(i_key))
          !Set state value for this key
          test_state= test_arrhenius%get_state(test_arrhenius%keys(i_key),state_id)
        end do

        test_temp = test_arrhenius%temp
        test_pressure = test_arrhenius%pressure

        !Todo get temp and press
        print*, test_state
        print*, test_temp

      case (2)

        !size_state_per_cell=arrhenius->size_state;
        size_state_per_cell=4
        size_state = size_state_per_cell * N_CELLS

        ! Set the initial state values
        allocate(test_state(size_state))

        !test_state=arrhenius->state;

        test_state(:) = 0.2
        !do i=1, size_state
        !   test_state(i) = test_state(i) + 0.01*i
        !end do
        print*, test_state

        !...

    end select




  end subroutine init_test




  subroutine deallocate_data()

    !Deallocate arrays
    deallocate(test_state)

    !Deallocate objects
    select case (id_test)
      case (1)
        deallocate(test_arrhenius)

    end select

  end subroutine deallocate_data

end program unit_test_driver
