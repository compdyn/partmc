! Copyright (C) 2019 Christian Guzman and Matt Dawson
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The unit test driver program.
program unit_test_driver

  use pmc_camp_core
  use pmc_camp_state
  use pmc_mpi
  use pmc_rand
  use pmc_unit_test_data
  use UNIT_TEST_MODULE_

  implicit none

  ! Number of grid cells to solve
  integer(kind=i_kind), parameter :: N_CELLS = 30

  ! State pointer for building array of states
  type :: state_ptr
    class(camp_state_t), pointer :: val
  end type state_ptr

  ! Unit test object
  class(unit_test_data_t), pointer :: unit_test
  ! CAMP core
  class(camp_core_t), pointer :: camp_core

  type(state_ptr) :: grid_cell_state(N_CELLS)
  integer(kind=i_kind) :: grid_cell_state_id(N_CELLS)
  integer(kind=i_kind) :: i_cell
  real(kind=dp) :: time
  logical :: passed

  ! initialize MPI
  call pmc_mpi_init()

  ! Seed the random number generator
  call pmc_srand(0,0)

  ! Instantiate the unit test
  unit_test => UNIT_TEST_TYPE_

  ! Initialize the model
  camp_core => camp_core_t(unit_test%input_file_name())
  call camp_core%initialize()

  !... finish model init ...

#ifdef PMC_USE_MPI
  !... core from primary process to run on another process ...
#endif

  ! Initialize the solver
  call camp_core%solver_initialize()

  ! Set the initial states (unit_test_data_t classes should provide a certain number of
  ! unique initial states that can be analyzed during solving)
  do i_cell=1, N_CELLS
    grid_cell_state(i_cell)%val => camp_core%new_state()
    grid_cell_state_id(i_cell) = pmc_rand_int(unit_test%num_unique_states())

    ! Set the initial state for this grid cell to grid_cell_state_id
    call unit_test%initialize_state(i_cell, camp_core, &
                                    grid_cell_state(i_cell)%val, &
                                    grid_cell_state_id(i_cell))
  end do

  ! Loop over time
  ! Loop over cells
  !   call camp_core%solve(grid_cell_state(i_cell), time_step)

  ! (or solve multi-cell)

  passed = unit_test%analyze_state(i_cell, camp_core, &
                                   grid_cell_state(i_cell)%val, &
                                   grid_cell_state_id(i_cell), time)

  ! output state to file

  ! deallocate pointers and arrays
  deallocate(unit_test)
  deallocate(camp_core)
  do i_cell = 1, N_CELLS
    deallocate(grid_cell_state(i_cell)%val)
  end do

  ! Finalize the test
  if (passed) then
    if (pmc_mpi_rank().eq.0) write(*,*) "Unit test - PASS"
    call pmc_mpi_finalize()
  else
    if (pmc_mpi_rank().eq.0) write(*,*) "Unit test - FAIL"
    call pmc_mpi_finalize()
    stop 3
  end if

end program unit_test_driver
