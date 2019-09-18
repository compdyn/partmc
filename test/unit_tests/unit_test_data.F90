! Copyright (C) 2019 Christian Guzman and Matt Dawson
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_unit_test_data module
module pmc_unit_test_data

  implicit none
  private

  public :: unit_test_data_t

  !> Unit test data type
  type, abstract :: unit_test_data_t

  contains
    !> Get the name of the input file for the test
    procedure(input_file_name), deferred :: input_file_name
    !> Get the number of unique original states available
    procedure(num_unique_states), deferred :: num_unique_states
    !> Initialize a camp_state_t object based on a given index
    procedure(initialize_state), deferred :: initialize_state
    !> Get the number of time steps to run for this test
    procedure(num_time_steps), deferred :: num_time_steps
    !> Get the size of the time step for this test (s)
    procedure(time_step_size), deferred :: time_step_size
    !> Analyze results in a camp_state_t object
    procedure(analyze_state), deferred :: analyze_state
  end type unit_test_data_t

interface

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Input file name
  function input_file_name(this)
    import :: unit_test_data_t

    !> Input file name
    character(len=:), allocatable :: input_file_name
    !> Unit test data
    class(unit_test_data_t), intent(in) :: this

  end function input_file_name

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Number of unique initial states available from the unit test
  function num_unique_states(this)
    use pmc_util,                                only : i_kind
    import :: unit_test_data_t

    !> Number of unique states
    integer(kind=i_kind) :: num_unique_states
    !> Unit test data
    class(unit_test_data_t), intent(in) :: this

  end function num_unique_states

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Initialize a camp_state_t object based on a given index
  subroutine initialize_state(this, grid_cell_id, camp_core, camp_state, &
      unique_state_id)
    use pmc_camp_core
    use pmc_camp_state
    use pmc_util,                                only : i_kind
    import :: unit_test_data_t

    !> Unit test data
    class(unit_test_data_t), intent(inout) :: this
    !> Grid cell id
    integer(kind=i_kind), intent(in) :: grid_cell_id
    !> CAMP core
    class(camp_core_t), intent(in) :: camp_core
    !> Grid cell state
    class(camp_state_t), intent(inout) :: camp_state
    !> Unique state id
    integer(kind=i_kind), intent(in) :: unique_state_id

  end subroutine initialize_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get the number of time steps for this test
  function num_time_steps(this)
    use pmc_util,                                only : i_kind
    import :: unit_test_data_t

    !> Number of time steps
    integer(kind=i_kind) :: num_time_steps
    !> Unit test data
    class(unit_test_data_t), intent(in) :: this

  end function num_time_steps

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get the size of the time step for this test
  function time_step_size(this)
    use pmc_util,                                only : dp
    import :: unit_test_data_t

    !> Time step size
    real(kind=dp) :: time_step_size
    !> Unit test data
    class(unit_test_data_t), intent(in) :: this

  end function time_step_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Analyze a camp_state_t object based on a given initial state id and the
  !> grid-cell environmental parameters
  !! This function should return true if the analysis passes, false otherwise
  function analyze_state(this, camp_core, camp_state, &
      unique_state_id, model_time_step) result (passed)
    use pmc_camp_core
    use pmc_camp_state
    use pmc_util,                                only : i_kind
    import :: unit_test_data_t

    !> Flag indicating whether the tests passed
    logical :: passed
    !> Unit test data
    class(unit_test_data_t), intent(inout) :: this
    !> CAMP core
    class(camp_core_t), intent(in) :: camp_core
    !> Grid cell state
    class(camp_state_t), intent(in) :: camp_state
    !> Unique state id
    integer(kind=i_kind), intent(in) :: unique_state_id
    !> Model time step most recently solved
    integer(kind=i_kind), intent(in) :: model_time_step

  end function analyze_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end interface

end module pmc_unit_test_data
