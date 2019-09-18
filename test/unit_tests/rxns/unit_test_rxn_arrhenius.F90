! Copyright (C) 2019 Christian Guzman and Matt Dawson
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_unit_test_rxn_arrhenius module
module pmc_unit_test_rxn_arrhenius

  use pmc_camp_core
  use pmc_camp_state
  use pmc_unit_test_data
  use pmc_util

  implicit none
  private

  public :: unit_test_rxn_arrhenius_t

  !> Number of available initial states
  integer(kind=i_kind), parameter :: NUM_INIT_STATES = 5
  !> Number of time steps
  integer(kind=i_kind), parameter :: NUM_TIME_STEPS_C = 100
  !> Time step size (s)
  real(kind=dp), parameter :: TIME_STEP_SIZE_C = 1.0


  !> The Arrhenius rxn unit test
  type, extends(unit_test_data_t) :: unit_test_rxn_arrhenius_t
    private
    !> Species ids
    integer(kind=i_kind) :: idx_A, idx_B, idx_C, idx_D
    !> Species initial concentrations
    real(kind=dp), dimension(NUM_INIT_STATES) :: conc_A, conc_B, conc_C, &
                                                 conc_D
    !> Flag inidicating local data is set up
    logical :: is_initialized = .false.
  contains
    !> Get the name of the input file for the test
    procedure :: input_file_name
    !> Get the number of unique original states available
    procedure :: num_unique_states
    !> Initialize a camp_state_t object based on a given index
    procedure :: initialize_state
    !> Number of time steps for this test
    procedure :: num_time_steps
    !> Time step size (s)
    procedure :: time_step_size
    !> Analyze results in a camp_state_t object
    procedure :: analyze_state
    !> Initialize the object data
    procedure, private :: initialize
  end type unit_test_rxn_arrhenius_t

  !> Constructor for the Arrhenius test
  interface unit_test_rxn_arrhenius_t
    procedure :: constructor
  end interface unit_test_rxn_arrhenius_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Constructor for the Arrhenius unit test
  function constructor() result(new_obj)

    !> A new unit test object
    type(unit_test_rxn_arrhenius_t), pointer :: new_obj

    allocate(new_obj)

    !... finish ...

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Input file name
  function input_file_name(this)

    !> Input file name
    character(len=:), allocatable :: input_file_name
    !> Unit test data
    class(unit_test_rxn_arrhenius_t), intent(in) :: this

    input_file_name = "arrhenius_config.json"

  end function input_file_name

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Number of unique initial states available from the unit test
  function num_unique_states(this)

    !> Number of unique states
    integer(kind=i_kind) :: num_unique_states
    !> Unit test data
    class(unit_test_rxn_arrhenius_t), intent(in) :: this

    num_unique_states = NUM_INIT_STATES

  end function num_unique_states

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Initialize a camp_state_t object based on a given index
  subroutine initialize_state(this, grid_cell_id, camp_core, camp_state, &
      unique_state_id)

    !> Unit test data
    class(unit_test_rxn_arrhenius_t), intent(inout) :: this
    !> Grid cell id
    integer(kind=i_kind), intent(in) :: grid_cell_id
    !> CAMP core
    class(camp_core_t), intent(in) :: camp_core
    !> Grid cell state
    class(camp_state_t), intent(inout) :: camp_state
    !> Unique state id
    integer(kind=i_kind), intent(in) :: unique_state_id

    camp_state%env_state%temp = 278.0d0
    camp_state%env_state%pressure = 101325.0d0
    call camp_state%update_env_state( )

    camp_state%state_var( : ) = 0.0d0

    ! ... finish ...

  end subroutine initialize_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Number of time steps for this test
  function num_time_steps(this)

    !> Number of time steps
    integer(kind=i_kind) :: num_time_steps
    !> Unit test data
    class(unit_test_rxn_arrhenius_t), intent(in) :: this

    num_time_steps = NUM_TIME_STEPS_C

  end function num_time_steps

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Time step size for this test (s)
  function time_step_size(this)

    !> Time step size (s)
    real(kind=dp) :: time_step_size
    !> Unit test data
    class(unit_test_rxn_arrhenius_t), intent(in) :: this

    time_step_size = TIME_STEP_SIZE_C

  end function time_step_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Analyze a camp_state_t object based on a given initial state id and the
  !> grid-cell environmental parameters
  !! This function should return true if the analysis passes, false otherwise
  function analyze_state(this, camp_core, camp_state, &
      unique_state_id, model_time_step) result (passed)

    !> Flag indicating whether the tests passed
    logical :: passed
    !> Unit test data
    class(unit_test_rxn_arrhenius_t), intent(inout) :: this
    !> CAMP core
    class(camp_core_t), intent(in) :: camp_core
    !> Grid cell state
    class(camp_state_t), intent(in) :: camp_state
    !> Unique state id
    integer(kind=i_kind), intent(in) :: unique_state_id
    !> Lastest time step solved
    integer(kind=i_kind), intent(in) :: model_time_step

    ! ... finish ...

    passed = .true.

  end function analyze_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Initialize the object data
  subroutine initialize(this, camp_core)

    use pmc_chem_spec_data

    !> Unit test data
    class(unit_test_rxn_arrhenius_t), intent(inout) :: this
    !> CAMP core
    class(camp_core_t), intent(in) :: camp_core

    type(chem_spec_data_t), pointer :: chem_spec_data

    if (this%is_initialized) return

    ! Get the chemical species data
    call assert(672415022, camp_core%get_chem_spec_data(chem_spec_data))

    ! Get the species ids

    !TODO finish

  end subroutine initialize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_unit_test_rxn_arrhenius
