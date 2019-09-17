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
  integer(kind=i_kind), parameter :: SIZE_STATE = 4

  !Test variables
  real(kind=dp) :: k1, k2, conc_D
  real(kind=dp) :: time_step !todo is need?

  !> The Arrhenius rxn unit test
  ! Varibles inside are used outside the module
  type :: unit_test_rxn_arrhenius_t
  !type, extends(unit_test_data_t) :: unit_test_rxn_arrhenius_t
    !Maybe add the extends for the number of cells? It should be on all tests anyway soo it saves code

    integer(kind=i_kind) :: size_state = SIZE_STATE
    real(kind=dp) :: state(SIZE_STATE)
    character(len=200), dimension(SIZE_STATE) :: keys

    integer(kind=dp) :: state_id(SIZE_STATE)

    !Environment variables
    real(kind=dp) :: temp, pressure

    !private

  contains
    !> Get the number of unique original states available
    procedure :: num_unique_states
    !> Initialize a camp_state_t object based on a given index
    procedure :: initialize_state
    !> Analyze results in a camp_state_t object
    procedure :: analyze_state
    !>Initialize state
    procedure:: get_state
    !> Check results correctness
    procedure :: check_results


  end type unit_test_rxn_arrhenius_t

  !> Constructor for the Arrhenius test
  interface unit_test_rxn_arrhenius_t
    procedure :: constructor
  end interface unit_test_rxn_arrhenius_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Constructor for the Arrhenius unit test
  function constructor() result(this)

    !> A new unit test object
    type(unit_test_rxn_arrhenius_t), pointer :: this

    real(kind=dp) :: conv

    allocate(this)

    !Data initialization
    time_step = 1.0

    this%temp = 272.5d0
    this%pressure = 101253.3d0
    conv = const%avagadro / const%univ_gas_const * 10.0d0**(-12.0d0) * &
            this%pressure / this%temp
    conc_D = 1.2d0 / conv
    k1 = conc_D * conv + 1476.0d0 * exp( -5.5d-21 / &
            (const%boltzmann * this%temp) ) * (this%temp/300.0d0)**(150.0d0) * &
            (1.0d0 + 0.15d0 * this%pressure) / 60.0d0
    k2 = 21.0d0 * exp( -4000.0d0/this%temp ) * (this%temp/315.0d0)**(11.0d0) * &
            (1.0d0 + 0.05d0 * this%pressure)

    ! Species names
    this%keys(1) = "A"
    this%keys(2) = "B"
    this%keys(3) = "C"
    this%keys(4) = "D"

  end function constructor


  function get_state(this, key, state_id) result(state)

    !> Unit test data
    class(unit_test_rxn_arrhenius_t), intent(inout) :: this

    character(len=*), intent(in) :: key
    integer(kind=i_kind), intent(in) :: state_id
    real(kind=dp) :: state(SIZE_STATE)

    !Set state
    if (key.eq.'A') then
      state(state_id) = 1.0
    end if
    if (key.eq.'B') then
      state(state_id) = 0.0
    end if
    if (key.eq.'C') then
      state(state_id) = 0.0
    end if
    if (key.eq.'D') then
      state(state_id) = conc_D
    end if
    
  end function get_state


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !! This function should return true if the checking passes, false otherwise
  function check_results(this) result(passed)

  !> Flag indicating whether the tests passed
  logical :: passed
  !> Unit test data
  class(unit_test_rxn_arrhenius_t), intent(inout) :: this

    !

  end function check_results








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

    ! ... finish ...

  end subroutine initialize_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Analyze a camp_state_t object based on a given initial state id and the
  !> grid-cell environmental parameters
  !! This function should return true if the analysis passes, false otherwise
  function analyze_state(this, grid_cell_id, camp_core, camp_state, &
      unique_state_id, model_time) result (passed)

    !> Flag indicating whether the tests passed
    logical :: passed
    !> Unit test data
    class(unit_test_rxn_arrhenius_t), intent(inout) :: this
    !> Grid cell id
    integer(kind=i_kind), intent(in) :: grid_cell_id
    !> CAMP core
    class(camp_core_t), intent(in) :: camp_core
    !> Grid cell state
    class(camp_state_t), intent(in) :: camp_state
    !> Unique state id
    integer(kind=i_kind), intent(in) :: unique_state_id
    !> Total integration time
    real(kind=dp), intent(in) :: model_time

    ! ... finish ...

    passed = .true.

  end function analyze_state



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_unit_test_rxn_arrhenius
