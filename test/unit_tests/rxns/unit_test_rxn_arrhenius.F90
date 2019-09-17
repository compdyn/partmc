! Copyright (C) 2019 Christian Guzman and Matt Dawson
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_unit_test_rxn_arrhenius module
program unit_test_rxn_arrhenius

  use pmc_camp_core
  use pmc_camp_state
  use pmc_util
  use pmc_unit_test_driver

  implicit none

  !> Number of available initial states
  integer(kind=i_kind), parameter :: NUM_INIT_STATES = 5
  integer(kind=i_kind), parameter :: SIZE_STATE = 4
  integer(kind=i_kind), parameter :: NUM_TIME_STEP = 20
  character(len=200), dimension(SIZE_STATE) :: keys
  character(len=:), allocatable :: input_file_path


  real(kind=dp) :: k1, k2, conc_D
  real(kind=dp) :: time_step = 1.0

  real(kind=dp) :: state(SIZE_STATE)
  real(kind=dp) :: state_analytic(SIZE_STATE)

  !Environment variables
  real(kind=dp) :: temp, pressure
  real(kind=dp) :: conv
  integer(kind=i_kind) :: i_time, time

  type(unit_test_driver_t), pointer :: test_driver

  !Data initialization
  input_file_path = '../input_files/arrhenius_config.json'
  time_step = 1.0

  temp = 272.5d0
  pressure = 101253.3d0
  conv = const%avagadro / const%univ_gas_const * 10.0d0**(-12.0d0) * &
          pressure / temp
  conc_D = 1.2d0 / conv
  k1 = conc_D * conv + 1476.0d0 * exp( -5.5d-21 / &
          (const%boltzmann * temp) ) * (temp/300.0d0)**(150.0d0) * &
          (1.0d0 + 0.15d0 * pressure) / 60.0d0
  k2 = 21.0d0 * exp( -4000.0d0/temp ) * (temp/315.0d0)**(11.0d0) * &
          (1.0d0 + 0.05d0 * pressure)

  !TODO: Define id_A->1, id_B -> 2... to easier reading
  ! Species names
  keys(1) = "A"
  keys(2) = "B"
  keys(3) = "C"
  keys(4) = "D"

  !Init concentrations
  state(1) = 1.0
  state(2) = 0.0
  state(3) = 0.0
  state(4) = conc_D

  test_driver => unit_test_driver_t(size_state)

  call test_driver%init_test(input_file_path, state, temp, pressure, keys)

  !todo: we can get the time_step from unit_test_driver
  do i_time=1, NUM_TIME_STEP

    ! Get the analytic conc
    time = i_time * time_step
    state_analytic(1) = state(1) * exp(-(k1)*time)
    state_analytic(2) = state(1) * (k1/(k2-k1)) * &
            (exp(-k1*time) - exp(-k2*time))
    state_analytic(3) = state(1) * &
            (1.0 + (k1*exp(-k2*time) - k2*exp(-k1*time))/(k2-k1))
    state_analytic(4) = state(4)

    call test_driver%solve_and_compare(time_step, state_analytic)

  end do

  !call deallocate_data()

  !check if correct and enjoy :)
















#if 0

  contains

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

#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program unit_test_rxn_arrhenius
