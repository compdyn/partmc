! Copyright (C) 2019 Christian Guzman and Matt Dawson
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_unit_test_rxn_aqueous_equilibrium module
module pmc_unit_test_rxn_aqueous_equilibrium
  !todo: fix multi-cell
  !TODO: Fix strange missing species F error
  use pmc_camp_core
  use pmc_camp_state
  use pmc_mpi
  use pmc_unit_test_data
  use pmc_util
  use pmc_aero_rep_factory
  use pmc_aero_rep_single_particle
  use pmc_aero_rep_modal_binned_mass

  implicit none
  private

  public :: unit_test_rxn_aqueous_equilibrium_t

  !> Number of available initial states
  integer(kind=i_kind), parameter :: NUM_INIT_STATES = 5
  !> Number of species to evaluate
  integer(kind=i_kind), parameter :: NUM_SPEC = 12
  !> Number of time steps
  integer(kind=i_kind), parameter :: NUM_TIME_STEPS_C = 100
  !> Time step size (s)
  real(kind=dp), parameter :: TIME_STEP_SIZE_C = 1.0

  !> The aqueous_equilibrium rxn unit test
  type, extends(unit_test_data_t) :: unit_test_rxn_aqueous_equilibrium_t
    private
    !> Scenario
    integer(kind=i_kind) :: scenario = UNIT_TEST_SCENARIO_
    !> Species ids
    integer(kind=i_kind) ::  idx_A, idx_B, idx_C, idx_BC_act, idx_D, idx_E, idx_F, idx_G, &
            idx_H, idx_H2O, i_sect_unused, i_sect_the_mode
    !> Species initial concentrations
    real(kind=dp), dimension(NUM_INIT_STATES) :: &
            conc_A = (/ 13.5d0, 14.5d0, 15.5d0, 16.5d0, 17.5d0 /)
    !> Species initial concentrations
    real(kind=dp), dimension(NUM_INIT_STATES) :: &
            conc_B = (/ 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0 /)
    !> Species initial concentrations
    real(kind=dp), dimension(NUM_INIT_STATES) :: &
            conc_C = (/ 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0 /)
    !> Species initial concentrations
    real(kind=dp), dimension(NUM_INIT_STATES) :: &
            conc_BC_act = (/ 0.5d0, 0.5d0, 0.5d0, 0.5d0, 0.5d0 /)
    !> Species initial concentrations
    real(kind=dp), dimension(NUM_INIT_STATES) :: &
            conc_D = (/ 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0 /)
    !> Species initial concentrations
    real(kind=dp), dimension(NUM_INIT_STATES) :: &
            conc_E = (/ 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0 /)
    !> Species initial concentrations
    real(kind=dp), dimension(NUM_INIT_STATES) :: &
            conc_F = (/ 8.0d0, 8.0d0, 8.0d0, 8.0d0, 8.0d0 /)
    !> Species initial concentrations
    real(kind=dp), dimension(NUM_INIT_STATES) :: &
            conc_G = (/ 12.0d0, 12.0d0, 12.0d0, 12.0d0, 12.0d0 /)
    !> Species initial concentrations
    real(kind=dp), dimension(NUM_INIT_STATES) :: &
            conc_H = (/ 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0 /)
    !> Species initial concentrations
    real(kind=dp), dimension(NUM_INIT_STATES) :: &
            conc_H2O = (/ 1.4d-2, 1.4d-2, 1.4d-2, 1.4d-2, 1.4d-2 /)
    !> Initial temperatures
    real(kind=dp), dimension(NUM_INIT_STATES) :: &
            temperature = (/ 272.5d0, 248.3d0, 301.3d0, 276.0d0, 245.1d0 /)
    !> Initial pressures
    real(kind=dp), dimension(NUM_INIT_STATES) :: &
            pressure = (/ 101253.3d0, 92834.4d0, 90328.4d0, 100495.4d0, 96943.4d0 /)
    !> Flag inidicating local data is set up
    logical :: is_initialized = .false.
    !> For setting rates
    type(aero_rep_factory_t) :: aero_rep_factory
    type(aero_rep_update_data_modal_binned_mass_GMD_t) :: update_data_GMD
    type(aero_rep_update_data_modal_binned_mass_GSD_t) :: update_data_GSD
  contains
    !> Initialize the object data
    procedure :: initialize
    !> Get the name of the input file for the test
    procedure :: input_file_name
    !> Get the name of a file to output results to
    procedure :: output_file_name
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
    !> Output the results for a given cell
    procedure :: output_results
    !> Calculate true concentrations and output or analyze model results
    procedure, private :: analyze_or_output
    !> Determine the number of bytes required to pack the object onto a buffer
    procedure :: pack_size
    !> Pack the object onto a buffer, advancing position
    procedure :: bin_pack
    !> Unpack an object from a buffer, advancing position
    procedure :: bin_unpack
  end type unit_test_rxn_aqueous_equilibrium_t

  !> Constructor for the aqueous_equilibrium test
  interface unit_test_rxn_aqueous_equilibrium_t
    procedure :: constructor
  end interface unit_test_rxn_aqueous_equilibrium_t

contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Constructor for the aqueous_equilibrium unit test
  function constructor() result(new_obj)

    !> A new unit test object
    type(unit_test_rxn_aqueous_equilibrium_t), pointer :: new_obj

    allocate(new_obj)

  end function constructor

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Initialize the object data
  subroutine initialize( this, camp_core )

    use pmc_aero_rep_data
    use pmc_rxn_data
    use pmc_rxn_factory
    use pmc_mechanism_data

    !> Unit test data
    class(unit_test_rxn_aqueous_equilibrium_t), intent(inout) :: this
    !> CAMP core
    class(camp_core_t), intent(in) :: camp_core
    !> rxn data
    class(rxn_data_t), pointer :: rxn

    character(len=:), allocatable :: key, str_val, rep_name, idx_prefix
    integer(kind=i_kind) ::  i_mech_rxn_rain, i_mech_rxn_cloud, i_rxn
    real(kind=dp) :: rate_rain, rate_cloud
    class(aero_rep_data_t), pointer :: aero_rep_ptr
    type(string_t), allocatable :: unique_names(:)

    !> For setting rates
    type(rxn_factory_t) :: rxn_factory

    call assert( 654980185, .not. this%is_initialized )

    ! Find the aerosol representation
    key = "my aero rep 2"
    call assert(260845179, camp_core%get_aero_rep(key, aero_rep_ptr))

    ! Set the aerosol representation id
    if (this%scenario.eq.2) then
      select type (aero_rep_ptr)
      type is (aero_rep_modal_binned_mass_t)
        call this%aero_rep_factory%initialize_update_data( aero_rep_ptr, &
                this%update_data_GMD)
        call this%aero_rep_factory%initialize_update_data( aero_rep_ptr, &
                this%update_data_GSD)
        call assert_msg(126380597, &
                aero_rep_ptr%get_section_id("unused mode", this%i_sect_unused), &
                "Could not get section id for the unused mode")
        call assert_msg(573748443, &
                aero_rep_ptr%get_section_id("the mode", this%i_sect_the_mode), &
                "Could not get section id for the unused mode")
      class default
        call die_msg(403591539, "Wrong aero rep type")
      end select

    end if

    ! Get species indices
    if (this%scenario.eq.1) then
      idx_prefix = ""
      print*, "scenario 1"
    else if (this%scenario.eq.2) then
      idx_prefix = "the mode."
      print*, "scenario 2"
    end if
    key = idx_prefix//"aqueous aerosol.A"
    this%idx_A = aero_rep_ptr%spec_state_id(key);
    key = idx_prefix//"aqueous aerosol.B"
    this%idx_B = aero_rep_ptr%spec_state_id(key);
    key = idx_prefix//"aqueous aerosol.C"
    this%idx_C = aero_rep_ptr%spec_state_id(key);
    key = idx_prefix//"aqueous aerosol.B-C"
    this%idx_BC_act = aero_rep_ptr%spec_state_id(key);
    key = idx_prefix//"aqueous aerosol.D"
    this%idx_D = aero_rep_ptr%spec_state_id(key);
    key = idx_prefix//"aqueous aerosol.E"
    this%idx_E = aero_rep_ptr%spec_state_id(key);
    key = idx_prefix//"organic aerosol.F"
    this%idx_F = aero_rep_ptr%spec_state_id(key);
    key = idx_prefix//"organic aerosol.G"
    this%idx_G = aero_rep_ptr%spec_state_id(key);
    key = idx_prefix//"organic aerosol.H"
    this%idx_H = aero_rep_ptr%spec_state_id(key);
    key = idx_prefix//"organic aerosol.H2O_aq"
    this%idx_H2O = aero_rep_ptr%spec_state_id(key);

    ! Make sure all the species were found
    call assert_msg( 204555021, this%idx_A.gt.0, "Missing species A" )
    call assert_msg( 883729391, this%idx_B.gt.0, "Missing species B" )
    call assert_msg( 883729392, this%idx_C.gt.0, "Missing species C" )
    call assert_msg( 883729393, this%idx_BC_act.gt.0, "Missing species BC_act" )
    call assert_msg( 883729394, this%idx_D.gt.0, "Missing species D" )
    call assert_msg( 883729395, this%idx_E.gt.0, "Missing species E" )
    call assert_msg( 883729396, this%idx_F.gt.0, "Missing species F" )
    call assert_msg( 883729397, this%idx_G.gt.0, "Missing species G" )
    call assert_msg( 883729398, this%idx_H.gt.0, "Missing species H" )
    call assert_msg( 883729399, this%idx_H2O.gt.0, "Missing species H2O" )

    ! Flag the test as being initialized
    this%is_initialized = .true.

  end subroutine initialize

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Input file name
  function input_file_name(this)

    !> Input file name
    character(len=:), allocatable :: input_file_name
    !> Unit test data
    class(unit_test_rxn_aqueous_equilibrium_t), intent(in) :: this

    ! Get the aqueous_equilibrium reaction mechanism json file
    if (this%scenario.eq.1) then
      input_file_name = 'rxn_aqueous_equilibrium_config.json'
    else if (this%scenario.eq.2) then
      input_file_name = 'rxn_aqueous_equilibrium_config_2.json'
    end if

  end function input_file_name

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Output file name
  function output_file_name(this)

    !> Output file name
    character(len=:), allocatable :: output_file_name
    !> Unit test data
    class(unit_test_rxn_aqueous_equilibrium_t), intent(in) :: this

    ! Get the aqueous_equilibrium reaction mechanism json file
    if (this%scenario.eq.1) then
      output_file_name = "rxn_aqueous_equilibrium_results.txt"
    else if (this%scenario.eq.2) then
      output_file_name = "rxn_aqueous_equilibrium_results_2.txt"
    end if

  end function output_file_name

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Number of unique initial states available from the unit test
  function num_unique_states(this)

    !> Number of unique states
    integer(kind=i_kind) :: num_unique_states
    !> Unit test data
    class(unit_test_rxn_aqueous_equilibrium_t), intent(in) :: this

    num_unique_states = NUM_INIT_STATES

  end function num_unique_states

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Initialize a camp_state_t object based on a given index
  subroutine initialize_state(this, grid_cell_id, camp_core, camp_state, &
          unique_state_id)

    !> Unit test data
    class(unit_test_rxn_aqueous_equilibrium_t), intent(inout) :: this
    !> Grid cell id
    integer(kind=i_kind), intent(in) :: grid_cell_id
    !> CAMP core
    class(camp_core_t), intent(in) :: camp_core
    !> Grid cell state
    class(camp_state_t), intent(inout) :: camp_state
    !> Unique state id
    integer(kind=i_kind), intent(in) :: unique_state_id

    real(kind=dp) :: rate_rain, rate_cloud

    ! Make sure the test data is initialized
    call assert(345723133, this%is_initialized)

    ! Update the GMD and GSD for the aerosol modes
    if (this%scenario.eq.2) then
      ! unused mode
      call this%update_data_GMD%set_GMD(this%i_sect_unused, 1.2d-6)
      call this%update_data_GSD%set_GSD(this%i_sect_the_mode, 1.2d0)
      call camp_core%update_aero_rep_data(this%update_data_GMD)
      call camp_core%update_aero_rep_data(this%update_data_GSD)
      ! the mode
      call this%update_data_GMD%set_GMD(this%i_sect_unused, 9.3d-7)
      call this%update_data_GSD%set_GSD(this%i_sect_the_mode, 0.9d0)
      call camp_core%update_aero_rep_data(this%update_data_GMD)
      call camp_core%update_aero_rep_data(this%update_data_GSD)
    end if

    ! Set the environmental conditions
    camp_state%env_state%temp     = this%temperature( unique_state_id )
    camp_state%env_state%pressure = this%pressure(    unique_state_id )
    call camp_state%update_env_state( )

    ! Set the species concentrations
    camp_state%state_var( this%idx_A ) = this%conc_A( unique_state_id )
    camp_state%state_var( this%idx_B ) = this%conc_B( unique_state_id )
    camp_state%state_var( this%idx_C ) = this%conc_C( unique_state_id )
    camp_state%state_var( this%idx_BC_act ) = this%conc_BC_act( unique_state_id )
    camp_state%state_var( this%idx_D ) = this%conc_D( unique_state_id )
    camp_state%state_var( this%idx_E ) = this%conc_E( unique_state_id )
    camp_state%state_var( this%idx_F ) = this%conc_F( unique_state_id )
    camp_state%state_var( this%idx_G ) = this%conc_G( unique_state_id )
    camp_state%state_var( this%idx_H ) = this%conc_H( unique_state_id )
    camp_state%state_var( this%idx_H2O ) = this%conc_H2O( unique_state_id )

  end subroutine initialize_state

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Number of time steps for this test
  function num_time_steps(this)

    !> Number of time steps
    integer(kind=i_kind) :: num_time_steps
    !> Unit test data
    class(unit_test_rxn_aqueous_equilibrium_t), intent(in) :: this

    num_time_steps = NUM_TIME_STEPS_C

  end function num_time_steps

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Time step size for this test (s)
  function time_step_size(this)

    !> Time step size (s)
    real(kind=dp) :: time_step_size
    !> Unit test data
    class(unit_test_rxn_aqueous_equilibrium_t), intent(in) :: this

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
    class(unit_test_rxn_aqueous_equilibrium_t), intent(inout) :: this
    !> CAMP core
    class(camp_core_t), intent(in) :: camp_core
    !> Grid cell state
    class(camp_state_t), intent(in) :: camp_state
    !> Unique state id
    integer(kind=i_kind), intent(in) :: unique_state_id
    !> Lastest time step solved
    integer(kind=i_kind), intent(in) :: model_time_step

    passed = this%analyze_or_output(camp_core, camp_state, unique_state_id, &
            model_time_step)

  end function analyze_state

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Output the results for a given grid cell
  subroutine output_results(this, camp_core, camp_state, unique_state_id, &
          model_time_step, output_file_unit)

    !> Unit test data
    class(unit_test_rxn_aqueous_equilibrium_t), intent(inout) :: this
    !> CAMP core
    class(camp_core_t), intent(in) :: camp_core
    !> Grid cell state
    class(camp_state_t), intent(in) :: camp_state
    !> Unique state id
    integer(kind=i_kind), intent(in) :: unique_state_id
    !> Lastest time step solved
    integer(kind=i_kind), intent(in) :: model_time_step
    !> Output file unit
    integer(kind=i_kind), intent(in) :: output_file_unit

    logical :: passed

    passed = this%analyze_or_output(camp_core, camp_state, unique_state_id, &
            model_time_step, output_file_unit)

  end subroutine output_results

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Private function to handle calculations for analysis or output
  function analyze_or_output(this, camp_core, camp_state, unique_state_id, &
          model_time_step, output_file_unit) result (passed)

    use pmc_constants

    !> Flag indicating whether the tests passed
    logical :: passed
    !> Unit test data
    class(unit_test_rxn_aqueous_equilibrium_t), intent(inout) :: this
    !> CAMP core
    class(camp_core_t), intent(in) :: camp_core
    !> Grid cell state
    class(camp_state_t), intent(in) :: camp_state
    !> Unique state id
    integer(kind=i_kind), intent(in) :: unique_state_id
    !> Lastest time step solved
    integer(kind=i_kind), intent(in) :: model_time_step
    !> Output file unit
    integer(kind=i_kind), intent(in), optional :: output_file_unit

    real(kind=dp) :: temperature, pressure
    real(kind=dp), dimension(NUM_SPEC) :: init_conc
    real(kind=dp), dimension(NUM_SPEC) :: true_conc
    real(kind=dp), dimension(NUM_SPEC) :: model_conc
    real(kind=dp) :: Keq_1, Keq_2, Keq_3, k1_forward, &
            k2_forward, k3_forward, k1_reverse, k2_reverse, k3_reverse, &
            total_init, equil_A, equil_B, equil_C, equil_D, equil_E, &
            equil_F, equil_G, equil_H, x, x0, time
    integer(kind=i_kind) :: i_spec

    passed = .true.

    ! Initial environmental state
    temperature    = this%temperature( unique_state_id )
    pressure       = this%pressure(    unique_state_id )

    ! Initial chemical state
    init_conc( 1 ) = this%conc_A(      unique_state_id )
    init_conc( 2 ) = this%conc_B(      unique_state_id )
    init_conc( 3 ) = this%conc_C(      unique_state_id )
    init_conc( 4 ) = this%conc_BC_act(      unique_state_id )
    init_conc( 5 ) = this%conc_D(      unique_state_id )
    init_conc( 6 ) = this%conc_E(      unique_state_id )
    init_conc( 7 ) = this%conc_F(      unique_state_id )
    init_conc( 8 ) = this%conc_G(      unique_state_id )
    init_conc( 9 ) = this%conc_H(      unique_state_id )
    init_conc( 10 ) = this%conc_H2O(      unique_state_id )

    ! Calculate the rate constants
    ! Henry's Law equilibrium constants (M/ppm)
    Keq_1 = 1.14d-2 * exp(2300.0d0 * (1.0d0/temperature - 1.0d0/298.0d0)) ! (M^2/M^2)
    Keq_2 = 12.3                                                   ! (M^3/M^2)
    Keq_3 = 2.35 * exp(1245.7d0 * (1.0d0/temperature - 1.0d0/298.0d0))    ! (M/M)

    ! Calculate the forward and reverse rate constants for reaction 1
    k1_reverse = 0.32d0 * init_conc( 4 )                           ! (1/M/s)
    k1_forward = Keq_1 * k1_reverse                                ! (1/M/s)

    ! Calculate the forward and reverse rate constants for reaction 2
    k2_reverse = 3.25e-3                                           ! (1/M/M/s)
    k2_forward = Keq_2 * k2_reverse                                ! (1/M/s)

    ! Calculate the forward and reverse rate constants for reaction 3
    k3_reverse = 1.56e-4                                           ! (1/s)
    k3_forward = Keq_3 * k3_reverse                                ! (1/s)

    ! Determine the equilibrium concentrations (ug/m3)
    !
    ! Reaction 1 (equil values in M)
    !
    ! K_eq = ([B][C])/([A]^2)
    ! [Total] = [A]i
    ! [B] = [C] = [Total] * (sqrt(1/Keq)-2)/(1/Keq-4)
    ! [A] = [Total] - [B] - [C]
    total_init = init_conc( 1 )/init_conc( 10 ) * 1000.0d0/48.0d0
    equil_B = (total_init * (sqrt(1.0d0/Keq_1)-2.0d0) / (1.0d0/Keq_1-4.0d0))
    equil_C = (total_init * (sqrt(1.0d0/Keq_1)-2.0d0) / (1.0d0/Keq_1-4.0d0))
    equil_A = (total_init * (1.0d0 - 2.0d0*(sqrt(1.0d0/Keq_1)-2.0d0) / &
            (1.0d0/Keq_1-4.0d0)))

    ! Reaction 2
    !
    ! K_eq = [F]^2/([D][E])
    ! [Total] = [F]i
    ! [D] = [E] = [Total] * (sqrt(Keq)-2)/(Keq-4)
    ! [F] = [Total] - [D] - [E]
    total_init = init_conc( 7 )/init_conc( 10 ) * 1000.0d0/28.0d0
    equil_D = (total_init * (sqrt(Keq_2)-2.0d0) / (Keq_2-4.0d0)) * &
            init_conc( 10 ) * 27.6d0 / 1000.0d0
    equil_E = (total_init * (sqrt(Keq_2)-2.0d0) / (Keq_2-4.0d0)) * &
            init_conc( 10 ) * 202.4d0 / 1000.0d0
    equil_F = (total_init * (1.0d0 - 2.0d0*(sqrt(Keq_2)-2.0d0) / &
            (Keq_2-4.0d0))) * init_conc( 10 ) * 28.0d0 / 1000.0d0

    ! Reaction 3
    !
    ! K_eq = [G]/[H]
    ! [Total] = [G]i
    ! [G] = [Total] / (1 + 1/Keq)
    ! [G] = [Total] / (Keq + 1)
    total_init = init_conc( 8 )/init_conc( 10 ) * 1000.0d0/35.67d0
    equil_H = (total_init / (1.0d0 + 1.0d0/Keq_3)) * &
            init_conc( 10 ) * 284.2 / 1000.0d0
    equil_G = (total_init / (Keq_3 + 1.0d0)) * &
            init_conc( 10 ) * 35.67d0 / 1000.0d0


    ! Get the current model time
    time = model_time_step * this%time_step_size( )

    ! Calculate the true current concentrations
    ! Two-reactant, two-product reactions
    ! FIXME This is not the analytic answer (waiting to use MatLab to solve)
    ! x = [A] - [A_eq]
    ! x0 = [A_init] - [A_eq]
    ! [A] = x + [A_eq] = x0exp(-t/tau) + [A_eq]
    ! 1/tau = k_f + k_b
    ! [A] = ([A_init] - [A_eq]) * exp(-t *(k_f + k_b)) + [A_eq]

    x0 = init_conc( 1 ) * 1000.0d0 / 48.0d0 / init_conc( 10 )
    x = 2.0 / ((1.0d0+2.0d0/x0)* &
            exp(-2.0d0*(k1_forward*equil_A+k1_reverse*equil_B)*time)-1.0d0)
    true_conc( 1 ) = (equil_A + x) * init_conc( 10 ) * 48.0d0 / 1000.0d0
    true_conc( 2 ) = (equil_B - x) * init_conc( 10 ) * 32.67d0 / 1000.0d0
    true_conc( 3 ) = (equil_C - x) * init_conc( 10 ) * 114.3d0 / 1000.0d0
    true_conc( 5 ) = (equil_C - x) * init_conc( 10 ) * 114.3d0 / 1000.0d0
    true_conc( 6 ) = (init_conc( 6 ) - equil_E) * &
            exp(-time * (k2_forward + k2_reverse)) + equil_E
    true_conc( 7 ) = (init_conc( 7 ) - equil_F) * &
            exp(-time * (k2_forward + k2_reverse)) + equil_F

    ! One-reactant, one-product reaction
    ! x = [G] - [G_eq]
    ! x0 = [G_init] - [G_eq]
    ! [G] = x + [G_eq] = x0 * exp(-t/tau) + [G_eq]
    ! 1/tau = k_f + f_r
    ! [G] = ([G_init] - [G_eq]) * exp(-t * (k_f + k_r)) + [G_eq]
    ! [H] = ([H_init] - [H_eq]) * exp(-t * (k_f + k_r)) + [H_eq]
    true_conc( 8 ) = (init_conc( 8 )- equil_G) * &
            exp(-time * (k3_forward + k3_reverse)) + equil_G
    true_conc( 9 ) = (init_conc( 9 ) - equil_H) * &
            exp(-time * (k3_forward + k3_reverse)) + equil_H

    ! BC_act & H20
    true_conc( 4 ) = init_conc( 4 )
    true_conc( 10 ) = init_conc( 10 )

    ! Get the model results
    model_conc( 1 ) = camp_state%state_var( this%idx_A )
    model_conc( 2 ) = camp_state%state_var( this%idx_B )
    model_conc( 3 ) = camp_state%state_var( this%idx_C )
    model_conc( 4 ) = camp_state%state_var( this%idx_BC_act )
    model_conc( 5 ) = camp_state%state_var( this%idx_D )
    model_conc( 6 ) = camp_state%state_var( this%idx_E )
    model_conc( 7 ) = camp_state%state_var( this%idx_F )
    model_conc( 8 ) = camp_state%state_var( this%idx_G )
    model_conc( 9 ) = camp_state%state_var( this%idx_H )
    model_conc( 10 ) = camp_state%state_var( this%idx_H2O )

    ! Analyze or output results
    if( present( output_file_unit ) ) then

      ! Output the results
      write( OUTPUT_FILE_UNIT, * ) time, &
              ' ', true_conc( 1 ), ' ', model_conc( 1 ), &
              ' ', true_conc( 2 ), ' ', model_conc( 2 ), &
              ' ', true_conc( 3 ), ' ', model_conc( 3 ), &
              ' ', true_conc( 4 ), ' ', model_conc( 4 ), &
              ' ', true_conc( 5 ), ' ', model_conc( 5 ), &
              ' ', true_conc( 6 ), ' ', model_conc( 6 ), &
              ' ', true_conc( 7 ), ' ', model_conc( 7 ), &
              ' ', true_conc( 8 ), ' ', model_conc( 8 ), &
              ' ', true_conc( 9 ), ' ', model_conc( 9 ), &
              ' ', true_conc( 10 ), ' ', model_conc( 10 )

    else

      ! Analyze the results (single-particle only)
      if (this%scenario.eq.1) then
        do i_spec = 1, NUM_SPEC
          call assert_msg( 334093929, &
                  almost_equal( model_conc( i_spec ), true_conc( i_spec ), 1.0d-2 ) &
                          .or.( model_conc( i_spec ) .lt. 1.0d-5 * init_conc( i_spec ) .and. &
                          true_conc( i_spec ) .lt. 1.0d-5 * init_conc( i_spec ) ), &
                  "time: "//trim( to_string( time ) )//"; scenario: "//&
                          trim(to_string( unique_state_id ) )//"; species: "// &
                          trim(to_string( i_spec ) )//"; mod: "// &
                          trim(to_string( model_conc( i_spec ) ) )//"; true: "// &
                          trim(to_string( true_conc( i_spec ) ) ) )
        end do
      end if
    end if

  end function analyze_or_output

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Determine the number of bytes required to pack the object onto a buffer
  integer(kind=i_kind) function pack_size(this, comm)

    !> Unit test data
    class(unit_test_rxn_aqueous_equilibrium_t), intent(in) :: this
    !> MPI communicator
    integer, intent(in), optional :: comm

#ifdef PMC_USE_MPI
    integer :: l_comm

    if (present(comm)) then
      l_comm = comm
    else
      l_comm = MPI_COMM_WORLD
    endif

    pack_size = &
            pmc_mpi_pack_size_integer(this%idx_A, l_comm) + &
                    pmc_mpi_pack_size_integer(this%idx_B, l_comm) + &
                    pmc_mpi_pack_size_integer(this%idx_C, l_comm) + &
                    pmc_mpi_pack_size_integer(this%idx_BC_act, l_comm) + &
                    pmc_mpi_pack_size_integer(this%idx_D, l_comm) + &
                    pmc_mpi_pack_size_integer(this%idx_E, l_comm) + &
                    pmc_mpi_pack_size_integer(this%idx_F, l_comm) + &
                    pmc_mpi_pack_size_integer(this%idx_G, l_comm) + &
                    pmc_mpi_pack_size_integer(this%idx_H, l_comm) + &
                    pmc_mpi_pack_size_integer(this%idx_H2O, l_comm) + &
                    pmc_mpi_pack_size_integer(this%i_sect_unused, l_comm) + &
                    pmc_mpi_pack_size_integer(this%i_sect_the_mode, l_comm) + &
                    this%update_data_GMD%pack_size() + &
                    this%update_data_GSD%pack_size()
#else
    pack_size = 0
#endif

  end function pack_size

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Pack the object onto a buffer, advancing position
  subroutine bin_pack(this, buffer, pos, comm)

    !> Unit test data
    class(unit_test_rxn_aqueous_equilibrium_t), intent(in) :: this
    !> Memory buffer
    character, intent(inout) :: buffer(:)
    !> Current buffer position
    integer, intent(inout) :: pos
    !> MPI communicator
    integer, intent(in), optional ::comm

#ifdef PMC_USE_MPI
    integer :: prev_position, l_comm

    if (present(comm)) then
      l_comm = comm
    else
      l_comm = MPI_COMM_WORLD
    endif

    call assert_msg(282594670, this%is_initialized, &
            "Trying to pack an uninitialized unit test on a buffer")

    prev_position = pos
    call pmc_mpi_pack_integer(buffer, pos, this%idx_A, l_comm)
    call pmc_mpi_pack_integer(buffer, pos, this%idx_B, l_comm)
    call pmc_mpi_pack_integer(buffer, pos, this%idx_C, l_comm)
    call pmc_mpi_pack_integer(buffer, pos, this%idx_BC_act, l_comm)
    call pmc_mpi_pack_integer(buffer, pos, this%idx_D, l_comm)
    call pmc_mpi_pack_integer(buffer, pos, this%idx_E, l_comm)
    call pmc_mpi_pack_integer(buffer, pos, this%idx_F, l_comm)
    call pmc_mpi_pack_integer(buffer, pos, this%idx_G, l_comm)
    call pmc_mpi_pack_integer(buffer, pos, this%idx_H, l_comm)
    call pmc_mpi_pack_integer(buffer, pos, this%idx_H2O, l_comm)
    call pmc_mpi_pack_integer(buffer, pos, this%i_sect_unused, l_comm)
    call pmc_mpi_pack_integer(buffer, pos, this%i_sect_the_mode, l_comm)
    call this%update_data_GMD%bin_pack(buffer, pos)
    call this%update_data_GSD%bin_pack(buffer, pos)
    call assert(897212942, &
            pos - prev_position <= this%pack_size(l_comm))
#endif

  end subroutine bin_pack

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Unpack an object from a buffer, advancing position
  subroutine bin_unpack(this, buffer, pos, comm)

    !> Unit test data
    class(unit_test_rxn_aqueous_equilibrium_t), intent(inout) :: this
    !> Memory buffer
    character, intent(inout) :: buffer(:)
    !> Current buffer position
    integer, intent(inout) :: pos
    !> MPI communicator
    integer, intent(in), optional :: comm

#ifdef PMC_USE_MPI
    integer :: prev_position, l_comm

    if (present(comm)) then
      l_comm = comm
    else
      l_comm = MPI_COMM_WORLD
    endif

    call assert_msg(768100732, .not. this%is_initialized, &
            "Trying to overwrite an initialized unit test object")

    prev_position = pos
    call pmc_mpi_unpack_integer(buffer, pos, this%idx_A, l_comm)
    call pmc_mpi_unpack_integer(buffer, pos, this%idx_B, l_comm)
    call pmc_mpi_unpack_integer(buffer, pos, this%idx_C, l_comm)
    call pmc_mpi_unpack_integer(buffer, pos, this%idx_BC_act, l_comm)
    call pmc_mpi_unpack_integer(buffer, pos, this%idx_D, l_comm)
    call pmc_mpi_unpack_integer(buffer, pos, this%idx_E, l_comm)
    call pmc_mpi_unpack_integer(buffer, pos, this%idx_F, l_comm)
    call pmc_mpi_unpack_integer(buffer, pos, this%idx_G, l_comm)
    call pmc_mpi_unpack_integer(buffer, pos, this%idx_H, l_comm)
    call pmc_mpi_unpack_integer(buffer, pos, this%idx_H2O, l_comm)
    call pmc_mpi_unpack_integer(buffer, pos, this%i_sect_unused, l_comm)
    call pmc_mpi_unpack_integer(buffer, pos, this%i_sect_the_mode, l_comm)
    call this%update_data_GMD%bin_unpack(buffer, pos)
    call this%update_data_GSD%bin_unpack(buffer, pos)
    call assert(466926084, &
            pos - prev_position <= this%pack_size(l_comm))

    this%is_initialized = .true.
#endif

  end subroutine bin_unpack

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_unit_test_rxn_aqueous_equilibrium
