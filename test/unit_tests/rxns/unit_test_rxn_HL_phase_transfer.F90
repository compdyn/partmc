! Copyright (C) 2019 Christian Guzman and Matt Dawson
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_unit_test_rxn_HL_phase_transfer module
module pmc_unit_test_rxn_HL_phase_transfer
  !todo: fix multi-cell
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

  public :: unit_test_rxn_HL_phase_transfer_t

  !> Number of available initial states
  integer(kind=i_kind), parameter :: NUM_INIT_STATES = 5
  !> Number of species to evaluate
  integer(kind=i_kind), parameter :: NUM_SPEC = 13
  !> Number of time steps
  integer(kind=i_kind), parameter :: NUM_TIME_STEPS_C = 100
  !> Time step size (s)
  real(kind=dp), parameter :: TIME_STEP_SIZE_C = 1.0

  !> The HL_phase_transfer rxn unit test
  type, extends(unit_test_data_t) :: unit_test_rxn_HL_phase_transfer_t
    private
    !> Scenario
    integer(kind=i_kind) :: scenario = UNIT_TEST_SCENARIO_
    !> Species ids
    integer(kind=i_kind) ::  idx_O3, idx_O3_aq, &
            idx_H2O2, idx_H2O2_aq, idx_H2O_aq, idx_HNO3, idx_HNO3_aq, &
            idx_NH3, idx_NH3_aq, idx_H2O, idx_Na_p, idx_Cl_m, idx_Ca_pp, i_sect_unused, i_sect_the_mode
    !> Species initial concentrations
    real(kind=dp), dimension(NUM_INIT_STATES) :: &
            conc_O3 = (/ 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0 /)
    !> Species initial concentrations
    real(kind=dp), dimension(NUM_INIT_STATES) :: &
            conc_O3_aq = (/ 1.0d-3, 1.0d-3, 1.0d-3, 1.0d-3, 1.0d-3 /)
    !> Species initial concentrations
    real(kind=dp), dimension(NUM_INIT_STATES) :: &
            conc_H2O2 = (/ 1.0d0, 2.0d0, 3.0d0, 4.0d0, 5.0d0 /)
    !> Species initial concentrations
    real(kind=dp), dimension(NUM_INIT_STATES) :: &
            conc_H2O2_aq = (/ 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0 /)
    !> Species initial concentrations
    real(kind=dp), dimension(NUM_INIT_STATES) :: &
            conc_HNO3 = (/ 1.0d0, 2.0d0, 3.0d0, 4.0d0, 5.0d0 /)
    !> Species initial concentrations
    real(kind=dp), dimension(NUM_INIT_STATES) :: &
            conc_HNO3_aq = (/ 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0 /)
    !> Species initial concentrations
    real(kind=dp), dimension(NUM_INIT_STATES) :: &
            conc_NH3 = (/ 1.0d0, 2.0d0, 3.0d0, 4.0d0, 5.0d0 /)
    !> Species initial concentrations
    real(kind=dp), dimension(NUM_INIT_STATES) :: &
            conc_NH3_aq = (/ 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0 /)
    !> Species initial concentrations
    real(kind=dp), dimension(NUM_INIT_STATES) :: &
            conc_Na_p = (/ 2.5d0, 2.5d0, 2.5d0, 2.5d0, 2.5d0 /)
    !> Species initial concentrations
    real(kind=dp), dimension(NUM_INIT_STATES) :: &
            conc_Cl_m = (/ 5.3d0, 5.3d0, 5.3d0, 5.3d0, 5.3d0 /)
    !> Species initial concentrations
    real(kind=dp), dimension(NUM_INIT_STATES) :: &
            conc_Ca_pp = (/ 1.3d0, 1.3d0, 1.3d0, 1.3d0, 1.3d0 /)
    !> Species initial concentrations
    real(kind=dp), dimension(NUM_INIT_STATES) :: &
            conc_H2O = (/ 3000.0d0, 3000.0d0, 3000.0d0, 3000.0d0, 3000.0d0 /)
    !> Species initial concentrations
    real(kind=dp), dimension(NUM_INIT_STATES) :: &
            conc_H2O_aq = (/ 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0 /)
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
    type(aero_rep_update_data_single_particle_radius_t) :: radius_update
    type(aero_rep_update_data_single_particle_number_t) :: number_update
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
  end type unit_test_rxn_HL_phase_transfer_t

  !> Constructor for the HL_phase_transfer test
  interface unit_test_rxn_HL_phase_transfer_t
    procedure :: constructor
  end interface unit_test_rxn_HL_phase_transfer_t

contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Constructor for the HL_phase_transfer unit test
  function constructor() result(new_obj)

    !> A new unit test object
    type(unit_test_rxn_HL_phase_transfer_t), pointer :: new_obj

    allocate(new_obj)

  end function constructor

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Initialize the object data
  subroutine initialize( this, camp_core )

    use pmc_aero_rep_data
    use pmc_rxn_data
    use pmc_rxn_factory
    use pmc_mechanism_data
    use pmc_chem_spec_data

    !> Unit test data
    class(unit_test_rxn_HL_phase_transfer_t), intent(inout) :: this
    !> CAMP core
    class(camp_core_t), intent(in) :: camp_core
    !> rxn data
    class(rxn_data_t), pointer :: rxn

    character(len=:), allocatable :: key, str_val, rep_name, idx_prefix
    integer(kind=i_kind) ::  i_mech_rxn_rain, i_mech_rxn_cloud, i_rxn
    real(kind=dp) :: rate_rain, rate_cloud
    class(aero_rep_data_t), pointer :: aero_rep_ptr
    type(chem_spec_data_t), pointer :: chem_spec_data
    type(string_t), allocatable :: unique_names(:)

    !> For setting rates
    type(rxn_factory_t) :: rxn_factory

    call assert( 654980185, .not. this%is_initialized )

    ! Find the aerosol representation
    key = "my aero rep 2"
    call assert(260845179, camp_core%get_aero_rep(key, aero_rep_ptr))

    ! Set the aerosol representation id
    if (this%scenario.eq.1) then
      select type (aero_rep_ptr)
      type is (aero_rep_single_particle_t)
        call this%aero_rep_factory%initialize_update_data( aero_rep_ptr, &
                this%radius_update)
        call this%aero_rep_factory%initialize_update_data( aero_rep_ptr, &
                this%number_update)
      class default
        call die_msg(403591539, "Wrong aero rep type")
      end select

    else if (this%scenario.eq.2) then
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
        call die_msg(290304323, "Wrong aero rep type")
      end select
    end if

    ! Get the chemical species data
    call assert(191714381, camp_core%get_chem_spec_data(chem_spec_data))

    ! Get species indices
    if (this%scenario.eq.1) then
      idx_prefix = ""
    else if (this%scenario.eq.2) then
      idx_prefix = "the mode."
    end if
    key = "O3"
    this%idx_O3 = chem_spec_data%gas_state_id(key);
    key = idx_prefix//"aqueous aerosol.O3_aq"
    this%idx_O3_aq = aero_rep_ptr%spec_state_id(key);
    key = "H2O2"
    this%idx_H2O2 = chem_spec_data%gas_state_id(key);
    key = idx_prefix//"aqueous aerosol.H2O2_aq"
    this%idx_H2O2_aq = aero_rep_ptr%spec_state_id(key);
    if (this%scenario.eq.2) then
      key = "HNO3"
      this%idx_HNO3 = chem_spec_data%gas_state_id(key);
      key = idx_prefix//"aqueous aerosol.HNO3_aq"
      this%idx_HNO3_aq = aero_rep_ptr%spec_state_id(key);
      key = "NH3"
      this%idx_NH3 = chem_spec_data%gas_state_id(key);
      key = idx_prefix//"aqueous aerosol.NH3_aq"
      this%idx_NH3_aq = aero_rep_ptr%spec_state_id(key);
      key = "H2O"
      this%idx_H2O = chem_spec_data%gas_state_id(key);
      key = idx_prefix//"aqueous aerosol.Na_p"
      this%idx_Na_p = aero_rep_ptr%spec_state_id(key);
      key = idx_prefix//"aqueous aerosol.Cl_m"
      this%idx_Cl_m = aero_rep_ptr%spec_state_id(key);
      key = idx_prefix//"aqueous aerosol.Ca_pp"
      this%idx_Ca_pp = aero_rep_ptr%spec_state_id(key);
    end if
    key = idx_prefix//"aqueous aerosol.H2O_aq"
    this%idx_H2O_aq = aero_rep_ptr%spec_state_id(key);

    ! Make sure the expected species are in the model
    call assert_msg( 204555021, this%idx_O3.gt.0, "Missing species O3" )
    call assert_msg( 883729398, this%idx_O3_aq.gt.0, "Missing species O3_aq" )
    call assert_msg( 883729398, this%idx_H2O2.gt.0, "Missing species H2O2" )
    call assert_msg( 883729398, this%idx_H2O2_aq.gt.0, "Missing species H2O2_aq" )
    if (this%scenario.eq.2) then
      call assert_msg( 883729398, this%idx_HNO3.gt.0, "Missing species HNO3" )
      call assert_msg( 883729398, this%idx_HNO3_aq.gt.0, "Missing species HNO3_aq" )
      call assert_msg( 883729398, this%idx_NH3.gt.0, "Missing species NH3" )
      call assert_msg( 883729398, this%idx_NH3_aq.gt.0, "Missing species NH3_aq" )
      call assert_msg( 883729398, this%idx_H2O.gt.0, "Missing species H2O")
    end if
    call assert_msg( 883729398, this%idx_H2O_aq.gt.0, "Missing species H2O")

    ! Flag the test as being initialized
    this%is_initialized = .true.

  end subroutine initialize

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Input file name
  function input_file_name(this)

    !> Input file name
    character(len=:), allocatable :: input_file_name
    !> Unit test data
    class(unit_test_rxn_HL_phase_transfer_t), intent(in) :: this

    ! Get the HL_phase_transfer reaction mechanism json file
    if (this%scenario.eq.1) then
      input_file_name = 'rxn_HL_phase_transfer_config.json'
    else if (this%scenario.eq.2) then
      input_file_name = 'rxn_HL_phase_transfer_config_2.json'
    end if

  end function input_file_name

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Output file name
  function output_file_name(this)

    !> Output file name
    character(len=:), allocatable :: output_file_name
    !> Unit test data
    class(unit_test_rxn_HL_phase_transfer_t), intent(in) :: this

    ! Get the HL_phase_transfer reaction mechanism json file
    if (this%scenario.eq.1) then
      output_file_name = "rxn_HL_phase_transfer_results.txt"
    else if (this%scenario.eq.2) then
      output_file_name = "rxn_HL_phase_transfer_results_2.txt"
    end if

  end function output_file_name

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Number of unique initial states available from the unit test
  function num_unique_states(this)

    !> Number of unique states
    integer(kind=i_kind) :: num_unique_states
    !> Unit test data
    class(unit_test_rxn_HL_phase_transfer_t), intent(in) :: this

    num_unique_states = NUM_INIT_STATES

  end function num_unique_states

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Initialize a camp_state_t object based on a given index
  subroutine initialize_state(this, grid_cell_id, camp_core, camp_state, &
          unique_state_id)

    !> Unit test data
    class(unit_test_rxn_HL_phase_transfer_t), intent(inout) :: this
    !> Grid cell id
    integer(kind=i_kind), intent(in) :: grid_cell_id
    !> CAMP core
    class(camp_core_t), intent(in) :: camp_core
    !> Grid cell state
    class(camp_state_t), intent(inout) :: camp_state
    !> Unique state id
    integer(kind=i_kind), intent(in) :: unique_state_id

    real(kind=dp) :: radius, number_conc

    ! Make sure the test data is initialized
    call assert(345723133, this%is_initialized)

    ! Calculate the radius and number concentration to use
    ! ( the real values for the modal representation cannot be calculated
    !   because the number concentrations change sligthly during the run
    !   but the Jacobian checker can be run as a check. )
    if (this%scenario.eq.1) then
      radius = 1.5e-5             ! radius (m)
      number_conc = 1.3e6         ! particle number concentration (#/cc)
    else if (this%scenario.eq.2) then
      ! radius (m)
      radius = 9.37e-7 / 2.0 * exp(9.0/2.0 * 0.9 * 0.9)
      ! number conc
      number_conc = 1.0 / (const%pi/6.0 * (9.37e-7)**3.0 * &
              exp(9.0/2.0 * 0.9 * 0.9))
      number_conc = number_conc * 1.0e-9 * (1.0e-3 + 1.4e-2)
    end if

    ! Update the aerosol representation (single-particle only)
    if (this%scenario.eq.1) then
      call this%radius_update%set_radius(radius)
      call this%number_update%set_number(number_conc)
      call camp_core%update_aero_rep_data(this%radius_update)
      call camp_core%update_aero_rep_data(this%number_update)
    end if

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

    if (this%scenario.eq.1) then
      this%conc_O3_aq( unique_state_id ) = 1.0e-3
      this%conc_H2O_aq( unique_state_id ) = 1.4e-2
    else if (this%scenario.eq.2) then
      this%conc_O3_aq( unique_state_id ) = 1.0
      this%conc_H2O_aq( unique_state_id ) = 0.0
    end if

    camp_state%state_var( this%idx_O3 ) = this%conc_O3( unique_state_id )
    camp_state%state_var( this%idx_O3_aq ) = this%conc_O3_aq( unique_state_id )
    camp_state%state_var( this%idx_H2O2 ) = this%conc_H2O2( unique_state_id )
    camp_state%state_var( this%idx_H2O2_aq ) = this%conc_H2O2_aq( unique_state_id )

    if (this%scenario.eq.2) then
      camp_state%state_var( this%idx_HNO3 ) = this%conc_HNO3( unique_state_id )
      camp_state%state_var( this%idx_HNO3_aq ) = this%conc_HNO3_aq( unique_state_id )
      camp_state%state_var( this%idx_NH3 ) = this%conc_NH3( unique_state_id )
      camp_state%state_var( this%idx_NH3_aq ) = this%conc_NH3_aq( unique_state_id )
      camp_state%state_var( this%idx_Na_p ) = this%conc_Na_p( unique_state_id )
      camp_state%state_var( this%idx_Cl_m ) = this%conc_Cl_m( unique_state_id )
      camp_state%state_var( this%idx_Ca_pp ) = this%conc_Ca_pp( unique_state_id )
      camp_state%state_var( this%idx_H2O ) = this%conc_H2O( unique_state_id )
    end if

    camp_state%state_var( this%idx_H2O_aq ) = this%conc_H2O_aq( unique_state_id )

  end subroutine initialize_state

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Number of time steps for this test
  function num_time_steps(this)

    !> Number of time steps
    integer(kind=i_kind) :: num_time_steps
    !> Unit test data
    class(unit_test_rxn_HL_phase_transfer_t), intent(in) :: this

    num_time_steps = NUM_TIME_STEPS_C

  end function num_time_steps

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Time step size for this test (s)
  function time_step_size(this)

    !> Time step size (s)
    real(kind=dp) :: time_step_size
    !> Unit test data
    class(unit_test_rxn_HL_phase_transfer_t), intent(in) :: this

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
    class(unit_test_rxn_HL_phase_transfer_t), intent(inout) :: this
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
    class(unit_test_rxn_HL_phase_transfer_t), intent(inout) :: this
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
    class(unit_test_rxn_HL_phase_transfer_t), intent(inout) :: this
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
    real(kind=dp) :: n_star, del_H, del_S, del_G, alpha, &
            crms, M_to_ppm, ugm3_to_ppm, K_eq_O3, K_eq_H2O2, k_O3_forward, &
            k_O3_backward, k_H2O2_forward, k_H2O2_backward, equil_O3, &
            equil_O3_aq, equil_H2O2, equil_H2O2_aq, time
    real(kind=dp), target :: radius, number_conc
    integer(kind=i_kind) :: i_spec

    passed = .true.

    init_conc(:) = 0.0
    true_conc(:) = 0.0
    model_conc(:)= 0.0

    ! Calculate the radius and number concentration to use
    ! ( the real values for the modal representation cannot be calculated
    !   because the number concentrations change sligthly during the run
    !   but the Jacobian checker can be run as a check. )
    if (this%scenario.eq.1) then
      radius = 1.5e-5             ! radius (m)
      number_conc = 1.3e6         ! particle number concentration (#/cc)
    else if (this%scenario.eq.2) then
      ! radius (m)
      radius = 9.37e-7 / 2.0 * exp(9.0/2.0 * 0.9 * 0.9)
      ! number conc
      number_conc = 1.0 / (const%pi/6.0 * (9.37e-7)**3.0 * &
              exp(9.0/2.0 * 0.9 * 0.9))
      number_conc = number_conc * 1.0e-9 * (1.0e-3 + 1.4e-2)
    end if

    ! Initial environmental state
    temperature    = this%temperature( unique_state_id )
    pressure       = this%pressure(    unique_state_id )

    ! Initial chemical state
    init_conc( 1 ) = this%conc_O3(      unique_state_id )
    init_conc( 2 ) = this%conc_O3_aq(      unique_state_id )
    init_conc( 3 ) = this%conc_H2O2(      unique_state_id )
    init_conc( 4 ) = this%conc_H2O2_aq(      unique_state_id )
    if (this%scenario.eq.2) then
      init_conc( 5 ) = this%conc_HNO3(      unique_state_id )
      init_conc( 6 ) = this%conc_HNO3_aq(      unique_state_id )
      init_conc( 7 ) = this%conc_NH3(      unique_state_id )
      init_conc( 8 ) = this%conc_NH3_aq(      unique_state_id )
      init_conc( 9 ) = this%conc_Na_p(      unique_state_id )
      init_conc( 10 ) = this%conc_Cl_m(      unique_state_id )
      init_conc( 11 ) = this%conc_Ca_pp(      unique_state_id )
      init_conc( 12 ) = this%conc_H2O(      unique_state_id )
    end if
    init_conc( 13 ) = this%conc_H2O_aq(      unique_state_id )

    ! Calculate the rate constants

    ! Henry's Law equilibrium constants (M/ppm)
    ! O3 HLC Equil Const (M/ppm)
    K_eq_O3 = 1.14d-2 * exp(2300.0d0 * (1.0d0/temperature - 1.0d0/298.0d0)) / 1.0d6
    ! H2O2 HLC Equil Const (M/ppm)
    K_eq_H2O2 = 1.025d5 * exp(6340.0d0 * (1.0d0/temperature - 1.0d0/298.0d0)) / 1.0d6

    ! Determine the M -> ppm conversion using the total aerosol water
    M_to_ppm = number_conc * 1.0d-3 * init_conc( 13 ) * &
            const%univ_gas_const * temperature / pressure

    ! O3 rate constants
    n_star = 1.89d0
    del_H = - 10.0d0 * (n_star-1.0d0) + &
            7.53*(n_star**(2.0d0/3.0d0)-1.0d0) - 1.0d0
    del_S = - 32.0d0 * (n_star-1.0d0) + &
            9.21*(n_star**(2.0d0/3.0d0)-1.0d0) - 1.3d0
    del_G = (del_H - temperature * del_S/1000.0d0) * 4184.0d0
    alpha = exp(-del_G/(const%univ_gas_const*temperature))
    alpha = alpha / (1.0d0 + alpha)
    crms = sqrt(8.0d0*const%univ_gas_const*temperature/(const%pi*48.0d0))
    k_O3_forward = number_conc * ((radius**2 / (3.0d0 * 1.48d-5) + &
            4.0d0 * radius / (3.0d0 * crms * alpha))**(-1))           ! (1/s)
    k_O3_backward = k_O3_forward / (K_eq_O3 * M_to_ppm)             ! (1/s)

    ! H2O2 rate constants
    n_star = 1.74d0
    del_H = - 10.0d0 * (n_star-1.0d0) + &
            7.53*(n_star**(2.0d0/3.0d0)-1.0d0) - 1.0d0
    del_S = - 32.0d0 * (n_star-1.0d0) + &
            9.21*(n_star**(2.0d0/3.0d0)-1.0d0) - 1.3d0
    del_G = (del_H - temperature * del_S/1000.0d0) * 4184.0d0
    alpha = exp(-del_G/(const%univ_gas_const*temperature))
    alpha = alpha / (1.0d0 + alpha)
    crms = sqrt(8.0d0*const%univ_gas_const*temperature/(const%pi*34.0d0))
    k_H2O2_forward = number_conc * ((radius**2 / (3.0d0 * 1.46d-5) + &
            4.0d0 * radius / (3.0d0 * crms * alpha))**(-1))           ! (1/s)
    k_H2O2_backward = k_H2O2_forward / (K_eq_H2O2 * M_to_ppm)       ! (1/s)

    ! Determine the equilibrium concentrations
    ! [A_gas] = [A_total] / (1 + 1/K_HL)
    ! [A_aero] = [A_total] / (K_HL + 1)
    ugm3_to_ppm = const%univ_gas_const * temperature / (48.0d0 * pressure)
    equil_O3 = (init_conc( 1 ) + &
            init_conc( 2 )*number_conc*ugm3_to_ppm) / &
            (K_eq_O3*M_to_ppm + 1.0d0)
    equil_O3_aq = (init_conc( 1 )/ugm3_to_ppm/number_conc + &
            init_conc( 2 )) / &
            (1.0d0 + 1.0d0/(K_eq_O3*M_to_ppm))

    ugm3_to_ppm = const%univ_gas_const * temperature / (34.0d0 * pressure)
    equil_H2O2 = (init_conc( 3 ) + &
            init_conc( 4 )*number_conc*ugm3_to_ppm) / &
            (K_eq_H2O2*M_to_ppm + 1.0d0)
    equil_H2O2_aq = (init_conc( 3 )/ugm3_to_ppm/number_conc + &
            init_conc( 4 )) / &
            (1.0d0 + 1.0d0/(K_eq_H2O2*M_to_ppm))

    ! Get the current model time
    time = model_time_step * this%time_step_size( )

    ! Calculate the true current concentrations
    ! x = [A_gas] - [A_eq_gas]
    ! x0 = [A_init_gas] - [A_eq_gas]
    ! [A_gas] = x + [A_eq_gas] = x0exp(-t/tau) + [A_eq_gas]
    ! 1/tau = k_f + k_b
    ! [A_gas] = ([A_init_gas] - [A_eq_gas])
    !     * exp(-t *(k_f + k_b)) + [A_eq_gas]
    ! [A_aero] = ([A_init_aero] - [A_eq_aero])
    !     * exp(-t * (k_f + k_b)) + [A_eq_aero]
    true_conc( 1 ) =(init_conc( 1 ) - equil_O3) * &
            exp(-time * (k_O3_forward + k_O3_backward)) + equil_O3
    true_conc( 2 ) = (init_conc( 2 ) - equil_O3_aq) * &
            exp(-time * (k_O3_forward + k_O3_backward)) + equil_O3_aq
    true_conc( 3 ) = (init_conc( 3 ) - equil_H2O2) * &
            exp(-time * (k_H2O2_forward + k_H2O2_backward)) + equil_H2O2
    true_conc( 4 ) = &
            (init_conc( 4 ) - equil_H2O2_aq) * &
            exp(-time * (k_H2O2_forward + k_H2O2_backward)) + equil_H2O2_aq


      true_conc( 5 ) = init_conc( 5 )
      true_conc( 6 ) = init_conc( 6 )
      true_conc( 7 ) = init_conc( 7 )
      true_conc( 8 ) = init_conc( 8 )
      true_conc( 9 ) = init_conc( 9 )
      true_conc( 10 ) = init_conc( 10 )
      true_conc( 11 ) = init_conc( 11 )
      true_conc( 12 ) = init_conc( 12 )
      true_conc( 13 ) = init_conc( 13 )


    ! Get the model results
    model_conc( 1 ) = camp_state%state_var( this%idx_O3 )
    model_conc( 2 ) = camp_state%state_var( this%idx_O3_aq )
    model_conc( 3 ) = camp_state%state_var( this%idx_H2O2 )
    model_conc( 4 ) = camp_state%state_var( this%idx_H2O2_aq )

    if (this%scenario.eq.1) then

      model_conc( 5 ) = true_conc( 5 )
      model_conc( 6 ) = true_conc( 6 )
      model_conc( 7 ) = true_conc( 7 )
      model_conc( 8 ) = true_conc( 8 )
      model_conc( 9 ) = true_conc( 9 )
      model_conc( 10 ) = true_conc( 10 )
      model_conc( 11 ) = true_conc( 11 )
      model_conc( 12 ) = true_conc( 12 )

    else if (this%scenario.eq.2) then
      model_conc( 5 ) = camp_state%state_var( this%idx_HNO3 )
      model_conc( 6 ) = camp_state%state_var( this%idx_HNO3_aq )
      model_conc( 7 ) = camp_state%state_var( this%idx_NH3 )
      model_conc( 8 ) = camp_state%state_var( this%idx_NH3_aq )
      model_conc( 9 ) = camp_state%state_var( this%idx_Na_p )
      model_conc( 10 ) = camp_state%state_var( this%idx_Cl_m )
      model_conc( 11 ) = camp_state%state_var( this%idx_Ca_pp )
      model_conc( 12 ) = camp_state%state_var( this%idx_H2O )
    end if

    model_conc( 13 ) = camp_state%state_var( this%idx_H2O_aq )

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
              ' ', true_conc( 10 ), ' ', model_conc( 10 ), &
              ' ', true_conc( 11 ), ' ', model_conc( 11 ), &
              ' ', true_conc( 12 ), ' ', model_conc( 12 ), &
              ' ', true_conc( 13 ), ' ', model_conc( 13 )

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
    class(unit_test_rxn_HL_phase_transfer_t), intent(in) :: this
    !> MPI communicator
    integer, intent(in), optional :: comm

#ifdef PMC_USE_MPI
    integer :: l_comm

    if (present(comm)) then
      l_comm = comm
    else
      l_comm = MPI_COMM_WORLD
    endif

    if (this%scenario.eq.1) then
      pack_size = &
              pmc_mpi_pack_size_integer(this%idx_O3, l_comm) + &
                      pmc_mpi_pack_size_integer(this%idx_O3_aq, l_comm) + &
                      pmc_mpi_pack_size_integer(this%idx_H2O2, l_comm) + &
                      pmc_mpi_pack_size_integer(this%idx_H2O2_aq, l_comm) + &
                      pmc_mpi_pack_size_integer(this%idx_H2O_aq, l_comm) + &
                      pmc_mpi_pack_size_integer(this%i_sect_unused, l_comm) + &
                      pmc_mpi_pack_size_integer(this%i_sect_the_mode, l_comm) + &
                      this%radius_update%pack_size() + &
                      this%number_update%pack_size() + &
                      this%update_data_GSD%pack_size() + &
                      this%update_data_GSD%pack_size()

    else if (this%scenario.eq.2) then
      pack_size = &
              pmc_mpi_pack_size_integer(this%idx_O3, l_comm) + &
                      pmc_mpi_pack_size_integer(this%idx_O3_aq, l_comm) + &
                      pmc_mpi_pack_size_integer(this%idx_H2O2, l_comm) + &
                      pmc_mpi_pack_size_integer(this%idx_H2O2_aq, l_comm) + &
                      pmc_mpi_pack_size_integer(this%idx_HNO3, l_comm) + &
                      pmc_mpi_pack_size_integer(this%idx_HNO3_aq, l_comm) + &
                      pmc_mpi_pack_size_integer(this%idx_NH3, l_comm) + &
                      pmc_mpi_pack_size_integer(this%idx_NH3_aq, l_comm) + &
                      pmc_mpi_pack_size_integer(this%idx_Na_p, l_comm) + &
                      pmc_mpi_pack_size_integer(this%idx_Cl_m, l_comm) + &
                      pmc_mpi_pack_size_integer(this%idx_Ca_pp, l_comm) + &
                      pmc_mpi_pack_size_integer(this%idx_H2O, l_comm) + &
                      pmc_mpi_pack_size_integer(this%idx_H2O_aq, l_comm) + &
                      pmc_mpi_pack_size_integer(this%i_sect_unused, l_comm) + &
                      pmc_mpi_pack_size_integer(this%i_sect_the_mode, l_comm) + &
                      this%radius_update%pack_size() + &
                      this%number_update%pack_size() + &
                      this%update_data_GSD%pack_size() + &
                      this%update_data_GSD%pack_size()
    end if


#else
    pack_size = 0
#endif

  end function pack_size

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Pack the object onto a buffer, advancing position
  subroutine bin_pack(this, buffer, pos, comm)

    !> Unit test data
    class(unit_test_rxn_HL_phase_transfer_t), intent(in) :: this
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

    if (this%scenario.eq.1) then

      prev_position = pos
      call pmc_mpi_pack_integer(buffer, pos, this%idx_O3, l_comm)
      call pmc_mpi_pack_integer(buffer, pos, this%idx_O3_aq, l_comm)
      call pmc_mpi_pack_integer(buffer, pos, this%idx_H2O2, l_comm)
      call pmc_mpi_pack_integer(buffer, pos, this%idx_H2O2_aq, l_comm)
      call pmc_mpi_pack_integer(buffer, pos, this%idx_H2O_aq, l_comm)
      call pmc_mpi_pack_integer(buffer, pos, this%i_sect_unused, l_comm)
      call pmc_mpi_pack_integer(buffer, pos, this%i_sect_the_mode, l_comm)
      call this%radius_update%bin_pack(buffer, pos)
      call this%number_update%bin_pack(buffer, pos)
      call this%update_data_GMD%bin_pack(buffer, pos)
      call this%update_data_GSD%bin_pack(buffer, pos)
      call assert(897212942, &
              pos - prev_position <= this%pack_size(l_comm))

    else if (this%scenario.eq.2) then

      prev_position = pos
      call pmc_mpi_pack_integer(buffer, pos, this%idx_O3, l_comm)
      call pmc_mpi_pack_integer(buffer, pos, this%idx_O3_aq, l_comm)
      call pmc_mpi_pack_integer(buffer, pos, this%idx_H2O2, l_comm)
      call pmc_mpi_pack_integer(buffer, pos, this%idx_H2O2_aq, l_comm)
      call pmc_mpi_pack_integer(buffer, pos, this%idx_HNO3, l_comm)
      call pmc_mpi_pack_integer(buffer, pos, this%idx_HNO3_aq, l_comm)
      call pmc_mpi_pack_integer(buffer, pos, this%idx_NH3, l_comm)
      call pmc_mpi_pack_integer(buffer, pos, this%idx_NH3_aq, l_comm)
      call pmc_mpi_pack_integer(buffer, pos, this%idx_Na_p, l_comm)
      call pmc_mpi_pack_integer(buffer, pos, this%idx_Cl_m, l_comm)
      call pmc_mpi_pack_integer(buffer, pos, this%idx_Ca_pp, l_comm)
      call pmc_mpi_pack_integer(buffer, pos, this%idx_H2O, l_comm)
      call pmc_mpi_pack_integer(buffer, pos, this%idx_H2O_aq, l_comm)
      call pmc_mpi_pack_integer(buffer, pos, this%i_sect_unused, l_comm)
      call pmc_mpi_pack_integer(buffer, pos, this%i_sect_the_mode, l_comm)
      call this%radius_update%bin_pack(buffer, pos)
      call this%number_update%bin_pack(buffer, pos)
      call this%update_data_GMD%bin_pack(buffer, pos)
      call this%update_data_GSD%bin_pack(buffer, pos)
      call assert(897212942, &
              pos - prev_position <= this%pack_size(l_comm))

    end if



#endif

  end subroutine bin_pack

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Unpack an object from a buffer, advancing position
  subroutine bin_unpack(this, buffer, pos, comm)

    !> Unit test data
    class(unit_test_rxn_HL_phase_transfer_t), intent(inout) :: this
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

    if (this%scenario.eq.1) then

    prev_position = pos
    call pmc_mpi_unpack_integer(buffer, pos, this%idx_O3, l_comm)
    call pmc_mpi_unpack_integer(buffer, pos, this%idx_O3_aq, l_comm)
    call pmc_mpi_unpack_integer(buffer, pos, this%idx_H2O2, l_comm)
    call pmc_mpi_unpack_integer(buffer, pos, this%idx_H2O2_aq, l_comm)
    call pmc_mpi_unpack_integer(buffer, pos, this%idx_H2O_aq, l_comm)
    call pmc_mpi_unpack_integer(buffer, pos, this%i_sect_unused, l_comm)
    call pmc_mpi_unpack_integer(buffer, pos, this%i_sect_the_mode, l_comm)
    call this%radius_update%bin_unpack(buffer, pos)
    call this%number_update%bin_unpack(buffer, pos)
    call this%update_data_GMD%bin_unpack(buffer, pos)
    call this%update_data_GSD%bin_unpack(buffer, pos)
    call assert(466926084, &
            pos - prev_position <= this%pack_size(l_comm))

    else if (this%scenario.eq.2) then

      prev_position = pos
      call pmc_mpi_unpack_integer(buffer, pos, this%idx_O3, l_comm)
      call pmc_mpi_unpack_integer(buffer, pos, this%idx_O3_aq, l_comm)
      call pmc_mpi_unpack_integer(buffer, pos, this%idx_H2O2, l_comm)
      call pmc_mpi_unpack_integer(buffer, pos, this%idx_H2O2_aq, l_comm)
      call pmc_mpi_unpack_integer(buffer, pos, this%idx_HNO3, l_comm)
      call pmc_mpi_unpack_integer(buffer, pos, this%idx_HNO3_aq, l_comm)
      call pmc_mpi_unpack_integer(buffer, pos, this%idx_NH3, l_comm)
      call pmc_mpi_unpack_integer(buffer, pos, this%idx_NH3_aq, l_comm)
      call pmc_mpi_unpack_integer(buffer, pos, this%idx_Na_p, l_comm)
      call pmc_mpi_unpack_integer(buffer, pos, this%idx_Cl_m, l_comm)
      call pmc_mpi_unpack_integer(buffer, pos, this%idx_Ca_pp, l_comm)
      call pmc_mpi_unpack_integer(buffer, pos, this%idx_H2O, l_comm)
      call pmc_mpi_unpack_integer(buffer, pos, this%idx_H2O_aq, l_comm)
      call pmc_mpi_unpack_integer(buffer, pos, this%i_sect_unused, l_comm)
      call pmc_mpi_unpack_integer(buffer, pos, this%i_sect_the_mode, l_comm)
      call this%radius_update%bin_unpack(buffer, pos)
      call this%number_update%bin_unpack(buffer, pos)
      call this%update_data_GMD%bin_unpack(buffer, pos)
      call this%update_data_GSD%bin_unpack(buffer, pos)
      call assert(466926084, &
              pos - prev_position <= this%pack_size(l_comm))

    end if

    this%is_initialized = .true.
#endif

  end subroutine bin_unpack

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_unit_test_rxn_HL_phase_transfer
