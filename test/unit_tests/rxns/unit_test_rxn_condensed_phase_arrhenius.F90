! Copyright (C) 2019 Christian Guzman and Matt Dawson
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_unit_test_rxn_condensed_phase_arrhenius module
module pmc_unit_test_rxn_condensed_phase_arrhenius
  !todo: fix multi-cell
  !todo: condensed_2 is so slow idk why, num_time_steps are reduced for this reason
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

  public :: unit_test_rxn_condensed_phase_arrhenius_t

  !> Number of available initial states
  integer(kind=i_kind), parameter :: NUM_INIT_STATES = 5
  !> Number of species to evaluate
  integer(kind=i_kind), parameter :: NUM_SPEC = 9
  !> Number of time steps
  integer(kind=i_kind), parameter :: NUM_TIME_STEPS_C = 10
  !> Time step size (s)
  real(kind=dp), parameter :: TIME_STEP_SIZE_C = 1.0

  !> The condensed_phase_arrhenius rxn unit test
  type, extends(unit_test_data_t) :: unit_test_rxn_condensed_phase_arrhenius_t
    private
    !> Scenario
    integer(kind=i_kind) :: scenario = UNIT_TEST_SCENARIO_
    !> Species ids
    integer(kind=i_kind) ::  idx_A_aq, idx_B_aq, idx_C_aq, idx_D_aq, idx_H2O, &
            idx_A_org, idx_B_org, idx_C_org, idx_D_org, i_sect_unused, i_sect_the_mode
    !> Species initial concentrations
    real(kind=dp), dimension(NUM_INIT_STATES) :: &
            conc_A_aq = (/ 1.0d0, 2.0d0, 3.0d0, 4.0d0, 5.0d0 /)
    !> Species initial concentrations
    real(kind=dp), dimension(NUM_INIT_STATES) :: &
            conc_B_aq = (/ 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0 /)
    !> Species initial concentrations
    real(kind=dp), dimension(NUM_INIT_STATES) :: &
            conc_C_aq = (/ 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0 /)
    !> Species initial concentrations
    real(kind=dp), dimension(NUM_INIT_STATES) :: &
            conc_D_aq = (/ 1.2d-2, 1.2d-2, 1.2d-2, 1.2d-2, 1.2d-2 /)
    !> Species initial concentrations
    real(kind=dp), dimension(NUM_INIT_STATES) :: &
            conc_H2O = (/ 2.3d0, 2.3d0, 2.3d0, 2.3d0, 2.3d0 /)
    !> Species initial concentrations
    real(kind=dp), dimension(NUM_INIT_STATES) :: &
            conc_A_org = (/ 1.0d0, 2.0d0, 3.0d0, 4.0d0, 5.0d0 /)
    !> Species initial concentrations
    real(kind=dp), dimension(NUM_INIT_STATES) :: &
            conc_B_org = (/ 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0 /)
    !> Species initial concentrations
    real(kind=dp), dimension(NUM_INIT_STATES) :: &
            conc_C_org = (/ 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0 /)
    !> Species initial concentrations
    real(kind=dp), dimension(NUM_INIT_STATES) :: &
            conc_D_org = (/ 1.2d-2, 1.2d-2, 1.2d-2, 1.2d-2, 1.2d-2 /)
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
  end type unit_test_rxn_condensed_phase_arrhenius_t

  !> Constructor for the condensed_phase_arrhenius test
  interface unit_test_rxn_condensed_phase_arrhenius_t
    procedure :: constructor
  end interface unit_test_rxn_condensed_phase_arrhenius_t

contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Constructor for the condensed_phase_arrhenius unit test
  function constructor() result(new_obj)

    !> A new unit test object
    type(unit_test_rxn_condensed_phase_arrhenius_t), pointer :: new_obj

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
    class(unit_test_rxn_condensed_phase_arrhenius_t), intent(inout) :: this
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

    if (this%scenario.eq.2) then

      ! Set the aerosol representation id
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
    else if (this%scenario.eq.2) then
      idx_prefix = "the mode."
    end if
    key = idx_prefix//"aqueous aerosol.A"
    this%idx_A_aq = aero_rep_ptr%spec_state_id(key);
    key = idx_prefix//"aqueous aerosol.B"
    this%idx_B_aq = aero_rep_ptr%spec_state_id(key);
    key = idx_prefix//"aqueous aerosol.C"
    this%idx_C_aq = aero_rep_ptr%spec_state_id(key);
    key = idx_prefix//"aqueous aerosol.D"
    this%idx_D_aq = aero_rep_ptr%spec_state_id(key);
    key = idx_prefix//"aqueous aerosol.H2O_aq"
    this%idx_H2O = aero_rep_ptr%spec_state_id(key);
    key = idx_prefix//"organic aerosol.A"
    this%idx_A_org = aero_rep_ptr%spec_state_id(key);
    key = idx_prefix//"organic aerosol.B"
    this%idx_B_org = aero_rep_ptr%spec_state_id(key);
    key = idx_prefix//"organic aerosol.C"
    this%idx_C_org = aero_rep_ptr%spec_state_id(key);
    key = idx_prefix//"organic aerosol.D"
    this%idx_D_org = aero_rep_ptr%spec_state_id(key);

    ! Make sure all the species were found
    call assert_msg( 204555021, this%idx_A_aq.gt.0, "Missing species A_aq" )
    call assert_msg( 883729398, this%idx_B_aq.gt.0, "Missing species B_aq" )
    call assert_msg( 883729398, this%idx_C_aq.gt.0, "Missing species C_aq" )
    call assert_msg( 883729398, this%idx_D_aq.gt.0, "Missing species D_aq" )
    call assert_msg( 883729398, this%idx_H2O.gt.0, "Missing species H2O" )
    call assert_msg( 883729398, this%idx_A_org.gt.0, "Missing species A_org" )
    call assert_msg( 883729398, this%idx_B_org.gt.0, "Missing species B_org" )
    call assert_msg( 883729398, this%idx_C_org.gt.0, "Missing species C_org" )
    call assert_msg( 883729398, this%idx_D_org.gt.0, "Missing species D_org" )

    ! Flag the test as being initialized
    this%is_initialized = .true.

  end subroutine initialize

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Input file name
  function input_file_name(this)

    !> Input file name
    character(len=:), allocatable :: input_file_name
    !> Unit test data
    class(unit_test_rxn_condensed_phase_arrhenius_t), intent(in) :: this

    ! Get the condensed_phase_arrhenius reaction mechanism json file
    if (this%scenario.eq.1) then
      input_file_name = 'rxn_condensed_phase_arrhenius_config.json'
    else if (this%scenario.eq.2) then
      input_file_name = 'rxn_condensed_phase_arrhenius_config_2.json'
    end if

  end function input_file_name

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Output file name
  function output_file_name(this)

    !> Output file name
    character(len=:), allocatable :: output_file_name
    !> Unit test data
    class(unit_test_rxn_condensed_phase_arrhenius_t), intent(in) :: this

    ! Get the condensed_phase_arrhenius reaction mechanism json file
    if (this%scenario.eq.1) then
      output_file_name = "rxn_condensed_phase_arrhenius_results.txt"
    else if (this%scenario.eq.2) then
      output_file_name = "rxn_condensed_phase_arrhenius_results_2.txt"
    end if

  end function output_file_name

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Number of unique initial states available from the unit test
  function num_unique_states(this)

    !> Number of unique states
    integer(kind=i_kind) :: num_unique_states
    !> Unit test data
    class(unit_test_rxn_condensed_phase_arrhenius_t), intent(in) :: this

    num_unique_states = NUM_INIT_STATES

  end function num_unique_states

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Initialize a camp_state_t object based on a given index
  subroutine initialize_state(this, grid_cell_id, camp_core, camp_state, &
          unique_state_id)

    !> Unit test data
    class(unit_test_rxn_condensed_phase_arrhenius_t), intent(inout) :: this
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
    camp_state%state_var( this%idx_A_aq ) = this%conc_A_aq( unique_state_id )
    camp_state%state_var( this%idx_B_aq ) = this%conc_B_aq( unique_state_id )
    camp_state%state_var( this%idx_C_aq ) = this%conc_C_aq( unique_state_id )
    camp_state%state_var( this%idx_D_aq ) = this%conc_D_aq( unique_state_id )
    camp_state%state_var( this%idx_H2O ) = this%conc_H2O( unique_state_id )
    camp_state%state_var( this%idx_A_org ) = this%conc_A_org( unique_state_id )
    camp_state%state_var( this%idx_B_org ) = this%conc_B_org( unique_state_id )
    camp_state%state_var( this%idx_C_org ) = this%conc_C_org( unique_state_id )
    camp_state%state_var( this%idx_D_org ) = this%conc_D_org( unique_state_id )

  end subroutine initialize_state

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Number of time steps for this test
  function num_time_steps(this)

    !> Number of time steps
    integer(kind=i_kind) :: num_time_steps
    !> Unit test data
    class(unit_test_rxn_condensed_phase_arrhenius_t), intent(in) :: this

    num_time_steps = NUM_TIME_STEPS_C

  end function num_time_steps

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Time step size for this test (s)
  function time_step_size(this)

    !> Time step size (s)
    real(kind=dp) :: time_step_size
    !> Unit test data
    class(unit_test_rxn_condensed_phase_arrhenius_t), intent(in) :: this

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
    class(unit_test_rxn_condensed_phase_arrhenius_t), intent(inout) :: this
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
    class(unit_test_rxn_condensed_phase_arrhenius_t), intent(inout) :: this
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
    class(unit_test_rxn_condensed_phase_arrhenius_t), intent(inout) :: this
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
    real(kind=dp) :: k1_aq, k2_aq, k1_org, k2_org, MW_A, MW_B, MW_C, MW_D, conc_d, conc_water, time
    integer(kind=i_kind) :: i_spec

    passed = .true.

    init_conc(:) = 0.0
    true_conc(:) = 0.0

    ! Initial environmental state
    temperature    = this%temperature( unique_state_id )
    pressure       = this%pressure(    unique_state_id )

    ! Initial chemical state
    init_conc( 1 ) = this%conc_A_aq(      unique_state_id )
    init_conc( 2 ) = this%conc_B_aq(      unique_state_id )
    init_conc( 3 ) = this%conc_C_aq(      unique_state_id )
    init_conc( 4 ) = this%conc_D_aq(      unique_state_id )
    init_conc( 5 ) = this%conc_H2O(      unique_state_id )
    init_conc( 6 ) = this%conc_A_org(      unique_state_id )
    init_conc( 7 ) = this%conc_B_org(      unique_state_id )
    init_conc( 8 ) = this%conc_C_org(      unique_state_id )
    init_conc( 9 ) = this%conc_D_org(      unique_state_id )

    ! Calculate the rate constants

    conc_D = this%conc_D_aq(      unique_state_id )
    conc_water = this%conc_H2O(      unique_state_id )
    MW_A = 0.1572
    MW_B = 0.0219
    MW_C = 0.2049
    MW_D = 0.0345
    k1_aq = conc_D / MW_D / conc_water + 1476.0d0 * &
            exp( -5.5d-21 / (const%boltzmann * temperature) ) * &
            (temperature/300.0d0)**(150.0d0) * (1.0d0 + 0.15d0 * pressure) / 60.0d0
    k2_aq = 21.0d0 * exp( -4000.0d0/temperature ) * (temperature/315d0)**(11.0d0) * &
            (1.0d0 + 0.05d0 * pressure)
    k1_org = conc_D / (MW_D*1.0e9) + 1476.0d0 * &
            exp( -5.5d-21 / (const%boltzmann * temperature) ) * &
            (temperature/300.0d0)**(150.0d0) * (1.0d0 + 0.15d0 * pressure) / 60.0d0
    k2_org = 21.0d0 * exp( -4000.0d0/temperature ) * (temperature/315d0)**(11.0d0) * &
            (1.0d0 + 0.05d0 * pressure)

    ! Get the current model time
    time = model_time_step * this%time_step_size( )

    ! Calculate the true current concentrations
    true_conc( 1 ) = init_conc( 1 ) * exp(-k1_aq*time)
    true_conc( 2 ) = init_conc( 1 ) * (k1_aq/(k2_aq-k1_aq)) * &
            (exp(-k1_aq*time) - exp(-k2_aq*time)) * MW_B / MW_A
    true_conc( 3 ) = init_conc( 1 ) * MW_C / MW_A * &
            (1.0d0 + (k1_aq*exp(-k2_aq*time) - &
            k2_aq*exp(-k1_aq*time))/(k2_aq-k1_aq))
    true_conc( 4 ) = init_conc( 4 )
    true_conc( 5 ) = init_conc( 5 )
    true_conc( 6 ) = init_conc( 6 ) * exp(-k1_org*time)
    true_conc( 7 ) = init_conc( 6 ) * &
            (k1_org/(k2_org-k1_org)) * &
            (exp(-k1_org*time) - exp(-k2_org*time)) * MW_B / MW_A
    true_conc( 8 ) = init_conc( 6 ) * MW_C / MW_A * &
            (1.0d0 + (k1_org*exp(-k2_org*time) - &
            k2_org*exp(-k1_org*time))/(k2_org-k1_org))
    true_conc( 9 ) = init_conc( 9 )

    ! Get the model results
    model_conc( 1 ) = camp_state%state_var( this%idx_A_aq )
    model_conc( 2 ) = camp_state%state_var( this%idx_B_aq )
    model_conc( 3 ) = camp_state%state_var( this%idx_C_aq )
    model_conc( 4 ) = camp_state%state_var( this%idx_D_aq )
    model_conc( 5 ) = camp_state%state_var( this%idx_H2O )
    model_conc( 6 ) = camp_state%state_var( this%idx_A_org )
    model_conc( 7 ) = camp_state%state_var( this%idx_B_org )
    model_conc( 8 ) = camp_state%state_var( this%idx_C_org )
    model_conc( 9 ) = camp_state%state_var( this%idx_D_org )

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
              ' ', true_conc( 9 ), ' ', model_conc( 9 )

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
    class(unit_test_rxn_condensed_phase_arrhenius_t), intent(in) :: this
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
            pmc_mpi_pack_size_integer(this%idx_A_aq, l_comm) + &
                    pmc_mpi_pack_size_integer(this%idx_B_aq, l_comm) + &
                    pmc_mpi_pack_size_integer(this%idx_C_aq, l_comm) + &
                    pmc_mpi_pack_size_integer(this%idx_D_aq, l_comm) + &
                    pmc_mpi_pack_size_integer(this%idx_H2O, l_comm) + &
                    pmc_mpi_pack_size_integer(this%idx_A_org, l_comm) + &
                    pmc_mpi_pack_size_integer(this%idx_B_org, l_comm) + &
                    pmc_mpi_pack_size_integer(this%idx_C_org, l_comm) + &
                    pmc_mpi_pack_size_integer(this%idx_D_org, l_comm) + &
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
    class(unit_test_rxn_condensed_phase_arrhenius_t), intent(in) :: this
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
    call pmc_mpi_pack_integer(buffer, pos, this%idx_A_aq, l_comm)
    call pmc_mpi_pack_integer(buffer, pos, this%idx_B_aq, l_comm)
    call pmc_mpi_pack_integer(buffer, pos, this%idx_C_aq, l_comm)
    call pmc_mpi_pack_integer(buffer, pos, this%idx_D_aq, l_comm)
    call pmc_mpi_pack_integer(buffer, pos, this%idx_H2O, l_comm)
    call pmc_mpi_pack_integer(buffer, pos, this%idx_B_org, l_comm)
    call pmc_mpi_pack_integer(buffer, pos, this%idx_B_org, l_comm)
    call pmc_mpi_pack_integer(buffer, pos, this%idx_C_org, l_comm)
    call pmc_mpi_pack_integer(buffer, pos, this%idx_D_org, l_comm)
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
    class(unit_test_rxn_condensed_phase_arrhenius_t), intent(inout) :: this
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
    call pmc_mpi_unpack_integer(buffer, pos, this%idx_A_aq, l_comm)
    call pmc_mpi_unpack_integer(buffer, pos, this%idx_B_aq, l_comm)
    call pmc_mpi_unpack_integer(buffer, pos, this%idx_C_aq, l_comm)
    call pmc_mpi_unpack_integer(buffer, pos, this%idx_D_aq, l_comm)
    call pmc_mpi_unpack_integer(buffer, pos, this%idx_H2O, l_comm)
    call pmc_mpi_unpack_integer(buffer, pos, this%idx_A_org, l_comm)
    call pmc_mpi_unpack_integer(buffer, pos, this%idx_B_org, l_comm)
    call pmc_mpi_unpack_integer(buffer, pos, this%idx_C_org, l_comm)
    call pmc_mpi_unpack_integer(buffer, pos, this%idx_D_org, l_comm)
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

end module pmc_unit_test_rxn_condensed_phase_arrhenius
