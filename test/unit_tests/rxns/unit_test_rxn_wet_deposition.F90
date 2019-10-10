! Copyright (C) 2019 Christian Guzman and Matt Dawson
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_unit_test_rxn_wet_deposition module
module pmc_unit_test_rxn_wet_deposition
  !TODO: Check multi-cell failing (update_rate_data)
  use pmc_camp_core
  use pmc_camp_state
  use pmc_mpi
  use pmc_unit_test_data
  use pmc_util
  use pmc_rxn_wet_deposition

  implicit none
  private

  public :: unit_test_rxn_wet_deposition_t

  !> Number of available initial states
  integer(kind=i_kind), parameter :: NUM_INIT_STATES = 5
  !> Number of species to evaluate
  integer(kind=i_kind), parameter :: NUM_SPEC = 8
  !> Number of time steps
  integer(kind=i_kind), parameter :: NUM_TIME_STEPS_C = 100
  !> Time step size (s)
  real(kind=dp), parameter :: TIME_STEP_SIZE_C = 1.0


  !> The wet_deposition rxn unit test
  type, extends(unit_test_data_t) :: unit_test_rxn_wet_deposition_t
    private
    !> Species ids
    integer(kind=i_kind) ::  idx_1RA, idx_1RB, idx_1CB, idx_1CC
    integer(kind=i_kind) ::  idx_2RA, idx_2RB, idx_2CB, idx_2CC
    !> Species A initial concentrations
    real(kind=dp), dimension(NUM_INIT_STATES) :: &
            conc_1RA = (/ 1.0d0, 2.0d0, 3.0d0, 4.0d0, 5.0d0 /)
    !> Species B initial concentrations
    real(kind=dp), dimension(NUM_INIT_STATES) :: &
            conc_1RB = (/ 2.2d0, 3.3d0, 4.4d0, 5.5d0, 6.6d0 /)
    !> Species C initial concentrations
    real(kind=dp), dimension(NUM_INIT_STATES) :: &
            conc_1CB = (/ 1.7d0, 2.7d0, 3.7d0, 4.7d0, 5.7d0 /)
    !> Species C initial concentrations
    real(kind=dp), dimension(NUM_INIT_STATES) :: &
            conc_1CC = (/ 2.4d0, 3.4d0, 4.4d0, 5.4d0, 6.4d0 /)
    !> Species C initial concentrations
    real(kind=dp), dimension(NUM_INIT_STATES) :: &
            conc_2RA = (/ 2.5d0, 3.5d0, 4.5d0, 5.5d0, 6.5d0 /)
    !> Species C initial concentrations
    real(kind=dp), dimension(NUM_INIT_STATES) :: &
            conc_2RB = (/ 0.7d0, 1.7d0, 2.7d0, 3.7d0, 4.7d0 /)
    !> Species C initial concentrations
    real(kind=dp), dimension(NUM_INIT_STATES) :: &
            conc_2CB = (/ 1.6d0, 2.6d0, 3.6d0, 4.6d0, 5.6d0 /)
    !> Species C initial concentrations
    real(kind=dp), dimension(NUM_INIT_STATES) :: &
            conc_2CC = (/ 1.9d0, 2.9d0, 3.9d0, 4.9d0, 5.9d0 /)
    !> Initial temperatures
    real(kind=dp), dimension(NUM_INIT_STATES) :: &
            temperature = (/ 272.5d0, 248.3d0, 301.3d0, 276.0d0, 245.1d0 /)
    !> Initial pressures
    real(kind=dp), dimension(NUM_INIT_STATES) :: &
            pressure = (/ 101253.3d0, 92834.4d0, 90328.4d0, 100495.4d0, 96943.4d0 /)
    !> Flag inidicating local data is set up
    logical :: is_initialized = .false.
    !> For setting rates
    type(rxn_update_data_wet_deposition_t) :: rate_update_rain, rate_update_cloud
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
  end type unit_test_rxn_wet_deposition_t

  !> Constructor for the wet_deposition test
  interface unit_test_rxn_wet_deposition_t
    procedure :: constructor
  end interface unit_test_rxn_wet_deposition_t

contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Constructor for the wet_deposition unit test
  function constructor() result(new_obj)

    !> A new unit test object
    type(unit_test_rxn_wet_deposition_t), pointer :: new_obj

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
    class(unit_test_rxn_wet_deposition_t), intent(inout) :: this
    !> CAMP core
    class(camp_core_t), intent(in) :: camp_core
    !> rxn data
    class(rxn_data_t), pointer :: rxn

    character(len=:), allocatable :: key, str_val, rep_name
    integer(kind=i_kind) ::  i_mech_rxn_rain, i_mech_rxn_cloud, i_rxn
    real(kind=dp) :: rate_rain, rate_cloud
    class(aero_rep_data_t), pointer :: aero_rep
    type(string_t), allocatable :: unique_names(:)

    !> For setting rates
    type(mechanism_data_t), pointer :: mechanism
    type(rxn_factory_t) :: rxn_factory

    call assert( 654980185, .not. this%is_initialized )

    ! Find the mechanism
    key = "wet deposition"
    call assert(260845179, camp_core%get_mechanism(key, mechanism))

    ! Find the A wet_deposition reaction
    key = "rxn id"
    i_mech_rxn_rain = 0
    i_mech_rxn_cloud = 0
    do i_rxn = 1, mechanism%size()
      rxn => mechanism%get_rxn(i_rxn)
      if (rxn%property_set%get_string(key, str_val)) then
        if (trim(str_val).eq."rxn rain") then
          i_mech_rxn_rain = i_rxn
          select type (rxn_loss => rxn)
          class is (rxn_wet_deposition_t)
            call rxn_factory%initialize_update_data(rxn_loss, &
                    this%rate_update_rain)
          end select
        end if
        if (trim(str_val).eq."rxn cloud") then
          i_mech_rxn_cloud = i_rxn
          select type (rxn_loss => rxn)
          class is (rxn_wet_deposition_t)
            call rxn_factory%initialize_update_data(rxn_loss, &
                    this%rate_update_cloud)
          end select
        end if
      end if
    end do
    call assert(262750713, i_mech_rxn_rain.eq.1)
    call assert(508177387, i_mech_rxn_cloud.eq.2)

    ! Find species in the first aerosol representation
    rep_name = "my first particle"
    call assert( 672415022, camp_core%get_aero_rep(rep_name, aero_rep ) )

    ! Get the species ids
    unique_names = aero_rep%unique_names(phase_name = "rain", spec_name = "A" )
    call assert(778960523, size(unique_names).eq.1)
    this%idx_1RA = aero_rep%spec_state_id(unique_names(1)%string);

    unique_names = aero_rep%unique_names(phase_name = "rain", spec_name = "B" )
    call assert(778960523, size(unique_names).eq.1)
    this%idx_1RB = aero_rep%spec_state_id(unique_names(1)%string);

    unique_names = aero_rep%unique_names(phase_name = "cloud", spec_name = "B" )
    call assert(778960523, size(unique_names).eq.1)
    this%idx_1CB = aero_rep%spec_state_id(unique_names(1)%string);

    unique_names = aero_rep%unique_names(phase_name = "cloud", spec_name = "C" )
    call assert(778960523, size(unique_names).eq.1)
    this%idx_1CC = aero_rep%spec_state_id(unique_names(1)%string);

    ! Find species in the second aerosol representation
    rep_name = "my second particle"
    call assert( 672415022, camp_core%get_aero_rep(rep_name, aero_rep ) )

    ! Get the species ids
    unique_names = aero_rep%unique_names(phase_name = "rain", spec_name = "A" )
    call assert(778960523, size(unique_names).eq.1)
    this%idx_2RA = aero_rep%spec_state_id(unique_names(1)%string);

    unique_names = aero_rep%unique_names(phase_name = "rain", spec_name = "B" )
    call assert(778960523, size(unique_names).eq.1)
    this%idx_2RB = aero_rep%spec_state_id(unique_names(1)%string);

    unique_names = aero_rep%unique_names(phase_name = "cloud", spec_name = "B" )
    call assert(778960523, size(unique_names).eq.1)
    this%idx_2CB = aero_rep%spec_state_id(unique_names(1)%string);

    unique_names = aero_rep%unique_names(phase_name = "cloud", spec_name = "C" )
    call assert(778960523, size(unique_names).eq.1)
    this%idx_2CC = aero_rep%spec_state_id(unique_names(1)%string);

    ! Make sure all the species were found
    call assert_msg( 204555021, this%idx_1RA.gt.0, "Missing species 1RA" )
    call assert_msg( 883729398, this%idx_1RB.gt.0, "Missing species 1RB" )
    call assert_msg( 883729398, this%idx_1CB.gt.0, "Missing species 1CB" )
    call assert_msg( 883729398, this%idx_1CC.gt.0, "Missing species 1CC" )
    call assert_msg( 883729398, this%idx_2RA.gt.0, "Missing species 2RA" )
    call assert_msg( 883729398, this%idx_2RB.gt.0, "Missing species 2RB" )
    call assert_msg( 883729398, this%idx_2CB.gt.0, "Missing species 2CB" )
    call assert_msg( 883729398, this%idx_2CC.gt.0, "Missing species 2CC" )

    ! Flag the test as being initialized
    this%is_initialized = .true.

  end subroutine initialize

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Input file name
  function input_file_name(this)

    !> Input file name
    character(len=:), allocatable :: input_file_name
    !> Unit test data
    class(unit_test_rxn_wet_deposition_t), intent(in) :: this

    input_file_name = "rxn_wet_deposition_config.json"

  end function input_file_name

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Output file name
  function output_file_name(this)

    !> Output file name
    character(len=:), allocatable :: output_file_name
    !> Unit test data
    class(unit_test_rxn_wet_deposition_t), intent(in) :: this

    output_file_name = "rxn_wet_deposition_results.txt"

  end function output_file_name

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Number of unique initial states available from the unit test
  function num_unique_states(this)

    !> Number of unique states
    integer(kind=i_kind) :: num_unique_states
    !> Unit test data
    class(unit_test_rxn_wet_deposition_t), intent(in) :: this

    num_unique_states = NUM_INIT_STATES

  end function num_unique_states

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Initialize a camp_state_t object based on a given index
  subroutine initialize_state(this, grid_cell_id, camp_core, camp_state, &
          unique_state_id)

    !> Unit test data
    class(unit_test_rxn_wet_deposition_t), intent(inout) :: this
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

    !Set the rxn rates

    rate_rain = 0.954d0
    rate_cloud = 1.0d-2

    call this%rate_update_rain%set_rate(rate_rain)
    call this%rate_update_cloud%set_rate(43912.5d0)
    call camp_core%update_rxn_data(this%rate_update_rain)
    call camp_core%update_rxn_data(this%rate_update_cloud)

    ! Test re-setting of the rxn B rate
    call this%rate_update_cloud%set_rate(rate_cloud)
    call camp_core%update_rxn_data(this%rate_update_cloud)

    ! Set the environmental conditions
    camp_state%env_state%temp     = this%temperature( unique_state_id )
    camp_state%env_state%pressure = this%pressure(    unique_state_id )
    call camp_state%update_env_state( )

    ! Set the species concentrations
    camp_state%state_var( this%idx_1RA ) = this%conc_1RA( unique_state_id )
    camp_state%state_var( this%idx_1RB ) = this%conc_1RB( unique_state_id )
    camp_state%state_var( this%idx_1CB ) = this%conc_1CB( unique_state_id )
    camp_state%state_var( this%idx_1CC ) = this%conc_1CC( unique_state_id )
    camp_state%state_var( this%idx_2RA ) = this%conc_2RA( unique_state_id )
    camp_state%state_var( this%idx_2RB ) = this%conc_2RB( unique_state_id )
    camp_state%state_var( this%idx_2CB ) = this%conc_2CB( unique_state_id )
    camp_state%state_var( this%idx_2CC ) = this%conc_2CC( unique_state_id )

  end subroutine initialize_state

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Number of time steps for this test
  function num_time_steps(this)

    !> Number of time steps
    integer(kind=i_kind) :: num_time_steps
    !> Unit test data
    class(unit_test_rxn_wet_deposition_t), intent(in) :: this

    num_time_steps = NUM_TIME_STEPS_C

  end function num_time_steps

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Time step size for this test (s)
  function time_step_size(this)

    !> Time step size (s)
    real(kind=dp) :: time_step_size
    !> Unit test data
    class(unit_test_rxn_wet_deposition_t), intent(in) :: this

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
    class(unit_test_rxn_wet_deposition_t), intent(inout) :: this
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
    class(unit_test_rxn_wet_deposition_t), intent(inout) :: this
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
    class(unit_test_rxn_wet_deposition_t), intent(inout) :: this
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
    real(kind=dp) :: k_rain, k_cloud, conv, rate_rain, rate_cloud, time
    integer(kind=i_kind) :: i_spec

    passed = .true.

    ! Initial environmental state
    temperature    = this%temperature( unique_state_id )
    pressure       = this%pressure(    unique_state_id )

    ! Initial chemical state
    init_conc( 1 ) = this%conc_1RA(      unique_state_id )
    init_conc( 2 ) = this%conc_1RB(      unique_state_id )
    init_conc( 3 ) = this%conc_1CB(      unique_state_id )
    init_conc( 4 ) = this%conc_1CC(      unique_state_id )
    init_conc( 5 ) = this%conc_2RA(      unique_state_id )
    init_conc( 6 ) = this%conc_2RB(      unique_state_id )
    init_conc( 7 ) = this%conc_2CB(      unique_state_id )
    init_conc( 8 ) = this%conc_2CC(      unique_state_id )

    ! Calculate the rate constants

    rate_rain = 0.954d0
    rate_cloud = 1.0d-2
    k_rain = rate_rain
    k_cloud = rate_cloud * 12.3d0

    ! Get the current model time
    time = model_time_step * this%time_step_size( )

    ! Calculate the true current concentrations
    true_conc( 1 ) = init_conc( 1 ) * exp(-(k_rain)*time)
    true_conc( 2 ) = init_conc( 2 ) * exp(-(k_rain)*time)
    true_conc( 3 ) = init_conc( 3 ) * exp(-(k_cloud)*time)
    true_conc( 4 ) = init_conc( 4 ) * exp(-(k_cloud)*time)
    true_conc( 5 ) = init_conc( 5 ) * exp(-(k_rain)*time)
    true_conc( 6 ) = init_conc( 6 ) * exp(-(k_rain)*time)
    true_conc( 7 ) = init_conc( 7 ) * exp(-(k_cloud)*time)
    true_conc( 8 ) = init_conc( 8 ) * exp(-(k_cloud)*time)

    ! Get the model results
    model_conc( 1 ) = camp_state%state_var( this%idx_1RA )
    model_conc( 2 ) = camp_state%state_var( this%idx_1RB )
    model_conc( 3 ) = camp_state%state_var( this%idx_1CB )
    model_conc( 4 ) = camp_state%state_var( this%idx_1CC )
    model_conc( 5 ) = camp_state%state_var( this%idx_2RA )
    model_conc( 6 ) = camp_state%state_var( this%idx_2RB )
    model_conc( 7 ) = camp_state%state_var( this%idx_2CB )
    model_conc( 8 ) = camp_state%state_var( this%idx_2CC )

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
              ' ', true_conc( 8 ), ' ', model_conc( 8 )

    else

      ! Check the model results
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

  end function analyze_or_output

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Determine the number of bytes required to pack the object onto a buffer
  integer(kind=i_kind) function pack_size(this, comm)

    !> Unit test data
    class(unit_test_rxn_wet_deposition_t), intent(in) :: this
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
            pmc_mpi_pack_size_integer(this%idx_1RA, l_comm) + &
                    pmc_mpi_pack_size_integer(this%idx_1RB, l_comm) + &
                    pmc_mpi_pack_size_integer(this%idx_1CB, l_comm) + &
                    pmc_mpi_pack_size_integer(this%idx_1CC, l_comm) + &
                    pmc_mpi_pack_size_integer(this%idx_2RA, l_comm) + &
                    pmc_mpi_pack_size_integer(this%idx_2RB, l_comm) + &
                    pmc_mpi_pack_size_integer(this%idx_2CB, l_comm) + &
                    pmc_mpi_pack_size_integer(this%idx_2CC, l_comm) + &
                    this%rate_update_rain%pack_size() + &
                    this%rate_update_cloud%pack_size()
#else
    pack_size = 0
#endif

  end function pack_size

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Pack the object onto a buffer, advancing position
  subroutine bin_pack(this, buffer, pos, comm)

    !> Unit test data
    class(unit_test_rxn_wet_deposition_t), intent(in) :: this
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
    call pmc_mpi_pack_integer(buffer, pos, this%idx_1RA, l_comm)
    call pmc_mpi_pack_integer(buffer, pos, this%idx_1RB, l_comm)
    call pmc_mpi_pack_integer(buffer, pos, this%idx_1CB, l_comm)
    call pmc_mpi_pack_integer(buffer, pos, this%idx_1CC, l_comm)
    call pmc_mpi_pack_integer(buffer, pos, this%idx_2RA, l_comm)
    call pmc_mpi_pack_integer(buffer, pos, this%idx_2RB, l_comm)
    call pmc_mpi_pack_integer(buffer, pos, this%idx_2CB, l_comm)
    call pmc_mpi_pack_integer(buffer, pos, this%idx_2CC, l_comm)
    call this%rate_update_rain%bin_pack(buffer, pos)
    call this%rate_update_cloud%bin_pack(buffer, pos)
    call assert(897212942, &
            pos - prev_position <= this%pack_size(l_comm))
#endif

  end subroutine bin_pack

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Unpack an object from a buffer, advancing position
  subroutine bin_unpack(this, buffer, pos, comm)

    !> Unit test data
    class(unit_test_rxn_wet_deposition_t), intent(inout) :: this
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
    call pmc_mpi_unpack_integer(buffer, pos, this%idx_1RA, l_comm)
    call pmc_mpi_unpack_integer(buffer, pos, this%idx_1RB, l_comm)
    call pmc_mpi_unpack_integer(buffer, pos, this%idx_1CB, l_comm)
    call pmc_mpi_unpack_integer(buffer, pos, this%idx_1CC, l_comm)
    call pmc_mpi_unpack_integer(buffer, pos, this%idx_2RA, l_comm)
    call pmc_mpi_unpack_integer(buffer, pos, this%idx_2RB, l_comm)
    call pmc_mpi_unpack_integer(buffer, pos, this%idx_2CB, l_comm)
    call pmc_mpi_unpack_integer(buffer, pos, this%idx_2CC, l_comm)
    call this%rate_update_rain%bin_unpack(buffer, pos)
    call this%rate_update_cloud%bin_unpack(buffer, pos)
    call assert(466926084, &
            pos - prev_position <= this%pack_size(l_comm))

    this%is_initialized = .true.
#endif

  end subroutine bin_unpack

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_unit_test_rxn_wet_deposition
