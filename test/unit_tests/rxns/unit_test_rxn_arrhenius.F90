! Copyright (C) 2019 Christian Guzman and Matt Dawson
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_unit_test_rxn_arrhenius module
module pmc_unit_test_rxn_arrhenius

  use pmc_camp_core
  use pmc_camp_state
  use pmc_mpi
  use pmc_unit_test_data
  use pmc_util

  implicit none
  private

  public :: unit_test_rxn_arrhenius_t

  !> Number of available initial states
  integer(kind=i_kind), parameter :: NUM_INIT_STATES = 5
  !> Number of species to evaluate
  integer(kind=i_kind), parameter :: NUM_SPEC = 4
  !> Number of time steps
  integer(kind=i_kind), parameter :: NUM_TIME_STEPS_C = 100
  !> Time step size (s)
  real(kind=dp), parameter :: TIME_STEP_SIZE_C = 1.0


  !> The Arrhenius rxn unit test
  type, extends(unit_test_data_t) :: unit_test_rxn_arrhenius_t
    private
    !> Species ids
    integer(kind=i_kind) :: idx_A, idx_B, idx_C, idx_D
    !> Species A initial concentrations
    real(kind=dp), dimension(NUM_INIT_STATES) :: &
      conc_A = (/ 1.0d0, 2.0d0, 3.0d0, 4.0d0, 5.0d0 /)
    !> Species B initial concentrations
    real(kind=dp), dimension(NUM_INIT_STATES) :: &
      conc_B = (/ 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0 /)
    !> Species C initial concentrations
    real(kind=dp), dimension(NUM_INIT_STATES) :: &
      conc_C = (/ 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0 /)
    !> Species D initial concentrations
    real(kind=dp), dimension(NUM_INIT_STATES) :: &
      conc_D = (/ 1.2d0, 17.0d0, 13.0d0, 9.0d0, 4.0d0 /)
    !> Initial temperatures
    real(kind=dp), dimension(NUM_INIT_STATES) :: &
      temperature = (/ 272.5d0, 248.3d0, 301.3d0, 276.0d0, 245.1d0 /)
    !> Initial pressures
    real(kind=dp), dimension(NUM_INIT_STATES) :: &
      pressure = (/ 101253.3d0, 92834.4d0, 90328.4d0, 100495.4d0, 96943.4d0 /)
    !> Flag inidicating local data is set up
    logical :: is_initialized = .false.
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

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Initialize the object data
  subroutine initialize( this, camp_core )

    use pmc_chem_spec_data

    !> Unit test data
    class(unit_test_rxn_arrhenius_t), intent(inout) :: this
    !> CAMP core
    class(camp_core_t), intent(in) :: camp_core

    type(chem_spec_data_t), pointer :: chem_spec_data

    call assert( 654980185, .not. this%is_initialized )

    ! Get the chemical species data
    call assert( 672415022, camp_core%get_chem_spec_data( chem_spec_data ) )

    ! Get the species ids
    this%idx_A = chem_spec_data%gas_state_id( "A" )
    this%idx_B = chem_spec_data%gas_state_id( "B" )
    this%idx_C = chem_spec_data%gas_state_id( "C" )
    this%idx_D = chem_spec_data%gas_state_id( "D" )

    ! Make sure all the species were found
    call assert_msg( 204555021, this%idx_A.gt.0, "Missing species A" )
    call assert_msg( 883729398, this%idx_B.gt.0, "Missing species B" )
    call assert_msg( 996047743, this%idx_C.gt.0, "Missing species C" )
    call assert_msg( 208366089, this%idx_D.gt.0, "Missing species D" )

    ! Flag the test as being initialized
    this%is_initialized = .true.

  end subroutine initialize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Input file name
  function input_file_name(this)

    !> Input file name
    character(len=:), allocatable :: input_file_name
    !> Unit test data
    class(unit_test_rxn_arrhenius_t), intent(in) :: this

    input_file_name = "rxn_arrhenius_config.json"

  end function input_file_name

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Output file name
  function output_file_name(this)

    !> Output file name
    character(len=:), allocatable :: output_file_name
    !> Unit test data
    class(unit_test_rxn_arrhenius_t), intent(in) :: this

    output_file_name = "rxn_arrhenius_results.txt"

  end function output_file_name

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

    real(kind=dp) :: conv

    ! Make sure the test data is initialized
    call assert(345723133, this%is_initialized)

    ! Get a conversion factor needed to set the D concentration
    conv = const%avagadro / const%univ_gas_const * 10.0d0**(-12.0d0) * &
            this%pressure( unique_state_id ) / &
            this%temperature( unique_state_id )

    ! Set the environmental conditions
    call camp_state%env_states(1)%set_temperature_K( &
      this%temperature( unique_state_id ) )
    call camp_state%env_states(1)%set_pressure_Pa(   &
      this%pressure(    unique_state_id ) )

    ! Set the species concentrations
    camp_state%state_var( this%idx_A ) = this%conc_A( unique_state_id )
    camp_state%state_var( this%idx_B ) = this%conc_B( unique_state_id )
    camp_state%state_var( this%idx_C ) = this%conc_C( unique_state_id )
    camp_state%state_var( this%idx_D ) = this%conc_D( unique_state_id ) / conv

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

    passed = this%analyze_or_output(camp_core, camp_state, unique_state_id, &
                                    model_time_step)

  end function analyze_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Output the results for a given grid cell
  subroutine output_results(this, camp_core, camp_state, unique_state_id, &
      model_time_step, output_file_unit)

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
    class(unit_test_rxn_arrhenius_t), intent(inout) :: this
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
    real(kind=dp) :: k1, k2, conv, time
    integer(kind=i_kind) :: i_spec

    passed = .true.

    ! Initial environmental state
    temperature    = this%temperature( unique_state_id )
    pressure       = this%pressure(    unique_state_id )

    ! Conversion factor
    conv = const%avagadro / const%univ_gas_const * 10.0d0**(-12.0d0) * &
            pressure / temperature

    ! Initial chemical state
    init_conc( 1 ) = this%conc_A(      unique_state_id )
    init_conc( 2 ) = this%conc_B(      unique_state_id )
    init_conc( 3 ) = this%conc_C(      unique_state_id )
    init_conc( 4 ) = this%conc_D(      unique_state_id ) / conv

    ! Calculate the rate constants
    k1 = init_conc( 4 ) * conv + 1476.0d0 * exp( -5.5d-21 / &
            (const%boltzmann * temperature) ) * &
            (temperature/300.0d0)**(150.0d0) * &
            (1.0d0 + 0.15d0 * pressure) / 60.0d0
    k2 = 21.0d0 * exp( -4000.0d0/temperature ) * &
            (temperature/315.0d0)**(11.0d0) * &
            (1.0d0 + 0.05d0 * pressure)

    ! Get the current model time
    time = model_time_step * this%time_step_size( )

    ! Calculate the true current concentrations
    true_conc( 1 ) = init_conc( 1 ) * exp(-(k1)*time)
    true_conc( 2 ) = init_conc( 1 ) * (k1/(k2-k1)) * (exp(-k1*time) - exp(-k2*time))
    true_conc( 3 ) = init_conc( 1 ) * &
               (1.0 + (k1*exp(-k2*time) - k2*exp(-k1*time))/(k2-k1))
    true_conc( 4 ) = init_conc( 4 )

    ! Get the model results
    model_conc( 1 ) = camp_state%state_var( this%idx_A )
    model_conc( 2 ) = camp_state%state_var( this%idx_B )
    model_conc( 3 ) = camp_state%state_var( this%idx_C )
    model_conc( 4 ) = camp_state%state_var( this%idx_D )

    ! Analyze or output results
    if( present( output_file_unit ) ) then

      ! Output the results
      write( OUTPUT_FILE_UNIT, * ) time, &
          ' ', true_conc( 1 ), ' ', model_conc( 1 ), &
          ' ', true_conc( 2 ), ' ', model_conc( 2 ), &
          ' ', true_conc( 3 ), ' ', model_conc( 3 ), &
          ' ', true_conc( 4 ), ' ', model_conc( 4 )

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
    class(unit_test_rxn_arrhenius_t), intent(in) :: this
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
      pmc_mpi_pack_size_integer(this%idx_D, l_comm)
#else
    pack_size = 0
#endif

  end function pack_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Pack the object onto a buffer, advancing position
  subroutine bin_pack(this, buffer, pos, comm)

    !> Unit test data
    class(unit_test_rxn_arrhenius_t), intent(in) :: this
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
    call pmc_mpi_pack_integer(buffer, pos, this%idx_D, l_comm)
    call assert(897212942, &
                pos - prev_position <= this%pack_size(l_comm))
#endif

  end subroutine bin_pack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Unpack an object from a buffer, advancing position
  subroutine bin_unpack(this, buffer, pos, comm)

    !> Unit test data
    class(unit_test_rxn_arrhenius_t), intent(inout) :: this
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
    call pmc_mpi_unpack_integer(buffer, pos, this%idx_D, l_comm)
    call assert(466926084, &
                pos - prev_position <= this%pack_size(l_comm))

    this%is_initialized = .true.
#endif

  end subroutine bin_unpack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_unit_test_rxn_arrhenius
