! Copyright (C) 2017 Matt Dawson
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_camp_box_model_data_t type and related functions

!> A simple box model for \ref camp_chem "CAMP" mechanisms
!! \todo{ Modify this to run full PartMC scenarios with CAMP enabled
!!        so that aerosol physical processes can be included. }
module pmc_camp_box_model_data

  use pmc_camp_core
  use pmc_camp_state
  use pmc_rxn_data
  use pmc_util

  implicit none
  private

  public :: camp_box_model_data_t

  ! New-line character
  character(len=*), parameter :: new_line = char(10)

  !> Property time profile
  type :: profile_t
    !> Name of the profile
    character(len=:), allocatable :: name
    !> Current time step in the profile
    integer(kind=i_kind) :: curr_time_step
    !> Index for the current position in the transition arrays
    integer(kind=i_kind) :: curr_transition
    !> Time step for the next transition, starting with 0 (initial conditions)
    integer(kind=i_kind), allocatable :: transition_time_step(:)
    !> New values for each transition, starting with the initial value
    real(kind=dp), allocatable :: transition_value(:)
  contains
    !> Reset the profile
    procedure :: reset
    !> Advance the profile by one time step
    procedure :: advance
    !> Get the current profile value
    procedure :: current_value
  end type profile_t

#ifdef PMC_USE_JSON
  !> Constructor for profile_t
  interface profile_t
    procedure :: profile_constructor
  end interface profile_t
#endif

  !> Reaction rate time profile
  type, extends(profile_t) :: rxn_profile_t
    !> Reaction update object
    class(rxn_update_data_t), pointer :: update_data
  contains
    !> Update the reaction rate
    procedure :: update_rxn
  end type rxn_profile_t

#ifdef PMC_USE_JSON
  !> Constructor for rxn_profile_t
  interface rxn_profile_t
    procedure :: rxn_profile_constructor
  end interface rxn_profile_t
#endif

  !> CAMP Box model
  type :: camp_box_model_data_t
  private
    !> CAMP core
    type(camp_core_t), pointer :: camp_core
    !> CAMP state
    type(camp_state_t), pointer :: camp_state
    !> Initial CAMP state
    type(camp_state_t), pointer :: initial_camp_state
    !> Number of time steps
    integer(kind=i_kind) :: num_steps
    !> Time step length [s]
    real(kind=dp) :: time_step__s
    !> Total integration time [s]
    real(kind=dp) :: total_time__s
    !> Temperature [K] profile
    type(profile_t) :: temperature__K
    !> Pressure [Pa] profile
    type(profile_t) :: pressure__Pa
    !> Reaction profiles
    type(rxn_profile_t), allocatable :: rxn_profiles(:)
  contains
    !> Run the box model
    procedure :: run
  end type camp_box_model_data_t

  !> Constructor for camp_box_model_data_t
  interface camp_box_model_data_t
    procedure :: constructor
  end interface camp_box_model_data_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Constructor for the CAMP box model
  function constructor( config_file ) result( new_obj )

#ifdef PMC_USE_JSON
    use json_module
#endif

    !> Pointer to a new box model object
    type(camp_box_model_data_t), pointer :: new_obj
    !> Box model configuration file
    character(len=*), intent(in) :: config_file
#ifdef PMC_USE_JSON
    type(json_core), target :: json
    type(json_file) :: j_file
    type(json_value), pointer :: j_obj, j_next

    logical(kind=json_lk) :: found, valid
    character(kind=json_ck, len=:), allocatable :: unicode_str_val
    character(kind=json_ck, len=:), allocatable :: json_err_msg
    character(kind=json_ck, len=:), allocatable :: spec_name
    integer(kind=json_ik) :: int_value
    real(kind=json_rk) :: real_value

    logical :: file_exists
    integer(kind=i_kind) :: num_rates, i_rate, spec_id

    ! Create the new box model
    allocate( new_obj )

    ! Create the CAMP core
    new_obj%camp_core => camp_core_t( config_file )
    call new_obj%camp_core%initialize( )
    call new_obj%camp_core%solver_initialize( )

    ! Create a CAMP state for the box model runs and initial conditions
    new_obj%camp_state         => new_obj%camp_core%new_state( )
    new_obj%initial_camp_state => new_obj%camp_core%new_state( )

    ! Load the configuration data from the json file
    call j_file%initialize( )
    call j_file%get_core( json )
    call assert_msg( 135902099, trim( config_file ).eq."", &
                     "Received empty string for file path" )
    inquire( file = trim( config_file ), exist = file_exists )
    call assert_msg( 181758805, file_exists, "Cannot find file: "// &
                    trim( config_file ) )
    call j_file%load_file( filename = trim( config_file ) )

    ! Get the number of time steps
    call j_file%get( "output time step [s]", real_value, found )
    call assert_msg( 912341318, found, &
                     "Missing 'output time step [s]' "// &
                     "in box model data" )
    new_obj%time_step__s = real_value

    ! Get the integration time
    call j_file%get( "total integration time [s]", real_value, found )
    call assert_msg( 442175882, found, &
                     "Missing 'total integration time [s]' "// &
                     "in box model data" )
    new_obj%num_steps = ceiling( real_value / new_obj%time_step__s )
    new_obj%total_time__s  = real_value

    ! Get the temperature profile
    call j_file%get( "temperature [K]", j_obj, found )
    call assert_msg( 212013359, found, &
                     "Missing 'temperature [K]' "// &
                     "in box model data" )
    new_obj%temperature__K = profile_t( json, j_obj )

    ! Get the pressure profile
    call j_file%get( "pressure [Pa]", j_obj, found )
    call assert_msg( 115322763, found, &
                     "Missing 'pressure [Pa]' "// &
                     "in box model data" )
    new_obj%pressure__Pa = profile_t( json, j_obj )

    ! Get the reaction rate profiles
    call j_file%get( "rates", j_obj, found )
    if( found ) then
      call j_file%info( "rates", n_children = num_rates )
      allocate( new_obj%rxn_profiles( num_rates ) )
      j_next => null( )
      i_rate = 0
      do while( associated( j_obj ) )
        i_rate = i_rate + 1
        new_obj%rxn_profiles( i_rate ) = &
          rxn_profile_t( new_obj%camp_core, json, j_obj )
        call json%get_next( j_obj, j_next )
        j_obj => j_next
      end do
      call assert( 763815888, i_rate .eq. num_rates )
    end if

    ! Get the initial species concentrations
    call j_file%get( "initial state", j_obj, found )
    if( found ) then
      j_next => null( )
      do while( associated( j_obj ) )
        call json%info( j_obj, name = spec_name )
        call json%get( j_obj, real_value )
        call assert_msg( 301872431, &
                         new_obj%camp_core%spec_state_id( spec_name, &
                                                          spec_id ), &
                         "Cannot find species '"//trim( spec_name )// &
                         "' in the chemical mechanism" )
        new_obj%camp_state%state_var( spec_id ) = real_value
      end do
    end if

    ! Set the initial environtmental conditions
    call new_obj%initial_camp_state%env_states(1)%set_temperature_K( &
      new_obj%temperature__K%current_value( ) )
    call new_obj%initial_camp_state%env_states(1)%set_pressure_Pa( &
      new_obj%pressure__Pa%current_value( ) )

    ! free the json file
    call j_file%destroy( )

#else
    new_obj => null( )
    call die_msg( 444795237, "JSON must be enabled for the CAMP box model" )
#endif

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Run the camp-chem box model
  subroutine run( this, output_file_unit )

    use pmc_solver_stats

    !> CAMP box model
    class(camp_box_model_data_t), intent(inout) :: this
    !> Output file unit
    integer(kind=i_kind), intent(in), optional :: output_file_unit

    type(solver_stats_t) :: solver_stats

    integer(kind=i_kind) :: f_unit
    integer(kind=i_kind) :: i_time, i_rate
    real(kind=dp) :: step_size__s, curr_time__s

    f_unit = 6
    if( present( output_file_unit ) ) f_unit = output_file_unit

    ! Initialize the model state
    this%camp_state%state_var( : ) = this%initial_camp_state%state_var( : )
    this%camp_state%env_var( : )   = this%initial_camp_state%env_var( : )

    ! Output file header and intitial conditions
    write(f_unit,*) 0.0, this%camp_state%env_var( : ), &
                    this%camp_state%state_var( : )

    ! Reset the profiles
    call this%temperature__K%reset( )
    call this%pressure__Pa%reset( )
    do i_rate = 1, size( this%rxn_profiles )
      call this%rxn_profiles( i_rate )%reset( )
    end do

    ! Solve the box model
    do i_time = 1, this%num_steps

      ! Set the current time and step size
      curr_time__s = i_time * this%time_step__s
      if( curr_time__s .gt. this%total_time__s ) then
        step_size__s = this%time_step__s - &
                       (curr_time__s - this%total_time__s )
      else
        step_size__s = this%time_step__s
      end if

      ! Update the model conditions
      call this%temperature__K%advance( )
      call this%camp_state%env_states(1)%set_temperature_K( &
        this%temperature__K%current_value( ) )
      call this%pressure__Pa%advance( )
      call this%camp_state%env_states(1)%set_pressure_Pa( &
        this%pressure__Pa%current_value( ) )
      do i_rate = 1, size( this%rxn_profiles )
        call this%rxn_profiles( i_rate )%advance( )
        call this%rxn_profiles( i_rate )%update_rxn( this%camp_core )
      end do

      ! Solve the chemistry
      call this%camp_core%solve( this%camp_state, step_size__s, &
                                 solver_stats = solver_stats )

      ! Output the model state
      write(f_unit,*) curr_time__s, this%camp_state%env_var( : ), &
                      this%camp_state%state_var( : )

    end do

  end subroutine run

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef PMC_USE_JSON
  !> Constructor for profile_t
  function profile_constructor( json, j_obj ) result( new_obj )

    use json_module

    !> New profile_t object
    type(profile_t) :: new_obj
    !> JSON core
    type(json_core), intent(inout) :: json
    !> JSON object to build profile from
    type(json_value), pointer, intent(inout)  :: j_obj

    type(json_value), pointer :: j_child, j_next
    integer(kind=json_ik) :: num_trans
    logical(kind=json_lk) :: found
    character(kind=json_ck, len=:), allocatable :: profile_name

    integer(kind=i_kind) :: i_trans

    ! Get the name of the profile
    call json%info( j_obj, name = profile_name )
    new_obj%name = profile_name

    ! Get the transitions
    call json%get( j_obj, "transitions", j_child, found )
    if( found ) then
      call json%info( j_obj, "transitions", n_children = num_trans )
      allocate( new_obj%transition_time_step( 0:num_trans ) )
      allocate( new_obj%transition_value(     0:num_trans ) )
      j_next => null( )
      call json%get( j_obj, "transitions(1)", j_child, found )
      i_trans = 0
      do while( associated( j_child ) )
        i_trans = i_trans + 1
        call json%get( j_child, "time step", &
                       new_obj%transition_time_step( i_trans ), found )
        call assert_msg( 325008338, found, &
                         "Missing 'time step' from element "// &
                         trim( to_string( i_trans ) )//" in transition '"// &
                         trim( new_obj%name )//"'" )
        call json%get( j_child, "value", &
                       new_obj%transition_value( i_trans ), found )
        call assert_msg( 483183389, found, &
                         "Missing 'value' from element "// &
                         trim( to_string( i_trans ) )//" in transition '"// &
                         trim( new_obj%name )//"'" )
        call json%get_next( j_child, j_next )
        j_child => j_next
      end do
      call assert( 249923619, i_trans .eq. num_trans )
    else
      allocate( new_obj%transition_time_step( 0:0 ) )
      allocate( new_obj%transition_value(     0:0 ) )
    end if

    ! Get the initial state
    call json%get( j_obj, "initial", j_child, found )
    call assert_msg( 382107479, found, &
                     "Missing 'initial' value in profile '"// &
                     trim( new_obj%name )//"'" )
    call json%get( j_child, new_obj%transition_value( 0 ) )
    new_obj%transition_time_step( 0 ) = 1

  end function profile_constructor
#endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Reset the profile to the initial state
  subroutine reset( this )

    !> Property profile
    class(profile_t), intent(inout) :: this

    this%curr_transition = 0
    this%curr_time_step  = 1

  end subroutine reset

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Advance the profile by one time step
  subroutine advance( this )

    !> Property profile
    class(profile_t), intent(inout) :: this

    this%curr_time_step = this%curr_time_step + 1

    if( this%curr_transition .eq. &
        size( this%transition_time_step ) - 1 ) return

    if( this%transition_time_step( this%curr_transition + 1 ) .eq. &
        this%curr_time_step ) &
      this%curr_transition = this%curr_transition + 1

  end subroutine advance

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get the current value of the profile
  function current_value( this )

    !> Current value of the profile
    real(kind=dp) :: current_value
    !> Property profile
    class(profile_t), intent(inout) :: this

    current_value = this%transition_value( this%curr_transition )

  end function current_value

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef PMC_USE_JSON
  !> Constructor for rxn_profile_t
  function rxn_profile_constructor( camp_core, json, j_obj ) result( new_obj )

    use pmc_mechanism_data
    use pmc_rxn_emission
    use pmc_rxn_factory
    use pmc_rxn_first_order_loss
    use pmc_rxn_photolysis
    use pmc_rxn_wet_deposition
    use json_module

    !> New reaction profile
    type(rxn_profile_t) :: new_obj
    !> CAMP core
    type(camp_core_t), intent(inout) :: camp_core
    !> JSON core
    type(json_core), intent(inout)  :: json
    !> JSON object to build the profile from
    type(json_value), pointer, intent(inout) :: j_obj

    type(mechanism_data_t), pointer :: mech
    class(rxn_data_t), pointer :: rxn
    type(rxn_factory_t) :: rxn_factory
    type(profile_t) :: base_profile

    integer(kind=i_kind) :: i_mech, i_rxn
    character(len=:), allocatable :: rxn_label
    logical :: found

    found = .false.

    ! Load the standard profile data
    base_profile                 = profile_t( json, j_obj )
    new_obj%name                 = base_profile%name
    new_obj%transition_time_step = base_profile%transition_time_step
    new_obj%transition_value     = base_profile%transition_value

    ! Find the reaction and initialize the update data object
    do i_mech = 1, size( camp_core%mechanism )
      mech => camp_core%mechanism( i_mech )%val
      do i_rxn = 1, mech%size( )
        rxn => mech%get_rxn( i_rxn )
        if( .not. rxn%property_set%get_string( "camp-box-model-id", &
                                               rxn_label ) ) continue
        if( .not. trim( rxn_label ) .eq. new_obj%name ) continue
        select type( rxn )
          class is( rxn_emission_t)
            allocate( rxn_update_data_emission_t::new_obj%update_data )
          class is( rxn_first_order_loss_t)
            allocate( rxn_update_data_first_order_loss_t::new_obj%update_data )
          class is( rxn_photolysis_t)
            allocate( rxn_update_data_photolysis_t::new_obj%update_data )
          class is( rxn_wet_deposition_t)
            allocate( rxn_update_data_wet_deposition_t::new_obj%update_data )
          class default
            call die_msg( 592455681, "Invalid reaction for rate updates" )
        end select
        call rxn_factory%initialize_update_data( rxn, new_obj%update_data )
        found = .true.
      end do
    end do

    call assert_msg( 800298506, found, "Could not find reaction label '"// &
                     trim( new_obj%name ) )

  end function rxn_profile_constructor
#endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Update a reaction with the current rate from the profile
  subroutine update_rxn( this, camp_core )

    use pmc_rxn_emission
    use pmc_rxn_first_order_loss
    use pmc_rxn_photolysis
    use pmc_rxn_wet_deposition

    !> Reaction rate profile
    class(rxn_profile_t) :: this
    !> CAMP core
    type(camp_core_t), intent(inout) :: camp_core

    class(rxn_update_data_t), pointer :: ud

    select type( ud => this%update_data )
      class is( rxn_update_data_emission_t )
        call ud%set_rate( this%current_value( ) )
      class is( rxn_update_data_first_order_loss_t )
        call ud%set_rate( this%current_value( ) )
      class is( rxn_update_data_photolysis_t )
        call ud%set_rate( this%current_value( ) )
      class is( rxn_update_data_wet_deposition_t )
        call ud%set_rate( this%current_value( ) )
      class default
        call die( 497670619 )
    end select

    call camp_core%update_rxn_data( this%update_data )

  end subroutine update_rxn

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_camp_box_model_data
