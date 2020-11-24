! Copyright (C) 2017 Matt Dawson
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_camp_box_model_data_t type and related functions

!> A simple box model for \ref camp_chem "CAMP" mechanisms
!! \todo{ Modify this to run full PartMC scenarios with CAMP enabled
!!        so that aerosol physical processes can be included. }
module pmc_camp_box_model_data

  use pmc_aero_rep_data
  use pmc_camp_core
  use pmc_camp_state
  use pmc_rxn_data
  use pmc_util

  implicit none
  private

  public :: camp_box_model_data_t

  ! New-line character
  character(len=*), parameter :: new_line = char(10)
  ! File unit to use for scipts output
  integer(kind=i_kind), parameter :: SCRIPTS_FILE_UNIT = 8

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
    !> Print the profile configuration
    procedure :: print => profile_do_print
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
    !> Finalize the reaction profile
    final :: rxn_profile_finalize
  end type rxn_profile_t

#ifdef PMC_USE_JSON
  !> Constructor for rxn_profile_t
  interface rxn_profile_t
    procedure :: rxn_profile_constructor
  end interface rxn_profile_t
#endif

  !> Aerosol representation time profile
  type, extends(profile_t) :: aero_rep_profile_t
    !> Index for update type
    integer(kind=i_kind) :: update_type
    !> Section id (modal/binned)
    integer(kind=i_kind) :: section_id = 0
    !> Aerosol representation update object
    class(aero_rep_update_data_t), pointer :: update_data
  contains
    !> Update the aerosol representation property of interest
    procedure :: update_aero_rep
    !> Finalize the profile
    final :: aero_rep_profile_finalize
  end type aero_rep_profile_t

#ifdef PMC_USE_JSON
  !> Constructor for aero_rep_profile_t
  interface aero_rep_profile_t
    procedure :: aero_rep_profile_constructor
  end interface aero_rep_profile_t
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
    !> State species names
    type(string_t), allocatable :: spec_names(:)
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
    !> Aerosol representation profiles
    type(aero_rep_profile_t), allocatable :: aero_rep_profiles(:)
  contains
    !> Run the box model
    procedure :: run
    !> Create a plotting script for the results
    procedure :: create_gnuplot_config_file
    !> Print the box model configuration
    procedure :: print => do_print
    !> Finalize the box model
    final :: finalize
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
    type(json_value), pointer :: j_obj, j_next, j_box_config

    logical(kind=json_lk) :: found, valid
    character(kind=json_ck, len=:), allocatable :: unicode_str_val
    character(kind=json_ck, len=:), allocatable :: json_err_msg
    character(kind=json_ck, len=:), allocatable :: spec_name
    integer(kind=json_ik) :: int_value
    real(kind=json_rk) :: real_value

    logical :: file_exists
    integer(kind=i_kind) :: num_rates, i_rate
    integer(kind=i_kind) :: num_aero_reps, i_aero_rep
    integer(kind=i_kind) :: spec_id

    ! Create the new box model
    allocate( new_obj )

    ! Create the CAMP core
    new_obj%camp_core => camp_core_t( config_file )
    call new_obj%camp_core%initialize( )

    ! Get the state species names
    new_obj%spec_names = new_obj%camp_core%unique_names( )

    ! Load the configuration data from the json file
    call j_file%initialize( )
    call j_file%get_core( json )
    json_err_msg = ""
    call assert_msg( 135902099, trim( config_file ).ne."", &
                     "Received empty string for file path" )
    inquire( file = trim( config_file ), exist = file_exists )
    call assert_msg( 181758805, file_exists, "Cannot find file: "// &
                    trim( config_file ) )
    call j_file%load_file( filename = trim( config_file ) )

    ! Get the box model configuration data
    call j_file%get( "camp-box-model", j_box_config, found )
    call assert_msg( 198537480, found, &
                     "Missing 'camp-box-model' configuration data "// &
                     "in configuration file.")
    call json%validate( j_box_config, valid, json_err_msg )
    call assert_msg( 992222859, valid, &
                     "Bad JSON format in 'camp-box-model' configuration "// &
                     "data: "//trim( json_err_msg ) )

    ! Get the number of time steps
    call json%get_child( j_box_config, "output time step [s]", j_obj, found )
    call assert_msg( 912341318, found, &
                     "Missing 'output time step [s]' "// &
                     "in box model data" )
    call json%get( j_obj, real_value )
    new_obj%time_step__s = real_value

    ! Get the integration time
    call json%get_child( j_box_config, "total integration time [s]", j_obj, &
                   found )
    call assert_msg( 442175882, found, &
                     "Missing 'total integration time [s]' "// &
                     "in box model data" )
    call json%get( j_obj, real_value )
    new_obj%num_steps = ceiling( real_value / new_obj%time_step__s )
    new_obj%total_time__s  = real_value

    ! Get the temperature profile
    call json%get_child( j_box_config, "temperature [K]", j_obj, found )
    call assert_msg( 212013359, found, &
                     "Missing 'temperature [K]' "// &
                     "in box model data" )
    new_obj%temperature__K = profile_t( json, j_obj )

    ! Get the pressure profile
    call json%get_child( j_box_config, "pressure [Pa]", j_obj, found )
    call assert_msg( 115322763, found, &
                     "Missing 'pressure [Pa]' "// &
                     "in box model data" )
    new_obj%pressure__Pa = profile_t( json, j_obj )

    ! Get the reaction rate profiles
    call json%get_child( j_box_config, "rates", j_obj, found )
    if( found ) then
      call json%info( j_obj, n_children = num_rates )
      allocate( new_obj%rxn_profiles( num_rates ) )
      j_next => null( )
      call json%get( j_box_config, "rates(1)", j_obj, found )
      call assert( 693229278, found )
      i_rate = 0
      do while( associated( j_obj ) )
        i_rate = i_rate + 1
        new_obj%rxn_profiles( i_rate ) = &
          rxn_profile_t( new_obj%camp_core, json, j_obj )
        call json%get_next( j_obj, j_next )
        j_obj => j_next
      end do
      call assert( 763815888, i_rate .eq. num_rates)
    end if

    ! Get the aerosol representation profiles
    call json%get_child( j_box_config, "aerosol representations", j_obj, &
                         found )
    if( found ) then
      call json%info( j_obj, n_children = num_aero_reps )
      allocate( new_obj%aero_rep_profiles( num_aero_reps ) )
      j_next => null( )
      call json%get( j_box_config, "aerosol representations(1)", j_obj, &
                     found )
      call assert( 877171198, found )
      i_aero_rep = 0
      do while( associated( j_obj ) )
        i_aero_rep = i_aero_rep + 1
        new_obj%aero_rep_profiles( i_aero_rep ) = &
          aero_rep_profile_t( new_obj%camp_core, json, j_obj )
        call json%get_next( j_obj, j_next )
        j_obj => j_next
      end do
      call assert( 472301285, i_aero_rep .eq. num_aero_reps )
    end if

    ! Initialize the solver now that all the update data objects are attached
    ! to reactions
    call new_obj%camp_core%solver_initialize( )

    ! Create a CAMP state for the box model runs and initial conditions
    new_obj%camp_state         => new_obj%camp_core%new_state( )
    new_obj%initial_camp_state => new_obj%camp_core%new_state( )

    ! Get the initial species concentrations
    call json%get( j_box_config, "initial state(1)", j_obj, found )
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
        new_obj%initial_camp_state%state_var( spec_id ) = real_value
        call json%get_next( j_obj, j_next )
        j_obj => j_next
      end do
    end if

    ! Set the initial environtmental conditions
    call new_obj%initial_camp_state%env_states(1)%set_temperature_K( &
      new_obj%temperature__K%current_value( ) )
    call new_obj%initial_camp_state%env_states(1)%set_pressure_Pa( &
      new_obj%pressure__Pa%current_value( ) )
    call new_obj%initial_camp_state%update_env_state( )

    ! free the json file
    call j_file%destroy( )
    call json%destroy( )

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

  !> Create a gnuplot configuration file for plotting box model results
  subroutine create_gnuplot_config_file( this, file_prefix )

    !> Box model data
    class(camp_box_model_data_t), intent(in) :: this
    !> Prefix for the results file (assumes a results file name of
    !> `'file_prefix'_results.txt`)
    character(len=*), intent(in) :: file_prefix

    character(len=:), allocatable :: spec_name, file_name
    integer(kind=i_kind) :: f_unit, i_char, i_spec

    f_unit = SCRIPTS_FILE_UNIT

    ! Create the gnuplot script
    file_name = file_prefix//".conf"
    open(unit=f_unit, file=file_name, status="replace", action="write")
    write(f_unit,*) "# GNU Plot configuration for CAMP box model output"
    write(f_unit,*) "# Run as: gnuplot "//file_name
    write(f_unit,*) "set terminal png truecolor"
    write(f_unit,*) "set autoscale"
    write(f_unit,*) "set xrange [ 0.0:", this%total_time__s, "]"

    ! output environmental conditions
    write(f_unit,*) "set output '"//file_prefix//"_temperature.png"
    write(f_unit,*) "plot\"
    write(f_unit,*) " '"//file_prefix//"_results.txt'\"
    write(f_unit,*) " using 1:2 title 'Temperature [K]'"
    write(f_unit,*) "set output '"//file_prefix//"_pressure.png"
    write(f_unit,*) "plot\"
    write(f_unit,*) " '"//file_prefix//"_results.txt'\"
    write(f_unit,*) " using 1:3 title 'Pressure [Pa]'"

    ! output species concentrations and parameter values
    do i_spec = 1, size( this%spec_names )
      spec_name = this%spec_names( i_spec )%string
      forall( i_char = 1 : len( spec_name ), &
              spec_name( i_char : i_char ) .eq. '/') &
        spec_name( i_char:i_char ) = '_'
      write(f_unit,*) "set output '"//file_prefix//"_"//spec_name//".png'"
      write(f_unit,*) "plot\"
      write(f_unit,*) " '"//file_prefix//"_results.txt'\"
      write(f_unit,*) " using 1:"//trim( to_string( i_spec + 3 ) )// &
                      " title '"//spec_name//"'"
    end do
    close(f_unit)

    deallocate( spec_name )
    deallocate( file_name )

  end subroutine create_gnuplot_config_file

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Print out the configuration of the box model
  subroutine do_print( this, file_unit )

    !> Box model data
    class(camp_box_model_data_t), intent(in) :: this
    !> Output file unit
    integer(kind=i_kind), intent(in), optional :: file_unit

    character(len=*), parameter :: fmt_blank_line = "('')"
    character(len=*), parameter :: fmt_state_hdr = &
      "(' | ',A50,' | ',A13,' |')"
    character(len=*), parameter :: fmt_state_data = &
      "(' | ',A50,' | ',ES13.3,' |')"
    integer(kind=i_kind) :: f_unit, i_rate, i_spec

    f_unit = 6
    if( present( file_unit ) ) f_unit = file_unit

    write(f_unit,*) "**********************"
    write(f_unit,*) "*** CAMP Box Model ***"
    write(f_unit,*) "**********************"
    write(f_unit,fmt_blank_line)
    write(f_unit,*) "total integration time", this%total_time__s, "s"
    write(f_unit,*) "time step size        ", this%time_step__s,  "s"
    write(f_unit,*) "number of time steps  ", this%num_steps
    write(f_unit,fmt_blank_line)
    write(f_unit,*) "** temperature profile **"
    call this%temperature__K%print( f_unit )
    write(f_unit,*) "** end temperature profile **"
    write(f_unit,fmt_blank_line)
    write(f_unit,*) "** pressure [Pa] profile **"
    call this%pressure__Pa%print( f_unit )
    write(f_unit,*) "** end pressure [Pa] profile **"
    write(f_unit,fmt_blank_line)
    write(f_unit,*) "****************************"
    write(f_unit,*) "** Reaction Rate Profiles **"
    write(f_unit,*) "****************************"
    write(f_unit,fmt_blank_line)
    do i_rate = 1, size( this%rxn_profiles )
      call this%rxn_profiles( i_rate )%print( f_unit )
    write(f_unit,fmt_blank_line)
    end do
    write(f_unit,*) "****************************"
    write(f_unit,*) "** Reaction Rate Profiles **"
    write(f_unit,*) "****************************"
    write(f_unit,fmt_blank_line)
    write(f_unit,*) "*******************"
    write(f_unit,*) "** Initial State **"
    write(f_unit,*) "*******************"
    write(f_unit,fmt_blank_line)
    write(f_unit,fmt_state_hdr) "Species", "Value"
    do i_spec = 1, size( this%spec_names )
      write(f_unit,fmt_state_data) this%spec_names( i_spec )%string, &
                                   this%initial_camp_state%state_var( i_spec )
    end do
    write(f_unit,fmt_blank_line)
    call this%camp_core%print( f_unit )
    write(f_unit,fmt_blank_line)
    write(f_unit,*) "***********************"
    write(f_unit,*) "** End Initial State **"
    write(f_unit,*) "***********************"
    write(f_unit,fmt_blank_line)
    write(f_unit,*) "**************************"
    write(f_unit,*) "*** End CAMP Box Model ***"
    write(f_unit,*) "**************************"

  end subroutine do_print

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Finalize the box model
  elemental subroutine finalize( this )

    !> CAMP box model
    type(camp_box_model_data_t), intent(inout) :: this

    if( associated( this%camp_core ) ) &
      deallocate( this%camp_core )
    if( associated( this%camp_state ) ) &
      deallocate( this%camp_state )
    if( associated( this%initial_camp_state) ) &
      deallocate( this%initial_camp_state )
    if( allocated( this%spec_names ) ) &
      deallocate( this%spec_names )
    if( allocated( this%rxn_profiles ) ) &
      deallocate( this%rxn_profiles )
    if( allocated( this%aero_rep_profiles ) ) &
      deallocate( this%aero_rep_profiles )

  end subroutine finalize

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

    type(json_value), pointer :: j_child, j_next, j_val
    integer(kind=json_ik) :: num_trans
    logical(kind=json_lk) :: found
    character(kind=json_ck, len=:), allocatable :: profile_name

    integer(kind=i_kind) :: i_trans

    ! Get the name of the profile
    call json%info( j_obj, name = profile_name )
    new_obj%name = profile_name

    ! Get the transitions
    call json%get_child( j_obj, "transitions", j_child, found )
    if( found ) then
      call json%info( j_obj, n_children = num_trans )
      allocate( new_obj%transition_time_step( 0:num_trans ) )
      allocate( new_obj%transition_value(     0:num_trans ) )
      j_next => null( )
      call json%get( j_obj, "transitions(1)", j_child, found )
      i_trans = 0
      do while( associated( j_child ) )
        i_trans = i_trans + 1

        ! get the time step
        call json%get_child( j_child, "time step", j_val, found )
        call assert_msg( 325008338, found, &
                         "Missing 'time step' from element "// &
                         trim( to_string( i_trans ) )//" in transition '"// &
                         trim( new_obj%name )//"'" )
        call json%get( j_val, new_obj%transition_time_step( i_trans ) )

        ! Get the value
        call json%get_child( j_child, "value", j_val, found )
        call assert_msg( 483183389, found, &
                         "Missing 'value' from element "// &
                         trim( to_string( i_trans ) )//" in transition '"// &
                         trim( new_obj%name )//"'" )
        call json%get( j_val, new_obj%transition_value( i_trans ) )
        call json%get_next( j_child, j_next )
        j_child => j_next
      end do
      call assert( 249923619, i_trans .eq. num_trans )
    else
      allocate( new_obj%transition_time_step( 0:0 ) )
      allocate( new_obj%transition_value(     0:0 ) )
    end if

    ! Get the initial state
    call json%get_child( j_obj, "initial", j_val, found )
    call assert_msg( 382107479, found, &
                     "Missing 'initial' value in profile '"// &
                     trim( new_obj%name )//"'" )
    call json%get( j_val, new_obj%transition_value( 0 ) )
    new_obj%transition_time_step( 0 ) = 1

    call new_obj%reset( )

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

  !> Print the profile configuration
  subroutine profile_do_print( this, file_unit )

    !> Property profile
    class(profile_t), intent(in) :: this
    !> Output file unit
    integer(kind=i_kind), intent(in), optional :: file_unit

    character(len=*), parameter :: fmt_trans_hdr = &
      "('   | ',A13,' | ',A13,' |')"
    character(len=*), parameter :: fmt_trans_data = &
      "('   | ',I13,' | ',ES13.3,' |')"
    integer(kind=i_kind) :: f_unit, i_trans

    f_unit = 6
    if( present( file_unit ) ) f_unit = file_unit

    write(f_unit,*) "Profile name: "//trim(this%name)
    write(f_unit,*) "  ** Transitions     **"
    write(f_unit,fmt_trans_hdr) "Time Step", "Value"
    do i_trans = 0, size( this%transition_time_step ) - 1
      write(f_unit,fmt_trans_data) this%transition_time_step( i_trans ), &
                                   this%transition_value( i_trans )
    end do
    write(f_unit,*) "  ** End Transitions **"

  end subroutine profile_do_print

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef PMC_USE_JSON
  !> Constructor for rxn_profile_t
  function rxn_profile_constructor( camp_core, json, j_obj ) result( new_obj )

    use pmc_mechanism_data
    use pmc_rxn_emission
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
                                               rxn_label ) ) cycle
        if( .not. trim( rxn_label ) .eq. new_obj%name ) cycle
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
        call camp_core%initialize_update_object( rxn, new_obj%update_data )
        found = .true.
        exit
      end do
      if( found ) exit
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

    call camp_core%update_data( this%update_data )

  end subroutine update_rxn

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Finalize the reaction profile
  elemental subroutine rxn_profile_finalize( this )

    !> Reaction profile
    type(rxn_profile_t), intent(inout) :: this

    if( associated( this%update_data ) ) deallocate( this%update_data )

  end subroutine rxn_profile_finalize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef PMC_USE_JSON
  !> Constructor for aero_rep_profile_t
  function aero_rep_profile_constructor( camp_core, json, j_obj ) &
      result( new_obj )

    use pmc_mechanism_data
    use pmc_aero_rep_modal_binned_mass
    use pmc_aero_rep_single_particle
    use json_module

    !> New aerosol representation profile
    type(aero_rep_profile_t) :: new_obj
    !> CAMP core
    type(camp_core_t), intent(inout) :: camp_core
    !> JSON core
    type(json_core), intent(inout)  :: json
    !> JSON object to build the profile from
    type(json_value), pointer, intent(inout) :: j_obj

    class(aero_rep_data_t), pointer :: aero_rep
    type(profile_t) :: base_profile

    character(kind=json_ck, len=:), allocatable :: update_type
    character(kind=json_ck, len=:), allocatable :: section_name
    logical(kind=json_lk) :: found

    integer(kind=i_kind) :: i_aero_rep
    character(len=:), allocatable :: aero_rep_label

    found = .false.

    ! Load the standard profile data
    base_profile                 = profile_t( json, j_obj )
    new_obj%name                 = base_profile%name
    new_obj%transition_time_step = base_profile%transition_time_step
    new_obj%transition_value     = base_profile%transition_value

    ! Get the type of value being updated
    call json%get( j_obj, "update type", update_type, found )
    call assert_msg( 591311579, found, &
                     "Missing update type for aerosol representation "// &
                     "profile '"//new_obj%name//"'" )

    ! Find the reaction and initialize the update data object
    do i_aero_rep = 1, size( camp_core%aero_rep )
      aero_rep => camp_core%aero_rep( i_aero_rep )%val
      if( .not. aero_rep%property_set%get_string( "camp-box-model-id", &
                                                   aero_rep_label ) ) cycle
      if( .not. trim( aero_rep_label ) .eq. new_obj%name ) cycle
      select type( aero_rep )
        class is( aero_rep_modal_binned_mass_t )

          ! Modal/binned updates must include a section id
          call json%get( j_obj, "section name", section_name, found )
          call assert_msg( 454366901, found, &
                           "Missing section name for modal/binned "// &
                           "aerosol representation '"//new_obj%name//"'" )
          call assert_msg( 112599854, aero_rep%get_section_id( section_name, &
                           new_obj%section_id ), &
                           "Cannot find aerosol section '"//section_name// &
                           "' in aerosol representation '"// &
                           new_obj%name//"'" )
          select case( update_type )
            case( "UPDATE_GMD" )
              allocate( aero_rep_update_data_modal_binned_mass_GMD_t::new_obj%update_data )
            case( "UPDATE_GSD" )
              allocate( aero_rep_update_data_modal_binned_mass_GSD_t::new_obj%update_data )
            case default
              call die_msg( 513320382, "Invalid update type for "// &
                            "modal/binned aerosol rep: "//update_type )
          end select
        class is( aero_rep_single_particle_t )
          select case( update_type )
            case( "UPDATE_NUMBER" )
              allocate( aero_rep_update_data_single_particle_number_t::new_obj%update_data )
            case default
              call die_msg( 379396160, "Invalid update type for "// &
                            "single particle aerosol rep: "//update_type )
          end select
        class default
          call die_msg( 203974949, "Invalid aerosol representation " )
      end select
      call camp_core%initialize_update_object( aero_rep, new_obj%update_data )
      found = .true.
      exit
    end do

    call assert_msg( 712992422, found, "Could not find aerosol "// &
                     "representation label '"// trim( new_obj%name ) )

  end function aero_rep_profile_constructor
#endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Update a reaction with the current rate from the profile
  subroutine update_aero_rep( this, camp_core )

    use pmc_aero_rep_modal_binned_mass
    use pmc_aero_rep_single_particle

    !> Aerosol representation profile
    class(aero_rep_profile_t) :: this
    !> CAMP core
    type(camp_core_t), intent(inout) :: camp_core

    class(aero_rep_update_data_t), pointer :: ud

    select type( ud => this%update_data )
      class is( aero_rep_update_data_modal_binned_mass_GMD_t )
        call assert( 832366630, this%section_id.gt.0 )
        call ud%set_GMD( this%section_id, this%current_value( ) )
      class is( aero_rep_update_data_modal_binned_mass_GSD_t )
        call assert( 264057359, this%section_id.gt.0 )
        call ud%set_GSD( this%section_id, this%current_value( ) )
      class is( aero_rep_update_data_single_particle_number_t )
        call ud%set_number__n_m3( 1, this%current_value( ) )
      class default
        call die( 623019372 )
    end select

    call camp_core%update_data( this%update_data )

  end subroutine update_aero_rep

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Finalize the aerosol representation profile
  elemental subroutine aero_rep_profile_finalize( this )

    !> Aerosol representation profile
    type(aero_rep_profile_t), intent(inout) :: this

    if( associated( this%update_data ) ) deallocate( this%update_data )

  end subroutine aero_rep_profile_finalize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_camp_box_model_data
