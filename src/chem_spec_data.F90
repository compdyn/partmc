! Copyright (C) 2017 Matt Dawson
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_chem_spec_data module.

!> \page phlex_species Phlexible Module for Chemistry: Chemical Species
!!
!! Chemical species in the \ref phlex_chem "phlex-chem" module are gas- or
!! aerosol-phase species that can participate in chemical or phase-transfer
!! reactions in the \ref phlex_mechanism "mechanism(s)". Each species must
!! have a unique name (the same name cannot be used for a gas-phase
!! species and an aerosol-phase species). Chemical species may be present
!! either in the gas-phase or any number of \ref phlex_aero_phase
!! "aerosol phases". For example, an aerosol-phase chemical species may be
!! present in an "organic" and an "aqueous" \ref phlex_aero_phase
!! "aerosol phase".
!!
!! Chemical species data include physical constants and species-specific
!! model parameters that are used during initialization to assemble reaction
!! and sub-model data for use during solving. Note that chemical species data
!! are  **only** available during initialization, and when using MPI are not
!! passed to child nodes. The primary node will, however, have access to data
!! in the \c pmc_phlex_core::phlex_core_t::chem_spec_data object for outputing
!! model data (e.g., species names).
!!
!! The input format for chemical species can be found \ref
!! input_format_species "here".

!> The chem_spec_data_t structure and associated subroutines.
module pmc_chem_spec_data

#ifdef PMC_USE_JSON
  use json_module
#endif
  use pmc_constants,                  only : phlex_real, phlex_int
  use pmc_property
  use pmc_util,                       only : die_msg, string_t, assert_msg

  use iso_c_binding

  implicit none
  private

  public :: chem_spec_data_t 

  !> State variable types (Must match values in phlex_solver.c)
  integer(kind=phlex_int), parameter, public :: CHEM_SPEC_UNKNOWN_TYPE = 0
  integer(kind=phlex_int), parameter, public :: CHEM_SPEC_VARIABLE = 1
  integer(kind=phlex_int), parameter, public :: CHEM_SPEC_CONSTANT = 2
  integer(kind=phlex_int), parameter, public :: CHEM_SPEC_PSSA = 3
  integer(kind=phlex_int), parameter, public :: CHEM_SPEC_ACTIVITY_COEFF = 4

  !> Species phase
  integer(kind=phlex_int), parameter, public :: CHEM_SPEC_UNKNOWN_PHASE = 0
  integer(kind=phlex_int), parameter, public :: CHEM_SPEC_GAS_PHASE = 1
  integer(kind=phlex_int), parameter, public :: CHEM_SPEC_AERO_PHASE = 2

  !> Reallocation increment
  integer(kind=phlex_int), parameter :: REALLOC_INC = 50
  !> Default absolute integration tolerance
  real(kind=phlex_real), parameter :: DEFAULT_ABS_TOL = 1.0e-14

  !> Chemical species data
  !!
  !! Time-invariant data related to a \ref phlex_species "chemical species"
  type chem_spec_data_t
    private
    !> Number of species
    integer(kind=phlex_int) :: num_spec = 0
    !> Species name
    type(string_t), pointer :: spec_name(:) => null()
    !> Species type
    integer(kind=phlex_int), pointer :: spec_type(:) => null()
    !> Species phase
    integer(kind=phlex_int), pointer :: spec_phase(:) => null()
    !> Species property set
    type(property_t), pointer :: property_set(:) => null()
  contains
    !> Load species from an input file
    procedure :: load
    !> Initialize the species set
    procedure :: initialize
    !> Get the number of species with specified conditions
    procedure :: size => get_size
    !> Check if a species name is in the set of chemical species
    procedure :: exists
    !> Get a list of species names
    procedure :: get_spec_names
    !> Get a species properties
    procedure :: get_property_set
    !> Get a species type
    procedure :: get_type
    !> Get a species phase
    procedure :: get_phase
    !> Get the absolute integration tolerance of a species
    procedure :: get_abs_tol
    !> Get a gas-phase species index in the \c
    !! pmc_phlex_state::phlex_state_t::state_var array.  Note that
    !! aerosol-phase species indices on the \c
    !! pmc_phlex_state::phlex_state_t::state_var array must be accessed from
    !! \c pmc_aero_rep_data::aero_rep_data_t::spec_state_id() for a
    !! particular \ref phlex_aero_rep "aerosol representation".
    procedure :: gas_state_id
    !> Get the name of a gas-phase species in the \c
    !! pmc_phlex_state::phlex_state_t::state_var array.  Note that
    !! aerosol-phase species names on the \c
    !! pmc_phlex_state::phlex_state_t::state_var array must be accessed from
    !! \c pmc_aero_rep_data::aero_rep_data_t::spec_state_id() for a
    !! particular \ref phlex_aero_rep "aerosol representation".
    procedure :: gas_state_name
    !> Print out the species data
    procedure :: print => do_print
    !> Finalize the chemical species data
    final :: finalize

    ! Private functions
    !> Ensure there is enough room in the species dataset to add a
    !! specified number of species
    procedure, private :: ensure_size
    !> Add a species
    procedure, private :: add
    !> Find a species index by name
    procedure, private :: find
  end type chem_spec_data_t

  ! Constructor for chem_spec_data_t
  interface chem_spec_data_t
    procedure :: constructor
  end interface chem_spec_data_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Constructor for chem_spec_data_t
  function constructor(init_size) result(new_obj)

    !> A new set of chemical species
    type(chem_spec_data_t), pointer :: new_obj
    !> Number of species to allocate space for initially
    integer(kind=phlex_int), intent(in), optional :: init_size

    integer(kind=phlex_int) :: alloc_size = REALLOC_INC

    if (present(init_size)) alloc_size = init_size
    allocate(new_obj)
    allocate(new_obj%spec_name(alloc_size))
    allocate(new_obj%spec_type(alloc_size))
    allocate(new_obj%spec_phase(alloc_size))
    allocate(new_obj%property_set(alloc_size))

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> \page input_format_species Input JSON Object Format: Chemical Species
  !!
  !! A \c json object containing information about a \ref phlex_species
  !! "chemical species" has the following format:
  !! \code{.json} 
  !! { "pmc-data" : [
  !!   {
  !!     "name" : "my species name",
  !!     "type" : "CHEM_SPEC",
  !!     "phase" : "SPEC_PHASE",
  !!     "type" : "SPEC_TYPE",
  !!     "some property" : 123.34,
  !!     "some other property" : true,
  !!     "nested properties" : {
  !!        "sub prop 1" : 12.43,
  !!        "sub prop other" : "some text"
  !!     },
  !!     ...
  !!   },
  !!   {
  !!     "name" : "my species name",
  !!     "type" : "CHEM_SPEC",
  !!     "phase" : "SPEC_PHASE",
  !!     "type" : "SPEC_TYPE",
  !!     "some property" : 123.34,
  !!     ...
  !!   },
  !!   ...
  !! ]}
  !! \endcode
  !! The key-value pair \b name is required and must contain the unique name
  !! used for this species in the \ref input_format_mechanism
  !! "mechanism object". (The same name cannot be used for a gas-phase species
  !! and an aerosol-phase species.) The key-value pair \b type is also
  !! required, and must be \b CHEM_SPEC.
  !!
  !! The key-value pair \b phase specifies the
  !! phase in which the species exists and can be \b GAS or \b AEROSOL. When
  !! the \b phase is not specified, it is assumed to be \b GAS. The \b type
  !! can be \b VARIABLE, \b CONSTANT or \b PSSA. When a \b type is not
  !! specified, it is assumed to be \b VARIABLE.
  !!
  !! All remaining data are optional and may include any valid \c json value,
  !! including nested objects. Multilple entries with the same species name
  !! will be merged into a single species, but duplicate property names for
  !! the same species will cause an error. However, nested objects with the
  !! same key name will be merged, if possible.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Load species from an input file
#ifdef PMC_USE_JSON
  subroutine load(this, json, j_obj)

    !> Species dataset
    class(chem_spec_data_t), intent(inout) :: this
    !> JSON core
    type(json_core), pointer, intent(in) :: json
    !> JSON object
    type(json_value), pointer, intent(in) :: j_obj

    type(json_value), pointer :: child, next
    character(kind=json_ck, len=:), allocatable :: key, unicode_str_val
    integer(kind=json_ik) :: var_type

    character(len=:), allocatable :: spec_name, str_val
    integer(kind=phlex_int) :: spec_type, spec_phase
    type(property_t), pointer :: property_set

    ! allocate space for the species property set
    property_set => property_t()

    ! initialize type and phase
    spec_type = CHEM_SPEC_UNKNOWN_TYPE
    spec_phase = CHEM_SPEC_UNKNOWN_PHASE
   
    ! cycle through the species properties to find the name, type and phase
    ! and load the remaining data into the species property set
    next => null()
    call json%get_child(j_obj, child)
    do while (associated(child))
      call json%info(child, name=key, var_type=var_type)
      
      ! species name
      if (key.eq."name") then
        if (var_type.ne.json_string) call die_msg(181339359, &
                "Received non-string species name")
        call json%get(child, unicode_str_val)
        spec_name = unicode_str_val
      
      ! variable type
      else if (key.eq."tracer type") then
        if (var_type.ne.json_string) call die_msg(249541360, &
                "Received non-string species type")
        call json%get(child, unicode_str_val)
        if (unicode_str_val.eq."VARIABLE") then
          spec_type = CHEM_SPEC_VARIABLE
        else if (unicode_str_val.eq."CONSTANT") then
          spec_type = CHEM_SPEC_CONSTANT
        else if (unicode_str_val.eq."PSSA") then
          spec_type = CHEM_SPEC_PSSA
        else if (unicode_str_val.eq."ION_PAIR") then
          spec_type = CHEM_SPEC_ACTIVITY_COEFF
          spec_phase = CHEM_SPEC_AERO_PHASE
        else
          str_val = unicode_str_val
          call die_msg(171550163, "Unknown chemical species type: "// &
                  str_val)
        end if
      
      ! species phase
      else if (key.eq."phase") then
        if (var_type.ne.json_string) call die_msg(273466657, &
                "Received non-string species phase")
        call json%get(child, unicode_str_val)
        if (unicode_str_val.eq."GAS") then
          call assert_msg(334314929, spec_type.ne.CHEM_SPEC_ACTIVITY_COEFF, &
                "Activity coefficients and ion pairs cannot be gas-phase "// &
                "species.")
          spec_phase = CHEM_SPEC_GAS_PHASE
        else if (unicode_str_val.eq."AEROSOL") then
          spec_phase = CHEM_SPEC_AERO_PHASE
        else
          str_val = unicode_str_val
          call die_msg(603704146, "Unknown chemical species phase: "// &
                  str_val)
        end if

      ! load remaining properties into the species property set
      else if (key.ne."type") then
        call property_set%load(json, child, .false.)
      end if

      call json%get_next(child, next)
      child => next
    end do

    ! Add or update the species data
    call this%add(spec_name, spec_type, spec_phase, property_set)

    ! deallocate the temporary property set
    deallocate(property_set)
#else
  subroutine load(this)

    !> Species dataset
    class(chem_spec_data_t), intent(inout) :: this

    call warn_msg(627397948, "No support for input files")
#endif
  end subroutine load

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Initialize the species set
  subroutine initialize(this)

    !> Species database
    class(chem_spec_data_t), intent(inout) :: this

    ! Species index
    integer(kind=phlex_int) :: i_spec

    do i_spec = 1, this%num_spec

      ! Set default value for type
      if (this%spec_type(i_spec).eq.CHEM_SPEC_UNKNOWN_TYPE) &
              this%spec_type(i_spec) = CHEM_SPEC_VARIABLE
      
      ! Set default value for phase
      if (this%spec_phase(i_spec).eq.CHEM_SPEC_UNKNOWN_PHASE) &
              this%spec_phase(i_spec) = CHEM_SPEC_GAS_PHASE

    end do

  end subroutine initialize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get the number of species with the given properties. If no properties
  !! are specified, return the total number of species. 
  integer(kind=phlex_int) function get_size(this, spec_type, spec_phase) &
                  result (num_spec)

    !> Species database
    class(chem_spec_data_t), intent(in) :: this
    !> State variable type for the species
    integer(kind=phlex_int), intent(in), optional :: spec_type
    !> Phase of the species
    integer(kind=phlex_int), intent(in), optional :: spec_phase

    ! species index
    integer(kind=phlex_int) :: i_spec

    ! add up the number of species
    if (present(spec_type).or.present(spec_phase)) then
      num_spec = 0
      do i_spec = 1, this%num_spec
        if (present(spec_type)) then
          if (spec_type.ne.this%spec_type(i_spec)) cycle
        end if
        if (present(spec_phase)) then
          if (spec_phase.ne.this%spec_phase(i_spec)) cycle
        end if
        num_spec = num_spec + 1
      end do
    else
      num_spec = this%num_spec
    end if

  end function get_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Check if a species name is in the set of chemical species
  logical function exists(this, spec_name) result (found)

    !> Species dataset
    class(chem_spec_data_t), intent(in) :: this
    !> Species name
    character(len=:), allocatable, intent(in) :: spec_name

    ! Index of species
    integer(kind=phlex_int) :: i_spec

    found = this%find(spec_name, i_spec)

  end function exists

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get a list of species names
  !!
  !! If a type or phase are specified, only species with the specified values
  !! are returned. Otherwise, all species names are returned.
  function get_spec_names(this, spec_type, spec_phase) result (spec_names)

    !> Species names
    type(string_t), allocatable :: spec_names(:)
    !> Species dataset
    class(chem_spec_data_t), intent(in) :: this
    !> State variable type for the species
    integer(kind=phlex_int), intent(in), optional :: spec_type
    !> Phase of the species
    integer(kind=phlex_int), intent(in), optional :: spec_phase

    ! species index and counter
    integer(kind=phlex_int) :: i_spec, num_spec

    ! add up the number of species
    if (present(spec_type).or.present(spec_phase)) then
      num_spec = 0
      do i_spec = 1, this%num_spec
        if (present(spec_type)) then
          if (spec_type.ne.this%spec_type(i_spec)) cycle
        end if
        if (present(spec_phase)) then
          if (spec_phase.ne.this%spec_phase(i_spec)) cycle
        end if
        num_spec = num_spec + 1
      end do
    else
      num_spec = this%num_spec
    end if

    ! Allocate space for the list
    allocate(spec_names(num_spec))

    ! Add species names to the list
    num_spec = 0
    do i_spec = 1, this%num_spec
      if (present(spec_type)) then
        if (spec_type.ne.this%spec_type(i_spec)) cycle
      end if
      if (present(spec_phase)) then
        if (spec_phase.ne.this%spec_phase(i_spec)) cycle
      end if
      num_spec = num_spec + 1
      spec_names(num_spec)%string = this%spec_name(i_spec)%string
    end do

  end function get_spec_names

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get a species property set. Returns true if the species is found, or
  !! false otherwise.
  logical function get_property_set(this, spec_name, property_set) &
                  result (found)

    !> Species dataset
    class(chem_spec_data_t), intent(in) :: this
    !> Species name to find properties of
    character(len=:), allocatable, intent(in) :: spec_name
    !> Pointer to species properties
    type(property_t), pointer, intent(out) :: property_set

    integer(kind=phlex_int) :: spec_id

    property_set => null()
    found = this%find(spec_name, spec_id)
    if (found) property_set => this%property_set(spec_id)

  end function get_property_set

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get a species type by species name. Returns true if the species is found
  !! or false otherwise.
  logical function get_type(this, spec_name, spec_type) result (found)

    !> Species dataset
    class(chem_spec_data_t), intent(in) :: this
    !> Species name to find properties of
    character(len=:), allocatable, intent(in) :: spec_name
    !> Species type
    integer(kind=phlex_int), intent(out) :: spec_type

    integer(kind=phlex_int) :: spec_id

    spec_type = CHEM_SPEC_UNKNOWN_TYPE
    found = this%find(spec_name, spec_id)
    if (found) spec_type = this%spec_type(spec_id)

  end function get_type

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get a species phase by name. Returns true if the species is found, or
  !! false otherwise.
  logical function get_phase(this, spec_name, spec_phase) &
                  result (found)

    !> Species dataset
    class(chem_spec_data_t), intent(in) :: this
    !> Species name to find properties of
    character(len=:), allocatable, intent(in) :: spec_name
    !> Species phase
    integer(kind=phlex_int), intent(out) :: spec_phase

    integer(kind=phlex_int) :: spec_id

    spec_phase = CHEM_SPEC_UNKNOWN_PHASE
    found = this%find(spec_name, spec_id)
    if (found) spec_phase = this%spec_phase(spec_id)

  end function get_phase

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get the absolute integration tolerance of a species by name. Returns 
  !! true if the species is found, or false otherwise.
  logical function get_abs_tol(this, spec_name, abs_tol) &
                  result (found)

    !> Species dataset
    class(chem_spec_data_t), intent(in) :: this
    !> Species name
    character(len=:), allocatable, intent(in) :: spec_name
    !> Absolute integration tolerance
    real(kind=phlex_real), intent(out) :: abs_tol

    character(len=:), allocatable :: key
    integer(kind=phlex_int) :: spec_id
    real(kind=phlex_real) :: val

    abs_tol = DEFAULT_ABS_TOL
    key = "absolute integration tolerance"
    found = this%find(spec_name, spec_id)
    if (found) then
      if (this%property_set(spec_id)%get_real(key, val)) then
        abs_tol = val
      end if
    end if

  end function get_abs_tol

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get a gas-phase species index in the \c
  !! pmc_phlex_state::phlex_state_t::state_var array.  Note that
  !! aerosol-phase species indices on the \c
  !! pmc_phlex_state::phlex_state_t::state_var array must be accessed from
  !! \c pmc_aero_rep_data::aero_rep_data_t::spec_state_id() for a
  !! particular \ref phlex_aero_rep "aerosol representation". Returns a valid
  !! state array index if the species is found, or 0 otherwise
  pure integer(kind=phlex_int) function gas_state_id(this, spec_name)

    !> Species dataset
    class(chem_spec_data_t), intent(in) :: this
    !> Species name
    character(len=:), allocatable, intent(in) :: spec_name

    integer(kind=phlex_int) :: i_spec
    
    gas_state_id = 0
    do i_spec = 1, this%num_spec
      if (this%spec_phase(i_spec).eq.CHEM_SPEC_GAS_PHASE) then
              gas_state_id = gas_state_id + 1
        if (trim(spec_name).eq.trim(this%spec_name(i_spec)%string)) return
      end if
    end do
    gas_state_id = 0

  end function gas_state_id

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get a gas-phase species name in the \c
  !! pmc_phlex_state::phlex_state_t::state_var array.  Note that
  !! aerosol-phase species names on the \c
  !! pmc_phlex_state::phlex_state_t::state_var array must be accessed from
  !! \c pmc_aero_rep_data::aero_rep_data_t::spec_state_id() for a
  !! particular \ref phlex_aero_rep "aerosol representation". Returns a valid
  !! state array index if the species is found, or 0 otherwise
  function gas_state_name(this, spec_id) result(spec_name)

    !> Species name
    character(len=:), allocatable :: spec_name
    !> Species dataset
    class(chem_spec_data_t), intent(in) :: this
    !> Species id
    integer(kind=phlex_int), intent(in) :: spec_id

    integer(kind=phlex_int) :: gas_state_id, i_spec
    
    gas_state_id = 0
    do i_spec = 1, this%num_spec
      if (this%spec_phase(i_spec).eq.CHEM_SPEC_GAS_PHASE) &
              gas_state_id = gas_state_id + 1
      if (gas_state_id.eq.spec_id) then
        spec_name = trim(this%spec_name(i_spec)%string)
        return
      end if
    end do
    spec_name = ""

  end function gas_state_name

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Print out the species data
  subroutine do_print(this, file_unit)

    !> Chemical species data
    class(chem_spec_data_t), intent(in) :: this
    !> File unit for output
    integer(kind=phlex_int), optional, intent(in) :: file_unit

    integer(kind=phlex_int) :: i_spec
    integer(kind=phlex_int) :: f_unit
    character(len=:), allocatable :: spec_phase, spec_type

    f_unit = 6
    if (present(file_unit)) f_unit = file_unit

    write(f_unit,*) "Number of species: ", this%num_spec
    do i_spec = 1, this%num_spec
      select case (this%spec_type(i_spec))
        case (CHEM_SPEC_UNKNOWN_TYPE)
          spec_type = "unknown"
        case (CHEM_SPEC_VARIABLE)
          spec_type = "variable"
        case (CHEM_SPEC_CONSTANT)
          spec_type = "constant"
        case (CHEM_SPEC_PSSA)
          spec_type = "PSSA"
        case (CHEM_SPEC_ACTIVITY_COEFF)
          spec_type = "binary activity coefficient"
        case default
          spec_type = "invalid"
      end select
      select case (this%spec_phase(i_spec))
        case (CHEM_SPEC_UNKNOWN_PHASE)
          spec_phase = "unknown"
        case (CHEM_SPEC_GAS_PHASE)
          spec_phase = "gas"
        case (CHEM_SPEC_AERO_PHASE)
          spec_phase = "aerosol"
        case default
          spec_phase = "invalid"
      end select

      write(f_unit,*) "  ", this%spec_name(i_spec)%string
      write(f_unit,*) "    phase: ", spec_phase, "; type: ", spec_type
      call this%property_set(i_spec)%print(f_unit)
    end do

  end subroutine do_print

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Finalize the chemical species data
  elemental subroutine finalize(this)

    !> Species dataset
    type(chem_spec_data_t), intent(inout) :: this

    deallocate(this%spec_name)
    deallocate(this%spec_type)
    deallocate(this%spec_phase)
    deallocate(this%property_set)

  end subroutine finalize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Ensure there is enough room in the species dataset to add a specified
  !! number of species
  subroutine ensure_size(this, num_spec)

    !> Species dataset
    class(chem_spec_data_t), intent(inout) :: this
    !> Number of new species to ensure space for
    integer(kind=phlex_int), intent(in) :: num_spec

    integer :: new_size
    type(string_t), pointer :: new_name(:)
    integer(kind=phlex_int), pointer :: new_type(:)
    integer(kind=phlex_int), pointer :: new_phase(:)
    type(property_t), pointer :: new_property_set(:)

    if (size(this%spec_name) .ge. this%num_spec + num_spec) return
    new_size = this%num_spec + num_spec + REALLOC_INC
    allocate(new_name(new_size))
    allocate(new_type(new_size))
    allocate(new_phase(new_size))
    allocate(new_property_set(new_size))
    new_name(1:this%num_spec) = this%spec_name(1:this%num_spec)
    new_type(1:this%num_spec) = this%spec_type(1:this%num_spec)
    new_phase(1:this%num_spec) = this%spec_phase(1:this%num_spec)
    call this%property_set(1:this%num_spec)%move(new_property_set( &
            1:this%num_spec))
    deallocate(this%spec_name)
    deallocate(this%spec_type)
    deallocate(this%spec_phase)
    deallocate(this%property_set)
    this%spec_name => new_name
    this%spec_type => new_type
    this%spec_phase => new_phase
    this%property_set => new_property_set

  end subroutine ensure_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Add a new chemical species
  subroutine add(this, spec_name, spec_type, spec_phase, property_set)

    !> Species dataset
    class(chem_spec_data_t), intent(inout) :: this
    !> Name of species to add
    character(len=:), allocatable, intent(in) :: spec_name
    !> State variable type
    integer(kind=phlex_int), intent(inout) :: spec_type
    !> Species phase
    integer(kind=phlex_int), intent(inout) :: spec_phase
    !> Property set for new species
    type(property_t), intent(inout), optional :: property_set

    integer(kind=phlex_int) :: i_spec

    ! if the species exists, append the new data
    if (this%find(spec_name, i_spec)) then
     
      ! Check for a type mismatch
      if (spec_type.eq.CHEM_SPEC_UNKNOWN_TYPE) &
              spec_type = this%spec_type(i_spec)
      call assert_msg(596247182, &
              this%spec_type(i_spec).eq.CHEM_SPEC_UNKNOWN_TYPE &
              .or.spec_type.eq.this%spec_type(i_spec), &
              "Type mismatch for species "//spec_name)

      ! Check for a phase mismatch
      if (spec_phase.eq.CHEM_SPEC_UNKNOWN_PHASE) &
              spec_phase = this%spec_phase(i_spec)
      call assert_msg(612991075, &
              this%spec_phase(i_spec).eq.CHEM_SPEC_UNKNOWN_PHASE &
              .or.spec_phase.eq.this%spec_phase(i_spec), &
              "Phase mismatch for species "//spec_name)

      ! Update the species properties
      this%spec_type(i_spec) = spec_type
      this%spec_phase(i_spec) = spec_phase
      call this%property_set(i_spec)%update(property_set)

    ! ... otherwise, create a new species
    else
      call this%ensure_size(1)
      this%num_spec = this%num_spec + 1
      this%spec_name(this%num_spec)%string = spec_name
      this%spec_type(this%num_spec) = spec_type
      this%spec_phase(this%num_spec) = spec_phase
      if (present(property_set)) then
        call property_set%move(this%property_set(this%num_spec))
      else
        this%property_set(this%num_spec) = property_t()
      end if
    end if

  end subroutine add
          
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get the index of a chemical species by name. Returns true if the species
  !! is found or false otherwise.
  logical function find(this, spec_name, spec_id)

    !> Species dataset
    class(chem_spec_data_t), intent(in) :: this
    !> Species name
    character(len=:), allocatable, intent(in) :: spec_name
    !> Species id
    integer(kind=phlex_int), intent(out) :: spec_id

    find = .true.
    do spec_id = 1, this%num_spec
      if (this%spec_name(spec_id)%string .eq. spec_name) then
        return
      end if
    end do
    find = .false.
    spec_id = 0

  end function find

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_chem_spec_data
