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
!! have a unique name (the same species name cannot be used for a gas-phase
!! species and an aerosol-phase species). Chemical species may be present
!! either in the gas-phase or any number of \ref phlex_aero_phase
!! "aerosol phases". For example, an aerosol-phase chemical species may be
!! present in an "organic" and an "aqueous" \ref phlex_aero_phase
!! "aerosol phase".
!!
!! Chemical species data include physical constants and species-specific
!! model parameters that are used during initialization to assemble reaction
!! data for use during the model run. Note that chemical species data are
!! **only** available during initialization, and when using MPI are not passed
!! to child nodes. The primary node will, however, have access to data in the
!! \c pmc_phlex_core::phlex_core_t::chem_spec_data object for outputing model
!! data (e.g., species names).
!!
!! The input format for chemical species can be found \ref
!! input_format_species "here".

!> The chem_spec_data_t structure and associated subroutines.
module pmc_chem_spec_data

  use pmc_constants,                  only : i_kind
  use pmc_mpi
  use pmc_util,                       only : die_msg, string_t
  use pmc_property
#ifdef PMC_USE_MPI
  use mpi
#endif
#ifdef PMC_USE_JSON
  use json_module
#endif

  implicit none
  private

  public :: chem_spec_data_t

  !> Reallocation increment
  integer(kind=i_kind), parameter :: REALLOC_INC = 50
  !> Default absolute integration tolerance
  real(kind=dp), parameter :: DEFAULT_ABS_TOL = 1.0e-20

  !> Unknown species type
  integer(kind=i_kind), parameter, public :: UNKNOWN_SPEC = 0
  !> Gas-phase species
  integer(kind=i_kind), parameter, public :: GAS_SPEC = 1
  !> Aerosol-phase species
  integer(kind=i_kind), parameter, public :: AERO_SPEC = 2

  !> Chemical species data
  !!
  !! Time-invariant data related to a \ref phlex_species "chemical species"
  type chem_spec_data_t
    private
    !> Number of species
    integer(kind=i_kind) :: num_spec = 0
    !> Species name
    type(string_t), pointer :: spec_name(:) => null()
    !> Species type
    integer(kind=i_kind), pointer :: spec_type(:) => null()
    !> Species property set
    type(property_t), pointer :: property_set(:) => null()
  contains
    !> Load species from an input file
    procedure :: load
    !> Get the size of the species database
    procedure :: size => get_size
    !> Get the number of species by type
    procedure :: num_spec_by_type
    !> Get an list of species names by type
    procedure :: spec_names_by_type
    !> Check if a species name is in the set of chemical species
    procedure :: exists
    !> Get a species properties
    procedure :: get_property_set
    !> Get a species type
    procedure :: get_type
    !> Get a gas-phase species index in the \c
    !! pmc_phlex_state::phlex_state_t::state_var array.  Note that
    !! aerosol-phase species indices on the \c
    !! pmc_phlex_state::phlex_state_t::state_var array must be accessed from
    !! \c pmc_aero_rep_data::aero_rep_data_t::species_state_id() for a
    !! particular \ref phlex_aero_rep "aerosol representation".
    procedure :: gas_state_id
    !> Get an array of absolute integration tolerances for the gas-phase 
    !! species on the \c pmc_phlex_state::phlex_state_t::state_var array
    procedure :: gas_abs_tol
    !> Get the absolute integration tolerance of a species by name
    procedure :: get_abs_tol

    ! Private functions
    !> Add a species
    procedure, private :: add
    !> Ensure there is enough room in the species dataset to add a
    !! specified number of species
    procedure, private :: ensure_size
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
    integer(i_kind), intent(in), optional :: init_size

    integer(i_kind) :: alloc_size = REALLOC_INC

    if (present(init_size)) alloc_size = init_size
    allocate(new_obj)
    allocate(new_obj%spec_name(alloc_size))
    allocate(new_obj%spec_type(alloc_size))
    allocate(new_obj%property_set(alloc_size))

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Ensure there is enough room in the species dataset to add a specified
  !! number of species
  subroutine ensure_size(this, num_spec)

    !> Species dataset
    class(chem_spec_data_t), intent(inout) :: this
    !> Number of new species to ensure space for
    integer(i_kind), intent(in) :: num_spec

    integer :: new_size
    type(string_t), pointer :: new_name(:)
    integer(kind=i_kind), pointer :: new_type(:)
    type(property_t), pointer :: new_property_set(:)

    if (size(this%spec_name) .ge. this%num_spec + num_spec) return
    new_size = this%num_spec + num_spec + REALLOC_INC
    allocate(new_name(new_size))
    allocate(new_type(new_size))
    allocate(new_property_set(new_size))
    new_name(1:this%num_spec) = this%spec_name(1:this%num_spec)
    new_type(1:this%num_spec) = this%spec_type(1:this%num_spec)
    call this%property_set(1:this%num_spec)%move(new_property_set(1:this%num_spec))
    deallocate(this%spec_name)
    deallocate(this%spec_type)
    deallocate(this%property_set)
    this%spec_name => new_name
    this%spec_type => new_type
    this%property_set => new_property_set

  end subroutine ensure_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Add a new chemical species
  subroutine add(this, spec_name, spec_type, property_set)

    !> Species dataset
    class(chem_spec_data_t), intent(inout) :: this
    !> Name of species to add
    character(len=:), allocatable, intent(in) :: spec_name
    !> Type of chemical species
    integer(kind=i_kind), intent(in) :: spec_type
    !> Property set for new species
    type(property_t), intent(inout), optional :: property_set

    integer(kind=i_kind) :: i_spec

    if (spec_type.ne.GAS_SPEC .and. spec_type.ne.AERO_SPEC) &
            call die_msg(788361092, "Received invalid type for species "// &
                spec_name//": "//to_string(spec_type))
    i_spec = this%find(spec_name)
    if (i_spec.gt.0) then
      if (spec_type.ne.this%spec_type(i_spec)) &
          call die_msg(612991075, "Type mismatch for species "//spec_name)
      call this%property_set(i_spec)%update(property_set)
    else
      call this%ensure_size(1)
      this%num_spec = this%num_spec + 1
      this%spec_name(this%num_spec) = string_t(spec_name)
      this%spec_type(this%num_spec) = spec_type
      if (present(property_set)) then
        this%property_set(this%num_spec) = property_set
      else
        this%property_set(this%num_spec) = property_t()
      end if
    end if

  end subroutine add
          
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> \page input_format_species Input JSON Object Format: Chemical Species
  !!
  !! A \c json object containing information about a \ref phlex_species
  !! "chemical species" of the form:
  !! \code{.json} 
  !! { "pmc-data" : [
  !!   {
  !!     "name" : "my species name",
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
  !!     "type" : "SPEC_TYPE",
  !!     "some property" : 123.34,
  !!     ...
  !!   },
  !!   ...
  !! ]}
  !! \endcode
  !! The key-value pair \b name is required and must contain the unique name
  !! used for this species in the \ref input_format_mechanism
  !! "mechanism object". (The same name may not be used for a gas species and
  !! an aerosol species.) The key-value pair \b type is also required and its
  !! value must be \b GAS_SPEC or \b AERO_SPEC. All remaining data are
  !! optional and may include any valid \c json value, including nested
  !! objects. Multilple entries with the same species name will be merged
  !! into a single species, but duplicate property names for the same species
  !! will cause an error.

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
    integer(kind=i_kind) :: spec_type = UNKNOWN_SPEC
    type(property_t), pointer :: property_set

    allocate(property_set)
    property_set = property_t()

    next => null()
    call json%get_child(j_obj, child)
    do while (associated(child))
      call json%info(child, name=key, var_type=var_type)
      if (key.eq."name") then
        if (var_type.ne.json_string) call die_msg(181339359, &
                "Received non-string species name")
        call json%get(child, unicode_str_val)
        spec_name = unicode_str_val
      else if (key.eq."type") then
        if (var_type.ne.json_string) call die_msg(249541360, &
                "Received non-string species type")
        call json%get(child, unicode_str_val)
        if (unicode_str_val.eq."GAS_SPEC") then
          spec_type = GAS_SPEC
        else if (unicode_str_val.eq."AERO_SPEC") then
          spec_type = AERO_SPEC
        else
          str_val = unicode_str_val
          call die_msg(171550163, "Unknown chemical species type: "// &
                  str_val)
        end if
      else
        call property_set%load(json, child, .false.)
      end if
      call json%get_next(child, next)
      child => next
    end do

    call this%add(spec_name, spec_type, property_set)
#else
  subroutine load(this)

    !> Species dataset
    class(chem_spec_data_t), intent(inout) :: this

    call warn_msg(627397948, "No support for input files")
#endif
  end subroutine load

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get the size of a species database
  integer(kind=i_kind) function get_size(this)

    !> Species database
    class(chem_spec_data_t), intent(in) :: this

    get_size = this%num_spec

  end function get_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get the number of species of a certain type
  integer(kind=i_kind) function num_spec_by_type(this, &
                  spec_type) result (num_spec)

    !> Species database
    class(chem_spec_data_t), intent(in) :: this
    !> Species type
    integer(kind=i_kind), intent(in) :: spec_type

    integer(kind=i_kind) :: i_spec

    num_spec = 0
    do i_spec = 1, this%num_spec
      if (this%spec_type(i_spec).eq.spec_type) num_spec = num_spec + 1
    end do

  end function num_spec_by_type

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get a list of species names by type
  function spec_names_by_type(this, spec_type) &
                  result (spec_names)

    !> List of species names
    type(string_t), allocatable :: spec_names(:)
    !> Species database
    class(chem_spec_data_t), intent(in) :: this
    !> Species type
    integer(kind=i_kind), intent(in) :: spec_type

    integer(kind=i_kind) :: i_spec, j_spec

    allocate(spec_names(this%num_spec_by_type(spec_type)))
    j_spec = 1
    do i_spec = 1, this%num_spec
      if (this%spec_type(i_spec).eq.spec_type) then
        spec_names(j_spec) = this%spec_name(i_spec) 
        j_spec = j_spec + 1
      end if
    end do

  end function spec_names_by_type

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Check if a species name is in the set of chemical species
  logical function exists(this, spec_name) result (found)

    !> Species dataset
    class(chem_spec_data_t), intent(in) :: this
    !> Species name
    character(len=:), allocatable, intent(in) :: spec_name

    found = this%find(spec_name).ne.0

  end function exists

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get the index of a chemical species by name. Return 0 if the species is
  !! not found.
  integer(kind=i_kind) function find(this, spec_name) &
      result (spec_id)

    !> Species dataset
    class(chem_spec_data_t), intent(in) :: this
    !> Species name
    character(len=:), allocatable, intent(in) :: spec_name

    integer(kind=i_kind) :: i_spec

    spec_id = 0
    do i_spec = 1, this%num_spec
      if (this%spec_name(i_spec)%string .eq. spec_name) then
        spec_id = i_spec
        return
      end if
    end do

  end function find

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get a species property set
  function get_property_set(this, spec_name) result (property_set)

    !> Pointer to species properties
    class(property_t), pointer :: property_set
    !> Species dataset
    class(chem_spec_data_t), intent(in) :: this
    !> Species name to find properties of
    character(len=:), allocatable, intent(in) :: spec_name

    integer(i_kind) :: spec_id

    property_set => null()
    spec_id = this%find(spec_name)
    if (spec_id.gt.0) property_set => this%property_set(spec_id)

  end function get_property_set

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get a species type
  function get_type(this, spec_name) result(spec_type)

    !> Species type
    integer(kind=i_kind) :: spec_type
    !> Species dataset
    class(chem_spec_data_t), intent(in) :: this
    !> Species name to find properties of
    character(len=:), allocatable, intent(in) :: spec_name

    integer(i_kind) :: spec_id

    spec_type = UNKNOWN_SPEC
    spec_id = this%find(spec_name)
    if (spec_id.gt.0) spec_type = this%spec_type(spec_id)

  end function get_type

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get a gas-phase species index in the \c
  !! pmc_phlex_state::phlex_state_t::state_var array.  Note that
  !! aerosol-phase species indices on the \c
  !! pmc_phlex_state::phlex_state_t::state_var array must be accessed from
  !! \c pmc_aero_rep_data::aero_rep_data_t::species_state_id() for a
  !! particular \ref phlex_aero_rep "aerosol representation".
  integer(kind=i_kind) function gas_state_id(this, spec_name)

    !> Species dataset
    class(chem_spec_data_t), intent(in) :: this
    !> Species name
    character(len=:), allocatable, intent(in) :: spec_name

    gas_state_id = this%find(spec_name)

  end function gas_state_id

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Add gas-phase species tolerances to the absolute integration tolerance
  !! array
  subroutine gas_abs_tol(this, abs_tol)

    !> Species dataset
    class(chem_spec_data_t), intent(in) :: this
    !> Array of species absolute integration tolerances. This array should be
    !! the same size as the pmc_phlex_state::phlex_state_t::state_var array
    real(kind=dp), allocatable :: abs_tol(:)

    integer(kind=i_kind) :: i_spec, state_id

    state_id = 1
    do i_spec = 1, this%num_spec
      if (this%spec_type(i_spec).eq.GAS_SPEC) then
        abs_tol(state_id) = this%get_abs_tol(this%spec_name(i_spec)%string)
        state_id = state_id + 1
      end if
    end do

  end subroutine gas_abs_tol

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get the absolute integration tolerance of a species by name
  function get_abs_tol(this, spec_name) result (abs_tol)

    !> Absolute integration tolerance
    real(kind=dp) :: abs_tol
    !> Species dataset
    class(chem_spec_data_t), intent(in) :: this
    !> Species name
    character(len=:), allocatable :: spec_name

    character(len=:), allocatable :: key
    integer(kind=i_kind) :: i_spec
    real(kind=dp) :: val

    key = "absolute integration tolerance"
    i_spec = this%find(spec_name)
    call assert_msg(956793138, i_spec.gt.0, "Could not find species "// &
            spec_name)
    if (this%property_set(i_spec)%get_real(key, val)) then
      abs_tol = val
    else
      abs_tol = DEFAULT_ABS_TOL
    end if

  end function get_abs_tol

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_chem_spec_data
