! Copyright (C) 2017 Matt Dawson
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

! TODO Determine if we can avoid sending this data through MPI.
! Everything needed to solve the chemistry should be in the state
! objects and the reaction data

!> \file
!> The pmc_chem_spec_data module.

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
  real(kind=dp), parameter :: DEFAULT_ABS_TOL = 1.0e-30

  !> Unknown species type
  integer(kind=i_kind), parameter, public :: UNKNOWN_SPEC = 0
  !> Gas-phase species
  integer(kind=i_kind), parameter, public :: GAS_SPEC = 1
  !> Aerosol-phase species
  integer(kind=i_kind), parameter, public :: AERO_SPEC = 2

  !> Chemical species data
  !!
  !! Time-invariant data related to a chemical species
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
    procedure :: load => chem_spec_data_load
    !> Get the size of the species database
    procedure :: size => chem_spec_data_size
    !> Check if a species name is in the set of chemical species
    procedure :: exists => chem_spec_data_exists
    !> Get a species properties
    procedure :: get_property_set => chem_spec_data_get_property_set
    !> Get a species type
    procedure :: get_type => chem_spec_data_get_type
    !> Get a species index in the state variable
    procedure :: state_id => chem_spec_data_state_id
    !> Get an array of absolute integration tolerances corresponding
    !! to the dimensions of the state array
    procedure :: get_abs_tolerances => chem_spec_data_get_abs_tolerances

    !> Private functions
    !> Add a species
    procedure, private :: add => chem_spec_data_add
    !> Ensure there is enough room in the species dataset to add a
    !! specified number of species
    procedure, private :: ensure_size => chem_spec_data_ensure_size
    !> Find a species index by name
    procedure, private :: find => chem_spec_data_find
  end type chem_spec_data_t

  !> Constructor for chem_spec_data_t
  interface chem_spec_data_t
    procedure :: chem_spec_data_constructor
  end interface chem_spec_data_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Constructor for chem_spec_data_t
  function chem_spec_data_constructor(init_size) result(new_obj)

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

  end function chem_spec_data_constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Ensure there is enough room in the species dataset to add a specified
  !! number of species
  subroutine chem_spec_data_ensure_size(this, num_spec)

    !> Species dataset
    class(chem_spec_data_t), intent(inout) :: this
    !> Number of new species to ensure space for
    integer(i_kind) :: num_spec

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

  end subroutine chem_spec_data_ensure_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Add a new chemical species
  subroutine chem_spec_data_add(this, spec_name, spec_type, property_set)

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

  end subroutine chem_spec_data_add
          
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Load species from an input file
#ifdef PMC_USE_JSON
  !! j_obj is expected to be a JSON object containing data related to a
  !! chemical species required for building the chemical mechanism. It should
  !! be of the form:
  !! 
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
  !!
  !! The key-value pair "name" is required and must contain the unique name
  !! used for this species in the mechanism. (The same name may not be used
  !! for a gas species and an aerosol species.) The key-value pair "type" is 
  !! also required and its value must be GAS_SPEC or AERO_SPEC. All remaining
  !! data is optional and may include any valid JSON value, including nested
  !! objects. Multilple entries with the same species name will be merged
  !! into a single species, but duplicate property names for the same species
  !! will cause an error.
  !! Species data may be inter-mixed with json objects of other types (e.g.,
  !! reactions), but the there must be exactly one top-level key-value pair 
  !! named "pmc-data" per input file whose value is an array of json objects 
  !! with valid PMC types.
  subroutine chem_spec_data_load(this, json, j_obj)

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
  subroutine chem_spec_data_load(this)

    !> Species dataset
    class(chem_spec_data_t), intent(inout) :: this

    call warn_msg(627397948, "No support for input files")
#endif
  end subroutine chem_spec_data_load

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get the size of a species database
  integer(kind=i_kind) function chem_spec_data_size(this)

    !> Species database
    class(chem_spec_data_t), intent(in) :: this

    chem_spec_data_size = this%num_spec

  end function chem_spec_data_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Check if a species name is in the set of chemical species
  logical function chem_spec_data_exists(this, spec_name) result (found)

    !> Species dataset
    class(chem_spec_data_t), intent(in) :: this
    !> Species name
    character(len=:), allocatable, intent(in) :: spec_name

    found = this%find(spec_name).ne.0

  end function chem_spec_data_exists

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get the index of a chemical species by name. Return 0 if the species is
  !! not found.
  integer(kind=i_kind) function chem_spec_data_find(this, spec_name) &
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

  end function chem_spec_data_find

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get a species property set
  function chem_spec_data_get_property_set(this, spec_name) result (property_set)

    !> Pointer to species properties
    class(property_t), pointer :: property_set
    !> Species dataset
    class(chem_spec_data_t), intent(in) :: this
    !> Species name to find properties of
    character(len=:), allocatable :: spec_name

    integer(i_kind) :: spec_id

    property_set => null()
    spec_id = this%find(spec_name)
    if (spec_id.gt.0) property_set => this%property_set(spec_id)

  end function chem_spec_data_get_property_set

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get a species type
  function chem_spec_data_get_type(this, spec_name) result(spec_type)

    !> Species type
    integer(kind=i_kind) :: spec_type
    !> Species dataset
    class(chem_spec_data_t), intent(in) :: this
    !> Species name to find properties of
    character(len=:), allocatable :: spec_name

    integer(i_kind) :: spec_id

    spec_type = UNKNOWN_SPEC
    spec_id = this%find(spec_name)
    if (spec_id.gt.0) spec_type = this%spec_type(spec_id)

  end function chem_spec_data_get_type

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Return the index of a gas-phase species in the state array 
  function chem_spec_data_state_id(this, spec_name) result (spec_id)

    !> Index of the species in the dataset
    integer(kind=i_kind) :: spec_id
    !> Species dataset
    class(chem_spec_data_t), intent(in) :: this
    !> Species name
    character(len=:), allocatable, intent(in) :: spec_name

    spec_id = this%find(spec_name)

  end function chem_spec_data_state_id

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get an array of absolute integration tolerances corresponding to the
  !! dimensions of the state array
  function chem_spec_data_get_abs_tolerances(this) result (abs_tol)

    !> Array of species absolute integration tolerances
    real(kind=dp), pointer :: abs_tol(:)
    !> Species dataset
    class(chem_spec_data_t), intent(in) :: this

    character(len=:), allocatable :: key
    integer(kind=i_kind) :: i_spec
    real(kind=dp) :: val

    allocate(abs_tol(this%num_spec))
    key = "absolute integration tolerance"
    do i_spec = 1, this%num_spec
      if (this%property_set(i_spec)%get_real(key, val)) then
        abs_tol(i_spec) = val
      else
        abs_tol(i_spec) = DEFAULT_ABS_TOL
      end if
    end do

  end function chem_spec_data_get_abs_tolerances

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_chem_spec_data
