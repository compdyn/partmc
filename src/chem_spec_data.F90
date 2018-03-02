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

  use iso_c_binding

  implicit none
  private

  public :: chem_spec_data_t 

  !> State variable types (Must match StateVarType c enum defined in phlex_solver.h)
  integer(kind=i_kind), parameter, public :: CHEM_SPEC_UNKNOWN_TYPE = 0
  integer(kind=i_kind), parameter, public :: CHEM_SPEC_VARIABLE = 1
  integer(kind=i_kind), parameter, public :: CHEM_SPEC_CONSTANT = 2
  integer(kind=i_kind), parameter, public :: CHEM_SPEC_PSSA = 3

  !> Species phase
  integer(kind=i_kind), parameter, public :: CHEM_SPEC_UNKNOWN_PHASE = 0
  integer(kind=i_kind), parameter, public :: CHEM_SPEC_GAS_PHASE = 1
  integer(kind=i_kind), parameter, public :: CHEM_SPEC_AERO_PHASE = 2

  !> Reallocation increment
  integer(kind=i_kind), parameter :: REALLOC_INC = 50
  !> Default absolute integration tolerance
  real(kind=dp), parameter :: DEFAULT_ABS_TOL = 1.0e-14

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
    !> Species phase
    integer(kind=i_kind), pointer :: spec_phase(:) => null()
    !> Species property set
    type(property_t), pointer :: property_set(:) => null()
  contains
    !> Load species from an input file
    procedure :: load
    !> Get the number of species with specified conditions
    procedure :: size => get_size
    !> Check if a species name is in the set of chemical species
    procedure :: exists
    !> Get the name of a species
    procedure :: get_name

    procedure :: get_property_set_by_name
    procedure :: get_property_set_by_id
    !> Get a species properties
    generic :: get_property_set => get_property_set_by_name, &
            get_property_set_by_id
    
    procedure :: get_type_by_name
    procedure :: get_type_by_id
    !> Get a species type
    generic :: get_type => get_type_by_name, get_type_by_id
    
    procedure :: get_phase_by_name
    procedure :: get_phase_by_id
    !> Get a species phase
    generic :: get_phase => get_phase_by_name, get_phase_by_id
    
    procedure :: get_abs_tol_by_name
    procedure :: get_abs_tol_by_id
    !> Get the absolute integration tolerance of a species
    generic :: get_abs_tol => get_abs_tol_by_name, get_abs_tol_by_id
    
    !> Get a gas-phase species index in the \c
    !! pmc_phlex_state::phlex_state_t::state_var array.  Note that
    !! aerosol-phase species indices on the \c
    !! pmc_phlex_state::phlex_state_t::state_var array must be accessed from
    !! \c pmc_aero_rep_data::aero_rep_data_t::species_state_id() for a
    !! particular \ref phlex_aero_rep "aerosol representation".
    procedure :: gas_state_id
    !> Print out the species data
    procedure :: print => do_print

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
    allocate(new_obj%spec_phase(alloc_size))
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
    integer(kind=i_kind), pointer :: new_phase(:)
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
    call this%property_set(1:this%num_spec)%move(new_property_set(1:this%num_spec))
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
    integer(kind=i_kind), intent(in) :: spec_type
    !> Species phase
    integer(kind=i_kind), intent(in) :: spec_phase
    !> Property set for new species
    type(property_t), intent(inout), optional :: property_set

    integer(kind=i_kind) :: i_spec

    ! check for a valid phase
    call assert_msg(821623290, spec_phase.ne.CHEM_SPEC_UNKNOWN_PHASE, &
            "Received unknown species type for :"//spec_name)
    
    ! if the species exists, append the new data
    if (this%find(spec_name, i_spec)) then
      call assert_msg(596247182, spec_type.eq.this%spec_type(i_spec), &
          "Type mismatch for species "//spec_name)
      call assert_msg(612991075, spec_phase.eq.this%spec_phase(i_spec), &
          "Phase mismatch for species "//spec_name)
      call this%property_set(i_spec)%update(property_set)

    ! ... otherwise, create a new species
    else
      call this%ensure_size(1)
      this%num_spec = this%num_spec + 1
      this%spec_name(this%num_spec) = string_t(spec_name)
      this%spec_type(this%num_spec) = spec_type
      this%spec_phase(this%num_spec) = spec_phase
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
  !! "mechanism object". (The same name may not be used for a gas species and
  !! an aerosol species.) The key-value pair \b phase is also required and its
  !! value must be \b GAS or \b AEROSOL. The \b type must be \b VARIABLE, 
  !! \b CONSTANT or \b PSSA, however if a \b type is not specified, it will
  !! be assumed to be \b VARIABLE.
  !!
  !! All remaining data are
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
    integer(kind=i_kind) :: spec_type = CHEM_SPEC_UNKNOWN_TYPE
    integer(kind=i_kind) :: spec_phase = CHEM_SPEC_UNKNOWN_PHASE
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
        if (unicode_str_val.eq."VARIABLE") then
          spec_type = CHEM_SPEC_VARIABLE
        else if (unicode_str_val.eq."CONSTANT") then
          spec_type = CHEM_SPEC_CONSTANT
        else if (unicode_str_val.eq."PSSA") then
          spec_type = CHEM_SPEC_PSSA
        else
          str_val = unicode_str_val
          call die_msg(171550163, "Unknown chemical species type: "// &
                  str_val)
        end if
      else if (key.eq."phase") then
        if (var_type.ne.json_string) call die_msg(273466657, &
                "Received non-string species phase")
        call json%get(child, unicode_str_val)
        if (unicode_str_val.eq."GAS") then
          spec_phase = CHEM_SPEC_GAS_PHASE
        else if (unicode_str_val.eq."AEROSOL") then
          spec_phase = CHEM_SPEC_AERO_PHASE
        else
          str_val = unicode_str_val
          call die_msg(603704146, "Unknown chemical species phase: "// &
                  str_val)
        end if
      else
        call property_set%load(json, child, .false.)
      end if
      call json%get_next(child, next)
      child => next
    end do

    ! Set default values for type and phase
    if (spec_type.eq.CHEM_SPEC_UNKNOWN_TYPE) spec_type = CHEM_SPEC_VARIABLE
    if (spec_phase.eq.CHEM_SPEC_UNKNOWN_PHASE) spec_phase = CHEM_SPEC_GAS_PHASE

    ! Add or update the species data
    call this%add(spec_name, spec_type, spec_phase, property_set)
#else
  subroutine load(this)

    !> Species dataset
    class(chem_spec_data_t), intent(inout) :: this

    call warn_msg(627397948, "No support for input files")
#endif
  end subroutine load

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get the number of species with the given properties. If no properties
  !! are specified, return the total number of species. 
  integer(kind=i_kind) function get_size(this, spec_type, spec_phase) &
                  result (num_spec)

    !> Species database
    class(chem_spec_data_t), intent(in) :: this
    !> State variable type for the species
    integer(kind=i_kind), intent(in), optional :: spec_type
    !> Phase of the species
    integer(kind=i_kind), intent(in), optional :: spec_phase

    ! species index
    integer(kind=i_kind) :: i_spec

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
    integer(kind=i_kind) :: i_spec

    found = this%find(spec_name, i_spec)

  end function exists

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get the name of a species. Returns true if the species is found, and 
  !! false otherwise.
  logical function get_name(this, spec_id, spec_name)

    !> Species dataset
    class(chem_spec_data_t), intent(in) :: this
    !> Species id
    integer(kind=i_kind), intent(in) :: spec_id
    !> Species name
    character(len=:), allocatable, intent(out) :: spec_name

    get_name = .false.
    if (spec_id.lt.0.or.spec_id.gt.this%num_spec) return
    spec_name = this%spec_name(spec_id)%string
    get_name = .true.

  end function get_name

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get the index of a chemical species by name. Returns true if the species
  !! is found or false otherwise.
  logical function find(this, spec_name, spec_id)

    !> Species dataset
    class(chem_spec_data_t), intent(in) :: this
    !> Species name
    character(len=:), allocatable, intent(in) :: spec_name
    !> Species id
    integer(kind=i_kind), intent(out) :: spec_id

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

  !> Get a species property set. Returns true if the species is found, or
  !! false otherwise.
  logical function get_property_set_by_name(this, spec_name, property_set) &
                  result (found)

    !> Species dataset
    class(chem_spec_data_t), intent(in) :: this
    !> Species name to find properties of
    character(len=:), allocatable, intent(in) :: spec_name
    !> Pointer to species properties
    type(property_t), pointer, intent(out) :: property_set

    integer(i_kind) :: spec_id

    property_set => null()
    found = this%find(spec_name, spec_id)
    if (found) property_set => this%property_set(spec_id)

  end function get_property_set_by_name

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get a species property set by species id. Returns true if the species is
  !! found, or false otherwise.
  logical function get_property_set_by_id(this, spec_id, property_set) &
                  result (found)

    !> Species dataset
    class(chem_spec_data_t), intent(in) :: this
    !> Species id to find properties of
    integer(kind=i_kind), intent(in) :: spec_id
    !> Pointer to species properties
    type(property_t), pointer, intent(out) :: property_set

    found = spec_id.gt.0 .and. spec_id.le.this%num_spec
    if (found) property_set => this%property_set(spec_id)

  end function get_property_set_by_id

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get a species type by species name. Returns true if the species is found
  !! or false otherwise.
  logical function get_type_by_name(this, spec_name, spec_type) result (found)

    !> Species dataset
    class(chem_spec_data_t), intent(in) :: this
    !> Species name to find properties of
    character(len=:), allocatable, intent(in) :: spec_name
    !> Species type
    integer(kind=i_kind), intent(out) :: spec_type

    integer(i_kind) :: spec_id

    spec_type = CHEM_SPEC_UNKNOWN_TYPE
    found = this%find(spec_name, spec_id)
    if (found) spec_type = this%spec_type(spec_id)

  end function get_type_by_name

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get a species type by id. Returns true if the species is found, or false
  !! otherwise.
  logical function get_type_by_id(this, spec_id, spec_type) result (found)

    !> Species dataset
    class(chem_spec_data_t), intent(in) :: this
    !> Species id to find properties of
    integer(kind=i_kind), intent(in) :: spec_id
    !> Species type
    integer(kind=i_kind), intent(out) :: spec_type

    spec_type = CHEM_SPEC_UNKNOWN_TYPE
    found = spec_id.gt.0 .and. spec_id.le.this%num_spec
    if (found) spec_type = this%spec_type(spec_id)

  end function get_type_by_id

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get a species phase by name. Returns true if the species is found, or
  !! false otherwise.
  logical function get_phase_by_name(this, spec_name, spec_phase) &
                  result (found)

    !> Species dataset
    class(chem_spec_data_t), intent(in) :: this
    !> Species name to find properties of
    character(len=:), allocatable, intent(in) :: spec_name
    !> Species phase
    integer(kind=i_kind), intent(out) :: spec_phase

    integer(i_kind) :: spec_id

    spec_phase = CHEM_SPEC_UNKNOWN_PHASE
    found = this%find(spec_name, spec_id)
    if (found) spec_phase = this%spec_phase(spec_id)

  end function get_phase_by_name

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get a species phase by id. Returns true if the species is found, or 
  !! false otherwise
  logical function get_phase_by_id(this, spec_id, spec_phase) result (found)

    !> Species dataset
    class(chem_spec_data_t), intent(in) :: this
    !> Species id to find properties of
    integer(kind=i_kind), intent(in) :: spec_id
    !> Species phase
    integer(kind=i_kind), intent(out) :: spec_phase

    spec_phase = CHEM_SPEC_UNKNOWN_PHASE
    found = spec_id.gt.0 .and. spec_id.le.this%num_spec
    if (found) spec_phase = this%spec_phase(spec_id)

  end function get_phase_by_id

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get the absolute integration tolerance of a species by name. Returns 
  !! true if the species is found, or false otherwise.
  logical function get_abs_tol_by_name(this, spec_name, abs_tol) &
                  result (found)

    !> Species dataset
    class(chem_spec_data_t), intent(in) :: this
    !> Species name
    character(len=:), allocatable :: spec_name
    !> Absolute integration tolerance
    real(kind=dp), intent(out) :: abs_tol

    character(len=:), allocatable :: key
    integer(kind=i_kind) :: spec_id
    real(kind=dp) :: val

    abs_tol = DEFAULT_ABS_TOL
    key = "absolute integration tolerance"
    found = this%find(spec_name, spec_id)
    if (found) then
      if (this%property_set(spec_id)%get_real(key, val)) then
        abs_tol = val
      end if
    end if

  end function get_abs_tol_by_name

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get the absolute integration tolerance of a species by id, Returns true
  !! if the species is found, or false otherwise.
  logical function get_abs_tol_by_id(this, spec_id, abs_tol) result (found)

    !> Species dataset
    class(chem_spec_data_t), intent(in) :: this
    !> Species name
    integer(kind=i_kind), intent(in) :: spec_id
    !> Absolute integration tolerance
    real(kind=dp), intent(out) :: abs_tol

    character(len=:), allocatable :: key
    real(kind=dp) :: val

    abs_tol = DEFAULT_ABS_TOL
    key = "absolute integration tolerance"
    found = spec_id.gt.0 .and. spec_id.le.this%num_spec
    if (found) then
      if (this%property_set(spec_id)%get_real(key, val)) then
        abs_tol = val
      end if
    end if

  end function get_abs_tol_by_id

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get a gas-phase species index in the \c
  !! pmc_phlex_state::phlex_state_t::state_var array.  Note that
  !! aerosol-phase species indices on the \c
  !! pmc_phlex_state::phlex_state_t::state_var array must be accessed from
  !! \c pmc_aero_rep_data::aero_rep_data_t::species_state_id() for a
  !! particular \ref phlex_aero_rep "aerosol representation". Returns a valid
  !! state array index if the species is found, or 0 otherwise
  integer(kind=i_kind) function gas_state_id(this, spec_name)

    !> Species dataset
    class(chem_spec_data_t), intent(in) :: this
    !> Species name
    character(len=:), allocatable, intent(in) :: spec_name

    integer(kind=i_kind) :: i_spec
    
    gas_state_id = 0
    do i_spec = 1, this%num_spec
      if (this%spec_phase(i_spec).eq.CHEM_SPEC_GAS_PHASE) &
              gas_state_id = gas_state_id + 1
      if (trim(spec_name).eq.trim(this%spec_name(i_spec)%string)) return
    end do
    gas_state_id = 0

  end function gas_state_id

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Print out the species data
  subroutine do_print(this, file_unit)

    !> Chemical species data
    class(chem_spec_data_t), intent(in) :: this
    !> File unit for output
    integer(kind=i_kind), optional, intent(in) :: file_unit

    integer(kind=i_kind) :: i_spec
    integer(kind=i_kind) :: f_unit = 6

    if (present(file_unit)) f_unit = file_unit

    write(f_unit,*) "Number of species: ", this%num_spec
    do i_spec = 1, this%num_spec
      write(f_unit,*) "  ", this%spec_name(i_spec)%string
      call this%property_set(i_spec)%print(f_unit)
    end do

  end subroutine do_print

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_chem_spec_data
