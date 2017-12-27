! Copyright (C) 2017 Matt Dawson
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_aero_phase_data module.

!> The abstract aero_phase_data_t structure and associated subroutines.
module pmc_aero_phase_data

  use pmc_constants,                  only : i_kind, dp
  use pmc_mpi
  use pmc_util,                       only : die_msg, string_t
  use pmc_property
  use pmc_chem_spec_data
  use pmc_phlex_state
#ifdef PMC_USE_MPI
  use mpi
#endif
#ifdef PMC_USE_JSON
  use json_module
#endif

  implicit none
  private

  public :: aero_phase_data_t

  !> Reallocation increment
  integer(kind=i_kind), parameter :: REALLOC_INC = 50

  !> Aerosol phase data type
  !!
  !! Aerosol phase information. The chemistry module uses a set of phases to
  !! model chemical reactions for aerosol species. Chemical reactions can take
  !! place within a phase or across an interface between phases (e.g., for 
  !! condensation/evaporation). Phases may or may not correspond to aerosol
  !! representations directly. For example, a binned aerosol representation
  !! may have one aerosol phase per bin, or an organic and aqueous phase in
  !! each bin, or three concentric aerosol phases (layers).
  type :: aero_phase_data_t
    private
    !> Name of the aerosol phase
    character(len=:), allocatable :: phase_name
    !> Number of species in the phase
    integer(kind=i_kind) :: num_spec = 0
    !> Species names. These are species that are present in the aerosol
    !! phase. These species must exist in the chem_spec_data_t variable
    !! during initialization.
    type(string_t), pointer :: spec_name(:) => null()
    !> Aerosol phase parameters. These will be available during 
    !! initialization, but not during integration. All information required
    !! by functions of the aerosol representation related to a phase must be 
    !! saved by the aero_rep_data_t-exdending type in the condensed data 
    !! arrays.
    type(property_t), pointer :: property_set => null()
    !> Condensed representaiton data. Theses arrays will be available during
    !! integration, and should contain any information required by the
    !! functions of the aerosol representation that cannot be obtained
    !! from the phlex_state_t object. (floating-point)
    real(kind=dp), allocatable :: condensed_data_real(:)
    !> Condensed reaction data (integers)
    integer(kind=i_kind), allocatable ::  condensed_data_int(:)
  contains
    !> Aerosol representation initialization
    procedure :: initialize
    !> Get the name of the aerosol phase
    procedure :: name => get_name
    !> Get property data associated with this phase
    procedure :: get_property_set
    !> Get a list of species names in this phase
    procedure :: get_species
    !> Get an aerosol species state id
    procedure :: state_id
    !> Get the total mass in a phase (ug/m^3)
    procedure :: total_mass
    !> Determine the number of bytes required to pack the given value
    procedure :: pack_size
    !> Packs the given value into the buffer, advancing position
    procedure :: bin_pack
    !> Unpacks the given value from the buffer, advancing position
    procedure :: bin_unpack
    !> Load data from an input file
    procedure :: load
    !> Get the number of species in the phase
    procedure :: size => get_size
    !> Print the aerosol phase data
    procedure :: print => do_print

    !> Private functions
    !> Add a species
    procedure, private :: add
    !> Ensure there is enough room in the species dataset to add a specified
    !! number of species
    procedure, private :: ensure_size
    !> Find a species index by name
    procedure, private :: find
  end type aero_phase_data_t

  !> Constructor for aero_phase_data_t
  interface aero_phase_data_t
    procedure :: constructor
  end interface aero_phase_data_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Constructor for aero_phase_data_t
  function constructor(phase_name, init_size) result(new_obj)

    !> A new set of aerosol-phase species
    type(aero_phase_data_t), pointer :: new_obj
    !> Name of the aerosol phase
    character(len=:), allocatable, intent(in), optional :: phase_name
    !> Number of species to allocate space for initially
    integer(kind=i_kind), intent(in), optional :: init_size

    integer(kind=i_kind) :: alloc_size = REALLOC_INC

    if (present(init_size)) alloc_size = init_size
    allocate(new_obj)
    if (present(phase_name)) new_obj%phase_name = phase_name
    allocate(new_obj%spec_name(alloc_size))

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Initialize the aerosol phase data, validating species names.
  subroutine initialize(this, chem_spec_data)

    !> Aerosol phase data
    class(aero_phase_data_t), intent(inout) :: this
    !> Chemical species data
    type(chem_spec_data_t), intent(in) :: chem_spec_data

    integer(kind=i_kind) :: i_spec

    do i_spec = 1, this%num_spec
      if (.not.chem_spec_data%exists(this%spec_name(i_spec)%string)) then
        call die_msg(589987734, "Aerosol phase species "// &
            trim(this%spec_name(i_spec)%string)//" missing in chem_spec_data.")
      end if
    end do

  end subroutine initialize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get the aerosol phase name
  function get_name(this) result (phase_name)

    !> The name of the aerosol phase
    character(len=:), allocatable :: phase_name
    !> Aerosol phase data
    class(aero_phase_data_t), intent(in) :: this

    phase_name = this%phase_name

  end function get_name

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get the aerosol phase property set
  function get_property_set(this) result (property_set)

    !> A pointer to the aerosol phase property set
    class(property_t), pointer :: property_set
    !> Aerosol phase data
    class(aero_phase_data_t), intent(in) :: this

    property_set => this%property_set

  end function get_property_set

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get the aerosol phase species list
  function get_species(this) result (species)

    !> A list of species in this phase
    type(string_t), allocatable :: species(:)
    !> Aerosol phase data
    class(aero_phase_data_t), intent(in) :: this

    integer(kind=i_kind) :: i_spec

    allocate(species(this%num_spec))
    do i_spec = 1, this%num_spec
      species(i_spec)%string = this%spec_name(i_spec)%string
    end do

  end function get_species

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get the index of a species in the state array for this phase
  integer(kind=i_kind) function state_id(this, spec_name) &
          result (id)

    !> Aerosol phase data
    class(aero_phase_data_t), intent(in) :: this
    !> Chemical species name
    character(len=:), allocatable, intent(in) :: spec_name

    id = this%find(spec_name)

  end function state_id

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Load species from an input file
#ifdef PMC_USE_JSON
  !! j_obj is expected to be a JSON object containing data related to an
  !! aerosol phase. It should be of the form:
  !!
  !! { "pmc-data" : [
  !!   {
  !!     "name" : "my aerosol phase"
  !!     "type" : "AERO_PHASE"
  !!     "species" : [
  !!       "a species",
  !!       "another species",
  !!       ...
  !!     ],
  !!     ...
  !!   },
  !!   ...
  !! ]}
  !!
  !! The key-value pair "name" is required and must contain the unique name
  !! used for this aerosol phase in the mechanism. A single phase may be 
  !! present in several aerosol groups (e.g., an aqueous phase is each bin
  !! of a binned aerosol representation), but the species associated with a
  !! particular phase will be constant throughout the model run.
  !!
  !! The key-value pair "type" is also required and its value must be
  !! AERO_PHASE. A list of species should be included in a key-value pair
  !! named "species" whose value is an array of species names. These names
  !! must correspond to species names loaded into the chem_spec_data_t 
  !! object. All other data is optional and may include any valid JSON value,
  !! including nested objects. Multiple entries with the same aerosol phase 
  !! name will be merged into a single phase, but duplicate property names for
  !! the same phase will cause an error.
  subroutine load(this, json, j_obj)

    !> Aerosol phase data
    class(aero_phase_data_t), intent(inout) :: this
    !> JSON core
    type(json_core), pointer, intent(in) :: json
    !> JSON object
    type(json_value), pointer, intent(in) :: j_obj

    type(json_value), pointer :: child, next, species
    character(kind=json_ck, len=:), allocatable :: key, unicode_str_val
    integer(kind=i_kind) :: var_type

    character(len=:), allocatable :: phase_name, str_val
    type(property_t), pointer :: property_set

    allocate(property_set)
    property_set = property_t()

    next => null()
    call json%get_child(j_obj, child)
    do while (associated(child))
      call json%info(child, name=key, var_type=var_type)
      if (key.eq."name") then
        if (var_type.ne.json_string) call die_msg(429142134, &
                "Received non-string aerosol phase name.")
        call json%get(child, unicode_str_val)
        this%phase_name = unicode_str_val
      else if (key.eq."species") then
        if (var_type.ne.json_array) call die_msg(293312378, &
                "Received non-array list of aerosol phase species: "//&
                to_string(var_type))
        call json%get_child(child, species)
        do while (associated(species))
          call json%info(species, var_type=var_type)
          if (var_type.ne.json_string) call die_msg(669858868, &
                  "Received non-string aerosol phase species name.")
          call json%get(species, unicode_str_val)
          str_val = unicode_str_val
          call this%add(str_val)
          call json%get_next(species, next)
          species => next
        end do
      else if (key.ne."type") then
        call property_set%load(json, child, .false.)
      end if 
      call json%get_next(child, next)
      child => next
    end do
  
    if (associated(this%property_set)) then
      call this%property_set%update(property_set)
    else 
      this%property_set => property_set
    end if
#else
  subroutine load(this)

    !> Aerosol phase data
    class(aero_phase_data_t), intent(in) :: this

    call warn_msg(236665532, "No support for input files.")
#endif
  end subroutine load
   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get the number of species in the phase
  integer(kind=i_kind) function get_size(this) result(num_spec)

    !> Aerosol phase data
    class(aero_phase_data_t), intent(in) :: this

    num_spec = this%num_spec

  end function get_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Print out the aerosol phase data
  subroutine do_print(this)

    !> Aerosol phase data
    class(aero_phase_data_t), intent(in) :: this

    integer(kind=i_kind) :: i_spec

    write(*,*) "Aerosol phase: ", this%phase_name
    write(*,*) "Number of species: ", this%num_spec
    write(*,*) "Species: ["
    do i_spec = 1, this%num_spec 
      write(*,*) this%spec_name(i_spec)%string
    end do
    write(*,*) "]"
    call this%property_set%print()
    write(*,*) "End aerosol phase: ", this%phase_name

  end subroutine do_print

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Add a new chemical species to the phase
  subroutine add(this, spec_name)

    !> Aerosol phase data
    class(aero_phase_data_t), intent(inout) :: this
    !> Name of the species to add
    character(len=:), allocatable, intent(in) :: spec_name

    integer(kind=i_kind) :: i_spec

    i_spec = this%find(spec_name)
    if (i_spec.ne.0) then
      call warn_msg(980242449, "Species "//spec_name//&
              " added more than once to phase "//this%name())
      return
    end if
    call this%ensure_size(1)
    this%num_spec = this%num_spec + 1
    this%spec_name(this%num_spec) = string_t(spec_name)

  end subroutine add

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Ensure there is enough room in the species dataset to add a specified 
  !! number of species
  subroutine ensure_size(this, num_spec)

    !> Aerosol phase data
    class(aero_phase_data_t), intent(inout) :: this
    !> Number of new species to ensure space for
    integer(kind=i_kind), intent(in) :: num_spec

    integer :: new_size
    type(string_t), pointer :: new_name(:)

    if (size(this%spec_name) .ge. this%num_spec + num_spec) return
    new_size = this%num_spec + num_spec + REALLOC_INC
    allocate(new_name(new_size))
    new_name(1:this%num_spec) = this%spec_name(1:this%num_spec)
    deallocate(this%spec_name)
    this%spec_name => new_name

  end subroutine ensure_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get the index of an aerosol-phase species by name. Return 0 if the species
  !! is not found
  integer(kind=i_kind) function find(this, spec_name) &
                  result (spec_id)

    !> Aerosol phase data
    class(aero_phase_data_t), intent(in) :: this
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

  !> Get the total mass in a phase based on the model state (ug/m^3)
  real(kind=dp) function total_mass(this, phlex_state, state_id)

    !> Aerosol phase data
    class(aero_phase_data_t), intent(in) :: this
    !> Current model state
    type(phlex_state_t), intent(in) :: phlex_state
    !> Beginning id in the state array for phase
    integer(kind=i_kind), intent(in) :: state_id

    integer :: i_spec

    total_mass = real(0.0, kind=dp)
    do i_spec = 0, this%num_spec-1
      total_mass = total_mass + phlex_state%state_var(state_id + i_spec)
    end do

  end function total_mass

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Determine the size of a binary required to pack the aerosol 
  !! representation data
  integer(kind=i_kind) function pack_size(this)

    !> Aerosol representation data
    class(aero_phase_data_t), intent(in) :: this
    
    pack_size = pmc_mpi_pack_size_integer(this%num_spec)

  end function pack_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Pack the given value to the buffer, advancing position
  subroutine bin_pack(this, buffer, pos)

    !> Aerosol representation data
    class(aero_phase_data_t), intent(in) :: this
    !> Memory buffer
    character, intent(inout) :: buffer(:)
    !> Current buffer position
    integer, intent(inout) :: pos

#ifdef PMC_USE_MPI
    integer :: prev_position

    prev_position = pos
    call pmc_mpi_pack_integer(buffer, pos, this%num_spec)
    call assert(561436372, &
         pos - prev_position <= this%pack_size())
#endif

  end subroutine bin_pack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Unpack the given value from the buffer, advancing position
  subroutine bin_unpack(this, buffer, pos)

    !> Aerosol representation data
    class(aero_phase_data_t), intent(out) :: this
    !> Memory buffer
    character, intent(inout) :: buffer(:)
    !> Current buffer position
    integer, intent(inout) :: pos

#ifdef PMC_USE_MPI
    integer :: prev_position

    prev_position = pos
    call pmc_mpi_unpack_integer(buffer, pos, this%num_spec)
    call assert(219217030, &
         pos - prev_position <= this%pack_size())
#endif

  end subroutine bin_unpack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_aero_phase_data
