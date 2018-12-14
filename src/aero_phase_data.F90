! Copyright (C) 2017-2018 Matt Dawson
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_aero_phase_data module.

!> \page phlex_aero_phase Phlexible Module for Chemistry: Aerosol Phase
!!
!! An \c pmc_aero_phase_data::aero_phase_data_t object describes a distinct
!! chemical phase within an aerosol. It is designed to allow the
!! implementation of the chemical and mass transfer processes to be
!! independent of the particular \ref phlex_aero_rep "aerosol representation"
!! used (e.g., bins, modes, single particles).
!!
!! A single \ref phlex_aero_phase "aerosol phase" may be present in several
!! \ref phlex_aero_rep "aerosol representations" (e.g., an aqueous phase in a
!! binned and a single-particle representation), but the \ref phlex_species
!! "chemical species" associated with a particular phase are constant
!! throughout the model run. Once loaded, \ref phlex_aero_phase
!! "aerosol phases" are made available to any \ref input_format_aero_rep
!! "aerosol representations" that want to implement them.
!! \ref phlex_aero_rep "Aerosol representations" are able to specify which
!! phases they implement and how many instances of that phase are present in
!! the \ref phlex_aero_rep "aerosol representation". For example, a binned
!! representation with 10 bins may implement 10 aqueous phases and 10 organic
!! phases, whereas a single  particle representation with a concentric shell
!! structure of 3 layers may implement 3 of each phase (assuming the chemistry
!! is solved for each particle individually).
!!
!! The set of \ref phlex_aero_phase "aerosol phases" is made available to the
!! \ref phlex_mechanism "mechanism(s)" during model intialization. Reactions
!! in the chemical mechanism are able to specify, by name, the phase
!! in which they take place, and which species in that phase are involved.
!! (How they decide this is up to the particular \ref  phlex_rxn
!! "reaction type".)  Any physical aerosol parameters, such as the surface
!! area between phases, the particle radius, or the number concentration,
!! required by a chemical reaction will be provided by the \ref phlex_aero_rep
!! "aerosol representation" at run time.
!!
!! The input format for an aerosol phase can be found \ref
!! input_format_aero_phase "here".

!> The abstract aero_phase_data_t structure and associated subroutines.
module pmc_aero_phase_data

#ifdef PMC_USE_JSON
  use json_module
#endif
#ifdef PMC_USE_MPI
  use mpi
#endif
  use pmc_chem_spec_data
  use pmc_constants,                  only : i_kind, dp
  use pmc_mpi
  use pmc_phlex_state
  use pmc_property
  use pmc_util,                       only : die_msg, string_t

  implicit none
  private

#define NUM_STATE_VAR_ this%condensed_data_int(1)
#define NUM_INT_PROP_ 1
#define NUM_REAL_PROP_ 0
#define SPEC_TYPE_(x) this%condensed_data_int(NUM_INT_PROP_+x)
#define MW_(x) this%condensed_data_real(NUM_REAL_PROP_+x)
#define DENSITY_(x) this%condensed_data_real(NUM_REAL_PROP_+NUM_STATE_VAR_+x)

  public :: aero_phase_data_t, aero_phase_data_ptr

  !> Reallocation increment
  integer(kind=i_kind), parameter :: REALLOC_INC = 50

  !> Aerosol phase data type
  !!
  !! \ref phlex_aero_phase "Aerosol phase" information.
  type :: aero_phase_data_t
    private
    !> Name of the aerosol phase
    character(len=:), allocatable :: phase_name
    !> Number of species in the phase
    integer(kind=i_kind) :: num_spec = 0
    !> Species names. These are species that are present in the aerosol
    !! phase. These species must exist in the \c
    !! pmc_phlex_core::phlex_core_t::chem_spec_data variable during
    !! initialization.
    type(string_t), pointer :: spec_name(:) => null()
    !> Aerosol phase parameters. These will be available during
    !! initialization, but not during solving.
    type(property_t), pointer :: property_set => null()
    !> Condensed phase data. Theses arrays will be available during
    !! solving, and should contain any information required by the
    !! functions of the aerosol phase that cannot be obtained
    !! from the \c pmc_phlex_state::phlex_state_t object. (floating-point)
    real(kind=dp), allocatable, public :: condensed_data_real(:)
    !> Condensed phase data. Theses arrays will be available during
    !! solving, and should contain any information required by the
    !! functions of the aerosol phase that cannot be obtained
    !! from the \c pmc_phlex_state::phlex_state_t object. (integer)
    integer(kind=i_kind), allocatable, public ::  condensed_data_int(:)
    !> Pointer to the set of chemical species data
    type(chem_spec_data_t), pointer :: chem_spec_data
  contains
    !> Load data from an input file
    procedure :: load
    !> Aerosol phase initialization
    procedure :: initialize
    !> Get the name of the aerosol phase
    procedure :: name => get_name
    !> Get the number of species in the phase
    procedure :: size => get_size
    !> Get property data associated with this phase
    procedure :: get_property_set
    !> Get a list of species names in this phase
    procedure :: get_species_names
    !> Get a species type by name
    procedure :: get_species_type
    !> Determine the number of bytes required to pack the given value
    procedure :: pack_size
    !> Packs the given value into the buffer, advancing position
    procedure :: bin_pack
    !> Unpacks the given value from the buffer, advancing position
    procedure :: bin_unpack
    !> Print the aerosol phase data
    procedure :: print => do_print
    !> Finalize the aerosol phase data
    final :: finalize

    ! Private functions
    !> Ensure there is enough room in the species dataset to add a specified
    !! number of species
    procedure, private :: ensure_size
    !> Add a species
    procedure, private :: add
    !> Find a species index by name
    procedure, private :: find
  end type aero_phase_data_t

  ! Constructor for aero_phase_data_t
  interface aero_phase_data_t
    procedure :: constructor
  end interface aero_phase_data_t

  !> Pointer type for building arrays
  type aero_phase_data_ptr
    type(aero_phase_data_t), pointer :: val => null()
  contains
    !> Dereference the pointer
    procedure :: dereference
    !> Finalize the pointer
    final :: ptr_finalize
  end type aero_phase_data_ptr

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
    if (present(phase_name)) then
      new_obj%phase_name = phase_name
    else
      new_obj%phase_name = ""
    endif
    allocate(new_obj%spec_name(alloc_size))

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> \page input_format_aero_phase Input JSON Object Format: Aerosol Phase
  !!
  !! A \c json object containing information about an \ref phlex_aero_phase
  !! "aerosol phase" has the following format:
  !! \code{.json}
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
  !! \endcode
  !! The key-value pair \b name is required and must contain the unique name
  !! used for this \ref phlex_aero_phase "aerosol phase" in the \ref
  !! input_format_mechanism "mechanism". The key-value pair \b type is also
  !! required and its value must be \b AERO_PHASE.
  !!
  !! A list of species names should be included in a key-value pair named
  !! \b species whose value is an array of species names. These names must
  !! correspond to \ref input_format_species "chemcical species" names.
  !! \ref input_format_species "Chemical species" included in the \b species
  !! array must have a \b phase of \b AEROSOL and must include key value pairs
  !! \b molecular \b weight \b [kg \b mol-1]
  !! (\f$\mbox{\si{\kilogram\per\mole}}\f$) and \b density \b [kg \m-3]
  !! (\f$\mbox{\si{\kilogram\per\cubic\metre}}\f$).
  !!
  !! All other data is optional and may include any valid \c json value.
  !! Multiple entries with the same aerosol phase \b name will be merged into
  !! a single phase, but duplicate property names for the same phase will
  !! cause an error.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Load species from an input file
#ifdef PMC_USE_JSON
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

    character(len=:), allocatable :: str_val
    type(property_t), pointer :: property_set

    ! allocate space for the phase property set
    property_set => property_t()

    ! cycle through the phase properties to find the name and species
    ! and load the remaining data into the phase property set
    next => null()
    call json%get_child(j_obj, child)
    do while (associated(child))
      call json%info(child, name=key, var_type=var_type)

      ! phase name
      if (key.eq."name") then
        if (var_type.ne.json_string) call die_msg(429142134, &
                "Received non-string aerosol phase name.")
        call json%get(child, unicode_str_val)
        this%phase_name = unicode_str_val

      ! chemical species in the phase
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

      ! load remaining properties into the phase property set
      else if (key.ne."type") then
        call property_set%load(json, child, .false., this%phase_name)
      end if

      call json%get_next(child, next)
      child => next
    end do

    ! save the property set
    if (associated(this%property_set)) then
      call this%property_set%update(property_set, this%phase_name)
      deallocate (property_set)
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

  !> Initialize the aerosol phase data, validating species names.
  subroutine initialize(this, chem_spec_data)

    !> Aerosol phase data
    class(aero_phase_data_t), intent(inout) :: this
    !> Chemical species data
    type(chem_spec_data_t), target, intent(in) :: chem_spec_data

    type(property_t), pointer :: spec_props
    integer(kind=i_kind) :: i_spec, i_spec_phase_type
    character(len=:), allocatable :: key_name

    ! Allocate space for the condensed data arrays
    allocate(this%condensed_data_int(NUM_INT_PROP_+this%num_spec))
    allocate(this%condensed_data_real(NUM_REAL_PROP_+2*this%num_spec))

    ! Set the number of species
    NUM_STATE_VAR_ = this%num_spec

    ! Find the aerosol-phase species, and save their needed properties
    do i_spec = 1, NUM_STATE_VAR_

      ! Get the species properties
      call assert_msg(140971956, chem_spec_data%get_property_set( &
              this%spec_name(i_spec)%string, spec_props), &
              "Missing property set for species '"// &
              this%spec_name(i_spec)%string// &
              "' in aerosol phase '"//this%phase_name//"'")

      ! Get the species type
      call assert_msg(129442398, chem_spec_data%get_type( &
              this%spec_name(i_spec)%string, SPEC_TYPE_(i_spec)), &
              "Missing type for species '"// &
              this%spec_name(i_spec)%string// &
              "' in aerosol phase '"//this%phase_name//"'")

      ! Make sure the species is an aerosol-phase species
      call assert_msg(136357145, chem_spec_data%get_phase( &
              this%spec_name(i_spec)%string, i_spec_phase_type ), &
              "Error getting phase for species "// &
              this%spec_name(i_spec)%string)
      call assert_msg(861388228, i_spec_phase_type.eq.CHEM_SPEC_AERO_PHASE, &
              "Trying to add non-aerosol phase species to aerosol phase "// &
              this%phase_name//"; species: "//this%spec_name(i_spec)%string)

      ! get the molecular weight and density of species
      ! present in the phase
      if (SPEC_TYPE_(i_spec).eq.CHEM_SPEC_VARIABLE .or. &
          SPEC_TYPE_(i_spec).eq.CHEM_SPEC_CONSTANT .or. &
          SPEC_TYPE_(i_spec).eq.CHEM_SPEC_PSSA) then

        ! Get the molecular weight
        key_name = "molecular weight [kg mol-1]"
        call assert_msg(512254139, &
                spec_props%get_real(key_name, MW_(i_spec)), &
                "Missing molecular weight for species '"// &
                this%spec_name(i_spec)%string// &
                "' in aerosol phase '"//this%phase_name//"'")

        ! Get the density
        key_name = "density [kg m-3]"
        call assert_msg(224966878, &
                spec_props%get_real(key_name, DENSITY_(i_spec)), &
                "Missing density for species '"// &
                this%spec_name(i_spec)%string// &
                "' in aerosol phase '"//this%phase_name//"'")

      ! activity coefficients do not need molecular weight or density
      else

        MW_(i_spec) = 0.0
        DENSITY_(i_spec) = 0.0

      end if

    end do

    ! Save a pointer to the chemical species data
    this%chem_spec_data => chem_spec_data

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

  !> Get the number of species in the phase
  integer(kind=i_kind) function get_size(this) result(num_spec)

    !> Aerosol phase data
    class(aero_phase_data_t), intent(in) :: this

    num_spec = this%num_spec

  end function get_size

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

  !> Get an aerosol phase species name
  function get_species_names(this) result (spec_names)

    !> Names of species in this phase
    type(string_t), allocatable :: spec_names(:)
    !> Aerosol phase data
    class(aero_phase_data_t), intent(in) :: this

    integer(kind=i_kind) :: i_spec

    allocate(spec_names(this%num_spec))
    do i_spec = 1, this%num_spec
      spec_names(i_spec)%string = this%spec_name(i_spec)%string
    end do

  end function get_species_names

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get an aerosol phase species type
  function get_species_type(this, spec_name) result (spec_type)

    !> The type of a species in this phase
    integer(kind=i_kind) :: spec_type
    !> Aerosol phase data
    class(aero_phase_data_t), intent(in) :: this
    !> Name of the species
    character(len=:), allocatable, intent(in) :: spec_name

    call assert_msg(163269315, this%find(spec_name).gt.0, &
            "Species '"//spec_name//"' is not in aerosol phase '"// &
            this%phase_name//"'.")
    call assert(255786656, this%chem_spec_data%get_type(spec_name, spec_type))

  end function get_species_type

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Determine the size of a binary required to pack the aerosol
  !! representation data
  integer(kind=i_kind) function pack_size(this, comm)

    !> Aerosol representation data
    class(aero_phase_data_t), intent(in) :: this
    !> MPI communicator
    integer, intent(in) :: comm

    pack_size = &
            pmc_mpi_pack_size_real_array(this%condensed_data_real, comm) + &
            pmc_mpi_pack_size_integer_array(this%condensed_data_int, comm)

  end function pack_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Pack the given value to the buffer, advancing position
  subroutine bin_pack(this, buffer, pos, comm)

    !> Aerosol representation data
    class(aero_phase_data_t), intent(in) :: this
    !> Memory buffer
    character, intent(inout) :: buffer(:)
    !> Current buffer position
    integer, intent(inout) :: pos
    !> MPI communicator
    integer, intent(in) :: comm

#ifdef PMC_USE_MPI
    integer :: prev_position

    prev_position = pos
    call pmc_mpi_pack_real_array(buffer, pos, this%condensed_data_real, comm)
    call pmc_mpi_pack_integer_array(buffer, pos, this%condensed_data_int,comm)
    call assert(561436372, &
         pos - prev_position <= this%pack_size(comm))
#endif

  end subroutine bin_pack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Unpack the given value from the buffer, advancing position
  subroutine bin_unpack(this, buffer, pos, comm)

    !> Aerosol representation data
    class(aero_phase_data_t), intent(out) :: this
    !> Memory buffer
    character, intent(inout) :: buffer(:)
    !> Current buffer position
    integer, intent(inout) :: pos
    !> MPI communicator
    integer, intent(in) :: comm

#ifdef PMC_USE_MPI
    integer :: prev_position

    prev_position = pos
    call pmc_mpi_unpack_real_array(buffer, pos, this%condensed_data_real,comm)
    call pmc_mpi_unpack_integer_array(buffer, pos, this%condensed_data_int,  &
                                                                         comm)
    call assert(219217030, &
         pos - prev_position <= this%pack_size(comm))
#endif

  end subroutine bin_unpack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Print out the aerosol phase data
  subroutine do_print(this, file_unit)

    !> Aerosol phase data
    class(aero_phase_data_t), intent(in) :: this
    !> File unit for output
    integer(kind=i_kind), optional :: file_unit

    integer(kind=i_kind) :: f_unit = 6
    integer(kind=i_kind) :: i_spec

    if (present(file_unit)) f_unit = file_unit
    write(f_unit,*) "Aerosol phase: ", this%phase_name
    write(f_unit,*) "Number of species: ", this%num_spec
    write(f_unit,*) "Species: ["
    do i_spec = 1, this%num_spec
      write(*,*) this%spec_name(i_spec)%string
    end do
    write(f_unit,*) "]"
    call this%property_set%print(f_unit)
    write(f_unit,*) "End aerosol phase: ", this%phase_name

  end subroutine do_print

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Finalize the aerosol phase data
  elemental subroutine finalize(this)

    !> Aerosol phase data
    type(aero_phase_data_t), intent(inout) :: this

    if (allocated(this%phase_name))    deallocate(this%phase_name)
    if (associated(this%spec_name))    deallocate(this%spec_name)
    if (associated(this%property_set)) deallocate(this%property_set)
    if (allocated(this%condensed_data_real)) &
                                       deallocate(this%condensed_data_real)
    if (allocated(this%condensed_data_int)) &
                                       deallocate(this%condensed_data_int)

  end subroutine finalize

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
    this%spec_name(this%num_spec)%string = spec_name

  end subroutine add

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get the index of an aerosol-phase species by name. Return 0 if the
  !! species is not found
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

  !> Dereference a pointer to aerosol phase data
  elemental subroutine dereference(this)

    !> Pointer to aerosol phase data
    class(aero_phase_data_ptr), intent(inout) :: this

    this%val => null()

  end subroutine dereference

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Finalize a pointer to aerosol phase data
  elemental subroutine ptr_finalize(this)

    !> Pointer to aerosol phase data
    type(aero_phase_data_ptr), intent(inout) :: this

    if (associated(this%val)) deallocate(this%val)

  end subroutine ptr_finalize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#undef NUM_STATE_VAR_
#undef NUM_INT_PROP_
#undef NUM_REAL_PROP_
#undef SPEC_TYPE_
#undef MW_
#undef DENSITY_

end module pmc_aero_phase_data
