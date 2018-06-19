! Copyright (C) 2017 Matt Dawson
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_aero_rep_data module.

!> \page phlex_aero_rep Phlexible Module for Chemistry: Aerosol Representation (general)
!!
!! An aerosol representation acts as an interface between an aerosol
!! micro-physics model and the \ref phlex_chem "phlex-chem" chemistry module.
!! Types that extend the abstract \c aero_rep_data_t type should be developed
!! for each type of aerosol representation used by an external model (e.g.,
!! binned, modal, single particle).
!!
!! The available aerosol representations are:
!!  - \subpage phlex_aero_rep_single_particle "Single Particle"
!!  - \subpage phlex_aero_rep_modal_binned_mass "Mass-only Binned/Modal"
!!
!! The general input format for an aerosol representation can be found
!! \subpage input_format_aero_rep "here".
!!
!! General instructions for adding a new aerosol representation can be found
!! \subpage phlex_aero_rep_add "here".

!> The abstract aero_rep_data_t structure and associated subroutines.
module pmc_aero_rep_data

#ifdef PMC_USE_JSON
  use json_module
#endif
#ifdef PMC_USE_MPI
  use mpi
#endif
  use pmc_aero_phase_data
  use pmc_chem_spec_data
  use pmc_constants,                  only : i_kind, dp
  use pmc_mpi
  use pmc_phlex_state
  use pmc_property
  use pmc_util,                       only : die_msg, string_t

  use iso_c_binding

  implicit none
  private

  public :: aero_rep_data_t, aero_rep_data_ptr, aero_rep_update_data_t

  !> Abstract aerosol representation data type
  !!
  !! Time-invariant data related to an aerosol representation. Derived types
  !! extending aero_rep_data_t should describe specific types of aerosol
  !! schemes (e.g., binned, modal, particle-resolved).
  !!
  !! See \ref phlex_aero_rep "Aerosol Representations" for details.
  type, abstract :: aero_rep_data_t
    private
    !> Name of the aerosol representation
    character(len=:), allocatable, public :: rep_name
    !> Aerosol phases associated with this aerosol scheme
    !!
    !! See \ref phlex_aero_phase "Aerosol Phases" for details.
    type(aero_phase_data_ptr), allocatable, public :: aero_phase(:)
    !> Aerosol representation parameters. These will be available during 
    !! initialization, but not during solving. All information required
    !! by functions of the aerosol representation  must be saved by the
    !! exdending type in the condensed data arrays.
    type(property_t), pointer, public :: property_set => null()
    !> Condensed representation data. Theses arrays will be available during
    !! solving, and should contain any information required by the
    !! functions of the aerosol representation that cannot be obtained
    !! from the pmc_phlex_state::phlex_state_t object. (floating-point)
    real(kind=dp), allocatable, public :: condensed_data_real(:)
    !> Condensed representation data. Theses arrays will be available during
    !! solving, and should contain any information required by the
    !! functions of the aerosol representation that cannot be obtained
    !! from the pmc_phlex_state::phlex_state_t object. (integer)
    integer(kind=i_kind), allocatable, public ::  condensed_data_int(:)
  contains
    !> Initialize the aerosol representation data, validating component data and
    !! loading any required information from the \c
    !! aero_rep_data_t::property_set. This routine should be called once for
    !! each aerosol representation at the beginning of a model run after all
    !! the input files have been read in. It ensures all data required during
    !! the model run are included in the condensed data arrays.
    procedure(initialize), deferred :: initialize
    !> Get the size of the section of the
    !! \c pmc_phlex_state::phlex_state_t::state_var array required for this
    !! aerosol representation
    procedure(get_size), deferred :: size
    !> Get a list of unique names for each element on the
    !! \c pmc_phlex_state::phlex_state_t::state_var array for this aerosol
    !! representation. The list may be restricted to a particular phase and/or
    !! aerosol species by including the phase_name and spec_name arguments.
    procedure(unique_names), deferred :: unique_names
    !> Get a species id on the \c pmc_phlex_state::phlex_state_t::state_var
    !! array by its unique name. These are unique ids for each element on the
    !! state array for this \ref phlex_aero_rep "aerosol representation" and
    !! are numbered:
    !!
    !!   \f[x_u \in x_f ... (x_f+n-1)\f]
    !!
    !! where \f$x_u\f$ is the id of the element corresponding to the species
    !! with unique name \f$u\f$ on the \c
    !! pmc_phlex_state::phlex_state_t::state_var array, \f$x_f\f$ is the index
    !! of the first element for this aerosol representation on the state array
    !! and \f$n\f$ is the total number of variables on the state array from
    !! this aerosol representation.
    procedure(spec_state_id), deferred :: spec_state_id
    !> Get the non-unique name of a species by its unique name
    procedure(spec_name), deferred :: spec_name
    !> Get the number of instances of an aerosol phase
    procedure(num_phase_instances), deferred :: num_phase_instances
    !> Load data from an input file
    procedure :: load
    !> Get the name of the aerosol representation
    procedure :: name => get_name
    !> Get ids for all instances of a phase in this aerosol representation for
    !! use during solving
    procedure :: phase_ids
    !> Determine the number of bytes required to pack the given value
    procedure :: pack_size
    !> Packs the given value into the buffer, advancing position
    procedure :: bin_pack
    !> Unpacks the given value from the buffer, advancing position
    procedure :: bin_unpack
    !> Print the aerosol representation data
    procedure :: print => do_print
  end type aero_rep_data_t

  !> Pointer to aero_rep_data_t extending types
  type :: aero_rep_data_ptr
    !> Pointer to an aerosol representation
    class(aero_rep_data_t), pointer :: val => null()
  contains
    !> Dereference the pointer
    procedure :: dereference
    !> Finalize the pointer
    final :: ptr_finalize
  end type aero_rep_data_ptr

  !> Update cookie
  type, abstract :: aero_rep_update_data_t
    !> Aerosol representation type
    integer(kind=c_int) :: aero_rep_type
    !> Update data
    type(c_ptr) :: update_data
  contains
    !> Get the aerosol representation type
    procedure :: get_type => aero_rep_update_data_get_type
    !> Get the update data
    procedure :: get_data => aero_rep_update_data_get_data
  end type aero_rep_update_data_t

interface
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Initialize the aerosol representation data, validating component data and
  !! loading any required information from the \c
  !! aero_rep_data_t::property_set. This routine should be called once for
  !! each aerosol representation at the beginning of a model run after all
  !! the input files have been read in. It ensures all data required during
  !! the model run are included in the condensed data arrays.
  subroutine initialize(this, aero_phase_set, spec_state_id)
    use pmc_util,                                     only : i_kind
    use pmc_chem_spec_data
    use pmc_aero_phase_data
    import :: aero_rep_data_t

    !> Aerosol representation data
    class(aero_rep_data_t), intent(inout) :: this
    !> The set of aerosol phases. Note that an aerosol representation may
    !! implement any number of instances of each phase.
    type(aero_phase_data_ptr), pointer, intent(in) :: aero_phase_set(:)
    !> Beginning state id for this aerosol representation in the
    !! \c pmc_phlex_state::phlex_state_t::state_var array
    integer(kind=i_kind), intent(in) :: spec_state_id

  end subroutine initialize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get the size of the section of the
  !! \c pmc_phlex_state::phlex_state_t::state_var array required for this
  !! aerosol representation
  function get_size(this) result (state_size)
    use pmc_util,                                     only : i_kind
    import :: aero_rep_data_t

    !> Size of the state array section
    integer(kind=i_kind) :: state_size
    !> Aerosol representation data
    class(aero_rep_data_t), intent(in) :: this

  end function get_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get a list of unique names for each element on the
  !! \c pmc_phlex_state::phlex_state_t::state_var array for this aerosol
  !! representation.
  function unique_names(this, phase_name, tracer_type, spec_name)
    use pmc_util,                                     only : string_t, i_kind
    import :: aero_rep_data_t

    !> List of unique names
    type(string_t), allocatable :: unique_names(:)
    !> Aerosol representation data
    class(aero_rep_data_t), intent(in) :: this
    !> Aerosol phase name
    character(len=:), allocatable, optional, intent(in) :: phase_name
    !> Tracer type
    integer(kind=i_kind), optional, intent(in) :: tracer_type
    !> Aerosol-phase species name
    character(len=:), allocatable, optional, intent(in) :: spec_name

  end function unique_names

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get a species id on the \c pmc_phlex_state::phlex_state_t::state_var
  !! array by unique name. These are unique ids for each element on the
  !! state array for this \ref phlex_aero_rep "aerosol representation" and
  !! are numbered:
  !!
  !!   \f[x_u \in x_f ... (x_f+n-1)\f]
  !!
  !! where \f$x_u\f$ is the id of the element corresponding to the species
  !! with unique name \f$u\f$ on the \c
  !! pmc_phlex_state::phlex_state_t::state_var array, \f$x_f\f$ is the index
  !! of the first element for this aerosol representation on the state array
  !! and \f$n\f$ is the total number of variables on the state array from
  !! this aerosol representation.
  !!
  !! If the name is not found, the return value is 0.
  function spec_state_id(this, unique_name) result (spec_id)
    use pmc_util,                                     only : i_kind
    import :: aero_rep_data_t

    !> Species state id
    integer(kind=i_kind) :: spec_id
    !> Aerosol representation data
    class(aero_rep_data_t), intent(in) :: this
    !> Unique name
    character(len=:), allocatable, intent(in) :: unique_name

  end function spec_state_id

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get the non-unique name of a chemical species by its unique name
  function spec_name(this, unique_name)
    use pmc_util,                                       only : i_kind
    import :: aero_rep_data_t

    !> Chemical species name
    character(len=:), allocatable :: spec_name 
    !> Aerosol representation data
    class(aero_rep_data_t), intent(in) :: this
    !> Unique name of the species  in this aerosol representation
    character(len=:), allocatable :: unique_name

  end function spec_name

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get the number of instances of a specified aerosol phase
  function num_phase_instances(this, phase_name)
    use pmc_util,                                       only : i_kind
    import :: aero_rep_data_t

    !> Number of instances of the aerosol phase
    integer(kind=i_kind) :: num_phase_instances
    !> Aerosol representation data
    class(aero_rep_data_t), intent(in) :: this
    !> Aerosol phase name
    character(len=:), allocatable, intent(in) :: phase_name

  end function num_phase_instances

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end interface

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> \page input_format_aero_rep Input JSON Object Format: Aerosol Representation (general)
  !!
  !! A \c json object containing information about an \ref phlex_aero_rep
  !! "aerosol representation" has the following format:
  !! \code{.json}
  !! { "pmc-data" : [
  !!   {
  !!     "name" : "my aero rep",
  !!     "type" : "AERO_REP_TYPE",
  !!     "some parameter" : 123.34,
  !!     "some other parameter" : true,
  !!     "nested parameters" : {
  !!       "sub param 1" : 12.43,
  !!       "sub param other" : "some text",
  !!       ...
  !!     },
  !!     ...
  !!   },
  !!   ...
  !! ]}
  !! \endcode
  !! Aerosol representations must have a unique \b name that will be used to
  !! identify the aerosol representation during initialization. The key-value
  !! pair \b type is also required and must correspond to a valid aerosol
  !! representation type. These include:
  !!
  !!   - \subpage phlex_aero_rep_single_particle "AERO_REP_SINGLE_PARTICLE"
  !!   - \subpage phlex_aero_rep_modal_binned_mass
  !!                    "AERO_REP_MODAL_BINNED_MASS" 
  !!
  !! All remaining data are optional and may include any valid \c json value. 
  !! However, extending types will have specific requirements for the
  !! remaining data. 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Load an aerosol representation from an input file
#ifdef PMC_USE_JSON
  subroutine load(this, json, j_obj)

    !> Aerosol representation data
    class(aero_rep_data_t), intent(inout) :: this
    !> JSON core
    type(json_core), pointer, intent(in) :: json
    !> JSON object
    type(json_value), pointer, intent(in) :: j_obj

    type(json_value), pointer :: child, next
    character(kind=json_ck, len=:), allocatable :: key, unicode_str_val
    integer(kind=json_ik) :: var_type
    logical :: found_name

    ! allocate space for the aerosol representation property set
    this%property_set => property_t()

    found_name = .false.

    ! cycle through the aerosol representation properties to find the name and
    ! load the remaining data into the aerosol representation property set
    next => null()
    call json%get_child(j_obj, child)
    do while (associated(child))
      call json%info(child, name=key, var_type=var_type)
      
      ! aerosol representation name
      if (key.eq."name") then
        call assert_msg(196193896, var_type.eq.json_string, &
                "Received non-string value for aerosol rep name")
        call json%get(child, unicode_str_val)
        this%rep_name = unicode_str_val
        found_name = .true.

      ! load remaining data tinto the property set
      else if (key.ne."type") then
        call this%property_set%load(json, child, .false.)
      end if

      call json%get_next(child, next)
      child => next
    end do
    call assert_msg(420903951, found_name, &
            "Received unnamed aerosol representation.")
#else
  subroutine load(this)

    !> Aerosol representation data
    class(aero_rep_data_t), intent(inout) :: this

    call warn_msg(433045149, "No support for input files")
#endif

  end subroutine load

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get the name of the aerosol representation
  function get_name(this)

    !> Aerosol representation name
    character(len=:), allocatable :: get_name
    !> Aerosol representation data
    class(aero_rep_data_t), intent(in) :: this

    get_name = this%rep_name

  end function get_name

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get a set of ids for all instances of a phase in this aerosol
  !! representation for use during solving
  function phase_ids(this, phase_name)

    !> List of phase ids
    integer(kind=i_kind), allocatable :: phase_ids(:)
    !> Aerosol representation data
    class(aero_rep_data_t), intent(in) :: this
    !> Aerosol phase name
    character(len=:), allocatable, intent(in) :: phase_name

    integer(kind=i_kind) :: num_instances, i_instance, i_phase

    num_instances = this%num_phase_instances(phase_name)
    allocate(phase_ids(num_instances))
    i_instance = 1
    do i_phase = 1, size(this%aero_phase)
      if (this%aero_phase(i_phase)%val%name().eq.phase_name) then
        phase_ids(i_instance) = i_phase
        i_instance = i_instance + 1
      end if
    end do

  end function phase_ids

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Determine the size of a binary required to pack the aerosol 
  !! representation data
  integer(kind=i_kind) function pack_size(this)

    !> Aerosol representation data
    class(aero_rep_data_t), intent(in) :: this
    
    pack_size = &
            pmc_mpi_pack_size_real_array(this%condensed_data_real) + &
            pmc_mpi_pack_size_integer_array(this%condensed_data_int) 

  end function pack_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Pack the given value to the buffer, advancing position
  subroutine bin_pack(this, buffer, pos)

    !> Aerosol representation data
    class(aero_rep_data_t), intent(in) :: this
    !> Memory buffer
    character, intent(inout) :: buffer(:)
    !> Current buffer position
    integer, intent(inout) :: pos

#ifdef PMC_USE_MPI
    integer :: prev_position

    prev_position = pos
    call pmc_mpi_pack_real_array(buffer, pos, this%condensed_data_real)
    call pmc_mpi_pack_integer_array(buffer, pos, this%condensed_data_int)
    call assert(257024095, &
         pos - prev_position <= this%pack_size())
#endif

  end subroutine bin_pack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Unpack the given value from the buffer, advancing position
  subroutine bin_unpack(this, buffer, pos)

    !> Aerosol representation data
    class(aero_rep_data_t), intent(out) :: this
    !> Memory buffer
    character, intent(inout) :: buffer(:)
    !> Current buffer position
    integer, intent(inout) :: pos

#ifdef PMC_USE_MPI
    integer :: prev_position

    prev_position = pos
    call pmc_mpi_unpack_real_array(buffer, pos, this%condensed_data_real)
    call pmc_mpi_unpack_integer_array(buffer, pos, this%condensed_data_int)
    call assert(954732699, &
         pos - prev_position <= this%pack_size())
#endif

  end subroutine bin_unpack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Print the aerosol representation data
  subroutine do_print(this, file_unit)

    !> Aerosol representation data
    class(aero_rep_data_t), intent(in) :: this
    !> File unit for output
    integer(kind=i_kind), optional :: file_unit

    integer(kind=i_kind) :: f_unit = 6

    if (present(file_unit)) f_unit = file_unit
    write(f_unit,*) "*** Aerosol Representation: ",trim(this%rep_name)," ***"
    if (associated(this%property_set)) call this%property_set%print(f_unit)
    if (allocated(this%condensed_data_int)) &
      write(f_unit,*) " *** condensed data int: ",this%condensed_data_int(:)
    if (allocated(this%condensed_data_real)) &
      write(f_unit,*) " *** condensed data real: ",this%condensed_data_real(:)

  end subroutine do_print

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Deference a pointer to an aerosol representation
  elemental subroutine dereference(this)

    !> Pointer to an aerosol representation
    class(aero_rep_data_ptr), intent(inout) :: this

    this%val => null()

  end subroutine dereference

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Finalize a pointer to an aerosol representation
  elemental subroutine ptr_finalize(this)

    !> Pointer to an aerosol representation
    type(aero_rep_data_ptr), intent(inout) :: this

    if (associated(this%val)) deallocate(this%val)

  end subroutine ptr_finalize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get the update data aerosol representation type
  function aero_rep_update_data_get_type(this) result (aero_rep_type)

    !> Aerosol representation type
    integer(kind=c_int) :: aero_rep_type
    !> Update data
    class(aero_rep_update_data_t), intent(in) :: this

    aero_rep_type = this%aero_rep_type

  end function aero_rep_update_data_get_type

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get the update data
  function aero_rep_update_data_get_data(this) result (update_data)

    !> Update data ptr
    type(c_ptr) :: update_data
    !> Update data
    class(aero_rep_update_data_t), intent(in) :: this

    update_data = this%update_data

  end function aero_rep_update_data_get_data

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_aero_rep_data
