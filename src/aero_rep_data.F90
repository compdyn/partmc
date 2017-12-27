! Copyright (C) 2017 Matt Dawson
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_aero_rep_data module.

!> The abstract aero_rep_data_t structure and associated subroutines.
module pmc_aero_rep_data

  use pmc_constants,                  only : i_kind, dp
  use pmc_mpi
  use pmc_util,                       only : die_msg, string_t
  use pmc_property
  use pmc_chem_spec_data
  use pmc_phlex_state
  use pmc_aero_phase_data
#ifdef PMC_USE_MPI
  use mpi
#endif
#ifdef PMC_USE_JSON
  use json_module
#endif

  implicit none
  private

  public :: aero_rep_data_t, aero_rep_data_ptr

  !> Abstract aerosol representation data type
  !!
  !! Time-invariant data related to an aerosol representation. Derived types
  !! extending aero_rep_data_t should describe specific types of aerosol
  !! schemes (e.g., binned, modal, particle-resolved).
  type, abstract :: aero_rep_data_t
    private
    !> Name of the aerosol representation
    character(len=:), allocatable :: rep_name
    !> Aerosol phases associated with this aerosol scheme
    !! FIXME This should be a private component, but extending types need
    !! access to it.
    class(aero_phase_data_ptr), allocatable, public :: aero_phase(:)
    !> Aerosol representation parameters. These will be available during 
    !! initialization, but not during integration. All information required
    !! by functions of the aerosol representation  must be saved by the
    !! exdending type in the condensed data arrays.
    type(property_t), pointer :: property_set => null()
    !> Condensed representaiton data. Theses arrays will be available during
    !! integration, and should contain any information required by the
    !! functions of the aerosol representation that cannot be obtained
    !! from the phlex_state_t object. (floating-point)
    !! FIXME This should be a private component, but extending types need
    !! access to it.
    real(kind=dp), allocatable, public :: condensed_data_real(:)
    !> Condensed reaction data (integers)
    !! FIXME This should be a private component, but extending types need
    !! access to it.
    integer(kind=i_kind), allocatable, public ::  condensed_data_int(:)
  contains
    !> Aerosol representation initialization
    procedure(pmc_aero_rep_data_initialize), deferred :: initialize
    !> Get the size of the state variable array required for this aerosol
    !! representation
    procedure(pmc_aero_rep_data_size), deferred :: size
    !> Get a list of unique species names for each element of the state
    !! variable array
    procedure(pmc_aero_rep_data_unique_names), deferred :: unique_names
    !> Get a species state id by unique name
    procedure(pmc_aero_rep_data_state_id_by_unique_name), deferred :: &
            state_id_by_unique_name
    !> Get an instance of the state variable for this aerosol representation
    procedure(pmc_aero_rep_data_new_state), deferred :: new_state
    !> Get aerosol species state id
    procedure(pmc_aero_rep_data_species_state_id), deferred :: &
            species_state_id
    !> Get surface area concentration (m^2/m^3)
    procedure(pmc_aero_rep_data_surface_area_conc), deferred :: &
            surface_area_conc
    !> Get the surface area concentration for a specific aerosol species
    !! (m^2/m^3)
    procedure(pmc_aero_rep_data_species_surface_area_conc), deferred :: &
            species_surface_area_conc
    !> Get the vapor pressure scaling for a particular species (unitless)
    procedure(pmc_aero_rep_data_vapor_pressure_scaling), deferred :: &
            vapor_pressure_scaling
    !> Get the name of the aerosol representation
    procedure :: name => pmc_aero_rep_data_name
    !> Get a phase id in this aerosol representation
    procedure :: phase_id => pmc_aero_rep_data_phase_id
    !> Determine the number of bytes required to pack the given value
    procedure :: pack_size => pmc_aero_rep_data_pack_size
    !> Packs the given value into the buffer, advancing position
    procedure :: bin_pack => pmc_aero_rep_data_bin_pack
    !> Unpacks the given value from the buffer, advancing position
    procedure :: bin_unpack => pmc_aero_rep_data_bin_unpack
    !> Load data from an input file
    procedure :: load => pmc_aero_rep_data_load
    !> Print the aerosol representation data
    procedure :: print => pmc_aero_rep_data_print
  end type aero_rep_data_t

  !> Pointer to aero_rep_data_t extending types
  type :: aero_rep_data_ptr
    !> Pointer to an aerosol representation
    class(aero_rep_data_t), pointer :: val
  end type aero_rep_data_ptr

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Initialize the aerosol representation data, validating component data and
  !! loading any required information from the  property_t object. This 
  !! routine should be called once for each aerosol representation
  !! at the beginning of a model run after all the input files have been
  !! read in. It ensures all data required during the model run are included
  !! in the condensed data arrays.
  interface pmc_aero_rep_data_initialize_if
    subroutine pmc_aero_rep_data_initialize(this, aero_phase_set, &
                  spec_state_id, aero_state_id, chem_spec_data)
      use pmc_util,                                     only : i_kind
      use pmc_chem_spec_data
      use pmc_aero_phase_data
      import :: aero_rep_data_t

      !> Aerosol representation data
      class(aero_rep_data_t), intent(inout) :: this
      !> The set of aerosol phases
      type(aero_phase_data_t), pointer, intent(in) :: aero_phase_set(:)
      !> Beginning state id for this aerosol representationin the model species
      !! state array
      integer(kind=i_kind), intent(in) :: spec_state_id
      !> Index for this representation in the model state aero_rep_state_t
      !! array
      integer(kind=i_kind), intent(in) :: aero_state_id
      !> Chemical species data
      type(chem_spec_data_t), intent(in) :: chem_spec_data

    end subroutine pmc_aero_rep_data_initialize
  end interface pmc_aero_rep_data_initialize_if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get the size of the state variable array required for this aerosol
  !! representation
  interface pmc_aero_rep_data_size_if
    function pmc_aero_rep_data_size(this) result (state_size)
      use pmc_util,                                     only : i_kind
      import :: aero_rep_data_t

      !> Size of the state array
      integer(kind=i_kind) :: state_size
      !> Aerosol representation data
      class(aero_rep_data_t), intent(in) :: this

    end function pmc_aero_rep_data_size
  end interface pmc_aero_rep_data_size_if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get a list of unique names for each element of the state variable array.
  interface pmc_aero_rep_data_unique_names_if
    function pmc_aero_rep_data_unique_names(this) result (unique_names)
      use pmc_util,                                     only : string_t
      import :: aero_rep_data_t

      !> List of unique names
      type(string_t), allocatable :: unique_names(:)
      !> Aerosol representation data
      class(aero_rep_data_t), intent(in) :: this

    end function pmc_aero_rep_data_unique_names
  end interface pmc_aero_rep_data_unique_names_if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get the species state id by the unique name
  interface pmc_aero_rep_data_state_id_by_unique_name_if
    function pmc_aero_rep_data_state_id_by_unique_name(this, &
                    unique_name) result (spec_id)
      use pmc_util,                                     only : i_kind
      import :: aero_rep_data_t

      !> Species state id
      integer(kind=i_kind) :: spec_id
      !> Aerosol representation data
      class(aero_rep_data_t), intent(in) :: this
      !> Unique name
      character(len=:), allocatable, intent(in) :: unique_name

    end function pmc_aero_rep_data_state_id_by_unique_name
  end interface pmc_aero_rep_data_state_id_by_unique_name_if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get an instance of the state type for this aerosol representation
  interface pmc_aero_rep_data_new_state_if
    function pmc_aero_rep_data_new_state(this) result (aero_rep_state)
      use pmc_aero_rep_state
      import :: aero_rep_data_t

      !> Aerosol representation state
      class(aero_rep_state_t), pointer :: aero_rep_state
      !> Aerosol representation data
      class(aero_rep_data_t), intent(in) :: this

    end function pmc_aero_rep_data_new_state
  end interface pmc_aero_rep_data_new_state_if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get an aerosol species state id(s). The returned array will have the
  !! index of the species in the specified phase in each aerosol group (e.g.,
  !! bin, mode). If the species is not present in a certain group-phase the
  !! index will be 0.
  !! 
  !! This function should only be called during initialization
  interface pmc_aero_rep_data_species_state_id_if
    function pmc_aero_rep_data_species_state_id(this, phase, species_name) &
                    result(spec_index)
      use pmc_util,                                     only : i_kind
      import :: aero_rep_data_t

      !> Species index array
      integer(kind=i_kind), allocatable :: spec_index(:)
      !> Aerosol representation data
      class(aero_rep_data_t), intent(in) :: this
      !> Aerosol phase
      character(len=:), allocatable, intent(in) :: phase
      !> Species name
      character(len=:), allocatable, intent(in) :: species_name

    end function pmc_aero_rep_data_species_state_id
  end interface pmc_aero_rep_data_species_state_id_if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get surface area concentration (m^2/m^3) between two phases. One phase
  !! may be 'gas' to indicate the gas-phase, or both phases may be aerosol
  !! phases.
  interface pmc_aero_rep_data_surface_area_conc_if 
    function pmc_aero_rep_data_surface_area_conc(this, i_phase1, i_phase2, &
                    phlex_state, jac_contrib) result(surface_area_conc)
      use pmc_util,                                     only : dp, i_kind
      use pmc_phlex_state
      import :: aero_rep_data_t

      !> Surface area concentration
      real(kind=dp) :: surface_area_conc
      !> Aerosol representation data
      class(aero_rep_data_t), intent(in) :: this
      !> Aerosol phase1
      integer(kind=i_kind), intent(in) :: i_phase1
      !> Aerosol phase2
      integer(kind=i_kind), intent(in) :: i_phase2
      !> Model state
      type(phlex_state_t), intent(in) :: phlex_state
      !> Contribution to Jacobian matrix. An array of the same size as the
      !! state array that, when present, will be filled with the partial
      !! derivatives of the result of this calculation with each state
      !! variable.
      real(kind=dp), allocatable, intent(out), optional :: jac_contrib(:)

    end function pmc_aero_rep_data_surface_area_conc
  end interface pmc_aero_rep_data_surface_area_conc_if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get the surface area concentration for a specific aerosol species
  !! (m^2/m^3). It is assumed the surface is between the gas-phase and an
  !! aerosol phase.
  interface pmc_aero_rep_data_species_surface_area_conc_if
    function pmc_aero_rep_data_species_surface_area_conc(this, i_phase, &
                    i_spec, phlex_state, jac_contrib) result(surface_area_conc)
      use pmc_util,                                     only : dp, i_kind
      use pmc_phlex_state
      import :: aero_rep_data_t

      !> Surface area concentration
      real(kind=dp) :: surface_area_conc
      !> Aerosol representation data
      class(aero_rep_data_t), intent(in) :: this
      !> Aerosol phase id
      integer(kind=i_kind), intent(in) :: i_phase
      !> Species id
      integer(kind=i_kind), intent(in) :: i_spec
      !> Model state
      type(phlex_state_t), intent(in) :: phlex_state
      !> Contribution to Jacobian matrix. An array of the same size as the
      !! state array that, when present, will be filled with the partial
      !! derivatives of the result of this calculation with each state
      !! variable.
      real(kind=dp), allocatable, intent(out), optional :: jac_contrib(:)

    end function pmc_aero_rep_data_species_surface_area_conc
  end interface pmc_aero_rep_data_species_surface_area_conc_if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            
  !> Get the vapor pressure scaling for a particular species (unitless)
  interface pmc_aero_rep_data_vapor_pressure_scaling_if  
    function pmc_aero_rep_data_vapor_pressure_scaling(this, i_spec, &
                    phlex_state, jac_contrib) result(vapor_pressure_scaling)
      use pmc_util,                                     only : dp, i_kind
      use pmc_phlex_state
      import :: aero_rep_data_t

      !> Vapor pressure scaling
      real(kind=dp) :: vapor_pressure_scaling
      !> Aerosol representation data
      class(aero_rep_data_t), intent(in) :: this
      !> Species name
      integer(kind=i_kind), intent(in) :: i_spec
      !> Model state
      type(phlex_state_t), intent(in) :: phlex_state
      !> Contribution to Jacobian matrix. An array of the same size as the
      !! state array that, when present, will be filled with the partial
      !! derivatives of the result of this calculation with each state
      !! variable.
      real(kind=dp), allocatable, intent(out), optional :: jac_contrib(:)

    end function pmc_aero_rep_data_vapor_pressure_scaling
  end interface pmc_aero_rep_data_vapor_pressure_scaling_if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get the name of the aerosol representation
  function pmc_aero_rep_data_name(this) result(name)

    !> Aerosol representation name
    character(len=:), allocatable :: name
    !> Aerosol representation data
    class(aero_rep_data_t), intent(in) :: this

    name = this%rep_name

  end function pmc_aero_rep_data_name

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get a phase id in this representation, for use in solver functions.
  !! Returns 0 if phase is not present.
  !!
  !! This function should only be called during initialization
  integer(kind=i_kind) function pmc_aero_rep_data_phase_id(this, phase_name) &
                  result (phase_id)

    !> Aerosol representation data
    class(aero_rep_data_t), intent(in) :: this
    !> Aerosol phase to find
    character(len=:), allocatable, intent(in) :: phase_name

    integer(kind=i_kind) :: i_phase

    phase_id = 0
    do i_phase = 1, size(this%aero_phase)
      if (this%aero_phase(i_phase)%val%name().eq.phase_name) then
        phase_id = i_phase
        return
      end if
    end do

  end function pmc_aero_rep_data_phase_id

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Determine the size of a binary required to pack the aerosol 
  !! representation data
  integer(kind=i_kind) function pmc_aero_rep_data_pack_size(this) &
                  result (pack_size)

    !> Aerosol representation data
    class(aero_rep_data_t), intent(in) :: this
    
    pack_size = &
            pmc_mpi_pack_size_real_array(this%condensed_data_real) + &
            pmc_mpi_pack_size_integer_array(this%condensed_data_int) 

  end function pmc_aero_rep_data_pack_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Pack the given value to the buffer, advancing position
  subroutine pmc_aero_rep_data_bin_pack(this, buffer, pos)

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

  end subroutine pmc_aero_rep_data_bin_pack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Unpack the given value from the buffer, advancing position
  subroutine pmc_aero_rep_data_bin_unpack(this, buffer, pos)

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

  end subroutine pmc_aero_rep_data_bin_unpack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Load an aerosol representation from an input file
#ifdef PMC_USE_JSON
  !! j_obj is expected to be a JSON object containing data related to an
  !! aerosol representation. It should be of the form:
  !! 
  !! { "pmc-data" : [
  !!   {
  !!     "name" : "my representation",
  !!     "type" : "AERO_REP_TYPE"
  !!     "some parameter" : 123.34,
  !!     "some other parameter" : true,
  !!     "nested parameters" : {
  !!        "sub param 1" : 12.43,
  !!        "sub param other" : "some text"
  !!   },
  !!   ...
  !! ]}
  !!
  !! Aerosol representations must have a unique name and have a "type" that
  !! corresponds to a valid aerosol representation type.
  !!
  !! All remaining data are optional and may include any valid JSON value, 
  !! including nested objects. However, extending types will have specific
  !! requirements for the remaining data. 
  subroutine pmc_aero_rep_data_load(this, json, j_obj)

    !> Aerosol representation data
    class(aero_rep_data_t), intent(out) :: this
    !> JSON core
    type(json_core), pointer, intent(in) :: json
    !> JSON object
    type(json_value), pointer, intent(in) :: j_obj

    type(json_value), pointer :: child, next, species
    character(kind=json_ck, len=:), allocatable :: key

    this%property_set => property_t()

    next => null()
    call json%get_child(j_obj, child)
    do while (associated(child))
      call json%info(child, name=key)
      if (key.ne."type" .and. key.ne."name") &
              call this%property_set%load(json, child, .false.)
      call json%get_next(child, next)
      child => next
    end do
#else
  subroutine pmc_aero_rep_data_load(this)

    !> Aerosol representation data
    class(aero_rep_data_t), intent(inout) :: this

    call warn_msg(433045149, "No support for input files")
#endif

  end subroutine pmc_aero_rep_data_load

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Print the aerosol representation data
  subroutine pmc_aero_rep_data_print(this)

    !> Aerosol representation data
    class(aero_rep_data_t), intent(in) :: this

    if (associated(this%property_set)) call this%property_set%print()

  end subroutine pmc_aero_rep_data_print

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_aero_rep_data
