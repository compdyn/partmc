! Copyright (C) 2017 Matt Dawson
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_rxn_data module.

!> The rxn_data_t structure and associated subroutines.
module pmc_rxn_data

  use pmc_constants,                  only : i_kind
  use pmc_mpi
  use pmc_util,                       only : die_msg, string_t
  use pmc_property
  use pmc_chem_spec_data
#ifdef PMC_USE_MPI
  use mpi
#endif
#ifdef PMC_USE_JSON
  use json_module
#endif

  implicit none
  private

  public :: rxn_data_t

  !> Abstract reaction data type
  !!
  !! Time-invariant data related to a chemical reaction. Types extending
  !! rxn_data_t should represent specific reaction types, with unique
  !! rate equations. Extending types should not have data members, as these
  !! will not be passed to the child nodes after initialization. Instead all
  !! data required by an extending type during a model run should be packed
  !! into the condensed_data arrays during initialization.
  type, abstract rxn_data_t
    private
    !> Reactant indices in the model_state_t instance
    integer(kind=i_kind), allocatable :: reactant_id(:)
    !> Product indices in the model_state_t instance
    integer(kind=i_kind), allocatable :: product_id(:)
    !> Product yields (unitless)
    real(kind=dp), allocatable :: product_yield(:)
    !> Reaction parameters. These will be available during initialization,
    !! but not during integration. All information required to calculate
    !! the time derivatives and jacobian matrix constributions must be
    !! saved by the exdending type.
    type(property_t), pointer :: property_set => null()
    !> Condensed reaction data. Theses arrays will be available during
    !! integration, and should contain any information required by the
    !! rate and Jacobian constribution functions that cannot be obtained
    !! from the model_state_t object.
    real(kind=dp), allocatable :: condensed_data_real(:)
    !> Condensed reaction data (integers)
    integer(kind=dp), allocatable ::  condensed_data_int(:)
  contains
    public
    !> Reaction initialization
    procedure(pmc_rxn_data_initialize), deferred :: initialize
    !> Rate calculation
    procedure(pmc_rxn_data_rate), deferred :: rate
    !> Jacobian matrix contribution
    procedure(pmc_rxn_data_jac_contrib), deferred :: jac_contrib
    !> Determine the number of bytes required to pack the given value
    procedure :: pack_size => pmc_rxn_data_pack_size
    !> Packs the given value into the buffer, advancing position
    procedure :: bin_pack => pmc_rxn_data_bin_pack
    !> Unpacks the given value from the buffer, advancing position
    procedure :: bin_unpack => pmc_rxn_data_bin_unpack
    !> Load data from an input file
    procedure :: load => pmc_rxn_data_load
  end type rxn_data_t

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Initialize the reaction data, validating component data and loading
  !! any required information from reactant, product and reaction 
  !! property_t objects. This routine should be called once for each reaction
  !! at the beginning of a model run after all the input files have been
  !! read in. It ensure all data required during the model run are included
  !! in the condensed data arrays.
  interface pmc_rxn_data_initialize_if
    subroutine pmc_rxn_data_initialize(this, chem_spec_data)
      import :: rxn_data_t

      !> Reaction data
      class(rxn_data_t), intent(inout) :: this
      !> Chemical species data
      type(chem_spec_data_t), intent(in) :: chem_spec_data

    end subroutine pmc_rxn_data_initialize
  end interface pmc_rxn_data_initialize_if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Calculate the contribution to the time derivative vector. The current 
  !! model state is provided for species concentrations, aerosol state and 
  !! environmental variables. All other parameters must have been saved to the
  !! reaction data instance during initialization.
  interface pmc_rxn_data_func_contrib_if
    subroutine pmc_rxn_data_func_contrib(this, model_state, func) 
      import :: rxn_data_t

      !> Reaction data
      class(rxn_data_t), intent(in) :: this
      !> Current model state
      type(model_state_t), intent(in) :: model_state
      !> Time derivative vector. This vector may include contributions from
      !! other reactions, so the contributions from this reaction should
      !! append, not overwrite, the values already in the vector
      real(kind=dp), allocatable, intent(inout) :: func(:)

    end function pmc_rxn_data_func_contrib
  end interface pmc_rxn_data_func_contrib_if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Calculate the contribution of this reaction to the Jacobian matrix.
  !! The current model state is provided for species concentrations and 
  !! aerosol state. All other parameters must have been saved to the reaction 
  !! data instance during initialization.
  interface pmc_rxn_data_jac_contrib_if
    subroutine pmc_rxn_data_jac_contrib(this, model_state, jac_matrix)
      import :: rxn_data_t

      !> Reaction data
      class(rxn_data_t), intent(in) :: this
      !> Current model state
      type(model_state_t), intent(in) :: model_state
      !> Jacobian matrix. This matrix may include contributions from other
      !! reactions, so the contributions from this reaction should append,
      !! not overwrite, the values already in the matrix.
      real(kind=dp), allocatable, intent(inout) :: jac_matrix(:,:)

    end subroutine pmc_rxn_data_jac_contrib
  end interface pmc_rxn_data_jac_contrib_if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Determine the size of a binary required to pack the reaction data
  integer(kind=i_kind) function pmc_rxn_data_pack_size(this) &
                  result (pack_size)

    !> Reaction data
    class(rxn_data_t), intent(in) :: this
    
    pack_size = &
            pmc_mpi_pack_size_integer_array(this%reactant_id) + &
            pmc_mpi_pack_size_integer_array(this%product_id) + &
            pmc_mpi_pack_size_real_array(this%product_yield) + &
            pmc_mpi_pack_size_real_array(this%condensed_data_real) + &
            pmc_mpi_pack_size_real_array(this%condensed_data_int)

  end function pmc_rxn_data_pack_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Pack the given value to the buffer, advancing position
  subroutine pmc_rxn_data_bin_pack(this, buffer, pos)

    !> Reaction data
    class(rxn_data_t), intent(in) :: this
    !> Memory buffer
    character, intent(inout) :: buffer(:)
    !> Current buffer position
    integer, intent(inout) :: pos

#ifdef PMC_USE_MPI
    integer :: prev_position

    prev_position = pos
    call pmc_mpi_pack_integer_array(buffer, pos, this%reactant_id)
    call pmc_mpi_pack_integer_array(buffer, pos, this%product_id)
    call pmc_mpi_pack_real_array(buffer, pos, this%product_yield)
    call pmc_mpi_pack_real_array(buffer, pos, this%condensed_data_real)
    call pmc_mpi_pack_real_array(buffer, pos, this%condensed_data_int)
    call assert(149359274, &
         pos - prev_position <= this%pack_size())
#endif

  end subroutine pmc_rxn_data_bin_pack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Unpack the given value from the buffer, advancing position
  subroutine pmc_rxn_data_bin_unpack(this, buffer, pos)

    !> Reaction data
    class(rxn_data_t), intent(out) :: this
    !> Memory buffer
    character, intent(in) :: buffer(:)
    !> Current buffer position
    integer, intent(inout) :: pos

#ifdef PMC_USE_MPI
    integer :: prev_position

    prev_position = pos
    call pmc_mpi_unpack_integer_array(buffer, pos, this%reactant_id)
    call pmc_mpi_unpack_integer_array(buffer, pos, this%product_id)
    call pmc_mpi_unpack_real_array(buffer, pos, this%product_yield)
    call pmc_mpi_unpack_real_array(buffer, pos, this%condensed_data_real)
    call pmc_mpi_unpack_real_array(buffer, pos, this%condensed_data_int)
    call assert(168345796, &
         pos - prev_position <= this%pack_size())
#endif

  end subroutine pmc_rxn_data_bin_unpack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Load reactions from an input file
#ifdef PMC_USE_JSON
  !! j_obj is expected to be a JSON object containing data related to a
  !! chemical reaction required for building the chemical mechanism. It should
  !! be of the form:
  !! 
  !! { "pmc-data" : [
  !!   {
  !!     "name" : "my mechanism",
  !!     "type" : "MECHANISM",
  !!     "reactions" : [
  !!     {
  !!       "rxn type" : "REACTION_TYPE",
  !!       "reactants" : {
  !!         "some species" : {},
  !!         "another species" : { "qty" : 2 }
  !!       },
  !!       "products" : {
  !!         "product species" : { "yield" : 0.42 },
  !!         "another prod species" : {}
  !!       },
  !!       "some parameter" : 123.34,
  !!       "some other parameter" : true,
  !!       "nested parameters" : {
  !!          "sub param 1" : 12.43,
  !!          "sub param other" : "some text"
  !!       },
  !!       ...
  !!     },
  !!     {
  !!       "rxn type" : "REACTION_TYPE",
  !!       "reactants" : {
  !!         "that one species" : { "qty" : 3 }
  !!       },
  !!       "some parameter" : 123.34,
  !!       ...
  !!     },
  !!     ...
  !!   ]},
  !!   ...
  !! ]}
  !!
  !! Reactions must be placed in an array with the key name "reactions" in a
  !! mechanism object. Reactions outside of mechanisms are not permitted.
  !! The key-value pair "rxn type" is required and its value must correspond
  !! to a valid reaction type. 
  !!
  !! All remaining data are optional and may include any valid JSON value, 
  !! including nested objects. However, extending types will have specific
  !! requirements for the remaining data. Additionally it is recommended to 
  !! use the above format for reactants and products when developing child
  !! reaction derived types, for consistency.
  !!
  !! Mechanism data may be inter-mixed with json objects of other types (e.g.,
  !! species), but the there must be exactly one top-level key-value pair 
  !! named "pmc-data" per input file whose value is an array of json objects 
  !! with valid PMC types.
  subroutine rxn_data_load(this, json, j_obj)

    !> Reaction data
    class(rxn_data_t), intent(out) :: this
    !> JSON core
    type(json_core), pointer, intent(in) :: json
    !> JSON object
    type(json_value), pointer, intent(in) :: j_obj

    type(json_value), pointer :: child, next, species
    character(kind=json_ck, len=:), allocatable :: key

    type(property_t), pointer :: property_set

    property_set = property_t()

    next => null()
    call json%get_child(j_obj, child)
    do while (associated(child))
      call json%info(child, name=key)
      if (key.ne."rxn type") call property_set%load(json, child, .false.)
      call json%get_next(child, next)
      child => next
    end do

    this%property_set => property_set
#else
  subroutine rxn_data_load(this)

    !> Reaction data
    class(rxn_data_t), intent(inout) :: this

    call warn_msg(332862889, "No support for input files")
#endif

  end subroutine rxn_data_load

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_rxn_data
