! Copyright (C) 2017 Matt Dawson
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_rxn_data module.

!> \page phlex_rxn Phlexible Module for Chemistry: Reactions
!!
!! Descriptions

!> The rxn_data_t structure and associated subroutines.
module pmc_rxn_data

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

  public :: rxn_data_t, rxn_data_ptr

  !> Gas-phase reaction
  integer(kind=i_kind), parameter, public :: GAS_RXN = 1
  !> Mixed-phase (gas and aerosol) reaction
  integer(kind=i_kind), parameter, public :: GAS_AERO_RXN = 2
  !> Aerosol-phase reaction
  integer(kind=i_kind), parameter, public :: AERO_RXN = 3

  !> Abstract reaction data type
  !!
  !! Time-invariant data related to a chemical reaction. Types extending
  !! rxn_data_t should represent specific reaction types, with unique
  !! rate equations. Extending types should not have data members, as these
  !! will not be passed to the child nodes after initialization. Instead all
  !! data required by an extending type during a model run should be packed
  !! into the condensed_data arrays during initialization.
  type, abstract :: rxn_data_t
    private
    !> Reaction phase
    !! FIXME This should be a private component, but extending types need
    !! access to it.
    integer(kind=i_kind), public :: rxn_phase
    !> Reaction parameters. These will be available during initialization,
    !! but not during integration. All information required to calculate
    !! the time derivatives and jacobian matrix constributions must be
    !! saved by the exdending type.
    !! FIXME This should be a private component, but extending types need
    !! access to it.
    type(property_t), pointer, public :: property_set => null()
    !> Condensed reaction data. Theses arrays will be available during
    !! integration, and should contain any information required by the
    !! rate and Jacobian constribution functions that cannot be obtained
    !! from the phlex_state_t object.
    !! FIXME This should be a private component, but extending types need
    !! access to it.
    real(kind=dp), allocatable, public :: condensed_data_real(:)
    !> Condensed reaction data (integers)
    !! FIXME This should be a private component, but extending types need
    !! access to it.
    integer(kind=i_kind), allocatable, public ::  condensed_data_int(:)
  contains
    !> Reaction initialization
    procedure(initialize), deferred :: initialize
    !> Rate contribution
    procedure(func_contrib), deferred :: func_contrib
    !> Jacobian matrix contribution
    procedure(jac_contrib), deferred :: jac_contrib
    !> Check the phase of the reaction against the phase being solved for.
    !! During GAS_RXN integrations, only GAS_RXN reactions are solved.
    !! During AERO_RXN integrations, only AERO_RXN and GAS_AERO_RXN
    !! reactions are solved. During GAS_AERO_RXN integrations, all 
    !! reactions are solved.
    procedure :: check_phase
    !> Determine the number of bytes required to pack the given value
    procedure :: pack_size
    !> Packs the given value into the buffer, advancing position
    procedure :: bin_pack
    !> Unpacks the given value from the buffer, advancing position
    procedure :: bin_unpack
    !> Load data from an input file
    procedure :: load
    !> Print the reaction data
    procedure :: print => do_print
  end type rxn_data_t

  !> Pointer type for building arrays of mixed reactions
  type :: rxn_data_ptr
    class(rxn_data_t), pointer :: val
  end type rxn_data_ptr

interface

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Initialize the reaction data, validating component data and loading
  !! any required information from reactant, product and reaction 
  !! property_t objects. This routine should be called once for each reaction
  !! at the beginning of a model run after all the input files have been
  !! read in. It ensure all data required during the model run are included
  !! in the condensed data arrays.
  subroutine initialize(this, chem_spec_data)
    import :: rxn_data_t, chem_spec_data_t

    !> Reaction data
    class(rxn_data_t), intent(inout) :: this
    !> Chemical species data
    type(chem_spec_data_t), intent(in) :: chem_spec_data

  end subroutine initialize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Calculate the contribution to the time derivative vector. The current 
  !! model state is provided for species concentrations, aerosol state and 
  !! environmental variables. All other parameters must have been saved to the
  !! reaction data instance during initialization.
  subroutine func_contrib(this, phlex_state, func) 
    import :: rxn_data_t, phlex_state_t, dp

    !> Reaction data
    class(rxn_data_t), intent(in) :: this
    !> Current model state
    type(phlex_state_t), intent(in) :: phlex_state
    !> Time derivative vector. This vector may include contributions from
    !! other reactions, so the contributions from this reaction should
    !! append, not overwrite, the values already in the vector
    real(kind=dp), pointer, intent(inout) :: func(:)

  end subroutine func_contrib

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Calculate the contribution of this reaction to the Jacobian matrix.
  !! The current model state is provided for species concentrations and 
  !! aerosol state. All other parameters must have been saved to the reaction 
  !! data instance during initialization.
  subroutine jac_contrib(this, phlex_state, jac_matrix)
    import :: rxn_data_t, phlex_state_t, dp

    !> Reaction data
    class(rxn_data_t), intent(in) :: this
    !> Current model state
    type(phlex_state_t), intent(in) :: phlex_state
    !> Jacobian matrix. This matrix may include contributions from other
    !! reactions, so the contributions from this reaction should append,
    !! not overwrite, the values already in the matrix.
    real(kind=dp), pointer, intent(inout) :: jac_matrix(:,:)

  end subroutine jac_contrib

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end interface

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Check that the phase of the reaction is being solved
  logical function check_phase(this, rxn_phase) result (valid_rxn)

    !> Reaction data
    class(rxn_data_t), intent(in) :: this
    !> Phase being solved
    integer(kind=i_kind), intent(in) :: rxn_phase

    if (rxn_phase.eq.GAS_AERO_RXN .or. &
        rxn_phase.eq.this%rxn_phase .or. &
        (rxn_phase.eq.AERO_RXN .and. this%rxn_phase.eq.GAS_AERO_RXN)) then
        valid_rxn = .true.
    else
      valid_rxn = .false.
    end if

  end function check_phase

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Determine the size of a binary required to pack the reaction data
  integer(kind=i_kind) function pack_size(this)

    !> Reaction data
    class(rxn_data_t), intent(in) :: this
    
    pack_size = &
            pmc_mpi_pack_size_real_array(this%condensed_data_real) + &
            pmc_mpi_pack_size_integer_array(this%condensed_data_int)

  end function pack_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Pack the given value to the buffer, advancing position
  subroutine bin_pack(this, buffer, pos)

    !> Reaction data
    class(rxn_data_t), intent(in) :: this
    !> Memory buffer
    character, intent(inout) :: buffer(:)
    !> Current buffer position
    integer, intent(inout) :: pos

#ifdef PMC_USE_MPI
    integer :: prev_position

    prev_position = pos
    call pmc_mpi_pack_real_array(buffer, pos, this%condensed_data_real)
    call pmc_mpi_pack_integer_array(buffer, pos, this%condensed_data_int)
    call assert(149359274, &
         pos - prev_position <= this%pack_size())
#endif

  end subroutine bin_pack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Unpack the given value from the buffer, advancing position
  subroutine bin_unpack(this, buffer, pos)

    !> Reaction data
    class(rxn_data_t), intent(out) :: this
    !> Memory buffer
    character, intent(inout) :: buffer(:)
    !> Current buffer position
    integer, intent(inout) :: pos

#ifdef PMC_USE_MPI
    integer :: prev_position

    prev_position = pos
    call pmc_mpi_unpack_real_array(buffer, pos, this%condensed_data_real)
    call pmc_mpi_unpack_integer_array(buffer, pos, this%condensed_data_int)
    call assert(168345796, &
         pos - prev_position <= this%pack_size())
#endif

  end subroutine bin_unpack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> \page input_format_rxn Input JSON Object Format: Reaction (general)
  !! 
  !! A \c json object containing information about a chemical reaction or 
  !! physical process in the gas phase, in an \ref phlex_aero_phase 
  !! "aerosol phase", or between two phases (phase-transfer). \ref phlex_rxn
  !! "Reactions" are used to build \ref phlex_mechanism "mechanisms" and are
  !! only found within an input \ref input_format_mechanism "mechanism object"
  !! in an array labelled \b reactions.
  !! \code{.json}
  !! { "pmc-data" : [
  !!   {
  !!     "name" : "my mechanism",
  !!     "type" : "MECHANISM",
  !!     "reactions" : [
  !!     {
  !!       "type" : "REACTION_TYPE",
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
  !!       "type" : "REACTION_TYPE",
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
  !! \endcode
  !! The key-value pair \b type is required and its value must correspond
  !! to a valid reaction type. Valid reaction types include:
  !!
  !!   - \subpage phlex_rxn_arrhenius "Arrhenius"
  !! 
  !! All remaining data are optional and may include any valid \c json value, 
  !! including nested objects. However, extending types (i.e. reactions) will
  !! have specific requirements for the remaining data. Additionally it is
  !! recommended to use the above format for reactants and products when
  !! developing derived types that extend \c rxn_data_t, and to use \b type
  !! values that match the name of the extending derived-type. For example, the 
  !! reaction type \c rxn_photolysis_t would have a \b type of \b PHOTOLYSIS.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef PMC_USE_JSON
  !> Load reactions from an input file
  subroutine load(this, json, j_obj)

    !> Reaction data
    class(rxn_data_t), intent(out) :: this
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
      if (key.ne."rxn type") call this%property_set%load(json, child, .false.)
      call json%get_next(child, next)
      child => next
    end do
#else
  subroutine load(this)

    !> Reaction data
    class(rxn_data_t), intent(inout) :: this

    call warn_msg(332862889, "No support for input files")
#endif

  end subroutine load

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Print the reaction data
  subroutine do_print(this)

    !> Reaction data
    class(rxn_data_t), intent(in) :: this

    if (associated(this%property_set)) call this%property_set%print()

  end subroutine do_print

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_rxn_data
