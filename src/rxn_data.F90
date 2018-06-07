! Copyright (C) 2017 Matt Dawson
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_rxn_data module.

!> \page phlex_rxn Phlexible Module for Chemistry: Reactions (general)
!!
!! A reaction represents a transformation of the model state due to a physical
!! or chemical process that occurs within a phase (gas or \ref
!! phlex_aero_phase "aerosol") or across the interface between two phases. In
!! the \ref phlex_chem "phlex-chem" model, reactions are grouped into \ref
!! phlex_mechanism "mechanisms", which are solved over time-steps specified by
!! the host model.
!!
!! The primary function of a reaction in the \ref phlex_chem "phlex-chem"
!! model is to provide the solver with contributions to the time derivative 
!! and Jacobian matrix for \ref phlex_species "chemical species"
!! concentrations based on the current model state described in a \c
!! pmc_phlex_state::phlex_state_t object.
!!
!! Specific reaction types extend the abstract \c pmc_rxn_data::rxn_data_t
!! type and generally accept a set of reactants and products whose names
!! correspond to \ref phlex_species "chemical species" names, as well as a
!! set of reaction parameters needed to describe a particular reaction.
!! During initialization, a reaction will have access to its set of parameters
!! as well as the parameters of any \ref phlex_species "species" and \ref
!! phlex_aero_rep "aerosol phase" in the \ref phlex_chem "phlex-chem" model,
!! however this information will not be available during a model run. The
!! information required by the reaction instance to calculate its contribution
!! to the time derivatve and Jacobian matrix must therefore be packed into
!! the condensed data arrays of the \c pmc_rxn_data::rxn_data_t object during
!! intialization.
!!  
!! Valid reaction types include:
!!
!!   - \subpage phlex_rxn_aqueous_equilibrium "Aqueous-phase Equilibrium"
!!   - \subpage phlex_rxn_arrhenius "Arrhenius"
!!   - \subpage phlex_rxn_CMAQ_H2O2 "CMAQ special reaction type for 2HO2 (+ H2O) -> H2O2"
!!   - \subpage phlex_rxn_CMAQ_OH_HNO3 "CMAQ special reaction type for OH + HNO3 -> NO3 + H2O"
!!   - \subpage phlex_rxn_condensed_phase_arrhenius "Condensed-Phase Arrhenius"
!!   - \subpage phlex_rxn_HL_phase_transfer "Henry's Law Phase Transfer"
!!   - \subpage phlex_rxn_PDFiTE_activity "PD-FiTE Activity"
!!   - \subpage phlex_rxn_photolysis "Photolysis"
!!   - \subpage phlex_rxn_SIMPOL_phase_transfer "SIMPOL.1 Phase Transfer"
!!   - \subpage phlex_rxn_troe "Troe (fall-off)"
!!   - \subpage phlex_rxn_ZSR_aerosol_water "ZSR Aerosol Water"
!!
!! The general input format for a reaction can be found
!! \subpage input_format_rxn "here".
!!
!! General instructions for adding a new reaction type can be found
!! \subpage phlex_rxn_add "here".

!> The rxn_data_t structure and associated subroutines.
module pmc_rxn_data

#ifdef PMC_USE_JSON
  use json_module
#endif
#ifdef PMC_USE_MPI
  use mpi
#endif
  use pmc_aero_rep_data
  use pmc_chem_spec_data
  use pmc_constants,                  only : i_kind, dp
  use pmc_mpi
  use pmc_phlex_state
  use pmc_property
  use pmc_util,                       only : die_msg, string_t

  use iso_c_binding

  implicit none
  private

  public :: rxn_data_t, rxn_data_ptr, rxn_update_data_t

  !> Gas-phase reaction
  integer(kind=i_kind), parameter, public :: GAS_RXN = 1
  !> Mixed-phase (gas and aerosol) reaction
  integer(kind=i_kind), parameter, public :: GAS_AERO_RXN = 2
  !> Aerosol-phase reaction
  integer(kind=i_kind), parameter, public :: AERO_RXN = 3

  !> Abstract reaction data type
  !!
  !! Time-invariant data related to a \ref phlex_rxn "reaction". Types
  !! extending \c rxn_data_t should represent specific reaction types, with
  !! unique rate equations. Extending types should not have data members, as
  !! these will not be passed to the child nodes after initialization. Instead
  !! all data required by an extending type during a model run should be
  !! packed into the condensed_data arrays during initialization.
  type, abstract :: rxn_data_t
    private
    !> Reaction phase
    integer(kind=i_kind), public :: rxn_phase
    !> Reaction parameters. These will be available during initialization,
    !! but not during integration. All information required to calculate
    !! the time derivatives and Jacobian matrix constributions must be
    !! saved by the exdending type.
    type(property_t), pointer, public :: property_set => null()
    !> Condensed reaction data. Theses arrays will be available during
    !! integration, and should contain any information required by the
    !! rate and Jacobian constribution functions that cannot be obtained
    !! from the \c pmc_phlex_state::phlex_state_t object. (floating-point)
    real(kind=dp), allocatable, public :: condensed_data_real(:)
    !> Condensed reaction data. Theses arrays will be available during
    !! integration, and should contain any information required by the
    !! rate and Jacobian constribution functions that cannot be obtained
    !! from the \c pmc_phlex_state::phlex_state_t object. (integer)
    integer(kind=i_kind), allocatable, public :: condensed_data_int(:)
  contains
    !> Reaction initialization. Takes species, phase and reaction parameters
    !! and packs required information into the condensed data arrays for use
    !! during the model run.
    !!
    !! This routine should be called once for each reaction
    !! at the beginning of a model run after all the input files have been
    !! read in.
    procedure(initialize), deferred :: initialize
    !> Load data from an input file
    procedure :: load
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
    !> Print the reaction data
    procedure :: print => do_print
  end type rxn_data_t

  !> Pointer type for building arrays of mixed reactions
  type :: rxn_data_ptr
    class(rxn_data_t), pointer :: val => null()
  contains
    !> Dereference the pointer
    procedure :: dereference
    !> Finalize the pointer
    final :: ptr_finalize
  end type rxn_data_ptr

  !> Update cookie
  type, abstract :: rxn_update_data_t
    !> Reaction type
    integer(kind=c_int) :: rxn_type
    !> Update data
    type(c_ptr) :: update_data
  contains
    !> Get the reaction type
    procedure :: get_type => rxn_update_data_get_type
    !> Get the update data
    procedure :: get_data => rxn_update_data_get_data
  end type rxn_update_data_t

interface

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Reaction initialization. Takes species, phase and reaction parameters
  !! and packs required information into the condensed data arrays for use
  !! during the model run.
  !!
  !! This routine should be called once for each reaction
  !! at the beginning of a model run after all the input files have been
  !! read in.
  subroutine initialize(this, chem_spec_data, aero_rep)
    import :: rxn_data_t, chem_spec_data_t, aero_rep_data_ptr

    !> Reaction data
    class(rxn_data_t), intent(inout) :: this
    !> Chemical species data
    type(chem_spec_data_t), intent(in) :: chem_spec_data
    !> Aerosol representations
    class(aero_rep_data_ptr), pointer, intent(in) :: aero_rep(:)

  end subroutine initialize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end interface

contains

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
  !!   - \subpage phlex_rxn_arrhenius "ARRHENIUS"
  !!   - \subpage phlex_rxn_aqueous_equilibrium "AQUEOUS_EQUILIBRIUM"
  !!   - \subpage phlex_rxn_CMAQ_H2O2 "CMAQ_H2O2"
  !!   - \subpage phlex_rxn_CMAQ_OH_HNO3 "CMAQ_OH_HNO3"
  !!   - \subpage phlex_rxn_condensed_phase_arrhenius "CONDENSED_PHASE_ARRHENIUS"
  !!   - \subpage phlex_rxn_HL_phase_transfer "HL_PHASE_TRANSFER"
  !!   - \subpage phlex_rxn_PDFiTE_activity "PDFITE_ACTIVITY"
  !!   - \subpage phlex_rxn_SIMPOL_phase_transfer "SIMPOL_PHASE_TRANSFER"
  !!   - \subpage phlex_rxn_photolysis "PHOTOLYSIS"
  !!   - \subpage phlex_rxn_troe "TROE"
  !!   - \subpage phlex_rxn_ZSR_aerosol_water "ZSR_AEROSOL_WATER"
  !! 
  !! All remaining data are optional and may include any valid \c json value, 
  !! including nested objects. However, extending types (i.e. reactions) will
  !! have specific requirements for the remaining data. Additionally it is
  !! recommended to use the above format for reactants and products when
  !! developing derived types that extend \c rxn_data_t, and to use \b type
  !! values that match the name of the extending derived-type. For example, the 
  !! reaction type \c rxn_photolysis_t would have a \b type of \b PHOTOLYSIS.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Load reactions from an input file
#ifdef PMC_USE_JSON
  subroutine load(this, json, j_obj)

    !> Reaction data
    class(rxn_data_t), intent(out) :: this
    !> JSON core
    type(json_core), pointer, intent(in) :: json
    !> JSON object
    type(json_value), pointer, intent(in) :: j_obj

    type(json_value), pointer :: child, next, species
    character(kind=json_ck, len=:), allocatable :: key

    ! allocate space for the reaction property set
    this%property_set => property_t()

    ! cycle through the reaction properties, loading them into the reaction
    ! property set
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

  !> Check the phase of the reaction against the phase being solved for.
  !! During GAS_RXN integrations, only GAS_RXN reactions are solved.
  !! During AERO_RXN integrations, only AERO_RXN and GAS_AERO_RXN
  !! reactions are solved. During GAS_AERO_RXN integrations, all 
  !! reactions are solved.
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

  !> Print the reaction data
  subroutine do_print(this, file_unit)

    !> Reaction data
    class(rxn_data_t), intent(in) :: this
    !> File unit for output
    integer(kind=i_kind), optional :: file_unit

    integer(kind=i_kind) :: f_unit = 6

    if (present(file_unit)) f_unit = file_unit
    write(f_unit,*) "*** Rxn ***"
    if (associated(this%property_set)) call this%property_set%print(f_unit)
    if (allocated(this%condensed_data_int)) &
      write (f_unit,*) "  *** condensed data int: ", &
            this%condensed_data_int(:)
    if (allocated(this%condensed_data_real)) &
      write (f_unit,*) "  *** condensed data real: ", &
            this%condensed_data_real(:)
  end subroutine do_print

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Dereference a pointer to a reaction
  elemental subroutine dereference(this)

    !> Pointer to a reaction
    class(rxn_data_ptr), intent(inout) :: this

    this%val => null()

  end subroutine dereference

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Finalize a pointer to a reaction
  elemental subroutine ptr_finalize(this)
    
    !> Pointer to a reaction
    type(rxn_data_ptr), intent(inout) :: this

    if (associated(this%val)) deallocate(this%val)

  end subroutine ptr_finalize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get the update data reaction type
  function rxn_update_data_get_type(this) result(rxn_type)

    !> Reaction type
    integer(kind=c_int) :: rxn_type
    !> Update data
    class(rxn_update_data_t), intent(in) :: this

    rxn_type = this%rxn_type

  end function rxn_update_data_get_type

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get the update data
  function rxn_update_data_get_data(this) result(update_data)

    !> Update data ptr
    type(c_ptr) :: update_data
    !> Update data
    class(rxn_update_data_t), intent(in) :: this

    update_data = this%update_data

  end function rxn_update_data_get_data

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_rxn_data
