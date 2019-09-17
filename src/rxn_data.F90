! Copyright (C) 2017-2018 Matt Dawson
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_rxn_data module.

!> \page camp_rxn CAMP: Reactions (general)
!!
!! A reaction represents a transformation of the model state due to a physical
!! or chemical process that occurs within a phase (gas or \ref
!! camp_aero_phase "aerosol") or across the interface between two phases. In
!! the \ref camp_chem "camp-chem" model, reactions are grouped into \ref
!! camp_mechanism "mechanisms", which are solved over time-steps specified by
!! the host model.
!!
!! The primary function of a reaction in the \ref camp_chem "camp-chem"
!! model is to provide the solver with contributions to the time derivative
!! and Jacobian matrix for \ref camp_species "chemical species"
!! concentrations based on the current model state described in a \c
!! pmc_camp_state::camp_state_t object.
!!
!! Specific reaction types extend the abstract \c pmc_rxn_data::rxn_data_t
!! type and generally accept a set of reactants and products whose names
!! correspond to \ref camp_species "chemical species" names, as well as a
!! set of reaction parameters needed to describe a particular reaction.
!! During initialization, a reaction will have access to its set of parameters
!! as well as the parameters of any \ref camp_species "species" and \ref
!! camp_aero_rep "aerosol phase" in the \ref camp_chem "camp-chem" model,
!! however this information will not be available during a model run. The
!! information required by the reaction instance to calculate its contribution
!! to the time derivatve and Jacobian matrix must therefore be packed into
!! the condensed data arrays of the \c pmc_rxn_data::rxn_data_t object during
!! intialization.
!!
!! Valid reaction types include:
!!
!!   - \subpage camp_rxn_aqueous_equilibrium "aqueous-phase equilibrium"
!!   - \subpage camp_rxn_arrhenius "Arrhenius"
!!   - \subpage camp_rxn_CMAQ_H2O2 "CMAQ special reaction type for 2HO2 (+ H2O) -> H2O2"
!!   - \subpage camp_rxn_CMAQ_OH_HNO3 "CMAQ special reaction type for OH + HNO3 -> NO3 + H2O"
!!   - \subpage camp_rxn_condensed_phase_arrhenius "condensed-chase Arrhenius"
!!   - \subpage camp_rxn_emission "emission"
!!   - \subpage camp_rxn_first_order_loss "first-order loss"
!!   - \subpage camp_rxn_HL_phase_transfer "Henry's Law phase transfer"
!!   - \subpage camp_rxn_photolysis "photolysis"
!!   - \subpage camp_rxn_SIMPOL_phase_transfer "SIMPOL.1 phase transfer"
!!   - \subpage camp_rxn_troe "Troe (fall-off)"
!!   - \subpage camp_rxn_wet_deposition "wet deposition"
!!
!! The general input format for a reaction can be found
!! \subpage input_format_rxn "here".
!!
!! General instructions for adding a new reaction type can be found
!! \subpage camp_rxn_add "here".

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
  use pmc_camp_state
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
  !! Time-invariant data related to a \ref camp_rxn "reaction". Types
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
    !! from the \c pmc_camp_state::camp_state_t object. (floating-point)
    real(kind=dp), allocatable, public :: condensed_data_real(:)
    !> Condensed reaction data. Theses arrays will be available during
    !! integration, and should contain any information required by the
    !! rate and Jacobian constribution functions that cannot be obtained
    !! from the \c pmc_camp_state::camp_state_t object. (integer)
    integer(kind=i_kind), allocatable, public :: condensed_data_int(:)
    !> Number of environment-dependent parameters
    !! These are parameters that need updated when environmental conditions
    !! change
    integer(kind=i_kind), public :: num_env_params = 0
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
    !> Index for this reaction in the solver data
    integer(kind=c_int)  :: rxn_solver_id = 0
    !> Grid cell to update
    integer(kind=c_int) :: cell_id = 1
    !> Update data
    type(c_ptr) :: update_data
  contains
    !> Get the reaction type
    procedure :: get_type => rxn_update_data_get_type
    !> Get the grid cell to update
    procedure :: get_cell_id => rxn_update_data_get_cell_id
    !> Get the update data
    procedure :: get_data => rxn_update_data_get_data
    !> Determine the number of bytes required to pack the given value
    procedure :: pack_size => rxn_update_data_pack_size
    !> Packs the given value into the buffer, advancing position
    procedure :: bin_pack => rxn_update_data_bin_pack
    !> Unpacks the given value from the buffer, advancing position
    procedure :: bin_unpack => rxn_update_data_bin_unpack
    !> Check whether the update data is attached to a reaction or set of
    !! reactions
    procedure(is_attached), deferred :: is_attached
    !> Extending type pack size (internal use only)
    procedure(internal_pack_size), deferred :: internal_pack_size
    !> Extending type bin pack (internal use only)
    procedure(internal_bin_pack), deferred :: internal_bin_pack
    !> Extending type bin unpack (internal use only)
    procedure(internal_bin_unpack), deferred :: internal_bin_unpack
    !> Print the update data
    procedure :: print => do_rxn_update_data_print
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
  subroutine initialize(this, chem_spec_data, aero_rep, n_cells)
    use pmc_util,                                only : i_kind
    import :: rxn_data_t, chem_spec_data_t, aero_rep_data_ptr

    !> Reaction data
    class(rxn_data_t), intent(inout) :: this
    !> Chemical species data
    type(chem_spec_data_t), intent(in) :: chem_spec_data
    !> Aerosol representations
    type(aero_rep_data_ptr), pointer, intent(in) :: aero_rep(:)
    !> Number of grid cells to solve simultaneously
    integer(kind=i_kind), intent(in) :: n_cells

  end subroutine initialize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Check whether the update data is attached to a reaction or set of
  !! reactions
  logical function is_attached(this)
      import :: rxn_update_data_t

      !> Reaction data
      class(rxn_update_data_t), intent(in) :: this

  end function is_attached

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Extending-type binary pack size (internal use only)
  integer(kind=i_kind) function internal_pack_size(this, comm)
    use pmc_util,                                only : i_kind
    import :: rxn_update_data_t

    !> Reaction data
    class(rxn_update_data_t), intent(in) :: this
    !> MPI communicator
    integer, intent(in) :: comm

  end function internal_pack_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Extending-type binary pack function (Internal use only)
  subroutine internal_bin_pack(this, buffer, pos, comm)
    import :: rxn_update_data_t

    !> Reaction data
    class(rxn_update_data_t), intent(in) :: this
    !> Memory buffer
    character, intent(inout) :: buffer(:)
    !> Current buffer position
    integer, intent(inout) :: pos
    !> MPI communicator
    integer, intent(in) :: comm

  end subroutine internal_bin_pack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Extending-type binary unpack function (Internal use only)
  subroutine internal_bin_unpack(this, buffer, pos, comm)
    import :: rxn_update_data_t

    !> Reaction data
    class(rxn_update_data_t), intent(inout) :: this
    !> Memory buffer
    character, intent(inout) :: buffer(:)
    !> Current buffer position
    integer, intent(inout) :: pos
    !> MPI communicator
    integer, intent(in) :: comm

  end subroutine internal_bin_unpack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end interface

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> \page input_format_rxn Input JSON Object Format: Reaction (general)
  !!
  !! A \c json object containing information about a chemical reaction or
  !! physical process in the gas phase, in an \ref camp_aero_phase
  !! "aerosol phase", or between two phases (phase-transfer). \ref camp_rxn
  !! "Reactions" are used to build \ref camp_mechanism "mechanisms" and are
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
  !!   - \subpage camp_rxn_arrhenius "ARRHENIUS"
  !!   - \subpage camp_rxn_aqueous_equilibrium "AQUEOUS_EQUILIBRIUM"
  !!   - \subpage camp_rxn_CMAQ_H2O2 "CMAQ_H2O2"
  !!   - \subpage camp_rxn_CMAQ_OH_HNO3 "CMAQ_OH_HNO3"
  !!   - \subpage camp_rxn_condensed_phase_arrhenius "CONDENSED_PHASE_ARRHENIUS"
  !!   - \subpage camp_rxn_HL_phase_transfer "HL_PHASE_TRANSFER"
  !!   - \subpage camp_rxn_SIMPOL_phase_transfer "SIMPOL_PHASE_TRANSFER"
  !!   - \subpage camp_rxn_photolysis "PHOTOLYSIS"
  !!   - \subpage camp_rxn_troe "TROE"
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
    class(rxn_data_t), intent(inout) :: this
    !> JSON core
    type(json_core), pointer, intent(in) :: json
    !> JSON object
    type(json_value), pointer, intent(in) :: j_obj

    type(json_value), pointer :: child, next
    character(kind=json_ck, len=:), allocatable :: key
    character(len=:), allocatable :: owner_name

    ! allocate space for the reaction property set
    this%property_set => property_t()

    ! No names currently for reactions, so use generic label
    owner_name = "reaction"

    ! cycle through the reaction properties, loading them into the reaction
    ! property set
    next => null()
    call json%get_child(j_obj, child)
    do while (associated(child))
      call json%info(child, name=key)
      if (key.ne."rxn type") call this%property_set%load(json, child, &
                                     .false., owner_name)

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
  integer(kind=i_kind) function pack_size(this, comm)

    !> Reaction data
    class(rxn_data_t), intent(in) :: this
    !> MPI communicator
    integer, intent(in) :: comm

    pack_size = &
            pmc_mpi_pack_size_integer(this%rxn_phase, comm) + &
            pmc_mpi_pack_size_real_array(this%condensed_data_real, comm) + &
            pmc_mpi_pack_size_integer_array(this%condensed_data_int, comm) + &
            pmc_mpi_pack_size_integer(this%num_env_params, comm)

  end function pack_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Pack the given value to the buffer, advancing position
  subroutine bin_pack(this, buffer, pos, comm)

    !> Reaction data
    class(rxn_data_t), intent(in) :: this
    !> Memory buffer
    character, intent(inout) :: buffer(:)
    !> Current buffer position
    integer, intent(inout) :: pos
    !> MPI communicator
    integer, intent(in) :: comm

#ifdef PMC_USE_MPI
    integer :: prev_position

    prev_position = pos
    call pmc_mpi_pack_integer(buffer, pos, this%rxn_phase, comm)
    call pmc_mpi_pack_real_array(buffer, pos, this%condensed_data_real, comm)
    call pmc_mpi_pack_integer_array(buffer, pos, this%condensed_data_int, comm)
    call pmc_mpi_pack_integer(buffer, pos, this%num_env_params, comm)
    call assert(149359274, &
         pos - prev_position <= this%pack_size(comm))
#endif

  end subroutine bin_pack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Unpack the given value from the buffer, advancing position
  subroutine bin_unpack(this, buffer, pos, comm)

    !> Reaction data
    class(rxn_data_t), intent(out) :: this
    !> Memory buffer
    character, intent(inout) :: buffer(:)
    !> Current buffer position
    integer, intent(inout) :: pos
    !> MPI communicator
    integer, intent(in) :: comm

#ifdef PMC_USE_MPI
    integer :: prev_position

    prev_position = pos
    call pmc_mpi_unpack_integer(buffer, pos, this%rxn_phase, comm)
    call pmc_mpi_unpack_real_array(buffer, pos, this%condensed_data_real,comm)
    call pmc_mpi_unpack_integer_array(buffer, pos, this%condensed_data_int,  &
                                                                         comm)
    call pmc_mpi_unpack_integer(buffer, pos, this%num_env_params, comm)
    call assert(168345796, &
         pos - prev_position <= this%pack_size(comm))
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

  !> Get the grid cell id to update
  function rxn_update_data_get_cell_id(this) result(cell_id)

    !> Grid cell id
    integer(kind=c_int) :: cell_id
    !> Update data
    class(rxn_update_data_t), intent(in) :: this

    cell_id = this%cell_id

  end function rxn_update_data_get_cell_id

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

  !> Determine the size of a binary required to pack the reaction data
  integer(kind=i_kind) function rxn_update_data_pack_size(this, comm) &
      result(pack_size)

    !> Reaction update data
    class(rxn_update_data_t), intent(in) :: this
    !> MPI communicator
    integer, intent(in), optional :: comm

#ifdef PMC_USE_MPI
    integer :: l_comm

    if (present(comm)) then
      l_comm = comm
    else
      l_comm = MPI_COMM_WORLD
    endif

    pack_size = &
      pmc_mpi_pack_size_integer(int(this%rxn_type, kind=i_kind), l_comm) +   &
      pmc_mpi_pack_size_integer(int(this%rxn_solver_id, kind=i_kind),        &
                                                                   l_comm) + &
      this%internal_pack_size(l_comm)
#else
    pack_size = 0
#endif

  end function rxn_update_data_pack_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Pack the given value to the buffer, advancing position
  subroutine rxn_update_data_bin_pack(this, buffer, pos, comm)

    !> Reaction update data
    class(rxn_update_data_t), intent(in) :: this
    !> Memory buffer
    character, intent(inout) :: buffer(:)
    !> Current buffer position
    integer, intent(inout) :: pos
    !> MPI communicator
    integer, intent(in), optional :: comm

#ifdef PMC_USE_MPI
    integer :: prev_position, l_comm

    if (present(comm)) then
      l_comm = comm
    else
      l_comm = MPI_COMM_WORLD
    endif

    prev_position = pos
    call pmc_mpi_pack_integer(buffer, pos, &
                              int(this%rxn_type, kind=i_kind), l_comm)
    call pmc_mpi_pack_integer(buffer, pos, &
                              int(this%rxn_solver_id, kind=i_kind), l_comm)
    call this%internal_bin_pack(buffer, pos, l_comm)
    call assert(713360087, &
         pos - prev_position <= this%pack_size(l_comm))
#endif

  end subroutine rxn_update_data_bin_pack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Unpack the given value from the buffer, advancing position
  subroutine rxn_update_data_bin_unpack(this, buffer, pos, comm)

    !> Reaction update data
    class(rxn_update_data_t), intent(out) :: this
    !> Memory buffer
    character, intent(inout) :: buffer(:)
    !> Current buffer position
    integer, intent(inout) :: pos
    !> MPI communicator
    integer, intent(in), optional :: comm

#ifdef PMC_USE_MPI
    integer :: prev_position, l_comm
    integer(kind=i_kind) :: temp_int

    if (present(comm)) then
      l_comm = comm
    else
      l_comm = MPI_COMM_WORLD
    endif

    prev_position = pos
    call pmc_mpi_unpack_integer(buffer, pos, temp_int, l_comm)
    this%rxn_type = int(temp_int, kind=c_int)
    call pmc_mpi_unpack_integer(buffer, pos, temp_int, l_comm)
    this%rxn_solver_id = int(temp_int, kind=c_int)
    call this%internal_bin_unpack(buffer, pos, l_comm)
    call assert(107364895, &
         pos - prev_position <= this%pack_size(l_comm))
#endif

  end subroutine rxn_update_data_bin_unpack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Print the update data
  subroutine do_rxn_update_data_print(this)

    !> Reaction update data
    class(rxn_update_data_t), intent(in) :: this

    write(*,*) "*** Reaction update data ***"
    write(*,*) "Rxn type", this%rxn_type
    write(*,*) "Rxn solver id", this%rxn_solver_id

  end subroutine do_rxn_update_data_print

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_rxn_data
