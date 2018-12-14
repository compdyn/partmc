! Copyright (C) 2017-2018 Matt Dawson
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_rxn_factory module.

!> \page phlex_rxn_add Phlexible Module for Chemistry: Adding a Reaction Type
!!
!! \b Note: these instructions are out-of-date. TODO update
!!
!! Adding a \ref phlex_rxn "reaction" to the \ref phlex_chem "phlex-chem"
!! module can be done in the following steps:
!!
!! ## Step 1. Create a new reaction module ##
!!   The module should be placed in the \c /src/rxns folder and extent the
!!   abstract \c pmc_rxn_data::rxn_data_t type, overriding all deferred
!!   functions, and providing a constructor that returns a pointer to a newly
!!   allocated instance of the new type:
!!
!! \code{.f90}
!! module rxn_foo
!!
!!   use ...
!!
!!   implicit none
!!   private
!!
!!   public :: rxn_foo_t
!!
!!   type, extends(rxn_data_t) :: rxn_foo_t
!!   contains
!!      ... (all deferred functions) ...
!!   end type rxn_foo_t
!!
!!   ! Constructor
!!   interface rxn_foo_t
!!     procedure :: constructor
!!   end interface rxn_foo_t
!!
!! contains
!!
!!   function constructor() result(new_obj)
!!     type(rxn_foo_t), pointer :: new_obj
!!     allocate(new_obj)
!!   end function constructor
!!
!!   ...
!!
!! end module pmc_rxn_foo
!! \endcode
!!
!! ## Step 2. Add the reaction to the \c pmc_rxn_factory module ##
!!
!! \code{.f90}
!! module pmc_rxn_factory
!!
!!  ...
!!
!!  ! Use all reaction modules
!!  ...
!!  use pmc_rxn_foo
!!
!!  ...
!!
!!  !> Identifiers for reaction types - used by binary packing/unpacking
!!  !! functions
!!  ...
!!  integer(kind=i_kind), parameter :: RXN_FOO = 32
!!
!!  ...
!!
!!  !> Create a new chemical reaction by type name
!!  function create(this, type_name) result (new_obj)
!!    ...
!!    select case (type_name)
!!      ...
!!      case ("FOO")
!!        new_obj => rxn_foo_t()
!!    ...
!!  end function create
!!
!!  ...
!!
!!  !> Pack the given value to the buffer, advancing position
!!  subroutine bin_pack(this, rxn, buffer, pos)
!!    ...
!!    select type (rxn)
!!      ...
!!      type is (rxn_foo_t)
!!        rxn_type = RXN_FOO
!!    ...
!!  end subroutine bin_pack
!!
!!  ...
!!
!!  !> Unpack the given value to the buffer, advancing position
!!  function bin_unpack(this, buffer, pos) result (rxn)
!!    ...
!!    select case (rxn_type)
!!      ...
!!      case (RXN_FOO)
!!        rxn => rxn_foo_t()
!!    ...
!!  end function bin_unpack
!!
!!  ...
!!
!! end module pmc_rxn_factory
!! \endcode
!!
!! # Step 4. Add the new module to the CMakeList file in the root directory. ##
!!
!! \code{.unparsed}
!! ...
!!
!! # partmc library
!!
!! set(REACTIONS
!!     ...
!!     src/rxns/pmc_foo.F90
!! )
!!
!! ...
!! \endcode
!!
!! ## Step 5. Add unit tests for the new \c rxn_foo_t type ##
!!
!! Unit testing should cover, at minimum, the initialization, time derivative
!! and Jacbian matrix functions, and in general 80% code coverage is
!! recommended. Some examples can be found in the \c /src/test folder.
!!
!! ## Step 6.  Update documentation ##
!!
!! TODO finish...
!!
!! ## Usage ##
!! The new \ref phlex_rxn "reaction type" is now ready to use. To include a
!! reaction of this type in a \ref phlex_mechanism "mechanism", add a \ref
!! input_format_rxn "reaction object" to a new or existing \ref
!! input_format_phlex_config "phlex-chem configuration file" as part of a
!! \ref input_format_mechanism "mechanism object". The reaction should have a
!! \b type corresponding to the newly created reaction type, along with any
!! required parameters:
!!
!! \code{.json}
!! { "pmc-data" : [
!!   {
!!     "name" : "my mechanism",
!!     "type" : "MECHANISM",
!!     "reactions" : [
!!       {
!!         "type" : "FOO",
!!         ...
!!       },
!!       ...
!!     ]
!!   },
!!   ...
!! ]}
!! \endcode
!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> The abstract rxn_factory_t structure and associated subroutines.
module pmc_rxn_factory

#ifdef PMC_USE_JSON
  use json_module
#endif
#ifdef PMC_USE_MPI
  use mpi
#endif
  use pmc_constants,                  only : i_kind, dp
  use pmc_mpi
  use pmc_rxn_data
  use pmc_util,                       only : die_msg, string_t, assert_msg, &
                                             warn_msg

  ! Use all reaction modules
  use pmc_rxn_aqueous_equilibrium
  use pmc_rxn_arrhenius
  use pmc_rxn_CMAQ_H2O2
  use pmc_rxn_CMAQ_OH_HNO3
  use pmc_rxn_condensed_phase_arrhenius
  use pmc_rxn_emission
  use pmc_rxn_first_order_loss
  use pmc_rxn_HL_phase_transfer
  use pmc_rxn_PDFiTE_activity
  use pmc_rxn_photolysis
  use pmc_rxn_SIMPOL_phase_transfer
  use pmc_rxn_troe
  use pmc_rxn_wet_deposition
  use pmc_rxn_ZSR_aerosol_water

  use iso_c_binding

  implicit none
  private

  public :: rxn_factory_t

  !> Identifiers for reaction types - used by binary packing/unpacking
  !! functions
  integer(kind=i_kind), parameter, public :: RXN_ARRHENIUS = 1
  integer(kind=i_kind), parameter, public :: RXN_TROE = 2
  integer(kind=i_kind), parameter, public :: RXN_CMAQ_H2O2 = 3
  integer(kind=i_kind), parameter, public :: RXN_CMAQ_OH_HNO3 = 4
  integer(kind=i_kind), parameter, public :: RXN_PHOTOLYSIS = 5
  integer(kind=i_kind), parameter, public :: RXN_HL_PHASE_TRANSFER = 6
  integer(kind=i_kind), parameter, public :: RXN_AQUEOUS_EQUILIBRIUM = 7
  integer(kind=i_kind), parameter, public :: RXN_ZSR_AEROSOL_WATER = 8
  integer(kind=i_kind), parameter, public :: RXN_PDFITE_ACTIVITY = 9
  integer(kind=i_kind), parameter, public :: RXN_SIMPOL_PHASE_TRANSFER = 10
  integer(kind=i_kind), parameter, public :: RXN_CONDENSED_PHASE_ARRHENIUS = 11
  integer(kind=i_kind), parameter, public :: RXN_FIRST_ORDER_LOSS = 12
  integer(kind=i_kind), parameter, public :: RXN_EMISSION = 13
  integer(kind=i_kind), parameter, public :: RXN_WET_DEPOSITION = 14

  !> Factory type for chemical reactions
  !!
  !! Provides new instances of types extending rxn_data_t by name or
  !! from input file data
  type :: rxn_factory_t
  contains
    !> Create a new chemical reaction by type name
    procedure :: create
    !> Create a new chemical reaction from input data
    procedure :: load
    !> Get the reaction type
    procedure :: get_type
    !> Get a new update data object
    procedure :: initialize_update_data
    !> Determine the number of bytes required to pack a given reaction
    procedure :: pack_size
    !> Pack a given reaction to the buffer, advancing the position
    procedure :: bin_pack
    !> Unpack a reaction from the buffer, advancing the position
    procedure :: bin_unpack
  end type rxn_factory_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Create a new chemical reaction by type name
  function create(this, type_name) result (new_obj)

    !> A chemical reaction
    class(rxn_data_t), pointer :: new_obj
    !> Aerosol representation factory
    class(rxn_factory_t), intent(in) :: this
    !> Name of the chemical reaction
    character(len=:), allocatable :: type_name

    new_obj => null()

    ! Create a new reaction instance based on the type name supplied
    select case (type_name)
      case ("ARRHENIUS")
        new_obj => rxn_arrhenius_t()
      case ("TROE")
        new_obj => rxn_troe_t()
      case ("CMAQ_H2O2")
        new_obj => rxn_CMAQ_H2O2_t()
      case ("CMAQ_OH_HNO3")
        new_obj => rxn_CMAQ_OH_HNO3_t()
      case ("PHOTOLYSIS")
        new_obj => rxn_photolysis_t()
      case ("HL_PHASE_TRANSFER")
        new_obj => rxn_HL_phase_transfer_t()
      case ("AQUEOUS_EQUILIBRIUM")
        new_obj => rxn_aqueous_equilibrium_t()
      case ("ZSR_AEROSOL_WATER")
        new_obj => rxn_ZSR_aerosol_water_t()
      case ("PDFITE_ACTIVITY")
        new_obj => rxn_PDFiTE_activity_t()
      case ("SIMPOL_PHASE_TRANSFER")
        new_obj => rxn_SIMPOL_phase_transfer_t()
      case ("CONDENSED_PHASE_ARRHENIUS")
        new_obj => rxn_condensed_phase_arrhenius_t()
      case ("FIRST_ORDER_LOSS")
        new_obj => rxn_first_order_loss_t()
      case ("EMISSION")
        new_obj => rxn_emission_t()
      case ("WET_DEPOSITION")
        new_obj => rxn_wet_deposition_t()
      case default
        call die_msg(367114278, "Unknown chemical reaction type: " &
                //type_name)
    end select

  end function create

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Load an aerosol represenation from input data
#ifdef PMC_USE_JSON
  function load(this, json, j_obj) result (new_obj)

    !> A chemical reaction
    class(rxn_data_t), pointer :: new_obj
    !> Chemical reaction factory
    class(rxn_factory_t), intent(in) :: this
    !> JSON core
    type(json_core), pointer, intent(in) :: json
    !> JSON object
    type(json_value), pointer, intent(in) :: j_obj

    character(kind=json_ck, len=:), allocatable :: unicode_type_name
    character(len=:), allocatable :: type_name
    logical(kind=json_lk) :: found

    new_obj => null()

    ! Get the reaction type
    call json%get(j_obj, "type", unicode_type_name, found)
    call assert_msg(137665576, found, 'Missing chemical reaction type.')
    type_name = unicode_type_name

    ! Create a new reaction instance of the type specified
    new_obj => this%create(type_name)

    ! Load reaction parameters from the json object
    call new_obj%load(json, j_obj)

#else
  !> Generic warning function when no input file support exists
  function load(this) result (new_obj)

    !> A chemical reaction
    class(rxn_data_t), pointer :: new_obj
    !> Chemical reaction factory
    class(rxn_factory_t), intent(in) :: this

    new_obj => null()

    call warn_msg(979827016, "No support for input files.")
#endif
  end function load

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get the reaction type as a RxnType
  integer(kind=i_kind) function get_type(this, rxn) result (rxn_type)

    !> Reaction factory
    class(rxn_factory_t), intent(in) :: this
    !> Reaction to get type of
    class(rxn_data_t), intent(in) :: rxn

    select type (rxn)
      type is (rxn_arrhenius_t)
        rxn_type = RXN_ARRHENIUS
      type is (rxn_troe_t)
        rxn_type = RXN_TROE
      type is (rxn_CMAQ_H2O2_t)
        rxn_type = RXN_CMAQ_H2O2
      type is (rxn_CMAQ_OH_HNO3_t)
        rxn_type = RXN_CMAQ_OH_HNO3
      type is (rxn_photolysis_t)
        rxn_type = RXN_PHOTOLYSIS
      type is (rxn_HL_phase_transfer_t)
        rxn_type = RXN_HL_PHASE_TRANSFER
      type is (rxn_aqueous_equilibrium_t)
        rxn_type = RXN_AQUEOUS_EQUILIBRIUM
      type is (rxn_ZSR_aerosol_water_t)
        rxn_type = RXN_ZSR_AEROSOL_WATER
      type is (rxn_PDFiTE_activity_t)
        rxn_type = RXN_PDFITE_ACTIVITY
      type is (rxn_SIMPOL_phase_transfer_t)
        rxn_type = RXN_SIMPOL_PHASE_TRANSFER
      type is (rxn_condensed_phase_arrhenius_t)
        rxn_type = RXN_CONDENSED_PHASE_ARRHENIUS
      type is (rxn_first_order_loss_t)
        rxn_type = RXN_FIRST_ORDER_LOSS
      type is (rxn_emission_t)
        rxn_type = RXN_EMISSION
      type is (rxn_wet_deposition_t)
        rxn_type = RXN_WET_DEPOSITION
      class default
        call die_msg(343941184, "Unknown reaction type.")
    end select

  end function get_type

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Initialize an update data object
  subroutine initialize_update_data(this, update_data)

    !> Reaction factory
    class(rxn_factory_t), intent(in) :: this
    !> Update data object
    class(rxn_update_data_t), intent(out) :: update_data

    select type (update_data)
      type is (rxn_update_data_wet_deposition_rate_t)
        call update_data%initialize(RXN_WET_DEPOSITION)
      type is (rxn_update_data_emission_rate_t)
        call update_data%initialize(RXN_EMISSION)
      type is (rxn_update_data_first_order_loss_rate_t)
        call update_data%initialize(RXN_FIRST_ORDER_LOSS)
      type is (rxn_update_data_photolysis_rate_t)
        call update_data%initialize(RXN_PHOTOLYSIS)
      class default
        call die_msg(239438576, "Internal error - update data type missing.")
    end select

  end subroutine initialize_update_data

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Determine the size of a binary required to pack a reaction
  integer(kind=i_kind) function pack_size(this, rxn, comm)

    !> Reaction factory
    class(rxn_factory_t) :: this
    !> Reaction to pack
    class(rxn_data_t), intent(in) :: rxn
    !> MPI communicator
    integer, intent(in) :: comm

    pack_size =  pmc_mpi_pack_size_integer(int(1, kind=i_kind), comm) + &
                 rxn%pack_size(comm)

  end function pack_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Pack the given value to the buffer, advancing position
  subroutine bin_pack(this, rxn, buffer, pos, comm)

    !> Reaction factory
    class(rxn_factory_t), intent(in) :: this
    !> Reaction to pack
    class(rxn_data_t), intent(in) :: rxn
    !> Memory buffer
    character, intent(inout) :: buffer(:)
    !> Current buffer position
    integer, intent(inout) :: pos
    !> MPI communicator
    integer, intent(in) :: comm

#ifdef PMC_USE_MPI
    integer :: rxn_type, i_rxn, prev_position

    prev_position = pos
    select type (rxn)
      type is (rxn_arrhenius_t)
        rxn_type = RXN_ARRHENIUS
      type is (rxn_troe_t)
        rxn_type = RXN_TROE
      type is (rxn_CMAQ_H2O2_t)
        rxn_type = RXN_CMAQ_H2O2
      type is (rxn_CMAQ_OH_HNO3_t)
        rxn_type = RXN_CMAQ_OH_HNO3
      type is (rxn_photolysis_t)
        rxn_type = RXN_PHOTOLYSIS
      type is (rxn_HL_phase_transfer_t)
        rxn_type = RXN_HL_PHASE_TRANSFER
      type is (rxn_aqueous_equilibrium_t)
        rxn_type = RXN_AQUEOUS_EQUILIBRIUM
      type is (rxn_ZSR_aerosol_water_t)
        rxn_type = RXN_ZSR_AEROSOL_WATER
      type is (rxn_PDFiTE_activity_t)
        rxn_type = RXN_PDFITE_ACTIVITY
      type is (rxn_SIMPOL_phase_transfer_t)
        rxn_type = RXN_SIMPOL_PHASE_TRANSFER
      type is (rxn_condensed_phase_arrhenius_t)
        rxn_type = RXN_CONDENSED_PHASE_ARRHENIUS
      type is (rxn_first_order_loss_t)
        rxn_type = RXN_FIRST_ORDER_LOSS
      type is (rxn_emission_t)
        rxn_type = RXN_EMISSION
      type is (rxn_wet_deposition_t)
        rxn_type = RXN_WET_DEPOSITION
      class default
        call die_msg(343941184, "Trying to pack reaction of unknown type.")
    end select
    call pmc_mpi_pack_integer(buffer, pos, rxn_type, comm)
    call rxn%bin_pack(buffer, pos, comm)
    call assert(194676336, &
         pos - prev_position <= this%pack_size(rxn, comm))
#endif

  end subroutine bin_pack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Unpack the given value to the buffer, advancing position
  function bin_unpack(this, buffer, pos, comm) result (rxn)

    !> Unpacked reaction
    class(rxn_data_t), pointer :: rxn
    !> Reaction factory
    class(rxn_factory_t), intent(in) :: this
    !> Memory buffer
    character, intent(inout) :: buffer(:)
    !> Current buffer position
    integer, intent(inout) :: pos
    !> MPI communicator
    integer, intent(in) :: comm

#ifdef PMC_USE_MPI
    integer :: rxn_type, i_rxn, prev_position

    prev_position = pos
    call pmc_mpi_unpack_integer(buffer, pos, rxn_type, comm)
    select case (rxn_type)
      case (RXN_ARRHENIUS)
        rxn => rxn_arrhenius_t()
      case (RXN_TROE)
        rxn => rxn_troe_t()
      case (RXN_CMAQ_H2O2)
        rxn => rxn_CMAQ_H2O2_t()
      case (RXN_CMAQ_OH_HNO3)
        rxn => rxn_CMAQ_OH_HNO3_t()
      case (RXN_PHOTOLYSIS)
        rxn => rxn_photolysis_t()
      case (RXN_HL_PHASE_TRANSFER)
        rxn => rxn_HL_phase_transfer_t()
      case (RXN_AQUEOUS_EQUILIBRIUM)
        rxn => rxn_aqueous_equilibrium_t()
      case (RXN_ZSR_AEROSOL_WATER)
        rxn => rxn_ZSR_aerosol_water_t()
      case (RXN_PDFITE_ACTIVITY)
        rxn => rxn_PDFiTE_activity_t()
      case (RXN_SIMPOL_PHASE_TRANSFER)
        rxn => rxn_SIMPOL_phase_transfer_t()
      case (RXN_CONDENSED_PHASE_ARRHENIUS)
        rxn => rxn_condensed_phase_arrhenius_t()
      case (RXN_FIRST_ORDER_LOSS)
        rxn => rxn_first_order_loss_t()
      case (RXN_EMISSION)
        rxn => rxn_emission_t()
      case (RXN_WET_DEPOSITION)
        rxn => rxn_wet_deposition_t()
      case default
        call die_msg(659290342, &
                "Trying to unpack reaction of unknown type:"// &
                to_string(rxn_type))
    end select
    call rxn%bin_unpack(buffer, pos, comm)
    call assert(880568259, &
         pos - prev_position <= this%pack_size(rxn, comm))
#endif

  end function bin_unpack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_rxn_factory
