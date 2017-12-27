! Copyright (C) 2017 Matt Dawson
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_rxn_factory module.

!> The abstract rxn_factory_t structure and associated subroutines.
module pmc_rxn_factory

  use pmc_constants,                  only : i_kind, dp
  use pmc_mpi
  use pmc_util,                       only : die_msg, string_t, assert_msg, &
                                             warn_msg
  use pmc_rxn_data
#ifdef PMC_USE_MPI
  use mpi
#endif
#ifdef PMC_USE_JSON
  use json_module
#endif

  ! Use all reaction modules
  use pmc_rxn_arrhenius

  implicit none
  private

  public :: rxn_factory_t

  !> Identifiers for reaction types - used by binary packing/unpacking 
  !! functions
  integer(kind=i_kind), parameter :: RXN_ARRHENIUS = 1

  !> Factory type for chemical reactions
  !!
  !! Provides new instances of types extending rxn_data_t by name or
  !! from input file data
  type :: rxn_factory_t
  contains
    !> Create a new chemical reaction by type name
    procedure :: create => pmc_rxn_factory_create
    !> Create a new chemical reaction from input data
    procedure :: load => pmc_rxn_factory_load
    !> Determine the number of bytes required to pack a given reaction
    procedure :: pack_size => pmc_rxn_factory_pack_size
    !> Pack a given reaction to the buffer, advancing the position
    procedure :: bin_pack => pmc_rxn_factory_bin_pack
    !> Unpack a reaction from the buffer, advancing the position
    procedure :: bin_unpack => pmc_rxn_factory_bin_unpack
  end type rxn_factory_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Create a new chemical reaction by type name
  function pmc_rxn_factory_create(this, type_name) result (new_obj)

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
      case default
        call die_msg(367114278, "Unknown chemical reaction type: " &
                //type_name) 
    end select

  end function pmc_rxn_factory_create

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Load an aerosol represenation from input data
#ifdef PMC_USE_JSON
  function pmc_rxn_factory_load(this, json, j_obj) result (new_obj)

    !> A chemical reaction
    class(rxn_data_t), pointer :: new_obj
    !> Chemical reaction factory
    class(rxn_factory_t), intent(in) :: this
    !> JSON core
    type(json_core), pointer, intent(in) :: json
    !> JSON object
    type(json_value), pointer, intent(in) :: j_obj

    type(json_value), pointer :: child
    character(kind=json_ck, len=:), allocatable :: key, unicode_type_name
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
  function pmc_rxn_factory_load(this) result (new_obj)

    !> A chemical reaction
    class(rxn_data_t), pointer :: new_obj
    !> Chemical reaction factory
    class(rxn_factory_t), intent(in) :: this

    new_obj => null()

    call warn_msg(979827016, "No support for input files.")
#endif
  end function pmc_rxn_factory_load

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Determine the size of a binary required to pack a reaction
  integer(kind=i_kind) function pmc_rxn_factory_pack_size(this, rxn) &
                  result (pack_size)

    !> Reaction factory
    class(rxn_factory_t) :: this
    !> Reaction to pack
    class(rxn_data_t), intent(in) :: rxn

    integer(kind=i_kind) :: i_rxn

    pack_size =  pmc_mpi_pack_size_integer(int(1, kind=i_kind)) + &
                 rxn%pack_size()

  end function pmc_rxn_factory_pack_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Pack the given value to the buffer, advancing position
  subroutine pmc_rxn_factory_bin_pack(this, rxn, buffer, pos)

    !> Reaction factory
    class(rxn_factory_t), intent(in) :: this
    !> Reaction to pack
    class(rxn_data_t), intent(in) :: rxn
    !> Memory buffer
    character, intent(inout) :: buffer(:)
    !> Current buffer position
    integer, intent(inout) :: pos

#ifdef PMC_USE_MPI
    integer :: rxn_type, i_rxn, prev_position

    prev_position = pos
    select type (rxn)
      type is (rxn_arrhenius_t)
        rxn_type = RXN_ARRHENIUS
      class default
        call die_msg(343941184, "Trying to pack reaction of unknown type.")
    end select
    call pmc_mpi_pack_integer(buffer, pos, rxn_type)
    call rxn%bin_pack(buffer, pos)
    call assert(194676336, &
         pos - prev_position <= this%pack_size(rxn))
#endif

  end subroutine pmc_rxn_factory_bin_pack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Unpack the given value to the buffer, advancing position
  function pmc_rxn_factory_bin_unpack(this, buffer, pos) result (rxn)

    !> Unpacked reaction
    class(rxn_data_t), pointer :: rxn
    !> Reaction factory
    class(rxn_factory_t), intent(in) :: this
    !> Memory buffer
    character, intent(inout) :: buffer(:)
    !> Current buffer position
    integer, intent(inout) :: pos

#ifdef PMC_USE_MPI
    integer :: rxn_type, i_rxn, prev_position

    prev_position = pos
    call pmc_mpi_unpack_integer(buffer, pos, rxn_type)
    select case (rxn_type)
      case (RXN_ARRHENIUS)
        rxn => rxn_arrhenius_t()
      case default
        call die_msg(659290342, "Trying to unpack reaction of unknown type:"// &
                to_string(rxn_type))
    end select
    call rxn%bin_unpack(buffer, pos)
    call assert(880568259, &
         pos - prev_position <= this%pack_size(rxn))
#endif

  end function pmc_rxn_factory_bin_unpack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_rxn_factory
