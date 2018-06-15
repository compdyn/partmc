! Copyright (C) 2017 Matt Dawson
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_sub_model_factory module.

!> \page phlex_sub_model_add Phlexible Module for Chemistry: Adding a Sub Model
!!
!! TODO write instructions
!!

!> The sub_model_factory_t type and associated subroutines
module pmc_sub_model_factory

#ifdef PMC_USE_JSON
  use json_module
#endif
#ifdef PMC_USE_MPI
  use mpi
#endif
  use pmc_constants,                    only : i_kind, dp
  use pmc_mpi
  use pmc_sub_model_data
  use pmc_util,                         only : die_msg, string_t, assert_msg, &
                                               warn_msg

  ! Use all sub-models
  use pmc_sub_model_UNIFAC

  implicit none
  private

  public :: sub_model_factory_t

  !> Identifiers for sub-models - used by binary packing/unpacking functions
  integer(kind=i_kind), parameter, public :: SUB_MODEL_UNIFAC = 1

  !> Factory type for sub-models
  !!
  !! Provides new instances of type extending sub_model_data_t by name or
  !! from input file data
  type :: sub_model_factory_t
  contains
    !> Create a new sub-model by type name
    procedure :: create
    !> Create a new aerosol representation from input data
    procedure :: load
    !> Get the aerosol representation type
    procedure :: get_type
    !> Get a new update data object
    procedure :: initialize_update_data
    !> Determine the number of bytes required to pack a given sub-model
    procedure :: pack_size
    !> Pack a given sub-model onto the buffer, advancing position
    procedure :: bin_pack
    !> Unpack a sub-model from the buffer, advancing position
    procedure :: bin_unpack
  end type sub_model_factory_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Create a new sub-model by type name
  function create(this, type_name) result (new_obj)

    !> A new sub-model
    class(sub_model_data_t), pointer :: new_obj
    !> Sub-model factory
    class(sub_model_factory_t), intent(in) :: this
    !> Type of the sub-model
    character(len=:), allocatable :: type_name

    new_obj => null()

    select case (type_name)
      case ("SUB_MODEL_UNIFAC")
        new_obj => sub_model_UNIFAC_t()
      case default
        call die_msg(293855421, "Unknown sub-model type: "//type_name)
    end select

  end function create

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Load a sub-model based on its type
#ifdef PMC_USE_JSON
  function load(this, json, j_obj) result (new_obj)

    !> A new sub-model
    class(sub_model_data_t), pointer :: new_obj
    !> Sub-model factory
    class(sub_model_factory_t), intent(in) :: this
    !> JSON core
    type(json_core), pointer, intent(in) :: json
    !> JSON object
    type(json_value), pointer, intent(in) :: j_obj

    type(json_value), pointer :: child
    character(kind=json_ck, len=:), allocatable :: key, unicode_type_name
    character(len=:), allocatable :: type_name
    logical(kind=json_lk) :: found

    new_obj => null()

    ! Get the sub-model type
    call json%get(j_obj, "type", unicode_type_name, found)
    call assert_msg(447218460, found, 'Missing sub-model type.')
    type_name = unicode_type_name
    
    ! Create a new sub-model instance of the type specified
    new_obj => this%create(type_name)    

    ! Load sub-model parameters from the json object
    call new_obj%load(json, j_obj)

#else
  !> Generic warning function when no input file support exists
  function load(this) result (new_obj)

    !> A new sub-model
    class(sub_model_data_t), pointer :: new_obj
    !> Sub-model factory
    class(sub_model_factory_t), intent(in) :: this

    new_obj => null()

    call warn_msg(545649418, "No support for input files.")
#endif
  end function load

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get the sub-model type as a constant
  integer(kind=i_kind) function get_type(this, sub_model) &
            result (sub_model_data_type)

    !> Sub-model factory
    class(sub_model_factory_t), intent(in) :: this
    !> Sub-model to get type of
    class(sub_model_data_t), intent(in) :: sub_model

    select type (sub_model)
      type is (sub_model_UNIFAC_t)
        sub_model_data_type = SUB_MODEL_UNIFAC
      class default
        call die_msg(695653684, "Unknown sub-model type")
    end select

  end function get_type

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Initialize an update data object
  subroutine initialize_update_data(this, update_data)

    !> Sub-model factory
    class(sub_model_factory_t), intent(in) :: this
    !> Update data object
    class(sub_model_update_data_t), intent(out) :: update_data

    select type (update_data)
      class default
        call die_msg(245232793, "Internal error - update data type missing")
    end select

  end subroutine initialize_update_data

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Determine the size of a binary required to pack a sub-model
  integer(kind=i_kind) function pack_size(this, sub_model)

    !> Sub-model factory
    class(sub_model_factory_t) :: this
    !> Sub-model to pack
    class(sub_model_data_t), intent(in) :: sub_model

    integer(kind=i_kind) :: i_sub_model

    pack_size =  pmc_mpi_pack_size_integer(int(1, kind=i_kind)) + &
                 sub_model%pack_size()

  end function pack_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Pack the given value to the buffer, advancing position
  subroutine bin_pack(this, sub_model, buffer, pos)

    !> Sub-model factory
    class(sub_model_factory_t), intent(in) :: this
    !> Sub-model to pack
    class(sub_model_data_t), intent(in) :: sub_model
    !> Memory buffer
    character, intent(inout) :: buffer(:)
    !> Current buffer position
    integer, intent(inout) :: pos

#ifdef PMC_USE_MPI
    integer :: sub_model_data_type, i_sub_model, prev_position

    prev_position = pos
    select type (sub_model)
      type is (sub_model_UNIFAC_t)
        sub_model_data_type = SUB_MODEL_UNIFAC
      class default
        call die_msg(850922257, "Trying to pack sub-model of unknown type.")
    end select
    call pmc_mpi_pack_integer(buffer, pos, sub_model_data_type)
    call sub_model%bin_pack(buffer, pos)
    call assert(340451545, &
         pos - prev_position <= this%pack_size(sub_model))
#endif

  end subroutine bin_pack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Unpack the given value to the buffer, advancing position
  function bin_unpack(this, buffer, pos) result (sub_model)

    !> Unpacked sub-model
    class(sub_model_data_t), pointer :: sub_model
    !> Sub-model factory
    class(sub_model_factory_t), intent(in) :: this
    !> Memory buffer
    character, intent(inout) :: buffer(:)
    !> Current buffer position
    integer, intent(inout) :: pos

#ifdef PMC_USE_MPI
    integer :: sub_model_data_type, i_sub_model, prev_position

    prev_position = pos
    call pmc_mpi_unpack_integer(buffer, pos, sub_model_data_type)
    select case (sub_model_data_type)
      case (SUB_MODEL_UNIFAC)
        sub_model => sub_model_UNIFAC_t()
      case default
        call die_msg(786366152, "Trying to unpack sub-model of unknown "// &
                "type: "//trim(to_string(sub_model_data_type)))
    end select
    call sub_model%bin_unpack(buffer, pos)
    call assert(209433801, &
         pos - prev_position <= this%pack_size(sub_model))
#endif

  end function bin_unpack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_sub_model_factory
