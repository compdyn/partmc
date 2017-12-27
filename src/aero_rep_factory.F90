! Copyright (C) 2017 Matt Dawson
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_aero_rep_factory module.

!> The abstract aero_rep_factory_t structure and associated subroutines.
module pmc_aero_rep_factory

  use pmc_constants,                  only : i_kind, dp
  use pmc_util,                       only : die_msg, string_t, assert_msg, &
                                             warn_msg
  use pmc_aero_rep_data
  use pmc_mpi
#ifdef PMC_USE_JSON
  use json_module
#endif
#ifdef PMC_USE_MPI
  use mpi
#endif

  ! Use all aerosol representation modules
  use pmc_aero_rep_single_particle

  implicit none
  private

  public :: aero_rep_factory_t

  !> Identifiers for aerosol representations - used by binary packing/unpacking 
  !! functions
  integer(kind=i_kind), parameter :: AERO_REP_SINGLE_PARTICLE = 1

  !> Factory type for aerosol representations
  !!
  !! Provides new instances of types extending aero_rep_data_t by name or
  !! from input file data
  type :: aero_rep_factory_t
  contains
    !> Create a new aerosol representation by type name
    procedure :: create
    !> Create a new aerosol representation from input data
    procedure :: load
    !> Determine the number of bytes required to pack a given aerosol representation
    procedure :: pack_size
    !> Pack a given aerosol representation to the buffer, advancing the position
    procedure :: bin_pack
    !> Unpack a aerosol representation from the buffer, advancing the position
    procedure :: bin_unpack
  end type aero_rep_factory_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Create a new aerosol representation by type name
  function create(this, type_name) result (new_obj)

    !> An aerosol representation
    class(aero_rep_data_t), pointer :: new_obj
    !> Aerosol representation factory
    class(aero_rep_factory_t), intent(in) :: this
    !> Name of the aerosol representation
    character(len=:), allocatable :: type_name

    new_obj => null()

    select case (type_name)
      case ("AERO_REP_SINGLE_PARTICLE")
        new_obj => aero_rep_single_particle_t()
      case default
        call die_msg(792930166, "Unknown aerosol representation type: " &
                //type_name) 
    end select

  end function create

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Load an aerosol represenation based on its type
#ifdef PMC_USE_JSON
  function load(this, json, j_obj) result (new_obj)

    !> An aerosol representation
    class(aero_rep_data_t), pointer :: new_obj
    !> Aerosol representation factory
    class(aero_rep_factory_t), intent(in) :: this
    !> JSON core
    type(json_core), pointer, intent(in) :: json
    !> JSON object
    type(json_value), pointer, intent(in) :: j_obj

    type(json_value), pointer :: child
    character(kind=json_ck, len=:), allocatable :: key, unicode_type_name
    character(len=:), allocatable :: type_name
    logical(kind=json_lk) :: found

    new_obj => null()

    call json%get(j_obj, "type", unicode_type_name, found)
    call assert_msg(904082395, found, 'Missing aerosol representation type.')
    type_name = unicode_type_name
    new_obj => this%create(type_name)    
    call new_obj%load(json, j_obj)
#else
  function load(this) result (new_obj)

    !> An aerosol representation
    class(aero_rep_data_t), pointer :: new_obj
    !> Aerosol representation factory
    class(aero_rep_factory_t), intent(in) :: this

    new_obj => null()

    call warn_msg(723960750, "No support for input files.")
#endif
  end function load

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Determine the size of a binary required to pack a aerosol representation
  integer(kind=i_kind) function pack_size(this, aero_rep)

    !> Aerosol representation factory
    class(aero_rep_factory_t) :: this
    !> Aerosol representation to pack
    class(aero_rep_data_t), intent(in) :: aero_rep

    integer(kind=i_kind) :: i_aero_rep

    pack_size =  pmc_mpi_pack_size_integer(int(1, kind=i_kind)) + &
                 aero_rep%pack_size()

  end function pack_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Pack the given value to the buffer, advancing position
  subroutine bin_pack(this, aero_rep, buffer, pos)

    !> Aerosol representation factory
    class(aero_rep_factory_t), intent(in) :: this
    !> Aerosol representation to pack
    class(aero_rep_data_t), intent(in) :: aero_rep
    !> Memory buffer
    character, intent(inout) :: buffer(:)
    !> Current buffer position
    integer, intent(inout) :: pos

#ifdef PMC_USE_MPI
    integer :: aero_rep_type, i_aero_rep, prev_position

    prev_position = pos
    select type (aero_rep)
      type is (aero_rep_single_particle_t)
        aero_rep_type = AERO_REP_SINGLE_PARTICLE
      class default
        call die_msg(278244560, "Trying to pack aerosol representation of unknown type.")
    end select
    call pmc_mpi_pack_integer(buffer, pos, aero_rep_type)
    call aero_rep%bin_pack(buffer, pos)
    call assert(897674844, &
         pos - prev_position <= this%pack_size(aero_rep))
#endif

  end subroutine bin_pack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Unpack the given value to the buffer, advancing position
  function bin_unpack(this, buffer, pos) result (aero_rep)

    !> Unpacked aerosol representation
    class(aero_rep_data_t), pointer :: aero_rep
    !> Aerosol representation factory
    class(aero_rep_factory_t), intent(in) :: this
    !> Memory buffer
    character, intent(inout) :: buffer(:)
    !> Current buffer position
    integer, intent(inout) :: pos

#ifdef PMC_USE_MPI
    integer :: aero_rep_type, i_aero_rep, prev_position

    prev_position = pos
    call pmc_mpi_unpack_integer(buffer, pos, aero_rep_type)
    select case (aero_rep_type)
      case (AERO_REP_SINGLE_PARTICLE)
        aero_rep => aero_rep_single_particle_t()
      case default
        call die_msg(106634417, "Trying to unpack aerosol representation of unknown type:"// &
                to_string(aero_rep_type))
    end select
    call aero_rep%bin_unpack(buffer, pos)
    call assert(948795857, &
         pos - prev_position <= this%pack_size(aero_rep))
#endif

  end function bin_unpack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_aero_rep_factory
