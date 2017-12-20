! Copyright (C) 2017 Matt Dawson
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_model_state module.

!> The model_state_t structure and associated subroutines.
module pmc_model_state

  use pmc_mpi
  use pmc_util,                       only : die_msg, string_t
  use pmc_env_state
  use pmc_aero_rep_state
  use pmc_rand
#ifdef PMC_USE_MPI
  use mpi
#endif

  implicit none
  private

  public :: model_state_t

  !> Model state
  !!
  !! Temporal state of the model
  type model_state_t
    !> State variable array. This array includes one entry for each
    !! variable whose state will be solved for during the mechanism
    !! integration.
    real(kind=dp), allocatable :: state_var(:)
    !> Aerosol representation states. This should contain time-varying
    !! information related to an aerosol representation. Each element
    !! in this array should correspond to an element with the same index
    !! in the model_data%aero_rep(:) array.
    type(aero_rep_state_ptr), allocatable :: aero_rep_state(:)
    !> Environmental conditions
    type(env_state_t) :: env_state
    !> Rxn phase being solved
    integer(kind=i_kind) :: rxn_phase = 0
    !> Identifier that can be used by expensive, repeated functions to
    !! determine whether the state has changed since the last time they
    !! were called
    character(len=PMC_UUID_LEN) :: uuid
  contains
    !> Reset the state identifier
    procedure :: reset_id => pmc_model_state_reset_id
    !> Get the unique state identifier
    procedure :: get_id => pmc_model_state_get_id
    !> Determine the size of a binary required to pack a given variable
    procedure :: pack_size => pmc_model_state_pack_size
    !> Pack the given value to the buffer, advancing position
    procedure :: bin_pack => pmc_model_state_bin_pack
    !> Unpack the given value from the buffer, advancing position
    procedure :: bin_unpack => pmc_model_state_bin_unpack
  end type model_state_t

  !> Constructor for model_state_t
  interface model_state_t
    procedure :: pmc_model_state_constructor
  end interface model_state_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Constructor for model_state_t
  function pmc_model_state_constructor() result (new_obj)

    !> New model state
    type(model_state_t), pointer :: new_obj

    allocate(new_obj)
    call pmc_srand(0,0)

  end function pmc_model_state_constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Reset the state identifier. This should be called every time the state 
  !! changes.
  subroutine pmc_model_state_reset_id(this)

    !> Model state
    class(model_state_t), intent(inout) :: this

    call uuid4_str(this%uuid)

  end subroutine pmc_model_state_reset_id

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get the unique state identifier. This will change every time the state
  !! changes, and can be used to avoid duplicating expensive calculations in
  !! the mechanisms
  character(len=PMC_UUID_LEN) function pmc_model_state_get_id(this)

    !> Model state
    class(model_state_t), intent(in) :: this

    pmc_model_state_get_id = this%uuid 

  end function pmc_model_state_get_id

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Determine the size of a binary required to pack a given variable
  integer(kind=i_kind) function pmc_model_state_pack_size(this) &
                  result (pack_size)

    !> Chemical species states
    class(model_state_t), intent(in) :: this

    pack_size = &
            this%env_state%pack_size()

  end function pmc_model_state_pack_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Pack the given value to the buffer, advancing position
  subroutine pmc_model_state_bin_pack(this, buffer, pos)

    !> Chemical species states
    class(model_state_t), intent(in) :: this
    !> Memory buffer
    character, intent(inout) :: buffer(:)
    !> Current buffer position
    integer, intent(inout) :: pos

#ifdef PMC_USE_MPI
     integer :: prev_position

     prev_position = pos
     call this%env_state%bin_pack(buffer, pos)
     call assert(204737650, &
             pos - prev_position <= this%pack_size())
#endif

   end subroutine pmc_model_state_bin_pack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Unpack the given value to the buffer, advancing position
  subroutine pmc_model_state_bin_unpack(this, buffer, pos)

    !> Chemical species states
    class(model_state_t), intent(inout) :: this
    !> Memory buffer
    character, intent(inout) :: buffer(:)
    !> Current buffer position
    integer, intent(inout) :: pos

#ifdef PMC_USE_MPI
     integer :: prev_position

     prev_position = pos
     call this%env_state%bin_pack(buffer, pos)
     call assert(820809161, &
             pos - prev_position <= this%pack_size())
#endif

   end subroutine pmc_model_state_bin_unpack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_model_state
