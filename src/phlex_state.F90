! Copyright (C) 2017 Matt Dawson
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_phlex_state module.

!> The phlex_state_t structure and associated subroutines.
module pmc_phlex_state

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

  public :: phlex_state_t

  !> Model state
  !!
  !! Temporal state of the model
  type phlex_state_t
    !> State variable array. This array includes one entry for each
    !! variable whose state will be solved for during the mechanism
    !! integration.
    real(kind=dp), allocatable :: state_var(:)
    !> Aerosol representation states. This should contain time-varying
    !! information related to an aerosol representation. Each element
    !! in this array should correspond to an element with the same index
    !! in the phlex_core%aero_rep(:) array.
    type(aero_rep_state_ptr), allocatable :: aero_rep_state(:)
    !> Environmental conditions
    type(env_state_t), pointer :: env_state
    !> Rxn phase being solved
    integer(kind=i_kind) :: rxn_phase = 0
    !> Identifier that can be used by expensive, repeated functions to
    !! determine whether the state has changed since the last time they
    !! were called
    character(len=PMC_UUID_LEN) :: uuid
  contains
    !> Reset the state identifier
    procedure :: reset_id
    !> Get the unique state identifier
    procedure :: get_id
    !> Determine the size of a binary required to pack a given variable
    procedure :: pack_size
    !> Pack the given value to the buffer, advancing position
    procedure :: bin_pack
    !> Unpack the given value from the buffer, advancing position
    procedure :: bin_unpack
  end type phlex_state_t

  ! Constructor for phlex_state_t
  interface phlex_state_t
    procedure :: constructor
  end interface phlex_state_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Constructor for phlex_state_t
  function constructor(env_state) result (new_obj)

    !> New model state
    type(phlex_state_t), pointer :: new_obj
    !> Environmental state
    type(env_state_t), target, intent(in), optional :: env_state

    allocate(new_obj)
    if (present(env_state)) then
      new_obj%env_state => env_state
    else
      allocate(new_obj%env_state)
    end if
    call pmc_srand(0,0)

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Reset the state identifier. This should be called every time the state 
  !! changes.
  subroutine reset_id(this)

    !> Model state
    class(phlex_state_t), intent(inout) :: this

    call uuid4_str(this%uuid)

  end subroutine reset_id

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get the unique state identifier. This will change every time the state
  !! changes, and can be used to avoid duplicating expensive calculations in
  !! the mechanisms
  character(len=PMC_UUID_LEN) function get_id(this)

    !> Model state
    class(phlex_state_t), intent(in) :: this

    get_id = this%uuid 

  end function get_id

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Determine the size of a binary required to pack a given variable
  integer(kind=i_kind) function pack_size(this)

    !> Chemical species states
    class(phlex_state_t), intent(in) :: this

    pack_size = &
            this%env_state%pack_size()

  end function pack_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Pack the given value to the buffer, advancing position
  subroutine bin_pack(this, buffer, pos)

    !> Chemical species states
    class(phlex_state_t), intent(in) :: this
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

   end subroutine bin_pack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Unpack the given value to the buffer, advancing position
  subroutine bin_unpack(this, buffer, pos)

    !> Chemical species states
    class(phlex_state_t), intent(inout) :: this
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

   end subroutine bin_unpack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_phlex_state
