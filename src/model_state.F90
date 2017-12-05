! Copyright (C) 2017 Matt Dawson
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_model_state module.

!> The model_state_t structure and associated subroutines.
module pmc_model_state

  use pmc_mpi
  use pmc_util,                       only : die_msg, string_t
  use pmc_chem_spec_state
  use pmc_env_state
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
    !> Chemical species state
    type(chem_spec_state_t) :: chem_spec_state
    !> Environmental conditions
    type(env_state_t) :: env_state
    ! TODO add aerosol representation states array
  contains
    !> Determine the size of a binary required to pack a given variable
    procedure :: pack_size => pmc_model_state_pack_size
    !> Pack the given value to the buffer, advancing position
    procedure :: bin_pack => pmc_model_state_bin_pack
    !> Unpack the given value from the buffer, advancing position
    procedure :: bin_unpack => pmc_model_state_bin_unpack
  end type model_state_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Determine the size of a binary required to pack a given variable
  integer(kind=i_kind) function pmc_model_state_pack_size(this) &
                  result (pack_size)

    !> Chemical species states
    class(model_state_t), intent(in) :: this

    pack_size = &
            this%chem_spec_state%pack_size() + &
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
     call this%chem_spec_state%bin_pack(buffer, pos)
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
     call this%chem_spec_state%bin_pack(buffer, pos)
     call this%env_state%bin_pack(buffer, pos)
     call assert(820809161, &
             pos - prev_position <= this%pack_size())
#endif

   end subroutine pmc_model_state_bin_unpack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_model_state
