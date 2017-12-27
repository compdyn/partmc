! Copyright (C) 2017 Matt Dawson
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_aero_rep_single_particle_state module.

!> The abstract aero_rep_single_particle_state_t structure and associated subroutines.
module pmc_aero_rep_single_particle_state

  use pmc_util,                                      only : dp, i_kind, &
                                                            string_t
  use pmc_aero_rep_state

  implicit none
  private

  public :: aero_rep_single_particle_state_t

  !> Single particle aerosol representation state
  !!
  !! State of a single particle
  type, extends(aero_rep_state_t) :: aero_rep_single_particle_state_t
  contains
    !> Determine the size of a binary required to pack a given variable
    procedure :: pack_size => pmc_aero_rep_single_particle_state_pack_size
    !> Pack the given value to the buffer, advancing position
    procedure :: bin_pack => pmc_aero_rep_single_particle_state_bin_pack
    !> Unpack the given value from the buffer, advancing position
    procedure :: bin_unpack => pmc_aero_rep_single_particle_state_bin_unpack
  end type aero_rep_single_particle_state_t

  !> Constructor for aero_rep_single_particle_state_t
  interface aero_rep_single_particle_state_t
    procedure :: pmc_aero_rep_single_particle_state_constructor
  end interface aero_rep_single_particle_state_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Constructor for aero_rep_single_particle_state_t
  function pmc_aero_rep_single_particle_state_constructor() result (new_obj)

    !> New aerosol state 
    type(aero_rep_single_particle_state_t), pointer :: new_obj

    allocate(new_obj)

  end function pmc_aero_rep_single_particle_state_constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Determine the size of a binary required to pack a given variable
  integer(kind=i_kind) function pmc_aero_rep_single_particle_state_pack_size(this) &
                  result (pack_size)

    !> Aerosol representation state
    class(aero_rep_single_particle_state_t), intent(in) :: this

    pack_size = 0

  end function pmc_aero_rep_single_particle_state_pack_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Pack the given value to the buffer, advancing position
  subroutine pmc_aero_rep_single_particle_state_bin_pack(this, buffer, pos)

    !> Aerosol representation state
    class(aero_rep_single_particle_state_t), intent(in) :: this
    !> Memory buffer
    character, allocatable, intent(inout) :: buffer(:)
    !> Current buffer position
    integer(kind=i_kind), intent(inout) :: pos

   end subroutine pmc_aero_rep_single_particle_state_bin_pack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Unpack the given value to the buffer, advancing position
  subroutine pmc_aero_rep_single_particle_state_bin_unpack(this, buffer, pos)

    !> Aerosol representation state
    class(aero_rep_single_particle_state_t), intent(inout) :: this
    !> Memory buffer
    character, allocatable, intent(inout) :: buffer(:)
    !> Current buffer position
    integer(kind=i_kind), intent(inout) :: pos

  end subroutine pmc_aero_rep_single_particle_state_bin_unpack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_aero_rep_single_particle_state
