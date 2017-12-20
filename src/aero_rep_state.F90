! Copyright (C) 2017 Matt Dawson
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_aero_rep_state module.

!> The abstract aero_rep_state_t structure and associated subroutines.
module pmc_aero_rep_state

  use pmc_constants,                  only : i_kind, dp
  use pmc_mpi
  use pmc_util,                       only : die_msg, string_t
#ifdef PMC_USE_MPI
  use mpi
#endif
#ifdef PMC_USE_JSON
  use json_module
#endif

  implicit none
  private

  public :: aero_rep_state_t, aero_rep_state_ptr

  !> Abstract aerosol representation state type
  !!
  !! State data related to an aerosol representation. Derived types extending 
  !! aero_rep_state_t should describe the state of specific types of aerosol
  !! schemes (e.g., binned, modal, particle-resolved).
  !!
  !! Depending on the scheme, message packing functions can be used to pass
  !! state data for the aerosol representation. Alternately, the aerosol
  !! state can be set by the host model at the beginning of each call to the
  !! chemistry module integration function
  type, abstract :: aero_rep_state_t
    private
  contains
    !> Determine the size of a binary required to pack a given variable
    procedure(pmc_aero_rep_state_pack_size), deferred :: pack_size 
    !> Pack the given value to the buffer, advancing position
    procedure(pmc_aero_rep_state_bin_pack), deferred :: bin_pack
    !> Unpack the given value from the buffer, advancing position
    procedure(pmc_aero_rep_state_bin_unpack), deferred :: bin_unpack
  end type aero_rep_state_t

  !> Pointer to an aero_rep_state_t extending type, for building arrays
  type aero_rep_state_ptr
    class(aero_rep_state_t), pointer :: val
  end type aero_rep_state_ptr

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Determine the size of a binary required to pack a given variable
  interface pmc_aero_rep_state_pack_size_if
    integer(kind=i_kind) function pmc_aero_rep_state_pack_size(this) &
                  result (pack_size)
      use pmc_util,                                     only: i_kind
      import :: aero_rep_state_t

      !> Aerosol representation state
      class(aero_rep_state_t), intent(in) :: this

    end function pmc_aero_rep_state_pack_size
  end interface pmc_aero_rep_state_pack_size_if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Pack the given value to the buffer, advancing position
  interface pmc_aero_rep_state_bin_pack_if
    subroutine pmc_aero_rep_state_bin_pack(this, buffer, pos)
      use pmc_util,                                     only: i_kind
      import :: aero_rep_state_t

      !> Aerosol representation state
      class(aero_rep_state_t), intent(in) :: this
      !> Memory buffer
      character, allocatable, intent(inout) :: buffer(:)
      !> Current buffer position
      integer(kind=i_kind), intent(inout) :: pos

     end subroutine pmc_aero_rep_state_bin_pack
   end interface pmc_aero_rep_state_bin_pack_if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Unpack the given value to the buffer, advancing position
  interface pmc_aero_rep_state_bin_unpack_if
    subroutine pmc_aero_rep_state_bin_unpack(this, buffer, pos)
      use pmc_util,                                     only: i_kind
      import :: aero_rep_state_t

      !> Aerosol representation state
      class(aero_rep_state_t), intent(inout) :: this
      !> Memory buffer
      character, allocatable, intent(inout) :: buffer(:)
      !> Current buffer position
      integer(kind=i_kind), intent(inout) :: pos

    end subroutine pmc_aero_rep_state_bin_unpack
  end interface pmc_aero_rep_state_bin_unpack_if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_aero_rep_state
