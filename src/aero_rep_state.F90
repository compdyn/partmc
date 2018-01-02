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
  !! \c aero_rep_state_t should describe the state, other than aerosol species
  !! concentrations, of specific types of aerosol schemes (e.g., binned,
  !! modal, particle-resolved). The only data that should be tracked using an
  !! \c aero_rep_state_t extending type are time-varying data that cannot be
  !! obtained or calculated from the species concentrations in the \c
  !! pmc_phlex_state::phlex_state_t::state_var array or the \c
  !! pmc_phlex_state::phlex_state_t::env_state object. In some cases (e.g.,
  !! particle-resolved), no information other than the species concentrations
  !! are required, so the \c aero_rep_state_t extending type does not include
  !! any member variables. 
  !!
  !! Depending on the scheme, message packing functions can be used to pass
  !! state data for the aerosol representation. Alternately, the aerosol
  !! state can be set by the host model at the prior to each call to
  !! \c pmc_phlex_core::phlex_core_t::solve()
  type, abstract :: aero_rep_state_t
    private
  contains
    !> Determine the size of a binary required to pack a given variable
    procedure(pack_size), deferred :: pack_size 
    !> Pack the given value to the buffer, advancing position
    procedure(bin_pack), deferred :: bin_pack
    !> Unpack the given value from the buffer, advancing position
    procedure(bin_unpack), deferred :: bin_unpack
  end type aero_rep_state_t

  !> Pointer to an aero_rep_state_t extending type, for building arrays
  type aero_rep_state_ptr
    class(aero_rep_state_t), pointer :: val
  end type aero_rep_state_ptr

interface

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Determine the size of a binary required to pack a given variable
  integer(kind=i_kind) function pack_size(this)
    use pmc_util,                                     only: i_kind
    import :: aero_rep_state_t

    !> Aerosol representation state
    class(aero_rep_state_t), intent(in) :: this

  end function pack_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Pack the given value to the buffer, advancing position
  subroutine bin_pack(this, buffer, pos)
    use pmc_util,                                     only: i_kind
    import :: aero_rep_state_t

    !> Aerosol representation state
    class(aero_rep_state_t), intent(in) :: this
    !> Memory buffer
    character, allocatable, intent(inout) :: buffer(:)
    !> Current buffer position
    integer(kind=i_kind), intent(inout) :: pos

   end subroutine bin_pack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Unpack the given value to the buffer, advancing position
  subroutine bin_unpack(this, buffer, pos)
    use pmc_util,                                     only: i_kind
    import :: aero_rep_state_t

    !> Aerosol representation state
    class(aero_rep_state_t), intent(inout) :: this
    !> Memory buffer
    character, allocatable, intent(inout) :: buffer(:)
    !> Current buffer position
    integer(kind=i_kind), intent(inout) :: pos

  end subroutine bin_unpack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end interface

end module pmc_aero_rep_state
