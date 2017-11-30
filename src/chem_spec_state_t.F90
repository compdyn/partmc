! Copyright (C) 2017 Matt Dawson
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_chem_spec_state module.

!> The chem_spec_state_t structure and associated subroutines.
module pmc_chem_spec_state

  use pmc_constants,                  only : i_kind
  use pmc_mpi
  use pmc_util,                       only : die_msg, string_t
#ifdef PMC_USE_MPI
  use mpi
#endif

  implicit none
  private

  public :: chem_spec_state_t

  !> Chemical species data
  !!
  !! Time-invariant data related to a chemical species
  type chem_spec_state_t
    !> Concentration TODO Determine how to handle units
    real(kind=dp), allocatable :: conc(:)
  contains
    !> Determine the size of a binary required to pack a given variable
    procedure :: pack_size => pmc_chem_spec_state_pack_size
    !> Pack the given value to the buffer, advancing position
    procedure :: bin_pack => pmc_chem_spec_state_bin_pack
    !> Unpack the given value from the buffer, advancing position
    procedure :: bin_unpack => pmc_chem_spec_state_bin_unpack
  end type chem_spec_state_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Determine the size of a binary required to pack a given variable
  integer(kind=i_kind) function pmc_chem_spec_state_pack_size(this) &
                  result (pack_size)

    !> Chemical species states
    class(chem_spec_state_t), intent(in) :: this

    pack_size = &
            pmc_mpi_pack_size_real_array(this%conc)

  end function pmc_chem_spec_state_pack_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Pack the given value to the buffer, advancing position
  subroutine pmc_chem_spec_state_bin_pack(this, buffer, pos)

    !> Chemical species states
    class(chem_spec_state_t), intent(in) :: this
    !> Memory buffer
    character, intent(inout) :: buffer(:)
    !> Current buffer position
    integer, intent(inout) :: pos

 #ifdef PMC_USE_MPI
     integer :: prev_position

     prev_position = pos
     call pmc_mpi_pack_real_array(buffer, pos, this%conc)
     call assert(361046754, &
             pos - prev_position <= this%pack_size())
 #endif

   end subroutine pmc_chem_spec_state_bin_pack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Unpack the given value to the buffer, advancing position
  subroutine pmc_chem_spec_state_bin_unpack(this, buffer, pos)

    !> Chemical species states
    class(chem_spec_state_t), intent(inout) :: this
    !> Memory buffer
    character, intent(in) :: buffer(:)
    !> Current buffer position
    integer, intent(inout) :: pos

 #ifdef PMC_USE_MPI
     integer :: prev_position

     prev_position = pos
     call pmc_mpi_unpack_real_array(buffer, pos, this%conc)
     call assert(184624599, &
             pos - prev_position <= this%pack_size())
 #endif

   end subroutine pmc_chem_spec_state_bin_unpack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_chem_spec_state
