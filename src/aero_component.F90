! Copyright (C) 2018 Jeffrey Curtis
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_aero_particle module.

!> The aero_particle_t structure and associated subroutines.
module pmc_aero_component

  use pmc_util
  use pmc_mpi
#ifdef PMC_USE_MPI
  use mpi
#endif

  !> Aerosol particle component data structure.
  !!
  !!

  type aero_component_t
     !> Source.
     integer :: source_id
     !> Time the component was created (s).
     real(kind=dp) :: create_time     
  end type aero_component_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Determines the number of bytes required to pack the given value.
  integer function pmc_mpi_pack_size_aero_component(val)

    !> Value to pack.
    type(aero_component_t) :: val

    pmc_mpi_pack_size_aero_component = &
         pmc_mpi_pack_size_integer(val%source_id) &
         + pmc_mpi_pack_size_real(val%create_time)

  end function pmc_mpi_pack_size_aero_component

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Packs the given value into the buffer, advancing position.
  subroutine pmc_mpi_pack_aero_component(buffer, position, val)

    !> Memory buffer.
    character, intent(inout) :: buffer(:)
    !> Current buffer position.
    integer, intent(inout) :: position
    !> Value to pack.
    type(aero_component_t), intent(in) :: val

#ifdef PMC_USE_MPI
    integer :: prev_position

    prev_position = position
    call pmc_mpi_pack_integer(buffer, position, val%source_id)
    call pmc_mpi_pack_real(buffer, position, val%create_time)
    call assert(297427248, position - prev_position &
         <= pmc_mpi_pack_size_aero_component(val))
#endif

  end subroutine pmc_mpi_pack_aero_component

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Unpacks the given value from the buffer, advancing position.
  subroutine pmc_mpi_unpack_aero_component(buffer, position, val)

    !> Memory buffer.
    character, intent(inout) :: buffer(:)
    !> Current buffer position.
    integer, intent(inout) :: position
    !> Value to pack.
    type(aero_component_t), intent(inout) :: val

#ifdef PMC_USE_MPI
    integer :: prev_position

    prev_position = position
    call pmc_mpi_unpack_integer(buffer, position, val%source_id)
    call pmc_mpi_unpack_real(buffer, position, val%create_time)
    call assert(297427277, position - prev_position &
         <= pmc_mpi_pack_size_aero_component(val))
#endif

  end subroutine pmc_mpi_unpack_aero_component

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_aero_component
