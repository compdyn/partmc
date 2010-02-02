! Copyright (C) 2010 Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_aero_weight module.

!> The aero_weight_t structure and associated subroutines.
module pmc_aero_weight

  use pmc_util
  use pmc_constants
  use pmc_rand
  use pmc_spec_file
  use pmc_mpi
#ifdef PMC_USE_MPI
  use mpi
#endif

  !> Type code for an undefined or invalid weighting.
  integer, parameter :: AERO_WEIGHT_TYPE_INVALID  = 0
  !> Type code for no (or flat) weighting.
  integer, parameter :: AERO_WEIGHT_TYPE_NONE     = 1
  !> Type code for power function weighting.
  integer, parameter :: AERO_WEIGHT_TYPE_POWER    = 2

  !> An aerosol size distribution weighting function.
  type aero_weight_t
     !> Weight type (given by module constants).
     integer :: type
     !> Reference radius at which the weight is 1.
     real(kind=dp) :: ref_radius
     !> Exponent for "power" weight.
     real(kind=dp) :: exponent
  end type aero_weight_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Allocates an aero_weight.
  subroutine aero_weight_allocate(aero_weight)

    !> Aerosol weight.
    type(aero_weight_t), intent(out) :: aero_weight

    call aero_weight_zero(aero_weight)

  end subroutine aero_weight_allocate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Allocates an aero_weight of the given size.
  subroutine aero_weight_allocate_size(aero_weight)

    !> Aerosol weight.
    type(aero_weight_t), intent(out) :: aero_weight

    call aero_weight_zero(aero_weight)

  end subroutine aero_weight_allocate_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Free all storage.
  subroutine aero_weight_deallocate(aero_weight)

    !> Aerosol weight.
    type(aero_weight_t), intent(in) :: aero_weight

  end subroutine aero_weight_deallocate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Zeros the contents of the \c aero_weight.
  subroutine aero_weight_zero(aero_weight)

    !> Aerosol weight.
    type(aero_weight_t), intent(out) :: aero_weight

    aero_weight%type = AERO_WEIGHT_TYPE_INVALID
    aero_weight%ref_radius = 0d0
    aero_weight%exponent = 0d0

  end subroutine aero_weight_zero

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Copy an aero_weight.
  subroutine aero_weight_copy(aero_weight_from, aero_weight_to)

    !> Aerosol weight original.
    type(aero_weight_t), intent(in) :: aero_weight_from
    !> Aerosol weight copy.
    type(aero_weight_t), intent(inout) :: aero_weight_to

    aero_weight_to%type = aero_weight_from%type
    aero_weight_to%ref_radius = aero_weight_from%ref_radius
    aero_weight_to%exponent = aero_weight_from%exponent

  end subroutine aero_weight_copy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Give the value of an aero_weight at a specific radius.
  real(kind=dp) function aero_weight_value(aero_weight, radius)

    !> Aerosol weight.
    type(aero_weight_t), intent(in) :: aero_weight
    !> Radius to compute weight at (m).
    real(kind=dp), intent(in) :: radius

    if (aero_weight%type == AERO_WEIGHT_TYPE_NONE) then
       aero_weight_value = 1d0
    elseif (aero_weight%type == AERO_WEIGHT_TYPE_POWER) then
       aero_weight_value &
            = (radius / aero_weight%ref_radius)**aero_weight%exponent
    else
       call die_msg(700421478, "unknown aero_weight type: " &
            // integer_to_string(aero_weight%type))
    end if

  end function aero_weight_value

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read an aero_weight from a spec file.
  subroutine spec_file_read_aero_weight(file, aero_weight)

    !> Spec file.
    type(spec_file_t), intent(inout) :: file
    !> Aerosol weight.
    type(aero_weight_t), intent(inout) :: aero_weight

    character(len=SPEC_LINE_MAX_VAR_LEN) :: weight_type

    call spec_file_read_string(file, 'weight', weight_type)
    if (trim(weight_type) == 'none') then
       aero_weight%type = AERO_WEIGHT_TYPE_NONE
    elseif (trim(weight_type) == 'power') then
       aero_weight%type = AERO_WEIGHT_TYPE_POWER
       call spec_file_read_real(file, 'ref_radius', aero_weight%ref_radius)
       call spec_file_read_real(file, 'exponent', aero_weight%exponent)
    else
       call spec_file_die_msg(456342050, file, "unknown weight_type: " &
            // trim(weight_type))
    end if

  end subroutine spec_file_read_aero_weight

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Determines the number of bytes required to pack the given value.
  integer function pmc_mpi_pack_size_aero_weight(val)

    !> Value to pack.
    type(aero_weight_t), intent(in) :: val

    pmc_mpi_pack_size_aero_weight = &
         pmc_mpi_pack_size_integer(val%type) &
         + pmc_mpi_pack_size_real(val%ref_radius) &
         + pmc_mpi_pack_size_real(val%exponent)

  end function pmc_mpi_pack_size_aero_weight

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Packs the given value into the buffer, advancing position.
  subroutine pmc_mpi_pack_aero_weight(buffer, position, val)

    !> Memory buffer.
    character, intent(inout) :: buffer(:)
    !> Current buffer position.
    integer, intent(inout) :: position
    !> Value to pack.
    type(aero_weight_t), intent(in) :: val

#ifdef PMC_USE_MPI
    integer :: prev_position

    prev_position = position
    call pmc_mpi_pack_integer(buffer, position, val%type)
    call pmc_mpi_pack_real(buffer, position, val%ref_radius)
    call pmc_mpi_pack_real(buffer, position, val%exponent)
    call assert(579699255, &
         position - prev_position == pmc_mpi_pack_size_aero_weight(val))
#endif

  end subroutine pmc_mpi_pack_aero_weight

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Unpacks the given value from the buffer, advancing position.
  subroutine pmc_mpi_unpack_aero_weight(buffer, position, val)

    !> Memory buffer.
    character, intent(inout) :: buffer(:)
    !> Current buffer position.
    integer, intent(inout) :: position
    !> Value to pack.
    type(aero_weight_t), intent(out) :: val

#ifdef PMC_USE_MPI
    integer :: prev_position

    prev_position = position
    call pmc_mpi_unpack_integer(buffer, position, val%type)
    call pmc_mpi_unpack_real(buffer, position, val%ref_radius)
    call pmc_mpi_unpack_real(buffer, position, val%exponent)
    call assert(874467577, &
         position - prev_position == pmc_mpi_pack_size_aero_weight(val))
#endif

  end subroutine pmc_mpi_unpack_aero_weight

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_aero_weight
