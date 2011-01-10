! Copyright (C) 2010-2011 Matthew West
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
  use pmc_netcdf
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
  !> Type code for MFA weighting.
  integer, parameter :: AERO_WEIGHT_TYPE_MFA      = 3

  !> An aerosol size distribution weighting function.
  type aero_weight_t
     !> Weight type (given by module constants).
     integer :: type
     !> Reference radius at which the weight is 1 (hard-coded at present).
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
    type(aero_weight_t), intent(inout) :: aero_weight

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
    elseif ((aero_weight%type == AERO_WEIGHT_TYPE_POWER) &
         .or. (aero_weight%type == AERO_WEIGHT_TYPE_MFA)) then
       aero_weight_value &
            = (radius / aero_weight%ref_radius)**aero_weight%exponent
    else
       call die_msg(700421478, "unknown aero_weight type: " &
            // trim(integer_to_string(aero_weight%type)))
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

    !> \page input_format_aero_weight Input File Format: Aerosol Weighting Function
    !!
    !! For efficiency the aerosol population can be weighted so that
    !! the true number distribution \f$n(D)\f$ is given by
    !! \f[ n(D) = w(D) c(D) \f]
    !! where \f$w(D)\f$ is a fixed weighting function, \f$c(D)\f$ is
    !! the computational (simulated) number distribution, and \f$D\f$
    !! is the diameter. Thus a large value of \f$w(D)\f$ means that
    !! relatively few computational particles are used at diameter
    !! \f$D\f$, while a small value of \f$w(D)\f$ means that
    !! relatively many computational particles will be used at that
    !! diameter.
    !!
    !! The aerosol weighting function is specified by the parameters:
    !!   - \b weight (string): the type of weighting function --- must
    !!     be one of: "none" for no weighting (\f$w(D) = 1\f$);
    !!     "power" for a power-law weighting (\f$w(D) = D^\alpha\f$),
    !!     or "mfa" for the mass flow algorithm weighting (\f$w(D) =
    !!     D^{-3}\f$ with dependent coagulation particle removal)
    !!   - if the \c weight is \c power then the next parameter is:
    !!     - \b exponent (real, dimensionless): the exponent
    !!       \f$\alpha\f$ in the power law relationship --- setting
    !!       the \c exponent to 0 is equivalent to no weighting, while
    !!       setting the exponent negative uses more computational
    !!       particles at larger diameters and setting the exponent
    !!       positive uses more computational particles at smaller
    !!       diameters; in practice exponents between 0 and -3 are
    !!       most useful
    !!
    !! See also:
    !!   - \ref spec_file_format --- the input file text format

    call spec_file_read_string(file, 'weight', weight_type)
    if (trim(weight_type) == 'none') then
       aero_weight%type = AERO_WEIGHT_TYPE_NONE
    elseif (trim(weight_type) == 'power') then
       aero_weight%type = AERO_WEIGHT_TYPE_POWER
       aero_weight%ref_radius = 1d0
       call spec_file_read_real(file, 'exponent', aero_weight%exponent)
    elseif (trim(weight_type) == 'mfa') then
       aero_weight%type = AERO_WEIGHT_TYPE_MFA
       aero_weight%ref_radius = 1d0
       aero_weight%exponent = -3d0
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
         position - prev_position <= pmc_mpi_pack_size_aero_weight(val))
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
    type(aero_weight_t), intent(inout) :: val

#ifdef PMC_USE_MPI
    integer :: prev_position

    prev_position = position
    call pmc_mpi_unpack_integer(buffer, position, val%type)
    call pmc_mpi_unpack_real(buffer, position, val%ref_radius)
    call pmc_mpi_unpack_real(buffer, position, val%exponent)
    call assert(874467577, &
         position - prev_position <= pmc_mpi_pack_size_aero_weight(val))
#endif

  end subroutine pmc_mpi_unpack_aero_weight

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Write full state.
  subroutine aero_weight_output_netcdf(aero_weight, ncid)

    !> Aero weight to write.
    type(aero_weight_t), intent(in) :: aero_weight
    !> NetCDF file ID, in data mode.
    integer, intent(in) :: ncid

    !> \page output_format_aero_weight Output File Format: Aerosol Weighting Function
    !!
    !! The aerosol weighting function NetCDF variables are:
    !!   - \b weight_type (dimensionless integer): the type of the
    !!     weighting function, with 0 = invalid weight, 1 = no weight
    !!     (\f$w(D) = 1\f$), 2 = power weight (\f$w(D) =
    !!     (D/D_0)^\alpha\f$), 3 = MFA weight (\f$w(D) =
    !!     (D/D_0)^{-3}\f$)
    !!   - \b weight_ref_radius (unit m): the reference radius
    !!     \f$R_0\f$ (corresponding to \f$D_0 = 2R_0\f$), set to
    !!     zero if \c weight_type does not need it
    !!   - \b weight_exponent (dimensionless): the exponent
    !!     \f$\alpha\f$ for the power \c weight_type, set to -3
    !!     for the MFA \c weight_type, and zero otherwise
    !!
    !! See also:
    !!   - \ref input_format_aero_weight --- the corresponding input format

    call pmc_nc_write_integer(ncid, aero_weight%type, "weight_type", &
         description="type of aerosol weighting function: 0 = invalid, " &
         // "1 = none (w(D) = 1), 2 = power (w(D) = (D/D_0)^alpha), " &
         // "3 = MFA (mass flow) (w(D) = (D/D_0)^(-3))")
    call pmc_nc_write_real(ncid, aero_weight%ref_radius, &
         "weight_ref_radius", unit="m", &
         description="reference radius R_0 corresponding to the " &
         // "reference diameter D_0 = 2 * R_0 if the weight_type " &
         // "requires it, otherwise zero")
    call pmc_nc_write_real(ncid, aero_weight%exponent, "weight_exponent", &
         unit="1", &
         description="exponent alpha for the power weight_type, " &
         // "set to -3 for MFA, and zero otherwise")

  end subroutine aero_weight_output_netcdf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read full state.
  subroutine aero_weight_input_netcdf(aero_weight, ncid)

    !> Environment state to read.
    type(aero_weight_t), intent(inout) :: aero_weight
    !> NetCDF file ID, in data mode.
    integer, intent(in) :: ncid

    call pmc_nc_read_integer(ncid, aero_weight%type, "weight_type")
    call pmc_nc_read_real(ncid, aero_weight%ref_radius, "weight_ref_radius")
    call pmc_nc_read_real(ncid, aero_weight%exponent, "weight_exponent")

  end subroutine aero_weight_input_netcdf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_aero_weight
