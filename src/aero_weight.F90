! Copyright (C) 2010-2015 Matthew West
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
  use pmc_aero_particle
  use pmc_aero_data
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
     !> Overall weight magnitude (m^{-3}).
     real(kind=dp) :: magnitude
     !> Exponent for "power" weight.
     real(kind=dp) :: exponent
  end type aero_weight_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Sets the \c aero_weight to a non-zero normalized value.
  elemental subroutine aero_weight_normalize(aero_weight)

    !> Aerosol weight.
    type(aero_weight_t), intent(inout) :: aero_weight

    aero_weight%magnitude = 1d0

  end subroutine aero_weight_normalize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Scale the weight by the given fraction, so <tt>new_weight =
  !> old_weight * factor</tt>.
  elemental subroutine aero_weight_scale(aero_weight, factor)

    !> Aerosol weight to scale.
    type(aero_weight_t), intent(inout) :: aero_weight
    !> Factor to scale by.
    real(kind=dp), intent(in) :: factor

    aero_weight%magnitude = aero_weight%magnitude * factor

  end subroutine aero_weight_scale

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Combine \c aero_weight_delta into \c aero_weight with a harmonic mean.
  elemental subroutine aero_weight_combine(aero_weight, aero_weight_delta)

    !> Aerosol weight to add into.
    type(aero_weight_t), intent(inout) :: aero_weight
    !> Aerosol weight to add from.
    type(aero_weight_t), intent(in) :: aero_weight_delta

    aero_weight%magnitude = harmonic_mean(aero_weight%magnitude, &
         aero_weight_delta%magnitude)

  end subroutine aero_weight_combine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Adjust source and destination weights to reflect moving \c
  !> sample_prop proportion of particles from \c aero_weight_from to
  !> \c aero_weight_to.
  elemental subroutine aero_weight_shift(aero_weight_from, &
       aero_weight_to, sample_prop, overwrite_to)

    !> Aerosol weight to shift from.
    type(aero_weight_t), intent(inout) :: aero_weight_from
    !> Aerosol weight to shift to.
    type(aero_weight_t), intent(inout) :: aero_weight_to
    !> Proportion of particles being transfered.
    real(kind=dp), intent(in) :: sample_prop
    !> Whether to overwrite the destination weight (default: no).
    logical, intent(in), optional :: overwrite_to

    real(kind=dp) :: magnitude_transfer
    logical :: use_overwrite_to

    magnitude_transfer = aero_weight_from%magnitude / sample_prop
    aero_weight_from%magnitude = harmonic_mean(aero_weight_from%magnitude, &
         - magnitude_transfer)
    use_overwrite_to = .false.
    if (present(overwrite_to)) then
       if (overwrite_to) then
          use_overwrite_to = .true.
       end if
    end if
    if (use_overwrite_to) then
       aero_weight_to%magnitude = magnitude_transfer
    else
       aero_weight_to%magnitude = harmonic_mean(aero_weight_to%magnitude, &
            magnitude_transfer)
    end if

  end subroutine aero_weight_shift

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Compute the number concentration at a given radius (m^{-3}).
  real(kind=dp) function aero_weight_num_conc_at_radius(aero_weight, radius)

    !> Aerosol weight.
    type(aero_weight_t), intent(in) :: aero_weight
    !> Radius to compute number concentration at (m).
    real(kind=dp), intent(in) :: radius

    aero_weight_num_conc_at_radius = 0d0
    if (aero_weight%type == AERO_WEIGHT_TYPE_NONE) then
       aero_weight_num_conc_at_radius = aero_weight%magnitude
    elseif ((aero_weight%type == AERO_WEIGHT_TYPE_POWER) &
         .or. (aero_weight%type == AERO_WEIGHT_TYPE_MFA)) then
       aero_weight_num_conc_at_radius = aero_weight%magnitude &
            * radius**aero_weight%exponent
    else
       call die_msg(700421478, "unknown aero_weight type: " &
            // trim(integer_to_string(aero_weight%type)))
    end if

  end function aero_weight_num_conc_at_radius

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Compute the number concentration for a particle (m^{-3}).
  real(kind=dp) function aero_weight_num_conc(aero_weight, aero_particle, &
       aero_data)

    !> Aerosol weight.
    type(aero_weight_t), intent(in) :: aero_weight
    !> Aerosol particle to compute number concentration for.
    type(aero_particle_t), intent(in) :: aero_particle
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data

    aero_weight_num_conc = aero_weight_num_conc_at_radius(aero_weight, &
         aero_particle_radius(aero_particle, aero_data))

  end function aero_weight_num_conc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Compute the radius at a given number concentration (m). Inverse
  !> of aero_weight_num_conc_at_radius().
  real(kind=dp) function aero_weight_radius_at_num_conc(aero_weight, num_conc)

    !> Aerosol weight.
    type(aero_weight_t), intent(in) :: aero_weight
    !> Number concentration to compute the radius at (m^{-3}).
    real(kind=dp), intent(in) :: num_conc

    aero_weight_radius_at_num_conc = 0d0
    if (aero_weight%type == AERO_WEIGHT_TYPE_NONE) then
       call die_msg(545688537, "cannot invert weight type 'none'")
    elseif ((aero_weight%type == AERO_WEIGHT_TYPE_POWER) &
         .or. (aero_weight%type == AERO_WEIGHT_TYPE_MFA)) then
       call assert_msg(902242996, aero_weight%exponent /= 0d0, &
            "cannot invert weight with zero exponent")
       aero_weight_radius_at_num_conc &
            = exp(log(num_conc / aero_weight%magnitude) / aero_weight%exponent)
    else
       call die_msg(307766845, "unknown aero_weight type: " &
            // trim(integer_to_string(aero_weight%type)))
    end if

  end function aero_weight_radius_at_num_conc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Ensures that a weight function exponent is valid.
  subroutine aero_weight_check_valid_exponent(aero_weight)

    !> Aerosol weight array.
    type(aero_weight_t), intent(in) :: aero_weight

    ! check we know about all the weight type
    call assert(269952052, (aero_weight%type == AERO_WEIGHT_TYPE_NONE) &
         .or. (aero_weight%type == AERO_WEIGHT_TYPE_POWER) &
         .or. (aero_weight%type == AERO_WEIGHT_TYPE_MFA))
    if (aero_weight%type == AERO_WEIGHT_TYPE_NONE) then
       call assert(853998284, aero_weight%exponent == 0d0)
    end if
    if (aero_weight%type == AERO_WEIGHT_TYPE_MFA) then
       call assert(829651126, aero_weight%exponent == -3d0)
    end if

  end subroutine aero_weight_check_valid_exponent

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Determine whether a weight function is monotone increasing,
  !> monotone decreasing, or neither.
  subroutine aero_weight_check_monotonicity(aero_weight, &
       monotone_increasing, monotone_decreasing)

    !> Aerosol weight array.
    type(aero_weight_t), intent(in) :: aero_weight
    !> Whether all weights are monotone increasing.
    logical, intent(out) :: monotone_increasing
    !> Whether all weights are monotone decreasing.
    logical, intent(out) :: monotone_decreasing

    ! check we know about all the weight types
    call assert(610698264, (aero_weight%type == AERO_WEIGHT_TYPE_NONE) &
         .or. (aero_weight%type == AERO_WEIGHT_TYPE_POWER) &
         .or. (aero_weight%type == AERO_WEIGHT_TYPE_MFA))
    monotone_increasing = .true.
    monotone_decreasing = .true.
    if (aero_weight%type == AERO_WEIGHT_TYPE_POWER) then
       if (aero_weight%exponent < 0d0) then
          monotone_increasing = .false.
       end if
       if (aero_weight%exponent > 0d0) then
          monotone_decreasing = .false.
       end if
    end if
    if (aero_weight%type == AERO_WEIGHT_TYPE_MFA) then
       monotone_increasing = .false.
    end if

  end subroutine aero_weight_check_monotonicity

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
    aero_weight%magnitude = 0d0
    if (trim(weight_type) == 'none') then
       aero_weight%type = AERO_WEIGHT_TYPE_NONE
    elseif (trim(weight_type) == 'power') then
       aero_weight%type = AERO_WEIGHT_TYPE_POWER
       call spec_file_read_real(file, 'exponent', aero_weight%exponent)
    elseif (trim(weight_type) == 'mfa') then
       aero_weight%type = AERO_WEIGHT_TYPE_MFA
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
         + pmc_mpi_pack_size_real(val%magnitude) &
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
    call pmc_mpi_pack_real(buffer, position, val%magnitude)
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
    call pmc_mpi_unpack_real(buffer, position, val%magnitude)
    call pmc_mpi_unpack_real(buffer, position, val%exponent)
    call assert(639056899, &
         position - prev_position <= pmc_mpi_pack_size_aero_weight(val))
#endif

  end subroutine pmc_mpi_unpack_aero_weight

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_aero_weight
