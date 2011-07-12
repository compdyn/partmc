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
  use pmc_aero_particle
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
     !> Computational volume (m^3).
     real(kind=dp) :: comp_vol
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
  elemental subroutine aero_weight_allocate(aero_weight)

    !> Aerosol weight.
    type(aero_weight_t), intent(out) :: aero_weight

    call aero_weight_zero(aero_weight)

  end subroutine aero_weight_allocate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Allocates an aero_weight of the given size.
  elemental subroutine aero_weight_allocate_size(aero_weight)

    !> Aerosol weight.
    type(aero_weight_t), intent(out) :: aero_weight

    call aero_weight_zero(aero_weight)

  end subroutine aero_weight_allocate_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Free all storage.
  elemental subroutine aero_weight_deallocate(aero_weight)

    !> Aerosol weight.
    type(aero_weight_t), intent(in) :: aero_weight

  end subroutine aero_weight_deallocate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Zeros the contents of the \c aero_weight.
  elemental subroutine aero_weight_zero(aero_weight)

    !> Aerosol weight.
    type(aero_weight_t), intent(inout) :: aero_weight

    aero_weight%comp_vol = 0d0
    aero_weight%type = AERO_WEIGHT_TYPE_INVALID
    aero_weight%ref_radius = 0d0
    aero_weight%exponent = 0d0

  end subroutine aero_weight_zero

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Copy an aero_weight.
  elemental subroutine aero_weight_copy(aero_weight_from, aero_weight_to)

    !> Aerosol weight original.
    type(aero_weight_t), intent(in) :: aero_weight_from
    !> Aerosol weight copy.
    type(aero_weight_t), intent(inout) :: aero_weight_to

    aero_weight_to%comp_vol = aero_weight_from%comp_vol
    aero_weight_to%type = aero_weight_from%type
    aero_weight_to%ref_radius = aero_weight_from%ref_radius
    aero_weight_to%exponent = aero_weight_from%exponent

  end subroutine aero_weight_copy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Copy an aero_weight array.
  subroutine aero_weight_array_copy(aero_weight_array_from, &
       aero_weight_array_to)

    !> Aerosol weight array original.
    type(aero_weight_t), intent(in) :: aero_weight_array_from(:)
    !> Aerosol weight array copy.
    type(aero_weight_t), allocatable, intent(inout) :: aero_weight_array_to(:)

    if (size(aero_weight_array_to) /= size(aero_weight_array_from)) then
       call aero_weight_deallocate(aero_weight_array_to)
       deallocate(aero_weight_array_to)
       allocate(aero_weight_array_to(size(aero_weight_array_from)))
       call aero_weight_allocate(aero_weight_array_to)
    end if
    call aero_weight_copy(aero_weight_array_from, aero_weight_array_to)

  end subroutine aero_weight_array_copy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Scale the computational volume by the given fraction, so
  !> <tt>new_comp_vol = old_comp_vol * fraction</tt>.
  elemental subroutine aero_weight_scale_comp_vol(aero_weight, fraction)

    !> Aerosol weight to halve.
    type(aero_weight_t), intent(inout) :: aero_weight
    !> Fraction to scale computational volume by.
    real(kind=dp), intent(in) :: fraction

    aero_weight%comp_vol = aero_weight%comp_vol * fraction

  end subroutine aero_weight_scale_comp_vol

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Add the computational volume of \c aero_weight_delta to \c aero_weight.
  elemental subroutine aero_weight_add_comp_vol(aero_weight, aero_weight_delta)

    !> Aerosol weight to add volume to.
    type(aero_weight_t), intent(inout) :: aero_weight
    !> Aerosol weight to add volume from.
    type(aero_weight_t), intent(in) :: aero_weight_delta

    aero_weight%comp_vol = aero_weight%comp_vol + aero_weight_delta%comp_vol

  end subroutine aero_weight_add_comp_vol

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Compute the number concentration at a given radius (m^{-3}).
  real(kind=dp) function aero_weight_num_conc_at_radius(aero_weight, radius)

    !> Aerosol weight.
    type(aero_weight_t), intent(in) :: aero_weight
    !> Radius to compute number concentration at (m).
    real(kind=dp), intent(in) :: radius

    if (aero_weight%type == AERO_WEIGHT_TYPE_NONE) then
       aero_weight_num_conc_at_radius = 1d0
    elseif ((aero_weight%type == AERO_WEIGHT_TYPE_POWER) &
         .or. (aero_weight%type == AERO_WEIGHT_TYPE_MFA)) then
       aero_weight_num_conc_at_radius &
            = (radius / aero_weight%ref_radius)**aero_weight%exponent
    else
       call die_msg(700421478, "unknown aero_weight type: " &
            // trim(integer_to_string(aero_weight%type)))
    end if
    aero_weight_num_conc_at_radius &
         = aero_weight_num_conc_at_radius / aero_weight%comp_vol

  end function aero_weight_num_conc_at_radius

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Compute the number concentration for a particle (m^{-3}).
  real(kind=dp) function aero_weight_num_conc(aero_weight, &
       aero_particle)

    !> Aerosol weight.
    type(aero_weight_t), intent(in) :: aero_weight
    !> Aerosol particle to compute number concentration for.
    type(aero_particle_t), intent(in) :: aero_particle

    aero_weight_num_conc = aero_weight_num_conc_at_radius(aero_weight, &
         aero_particle_radius(aero_particle))

  end function aero_weight_num_conc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Compute the number concentration for a particle (m^{-3}).
  real(kind=dp) function aero_weight_array_single_num_conc(aero_weight_array, &
       aero_particle)

    !> Aerosol weight array.
    type(aero_weight_t), intent(in) :: aero_weight_array(:)
    !> Aerosol particle to compute number concentration for.
    type(aero_particle_t), intent(in) :: aero_particle

    aero_weight_array_single_num_conc = aero_weight_num_conc( &
         aero_weight_array(aero_particle%weight_group), aero_particle)

  end function aero_weight_array_single_num_conc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Compute the total number concentration at a given radius (m^3).
  real(kind=dp) function aero_weight_array_num_conc_at_radius( &
       aero_weight_array, radius)

    !> Aerosol weight array.
    type(aero_weight_t), intent(in) :: aero_weight_array(:)
    !> Radius to compute number concentration at (m).
    real(kind=dp), intent(in) :: radius

    integer :: i_group
    real(kind=dp) :: num_conc(size(aero_weight_array))

    do i_group = 1,size(aero_weight_array)
       num_conc(i_group) &
            = aero_weight_num_conc_at_radius(aero_weight_array(i_group), &
            radius)
    end do
    ! harmonic mean (same as summing the computational volumes)
    aero_weight_array_num_conc_at_radius = 1d0 / sum(1d0 / num_conc)

  end function aero_weight_array_num_conc_at_radius

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Compute the number concentration for a particle (m^{-3}).
  real(kind=dp) function aero_weight_array_num_conc(aero_weight_array, &
       aero_particle)

    !> Aerosol weight array.
    type(aero_weight_t), intent(in) :: aero_weight_array(:)
    !> Aerosol particle to compute number concentration for.
    type(aero_particle_t), intent(in) :: aero_particle

    aero_weight_array_num_conc = aero_weight_array_num_conc_at_radius( &
         aero_weight_array, aero_particle_radius(aero_particle))

  end function aero_weight_array_num_conc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Compute the radius at a given number concentration (m). Inverse
  !> of aero_weight_num_conc_at_radius().
  real(kind=dp) function aero_weight_radius_at_num_conc(aero_weight, num_conc)

    !> Aerosol weight.
    type(aero_weight_t), intent(in) :: aero_weight
    !> Number concentration to compute the radius at (m^{-3}).
    real(kind=dp), intent(in) :: num_conc

    if (aero_weight%type == AERO_WEIGHT_TYPE_NONE) then
       call die_msg(902242996, "cannot invert weight type 'none'")
    elseif ((aero_weight%type == AERO_WEIGHT_TYPE_POWER) &
         .or. (aero_weight%type == AERO_WEIGHT_TYPE_MFA)) then
       call assert_msg(902242996, aero_weight%exponent /= 0d0, &
            "cannot invert weight with zero exponent")
       aero_weight_radius_at_num_conc = exp(log(aero_weight%ref_radius) &
            + log(num_conc * aero_weight%comp_vol) / aero_weight%exponent)
    else
       call die_msg(307766845, "unknown aero_weight type: " &
            // trim(integer_to_string(aero_weight%type)))
    end if

  end function aero_weight_radius_at_num_conc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Check whether a given aero_weight array is flat in total.
  logical function aero_weight_array_check_flat(aero_weight_array)

    !> Aerosol weight array.
    type(aero_weight_t), intent(in) :: aero_weight_array(:)

    integer :: i_group

    ! check we know about all the weight types
    do i_group = 1,size(aero_weight_array)
       call assert(269952052, &
            (aero_weight_array(i_group)%type == AERO_WEIGHT_TYPE_NONE) &
            .or. (aero_weight_array(i_group)%type == AERO_WEIGHT_TYPE_POWER) &
            .or. (aero_weight_array(i_group)%type == AERO_WEIGHT_TYPE_MFA))
       if (aero_weight_array(i_group)%type == AERO_WEIGHT_TYPE_NONE) then
          call assert(829651126, aero_weight_array(i_group)%exponent == 0d0)
       end if
       if (aero_weight_array(i_group)%type == AERO_WEIGHT_TYPE_MFA) then
          call assert(829651126, aero_weight_array(i_group)%exponent == -3d0)
       end if
    end do

    if (abs(sum(aero_weight_array%exponent)) &
         < 1d-20 * sum(abs(aero_weight_array%exponent))) then
       aero_weight_array_check_flat = .true.
    else
       aero_weight_array_check_flat = .false.
    end if

  end function aero_weight_array_check_flat

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Determine whether all weight functions in an array are monotone
  !> increasing, monotone decreasing, or neither.
  subroutine aero_weight_array_check_monotonicity(aero_weight_array, &
       monotone_increasing, monotone_decreasing)

    !> Aerosol weight array.
    type(aero_weight_t), intent(in) :: aero_weight_array(:)
    !> Whether all weights are monotone increasing.
    logical, intent(out) :: monotone_increasing
    !> Whether all weights are monotone decreasing.
    logical, intent(out) :: monotone_decreasing

    integer :: i_group

    monotone_increasing = .true.
    monotone_decreasing = .true.
    do i_group = 1,size(aero_weight_array)
       ! check we know about all the weight types
       call assert(610698264, &
            (aero_weight_array(i_group)%type == AERO_WEIGHT_TYPE_NONE) &
            .or. (aero_weight_array(i_group)%type == AERO_WEIGHT_TYPE_POWER) &
            .or. (aero_weight_array(i_group)%type == AERO_WEIGHT_TYPE_MFA))
       if (aero_weight_array(i_group)%type == AERO_WEIGHT_TYPE_POWER) then
          if (aero_weight_array(i_group)%exponent < 0d0) then
             monotone_increasing = .false.
          end if
          if (aero_weight_array(i_group)%exponent > 0d0) then
             monotone_decreasing = .false.
          end if
       end if
       if (aero_weight_array(i_group)%type == AERO_WEIGHT_TYPE_MFA) then
          monotone_increasing = .false.
       end if
    end do

  end subroutine aero_weight_array_check_monotonicity

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Compute the maximum and minimum number concentrations between the
  !> given radii.
  subroutine aero_weight_array_minmax_num_conc(aero_weight_array, radius_1, &
       radius_2, num_conc_min, num_conc_max)

    !> Aerosol weight.
    type(aero_weight_t), intent(in) :: aero_weight_array(:)
    !> First radius.
    real(kind=dp), intent(in) :: radius_1
    !> Second radius.
    real(kind=dp), intent(in) :: radius_2
    !> Minimum number concentration.
    real(kind=dp), intent(out) :: num_conc_min
    !> Maximum number concentration.
    real(kind=dp), intent(out) :: num_conc_max

    real(kind=dp) :: num_conc_1, num_conc_2
    logical :: monotone_increasing, monotone_decreasing

    call aero_weight_array_check_monotonicity(aero_weight_array, &
         monotone_increasing, monotone_decreasing)
    call assert(857727714, monotone_increasing .or. monotone_decreasing)

    num_conc_1 = aero_weight_array_num_conc_at_radius(aero_weight_array, &
         radius_1)
    num_conc_2 = aero_weight_array_num_conc_at_radius(aero_weight_array, &
         radius_2)
    num_conc_min = min(num_conc_1, num_conc_2)
    num_conc_max = max(num_conc_1, num_conc_2)

  end subroutine aero_weight_array_minmax_num_conc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Choose a random group at the given radius, with probability
  !> proportional to group volume at that radius.
  integer function aero_weight_array_rand_group(aero_weight_array, radius)

    !> Aerosol weight array.
    type(aero_weight_t), intent(in) :: aero_weight_array(:)
    !> Radius to sample group at (m).
    real(kind=dp), intent(in) :: radius

    real(kind=dp) :: comp_vols(size(aero_weight_array))
    integer :: i

    do i = 1,size(aero_weight_array)
       comp_vols(i) &
            = 1d0 / aero_weight_num_conc_at_radius(aero_weight_array(i), &
            radius)
    end do
    aero_weight_array_rand_group = sample_cts_pdf(comp_vols)

  end function aero_weight_array_rand_group

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
         pmc_mpi_pack_size_real(val%comp_vol) &
         + pmc_mpi_pack_size_integer(val%type) &
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
    call pmc_mpi_pack_real(buffer, position, val%comp_vol)
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
    call pmc_mpi_unpack_real(buffer, position, val%comp_vol)
    call pmc_mpi_unpack_integer(buffer, position, val%type)
    call pmc_mpi_unpack_real(buffer, position, val%ref_radius)
    call pmc_mpi_unpack_real(buffer, position, val%exponent)
    call assert(874467577, &
         position - prev_position <= pmc_mpi_pack_size_aero_weight(val))
#endif

  end subroutine pmc_mpi_unpack_aero_weight

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Write the \c aero_weight dimension to the given NetCDF file if it
  !> is not already present and in any case return the associated
  !> dimid.
  subroutine aero_weight_netcdf_dim_aero_weight(aero_weight_array, ncid, &
       dimid_aero_weight)

    !> Aero_weight structure array.
    type(aero_weight_t), intent(in) :: aero_weight_array(:)
    !> NetCDF file ID, in data mode.
    integer, intent(in) :: ncid
    !> Dimid of the species dimension.
    integer, intent(out) :: dimid_aero_weight

    integer :: status, i_weight, n_weight
    integer :: varid_aero_weight
    integer :: aero_weight_centers(size(aero_weight_array))

    ! try to get the dimension ID
    status = nf90_inq_dimid(ncid, "aero_weight", dimid_aero_weight)
    if (status == NF90_NOERR) return
    if (status /= NF90_EBADDIM) call pmc_nc_check(status)

    ! dimension not defined, so define now define it
    call pmc_nc_check(nf90_redef(ncid))

    n_weight = size(aero_weight_array)

    call pmc_nc_check(nf90_def_dim(ncid, "aero_weight", n_weight, &
         dimid_aero_weight))
    call pmc_nc_check(nf90_def_var(ncid, "aero_weight", NF90_INT, &
         dimid_aero_weight, varid_aero_weight))
    call pmc_nc_check(nf90_put_att(ncid, varid_aero_weight, "description", &
         "dummy dimension variable (no useful value)"))

    call pmc_nc_check(nf90_enddef(ncid))

    do i_weight = 1,n_weight
       aero_weight_centers(i_weight) = i_weight
    end do
    call pmc_nc_check(nf90_put_var(ncid, varid_aero_weight, &
         aero_weight_centers))

  end subroutine aero_weight_netcdf_dim_aero_weight

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Write full state.
  subroutine aero_weight_array_output_netcdf(aero_weight_array, ncid)

    !> Aero weight to write.
    type(aero_weight_t), intent(in) :: aero_weight_array(:)
    !> NetCDF file ID, in data mode.
    integer, intent(in) :: ncid

    integer :: dimid_aero_weight

    !> \page output_format_aero_weight_array Output File Format: Aerosol Weighting Functions
    !!
    !! The aerosol weighting function NetCDF dimensions are:
    !!   - \b aero_weight: number of aerosol weighting functions
    !!
    !! The aerosol weighting function NetCDF variables are:
    !!   - \b weight_comp_vol (unit m^3, dim \c aero_weight): the
    !!     computational volume associated with each weighting
    !!     function
    !!   - \b weight_type (no unit, dim \c aero_weight): the type of
    !!     each weighting function, with 0 = invalid weight, 1 = no
    !!     weight (\f$w(D) = 1\f$), 2 = power weight (\f$w(D) =
    !!     (D/D_0)^\alpha\f$), 3 = MFA weight (\f$w(D) =
    !!     (D/D_0)^{-3}\f$)
    !!   - \b weight_exponent (no unit, dim \c aero_weight): for each
    !!     weighting function, specifies the exponent \f$\alpha\f$ for
    !!     the power \c weight_type, the value -3 for the MFA \c
    !!     weight_type, and zero for any other \c weight_type

    call aero_weight_netcdf_dim_aero_weight(aero_weight_array, ncid, &
         dimid_aero_weight)

    call pmc_nc_write_real_1d(ncid, aero_weight_array%comp_vol, &
         "weight_comp_vol", (/ dimid_aero_weight /), unit="m^3", &
         description="computational volume for each weighting function")
    call pmc_nc_write_integer_1d(ncid, aero_weight_array%type, "weight_type", &
         (/ dimid_aero_weight /), &
         description="type of each aerosol weighting function: 0 = invalid, " &
         // "1 = none (w(D) = 1), 2 = power (w(D) = (D/D_0)^alpha), " &
         // "3 = MFA (mass flow) (w(D) = (D/D_0)^(-3))")
    call pmc_nc_write_real_1d(ncid, aero_weight_array%exponent, &
         "weight_exponent", (/ dimid_aero_weight /), unit="1", &
         description="exponent alpha for the power weight_type, " &
         // "set to -3 for MFA, and zero otherwise")

  end subroutine aero_weight_array_output_netcdf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read full state.
  subroutine aero_weight_array_input_netcdf(aero_weight_array, ncid)

    !> Environment state to read.
    type(aero_weight_t), intent(inout), allocatable :: aero_weight_array(:)
    !> NetCDF file ID, in data mode.
    integer, intent(in) :: ncid

    integer :: dimid_aero_weight, n_weight
    character(len=1000) :: name
    real(kind=dp), allocatable :: comp_vol(:), exponent(:)
    integer, allocatable :: type(:)

    call pmc_nc_check(nf90_inq_dimid(ncid, "aero_weight", &
         dimid_aero_weight))
    call pmc_nc_check(nf90_Inquire_Dimension(ncid, &
         dimid_aero_weight, name, n_weight))
    call assert(719221386, n_weight < 1000)

    allocate(comp_vol(n_weight))
    allocate(type(n_weight))
    allocate(exponent(n_weight))

    call pmc_nc_read_real_1d(ncid, comp_vol, "weight_comp_vol")
    call pmc_nc_read_integer_1d(ncid, type, "weight_type")
    call pmc_nc_read_real_1d(ncid, exponent, "weight_exponent")

    call assert(309191498, size(comp_vol) == size(type))
    call assert(588649520, size(comp_vol) == size(exponent))

    call aero_weight_deallocate(aero_weight_array)
    deallocate(aero_weight_array)
    allocate(aero_weight_array(size(comp_vol)))
    call aero_weight_allocate(aero_weight_array)

    aero_weight_array%comp_vol = comp_vol
    aero_weight_array%type = type
    aero_weight_array%ref_radius = 1d0
    aero_weight_array%exponent = exponent

    deallocate(comp_vol)
    deallocate(type)
    deallocate(exponent)

  end subroutine aero_weight_array_input_netcdf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_aero_weight
