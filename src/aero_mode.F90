! Copyright (C) 2005-2016 Nicole Riemer and Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_aero_mode module.

!> The aero_mode_t structure and associated subroutines.
module pmc_aero_mode

  use pmc_bin_grid
  use pmc_util
  use pmc_constants
  use pmc_spec_file
  use pmc_aero_data
  use pmc_aero_weight
  use pmc_mpi
  use pmc_rand
#ifdef PMC_USE_MPI
  use mpi
#endif

  !> Maximum length of a mode name.
  integer, parameter :: AERO_MODE_NAME_LEN = 300
  !> Maximum length of a mode type.
  integer, parameter :: AERO_MODE_TYPE_LEN = 20

  !> Type code for an undefined or invalid mode.
  integer, parameter :: AERO_MODE_TYPE_INVALID    = 0
  !> Type code for a log-normal mode.
  integer, parameter :: AERO_MODE_TYPE_LOG_NORMAL = 1
  !> Type code for an exponential mode.
  integer, parameter :: AERO_MODE_TYPE_EXP        = 2
  !> Type code for a mono-disperse mode.
  integer, parameter :: AERO_MODE_TYPE_MONO       = 3
  !> Type code for a sampled mode.
  integer, parameter :: AERO_MODE_TYPE_SAMPLED    = 4

  !> Type code for an undefined for invalid diameter type.
  integer, parameter :: AERO_MODE_DIAM_TYPE_INVALID   = 0
  !> Type code for geometric diameter.
  integer, parameter :: AERO_MODE_DIAM_TYPE_GEOMETRIC = 1
  !> Type code for mobility equivalent diameter.
  integer, parameter :: AERO_MODE_DIAM_TYPE_MOBILITY  = 2

  !> An aerosol size distribution mode.
  !!
  !! Each mode is assumed to be fully internally mixed so that every
  !! particle has the same composition. The composition is stored in
  !! \c vol_frac, while the other parameters define the size
  !! distribution (with \c type defining the type of size distribution
  !! function). See \ref input_format_mass_frac for descriptions of
  !! the parameters relevant to each mode type.
  type aero_mode_t
     !> Mode name, used to track particle sources.
     character(len=AERO_MODE_NAME_LEN) :: name
     !> Mode type (given by module constants).
     integer :: type
     !> Characteristic radius, with meaning dependent on mode type (m).
     real(kind=dp) :: char_radius
     !> Log base 10 of geometric standard deviation of radius, (m).
     real(kind=dp) :: log10_std_dev_radius
     !> Sample bin radii [length <tt>(N + 1)</tt>] (m).
     real(kind=dp), allocatable :: sample_radius(:)
     !> Sample bin number concentrations [length <tt>N</tt>] (m^{-3}).
     real(kind=dp), allocatable :: sample_num_conc(:)
     !> Total number concentration of mode (#/m^3).
     real(kind=dp) :: num_conc
     !> Species fractions by volume [length \c aero_data_n_spec(aero_data)] (1).
     real(kind=dp), allocatable :: vol_frac(:)
     !> Species fraction standard deviation
     !> [length \c aero_data_n_spec(aero_data)] (1).
     real(kind=dp), allocatable :: vol_frac_std(:)
     !> Source number.
     integer :: source
  end type aero_mode_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Return a string representation of a kernel type.
  character(len=AERO_MODE_TYPE_LEN) function aero_mode_type_to_string(type)

    !> Aero mode type.
    integer, intent(in) :: type

    if (type == AERO_MODE_TYPE_INVALID) then
       aero_mode_type_to_string = "invalid"
    elseif (type == AERO_MODE_TYPE_LOG_NORMAL) then
       aero_mode_type_to_string = "log_normal"
    elseif (type == AERO_MODE_TYPE_EXP) then
       aero_mode_type_to_string = "exp"
    elseif (type == AERO_MODE_TYPE_MONO) then
       aero_mode_type_to_string = "mono"
    elseif (type == AERO_MODE_TYPE_SAMPLED) then
       aero_mode_type_to_string = "sampled"
    else
       aero_mode_type_to_string = "unknown"
    end if

  end function aero_mode_type_to_string

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the total number concentration of a mode. (#/m^3)
  real(kind=dp) function aero_mode_total_num_conc(aero_mode)

    !> Aerosol mode.
    type(aero_mode_t), intent(in) :: aero_mode

    aero_mode_total_num_conc = 0d0
    if ((aero_mode%type == AERO_MODE_TYPE_LOG_NORMAL) &
         .or. (aero_mode%type == AERO_MODE_TYPE_EXP) &
       .or. (aero_mode%type == AERO_MODE_TYPE_MONO)) then
       aero_mode_total_num_conc = aero_mode%num_conc
    elseif (aero_mode%type == AERO_MODE_TYPE_SAMPLED) then
       aero_mode_total_num_conc = sum(aero_mode%sample_num_conc)
    else
       call die_msg(719625922, "unknown aero_mode type: " &
            // trim(integer_to_string(aero_mode%type)))
    end if

  end function aero_mode_total_num_conc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Compute a log-normal distribution.
  subroutine num_conc_log_normal(total_num_conc, geom_mean_radius, &
       log10_sigma_g, bin_grid, num_conc)

    !> Total number concentration of the mode (m^{-3}).
    real(kind=dp), intent(in) :: total_num_conc
    !> Geometric mean radius (m).
    real(kind=dp), intent(in) :: geom_mean_radius
    !> log_10(geom. std. dev.) (1).
    real(kind=dp), intent(in) :: log10_sigma_g
    !> Bin grid.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Number concentration (#(ln(r))d(ln(r))).
    real(kind=dp), intent(out) :: num_conc(bin_grid_size(bin_grid))

    integer :: k

    do k = 1,bin_grid_size(bin_grid)
       num_conc(k) = total_num_conc / (sqrt(2d0 * const%pi) &
            * log10_sigma_g) * dexp(-(dlog10(bin_grid%centers(k)) &
            - dlog10(geom_mean_radius))**2d0 &
            / (2d0 * log10_sigma_g**2d0)) / dlog(10d0)
    end do

    ! The formula above was originally for a distribution in
    ! log_10(r), while we are using log_e(r) for our bin grid. The
    ! division by dlog(10) at the end corrects for this.

    ! Remember that log_e(r) = log_10(r) * log_e(10).

  end subroutine num_conc_log_normal

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Compute a log-normal distribution in volume.
  subroutine vol_conc_log_normal(total_num_conc, &
       geom_mean_radius, log10_sigma_g, bin_grid, aero_data, vol_conc)

    !> Total number concentration of the mode (m^{-3}).
    real(kind=dp), intent(in) :: total_num_conc
    !> Geometric mean radius (m).
    real(kind=dp), intent(in) :: geom_mean_radius
    !> log_10(geom. std. dev.) (1).
    real(kind=dp), intent(in) :: log10_sigma_g
    !> Bin grid.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Volume concentration (V(ln(r))d(ln(r))).
    real(kind=dp), intent(out) :: vol_conc(bin_grid_size(bin_grid))

    real(kind=dp) :: num_conc(bin_grid_size(bin_grid))

    call num_conc_log_normal(total_num_conc, geom_mean_radius, &
         log10_sigma_g, bin_grid, num_conc)

    vol_conc = num_conc * aero_data_rad2vol(aero_data, bin_grid%centers)

  end subroutine vol_conc_log_normal

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Exponential distribution in volume.
  !!
  !! \f[ n(v) = \frac{1}{\rm mean-vol} \exp(- v / {\rm mean-vol}) \f]
  !! Normalized so that sum(num_conc(k) * log_width) = 1
  subroutine num_conc_exp(total_num_conc, radius_at_mean_vol, &
       bin_grid, aero_data, num_conc)

    !> Total number concentration of the mode (m^{-3}).
    real(kind=dp), intent(in) :: total_num_conc
    !> Radius at mean volume (m).
    real(kind=dp), intent(in) :: radius_at_mean_vol
    !> Bin grid.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Number concentration (#(ln(r))d(ln(r))).
    real(kind=dp), intent(out) :: num_conc(bin_grid_size(bin_grid))

    integer :: k
    real(kind=dp) :: mean_vol, num_conc_vol

    mean_vol = aero_data_rad2vol(aero_data, radius_at_mean_vol)
    do k = 1,bin_grid_size(bin_grid)
       num_conc_vol = total_num_conc / mean_vol &
            * exp(-(aero_data_rad2vol(aero_data, bin_grid%centers(k)) &
            / mean_vol))
       call vol_to_lnr(bin_grid%centers(k), num_conc_vol, num_conc(k))
    end do

  end subroutine num_conc_exp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Exponential distribution in volume.
  subroutine vol_conc_exp(total_num_conc, radius_at_mean_vol, &
       bin_grid, aero_data, vol_conc)

    !> Total number concentration of the mode (m^{-3}).
    real(kind=dp), intent(in) :: total_num_conc
    !> Radius at mean volume (m).
    real(kind=dp), intent(in) :: radius_at_mean_vol
    !> Bin grid.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Volume concentration (V(ln(r))d(ln(r))).
    real(kind=dp), intent(out) :: vol_conc(bin_grid_size(bin_grid))

    real(kind=dp) :: num_conc(bin_grid_size(bin_grid))

    call num_conc_exp(total_num_conc, radius_at_mean_vol, &
         bin_grid, aero_data, num_conc)
    vol_conc = num_conc * aero_data_rad2vol(aero_data, bin_grid%centers)

  end subroutine vol_conc_exp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Mono-disperse distribution.
  !> Normalized so that sum(num_conc(k) * log_width) = 1
  subroutine num_conc_mono(total_num_conc, radius, bin_grid, num_conc)

    !> Total number concentration of the mode (m^{-3}).
    real(kind=dp), intent(in) :: total_num_conc
    !> Radius of each particle (m^3).
    real(kind=dp), intent(in) :: radius
    !> Bin grid.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Number concentration (#(ln(r))d(ln(r))).
    real(kind=dp), intent(out) :: num_conc(bin_grid_size(bin_grid))

    integer :: k

    num_conc = 0d0
    k = bin_grid_find(bin_grid, radius)
    if ((k < 1) .or. (k > bin_grid_size(bin_grid))) then
       call warn_msg(825666877, "monodisperse radius outside of bin_grid")
    else
       num_conc(k) = total_num_conc / bin_grid%widths(k)
    end if

  end subroutine num_conc_mono

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Mono-disperse distribution in volume.
  subroutine vol_conc_mono(total_num_conc, radius, &
       bin_grid, aero_data, vol_conc)

    !> Total number concentration of the mode (m^{-3}).
    real(kind=dp), intent(in) :: total_num_conc
    !> Radius of each particle (m^3).
    real(kind=dp), intent(in) :: radius
    !> Bin grid.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Volume concentration (V(ln(r))d(ln(r))).
    real(kind=dp), intent(out) :: vol_conc(bin_grid_size(bin_grid))

    integer :: k

    vol_conc = 0d0
    k = bin_grid_find(bin_grid, radius)
    if ((k < 1) .or. (k > bin_grid_size(bin_grid))) then
       call warn_msg(420930707, "monodisperse radius outside of bin_grid")
    else
       vol_conc(k) = total_num_conc / bin_grid%widths(k) &
            * aero_data_rad2vol(aero_data, radius)
    end if

  end subroutine vol_conc_mono

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Sampled distribution, not normalized.
  subroutine num_conc_sampled(sample_radius, sample_num_conc, bin_grid, &
       num_conc)

    !> Sampled radius bin edges (m).
    real(kind=dp), intent(in) :: sample_radius(:)
    !> Sampled number concentrations (m^{-3}).
    real(kind=dp), intent(in) :: sample_num_conc(:)
    !> Bin grid.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Number concentration (#(ln(r))d(ln(r))).
    real(kind=dp), intent(out) :: num_conc(bin_grid_size(bin_grid))

    integer :: i_sample, n_sample, i_lower, i_upper, i_bin
    real(kind=dp) :: r_lower, r_upper
    real(kind=dp) :: r_bin_lower, r_bin_upper, r1, r2, ratio

    n_sample = size(sample_num_conc)
    call assert(188766208, size(sample_radius) == n_sample + 1)
    call assert(295384037, n_sample >= 1)

    num_conc = 0d0
    do i_sample = 1,n_sample
       r_lower = sample_radius(i_sample)
       r_upper = sample_radius(i_sample + 1)
       i_lower = bin_grid_find(bin_grid, r_lower)
       i_upper = bin_grid_find(bin_grid, r_upper)
       if (i_upper < 1) cycle
       if (i_lower > bin_grid_size(bin_grid)) cycle
       i_lower = max(1, i_lower)
       i_upper = min(bin_grid_size(bin_grid), i_upper)
       do i_bin = i_lower,i_upper
          r_bin_lower = bin_grid%edges(i_bin)
          r_bin_upper = bin_grid%edges(i_bin + 1)
          r1 = max(r_lower, r_bin_lower)
          r2 = min(r_upper, r_bin_upper)
          ratio = (log(r2) - log(r1)) / (log(r_upper) - log(r_lower))
          num_conc(i_bin) = num_conc(i_bin) + ratio &
               * sample_num_conc(i_sample) / bin_grid%widths(i_bin)
       end do
    end do

  end subroutine num_conc_sampled

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Sampled distribution in volume.
  subroutine vol_conc_sampled(sample_radius, sample_num_conc, &
       bin_grid, aero_data, vol_conc)

    !> Sampled radius bin edges (m).
    real(kind=dp), intent(in) :: sample_radius(:)
    !> Sampled number concentrations (m^{-3}).
    real(kind=dp), intent(in) :: sample_num_conc(:)
    !> Bin grid.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Volume concentration (V(ln(r))d(ln(r))).
    real(kind=dp), intent(out) :: vol_conc(bin_grid_size(bin_grid))

    real(kind=dp) :: num_conc(bin_grid_size(bin_grid))

    call num_conc_sampled(sample_radius, sample_num_conc, bin_grid, num_conc)
    vol_conc = num_conc * aero_data_rad2vol(aero_data, bin_grid%centers)

  end subroutine vol_conc_sampled

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Return the binned number concentration for an aero_mode.
  subroutine aero_mode_num_conc(aero_mode, bin_grid, aero_data, &
       num_conc)

    !> Aero mode for which to compute number concentration.
    type(aero_mode_t), intent(in) :: aero_mode
    !> Bin grid.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Number concentration (#(ln(r))d(ln(r))).
    real(kind=dp), intent(out) :: num_conc(bin_grid_size(bin_grid))

    if (aero_mode%type == AERO_MODE_TYPE_LOG_NORMAL) then
       call num_conc_log_normal(aero_mode%num_conc, aero_mode%char_radius, &
            aero_mode%log10_std_dev_radius, bin_grid, num_conc)
    elseif (aero_mode%type == AERO_MODE_TYPE_EXP) then
       call num_conc_exp(aero_mode%num_conc, &
            aero_mode%char_radius, bin_grid, aero_data, num_conc)
    elseif (aero_mode%type == AERO_MODE_TYPE_MONO) then
       call num_conc_mono(aero_mode%num_conc, aero_mode%char_radius, &
            bin_grid, num_conc)
    elseif (aero_mode%type == AERO_MODE_TYPE_SAMPLED) then
       call num_conc_sampled(aero_mode%sample_radius, &
            aero_mode%sample_num_conc, bin_grid, num_conc)
    else
       call die_msg(223903246, "unknown aero_mode type: " &
            // trim(integer_to_string(aero_mode%type)))
    end if

  end subroutine aero_mode_num_conc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Return the binned per-species volume concentration for an
  !> aero_mode.
  subroutine aero_mode_vol_conc(aero_mode, bin_grid, aero_data, &
       vol_conc)

    !> Aero mode for which to compute volume concentration.
    type(aero_mode_t), intent(in) :: aero_mode
    !> Bin grid.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Volume concentration (V(ln(r))d(ln(r))).
    real(kind=dp), intent(out) :: vol_conc(bin_grid_size(bin_grid), &
         aero_data_n_spec(aero_data))

    integer :: i_spec
    real(kind=dp) :: vol_conc_total(bin_grid_size(bin_grid))

    if (aero_mode%type == AERO_MODE_TYPE_LOG_NORMAL) then
       call vol_conc_log_normal(aero_mode%num_conc, &
            aero_mode%char_radius, aero_mode%log10_std_dev_radius, &
            bin_grid, aero_data, vol_conc_total)
    elseif (aero_mode%type == AERO_MODE_TYPE_EXP) then
       call vol_conc_exp(aero_mode%num_conc, &
            aero_mode%char_radius, bin_grid, aero_data, vol_conc_total)
    elseif (aero_mode%type == AERO_MODE_TYPE_MONO) then
       call vol_conc_mono(aero_mode%num_conc, &
            aero_mode%char_radius, bin_grid, aero_data, vol_conc_total)
    elseif (aero_mode%type == AERO_MODE_TYPE_SAMPLED) then
       call vol_conc_sampled(aero_mode%sample_radius, &
            aero_mode%sample_num_conc, bin_grid, aero_data, vol_conc_total)
    else
       call die_msg(314169653, "Unknown aero_mode type: " &
            // trim(integer_to_string(aero_mode%type)))
    end if
    call assert_msg(756593082, sum(aero_mode%vol_frac_std) == 0d0, &
         "cannot convert species fractions with non-zero standard deviation " &
         // "to binned distributions")
    do i_spec = 1,aero_data_n_spec(aero_data)
       vol_conc(:,i_spec) = vol_conc_total * aero_mode%vol_frac(i_spec)
    end do

  end subroutine aero_mode_vol_conc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Compute weighted sampled number concentrations.
  subroutine aero_mode_weighted_sampled_num_conc(aero_mode, aero_weight, &
       weighted_num_conc)

    !> Aerosol mode.
    type(aero_mode_t), intent(in) :: aero_mode
    !> Aerosol weight.
    type(aero_weight_t), intent(in) :: aero_weight
    !> Weighted number concentration
    real(kind=dp), intent(out) :: weighted_num_conc(:)

    integer :: i_sample
    real(kind=dp) :: x0, x1

    call assert(256667423, aero_mode%type == AERO_MODE_TYPE_SAMPLED)
    call assert(878731017, &
         size(weighted_num_conc) == size(aero_mode%sample_num_conc))

    if (aero_weight%type == AERO_WEIGHT_TYPE_NONE) then
       weighted_num_conc = aero_mode%sample_num_conc
    elseif ((aero_weight%type == AERO_WEIGHT_TYPE_POWER) &
         .or. (aero_weight%type == AERO_WEIGHT_TYPE_MFA)) then
       do i_sample = 1,size(aero_mode%sample_num_conc)
          x0 = log(aero_mode%sample_radius(i_sample))
          x1 = log(aero_mode%sample_radius(i_sample + 1))
          weighted_num_conc(i_sample) = aero_mode%sample_num_conc(i_sample) &
               / aero_weight%exponent * (exp(- aero_weight%exponent * x0) &
               - exp(- aero_weight%exponent * x1)) / (x1 - x0)
       end do
    else
       call die_msg(576124393, "unknown aero_weight type: " &
            // trim(integer_to_string(aero_weight%type)))
    end if

  end subroutine aero_mode_weighted_sampled_num_conc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Return the total number of computational particles for an \c aero_mode.
  real(kind=dp) function aero_mode_number(aero_mode, aero_weight)

    !> Aero_mode to sample radius from.
    type(aero_mode_t), intent(in) :: aero_mode
    !> Aero weight.
    type(aero_weight_t), intent(in) :: aero_weight

    real(kind=dp) :: x_mean_prime
    real(kind=dp), allocatable :: weighted_num_conc(:)

    aero_mode_number = 0d0
    if (aero_mode%type == AERO_MODE_TYPE_LOG_NORMAL) then
       if (aero_weight%type == AERO_WEIGHT_TYPE_NONE) then
          aero_mode_number = aero_mode%num_conc / aero_weight%magnitude
       elseif ((aero_weight%type == AERO_WEIGHT_TYPE_POWER) &
            .or. (aero_weight%type == AERO_WEIGHT_TYPE_MFA)) then
          x_mean_prime = log10(aero_mode%char_radius) &
               - aero_weight%exponent * aero_mode%log10_std_dev_radius**2 &
               * log(10d0)
          aero_mode_number = aero_mode%num_conc / aero_weight%magnitude &
               * exp((x_mean_prime**2 - log10(aero_mode%char_radius)**2) &
               / (2d0 * aero_mode%log10_std_dev_radius**2))
       else
          call die_msg(466668240, "unknown aero_weight type: " &
               // trim(integer_to_string(aero_weight%type)))
       end if
    elseif (aero_mode%type == AERO_MODE_TYPE_EXP) then
       if (aero_weight%type == AERO_WEIGHT_TYPE_NONE) then
          aero_mode_number = aero_mode%num_conc / aero_weight%magnitude
       else
          call die_msg(822252601, &
               "cannot use exponential modes with weighting")
       end if
    elseif (aero_mode%type == AERO_MODE_TYPE_MONO) then
       aero_mode_number = aero_mode%num_conc &
            / aero_weight_num_conc_at_radius(aero_weight, &
            aero_mode%char_radius)
    elseif (aero_mode%type == AERO_MODE_TYPE_SAMPLED) then
       allocate(weighted_num_conc(size(aero_mode%sample_num_conc)))
       call aero_mode_weighted_sampled_num_conc(aero_mode, aero_weight, &
            weighted_num_conc)
       aero_mode_number = sum(weighted_num_conc) / aero_weight%magnitude
       deallocate(weighted_num_conc)
    else
       call die_msg(901140225, "unknown aero_mode type: " &
            // trim(integer_to_string(aero_mode%type)))
    end if

  end function aero_mode_number

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Return a radius randomly sampled from the mode distribution.
  subroutine aero_mode_sample_radius(aero_mode, aero_data, aero_weight, &
       radius)

    !> Aero_mode to sample radius from.
    type(aero_mode_t), intent(in) :: aero_mode
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Aero weight.
    type(aero_weight_t), intent(in) :: aero_weight
    !> Sampled radius (m).
    real(kind=dp), intent(out) :: radius

    real(kind=dp) :: x_mean_prime, x0, x1, x, r, inv_nc0, inv_nc1, inv_nc
    integer :: i_sample
    real(kind=dp), allocatable :: weighted_num_conc(:)

    if (aero_mode%type == AERO_MODE_TYPE_LOG_NORMAL) then
       if (aero_weight%type == AERO_WEIGHT_TYPE_NONE) then
          x_mean_prime = log10(aero_mode%char_radius)
       elseif ((aero_weight%type == AERO_WEIGHT_TYPE_POWER) &
            .or. (aero_weight%type == AERO_WEIGHT_TYPE_MFA)) then
          x_mean_prime = log10(aero_mode%char_radius) &
               - aero_weight%exponent * aero_mode%log10_std_dev_radius**2 &
               * log(10d0)
       else
          call die_msg(517376844, "unknown aero_weight type: " &
               // trim(integer_to_string(aero_weight%type)))
       end if
       radius = 10d0**rand_normal(x_mean_prime, &
            aero_mode%log10_std_dev_radius)
    elseif (aero_mode%type == AERO_MODE_TYPE_SAMPLED) then
       allocate(weighted_num_conc(size(aero_mode%sample_num_conc)))
       call aero_mode_weighted_sampled_num_conc(aero_mode, aero_weight, &
            weighted_num_conc)
       i_sample = sample_cts_pdf(weighted_num_conc)
       deallocate(weighted_num_conc)
       if ((aero_weight%type == AERO_WEIGHT_TYPE_NONE) &
            .or. (((aero_weight%type == AERO_WEIGHT_TYPE_POWER) &
            .or. (aero_weight%type == AERO_WEIGHT_TYPE_MFA)) &
            .and. (aero_weight%exponent == 0d0))) then
          x0 = log(aero_mode%sample_radius(i_sample))
          x1 = log(aero_mode%sample_radius(i_sample + 1))
          r = pmc_random()
          x = (1d0 - r) * x0 + r * x1
          radius = exp(x)
       elseif ((aero_weight%type == AERO_WEIGHT_TYPE_POWER) &
            .or. (aero_weight%type == AERO_WEIGHT_TYPE_MFA)) then
          inv_nc0 = 1d0 / aero_weight_num_conc_at_radius(aero_weight, &
               aero_mode%sample_radius(i_sample))
          inv_nc1 = 1d0 / aero_weight_num_conc_at_radius(aero_weight, &
               aero_mode%sample_radius(i_sample + 1))
          r = pmc_random()
          inv_nc = (1d0 - r) * inv_nc0 + r * inv_nc1
          radius = aero_weight_radius_at_num_conc(aero_weight, 1d0 / inv_nc)
       else
          call die_msg(769131141, "unknown aero_weight type: " &
               // trim(integer_to_string(aero_weight%type)))
       end if
    elseif (aero_mode%type == AERO_MODE_TYPE_EXP) then
       if (aero_weight%type == AERO_WEIGHT_TYPE_NONE) then
          radius = aero_data_vol2rad(aero_data, -aero_data_rad2vol(aero_data, &
               aero_mode%char_radius) * log(pmc_random()))
       elseif ((aero_weight%type == AERO_WEIGHT_TYPE_POWER) &
            .or. (aero_weight%type == AERO_WEIGHT_TYPE_MFA)) then
          call die_msg(678481276, &
               "cannot use exponential modes with weighting")
       else
          call die_msg(301787712, "unknown aero_weight type: " &
               // trim(integer_to_string(aero_weight%type)))
       end if
    elseif (aero_mode%type == AERO_MODE_TYPE_MONO) then
       radius = aero_mode%char_radius
    else
       call die_msg(749122931, "Unknown aero_mode type: " &
            // trim(integer_to_string(aero_mode%type)))
    end if

  end subroutine aero_mode_sample_radius

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Return an array of volumes randomly sampled from the volume fractions.
  subroutine aero_mode_sample_vols(aero_mode, total_vol, vols)

    !> Aero_mode to sample from.
    type(aero_mode_t), intent(in) :: aero_mode
    !> Total volume (m^3).
    real(kind=dp), intent(in) :: total_vol
    !> Sampled volumes (m^3).
    real(kind=dp), intent(out) :: vols(size(aero_mode%vol_frac))

    integer, parameter :: AERO_MODE_MAX_SAMPLE_LOOPS = 1000000

    integer :: i_sample
    real(kind=dp) :: offset

    ! sampling of volume fractions, normalized to sum to 1 by
    ! projecting out the mean direction, and accepted if all
    ! non-negative, otherwise rejected and repeated (accept-reject)
    do i_sample = 1,AERO_MODE_MAX_SAMPLE_LOOPS
       call rand_normal_array_1d(aero_mode%vol_frac, aero_mode%vol_frac_std, &
            vols)
       ! add the correct amount of (offset * vol_frac) to vols to ensure
       ! that sum(vols) is 1
       offset = 1d0 - sum(vols)
       vols = vols + offset * aero_mode%vol_frac
       if (minval(vols) >= 0d0) exit
    end do
    if (i_sample == AERO_MODE_MAX_SAMPLE_LOOPS) then
       call die_msg(549015143, "Unable to sample non-negative volumes for " &
            // "mode: " // trim(aero_mode%name))
    end if
    vols = vols / sum(vols) * total_vol

  end subroutine aero_mode_sample_vols

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read volume fractions from a data file.
  subroutine spec_file_read_vol_frac(file, aero_data, vol_frac, vol_frac_std)

    !> Spec file to read mass fractions from.
    type(spec_file_t), intent(inout) :: file
    !> Aero_data data.
    type(aero_data_t), intent(in) :: aero_data
    !> Aerosol species volume fractions.
    real(kind=dp), allocatable, intent(inout) :: vol_frac(:)
    !> Aerosol species volume fraction standard deviations.
    real(kind=dp), allocatable, intent(inout) :: vol_frac_std(:)

    integer :: n_species, species, i
    character(len=SPEC_LINE_MAX_VAR_LEN), allocatable :: species_name(:)
    real(kind=dp), allocatable :: species_data(:,:)
    real(kind=dp) :: tot_vol_frac

    !> \page input_format_mass_frac Input File Format: Aerosol Mass Fractions
    !!
    !! An aerosol mass fractions file must consist of one line per
    !! aerosol species, with each line having the species name
    !! followed by the species mass fraction in each aerosol
    !! particle. The valid species names are those specfied by the
    !! \ref input_format_aero_data file, but not all species have to
    !! be listed. Any missing species will have proportions of
    !! zero. If the proportions do not sum to 1 then they will be
    !! normalized before use. For example, a mass fractions file file
    !! could contain:
    !! <pre>
    !! # species   proportion
    !! OC          0.3
    !! BC          0.7
    !! </pre>
    !! indicating that the particles are 30% organic carbon and 70%
    !! black carbon.
    !!
    !! Optionally, the standard deviation can also be provided for
    !! each species as a second number on each line. For example,
    !! <pre>
    !! # species   proportion std_dev
    !! OC          0.3        0.1
    !! BC          0.7        0.2
    !! </pre>
    !! indicates that the particles are on average 30% OC and 70% BC,
    !! but may vary to have particles with 20% OC and 80% BC, or 40%
    !! OC and 60% BC, for example. The standard deviations will be
    !! normalized by the sum of the proportions.
    !!
    !! Either all species in a given file must have standard
    !! deviations or none of them can.
    !!
    !! See also:
    !!   - \ref spec_file_format --- the input file text format
    !!   - \ref input_format_aero_dist --- the format for a complete
    !!     aerosol distribution with several modes
    !!   - \ref input_format_aero_mode --- the format for each mode
    !!     of an aerosol distribution

    ! read the aerosol data from the specified file
    call spec_file_read_real_named_array(file, 0, species_name, &
         species_data)

    ! check the data size
    n_species = size(species_data, 1)
    if (n_species < 1) then
       call die_msg(628123166, 'file ' // trim(file%name) &
            // ' must contain at least one line of data')
    end if
    if ((size(species_data, 2) /= 1) .and. (size(species_data, 2) /= 2)) then
       call die_msg(427666881, 'each line in file ' // trim(file%name) &
            // ' must contain exactly one or two data values')
    end if

    ! copy over the data
    if (allocated(vol_frac)) deallocate(vol_frac)
    if (allocated(vol_frac_std)) deallocate(vol_frac_std)
    allocate(vol_frac(aero_data_n_spec(aero_data)))
    allocate(vol_frac_std(aero_data_n_spec(aero_data)))
    vol_frac = 0d0
    vol_frac_std = 0d0
    do i = 1,n_species
       species = aero_data_spec_by_name(aero_data, species_name(i))
       if (species == 0) then
          call die_msg(775942501, 'unknown species ' // trim(species_name(i)) &
               // ' in file ' // trim(file%name))
       end if
       vol_frac(species) = species_data(i, 1)
       if (size(species_data, 2) == 2) then
          vol_frac_std(species) = species_data(i, 2)
       end if
    end do

    ! convert mass fractions to volume fractions
    vol_frac = vol_frac / aero_data%density
    vol_frac_std = vol_frac_std / aero_data%density

    ! normalize
    tot_vol_frac = sum(vol_frac)
    if ((minval(vol_frac) < 0d0) .or. (tot_vol_frac <= 0d0)) then
       call die_msg(356648030, 'fractions in ' // trim(file%name) &
            // ' are not positive')
    end if
    if (minval(vol_frac_std) < 0d0) then
       call die_msg(676576501, 'standard deviations in ' // trim(file%name) &
            // ' are not positive')
    end if
    vol_frac = vol_frac / tot_vol_frac
    vol_frac_std = vol_frac_std / tot_vol_frac

  end subroutine spec_file_read_vol_frac

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read a size distribution from a data file.
  subroutine spec_file_read_size_dist(file, sample_radius, sample_num_conc)

    !> Spec file to read size distribution from.
    type(spec_file_t), intent(inout) :: file
    !> Sample radius values (m).
    real(kind=dp), allocatable, intent(inout) :: sample_radius(:)
    !> Sample number concentrations (m^{-3}).
    real(kind=dp), allocatable, intent(inout) :: sample_num_conc(:)

    character(len=SPEC_LINE_MAX_VAR_LEN), allocatable :: names(:)
    real(kind=dp), allocatable :: data(:,:)
    integer :: n_sample, i_sample

    !> \page input_format_size_dist Input File Format: Size Distribution
    !!
    !! A size distribution file must consist of two lines:
    !! - the first line must begin with \c diam and be followed by
    !!   \f$N + 1\f$ space-separated real scalars, giving the diameters
    !!   \f$D_1,\ldots,D_{N+1}\f$ of bin edges (m) --- these must be
    !!   in increasing order, so \f$D_i < D_{i+1}\f$
    !! - the second line must begin with \c num_conc and be followed
    !!   by \f$N\f$ space-separated real scalars, giving the number
    !!   concenrations \f$C_1,\ldots,C_N\f$ in each bin (#/m^3) ---
    !!   \f$C_i\f$ is the total number concentrations of particles
    !!   with diameters in \f$[D_i, D_{i+1}]\f$
    !!
    !! The resulting size distribution is taken to be piecewise
    !! constant in log-diameter coordinates.
    !!
    !! Example: a size distribution could be:
    !! <pre>
    !! diam 1e-7 1e-6 1e-5  # bin edge diameters (m)
    !! num_conc 1e9 1e8     # bin number concentrations (m^{-3})
    !! </pre>
    !! This distribution has 1e9 particles per cubic meter with
    !! diameters between 0.1 micron and 1 micron, and 1e8 particles
    !! per cubic meter with diameters between 1 micron and 10 micron.
    !!
    !! See also:
    !!   - \ref spec_file_format --- the input file text format
    !!   - \ref input_format_aero_dist --- the format for a complete
    !!     aerosol distribution with several modes
    !!   - \ref input_format_aero_mode --- the format for each mode
    !!     of an aerosol distribution

    ! read the data from the file
    call spec_file_read_real_named_array(file, 1, names, data)
    call spec_file_assert_msg(311818741, file, size(names) == 1, &
         'must contain a line starting with "diam"')
    call spec_file_check_name(file, 'diam', names(1))
    n_sample = size(data,2) - 1
    call spec_file_assert_msg(669011124, file, n_sample >= 1, &
         'must have at least two diam values')

    if (allocated(sample_radius)) deallocate(sample_radius)
    allocate(sample_radius(n_sample + 1))
    sample_radius = diam2rad(data(1,:))
    do i_sample = 1,n_sample
       call spec_file_assert_msg(528089871, file, &
            sample_radius(i_sample) < sample_radius(i_sample + 1), &
            'diam values must be strictly increasing')
    end do

    call spec_file_read_real_named_array(file, 1, names, data)
    call spec_file_assert_msg(801676496, file, size(names) == 1, &
         'must contain a line starting with "num_conc"')
    call spec_file_check_name(file, 'num_conc', names(1))

    call spec_file_assert_msg(721029144, file, size(data, 2) == n_sample, &
         'must have one fewer num_conc than diam values')

    if (allocated(sample_num_conc)) deallocate(sample_num_conc)
    allocate(sample_num_conc(n_sample))
    sample_num_conc = data(1,:)
    do i_sample = 1,n_sample
       call spec_file_assert_msg(356490397, file, &
            sample_num_conc(i_sample) >= 0d0, &
            'num_conc values must be non-negative')
    end do

  end subroutine spec_file_read_size_dist

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read one mode of an aerosol distribution (number concentration,
  !> volume fractions, and mode shape).
  subroutine spec_file_read_aero_mode(file, aero_data, aero_mode, eof)

    !> Spec file.
    type(spec_file_t), intent(inout) :: file
    !> Aero_data data.
    type(aero_data_t), intent(inout) :: aero_data
    !> Aerosol mode.
    type(aero_mode_t), intent(inout) :: aero_mode
    !> If eof instead of reading data.
    logical :: eof

    character(len=SPEC_LINE_MAX_VAR_LEN) :: tmp_str, mode_type, diam_type_str
    character(len=SPEC_LINE_MAX_VAR_LEN) :: mass_frac_filename
    character(len=SPEC_LINE_MAX_VAR_LEN) :: size_dist_filename
    type(spec_line_t) :: line
    type(spec_file_t) :: mass_frac_file, size_dist_file
    real(kind=dp) :: diam, temp, pressure
    integer :: diam_type, i_radius
    real(kind=dp) :: log10_r0_mob, log10_r1_mob, log10_r2_mob, &
         r0_mob, r2_mob, r0_geom, r2_geom, log10_r0_geom, log10_r2_geom

    ! note that doxygen's automatic list creation breaks on the list
    ! below for some reason, so we just use html format

    !> \page input_format_aero_mode Input File Format: Aerosol Distribution Mode
    !!
    !! An aerosol distribution mode has the parameters:
    !! <ul>
    !! <li> \b mode_name (string): the name of the mode (for
    !!      informational purposes only)
    !! <li> \b mass_frac (string): name of file from which to read the
    !!      species mass fractions --- the file format should
    !!      be \subpage input_format_mass_frac
    !! <li> \b diam_type (string): the type of diameter for the mode
    !!      --- must be one of: \c geometric for geometric diameter;
    !!      or \c mobility for mobility equivalent diameter
    !! <li> if \c diam_type is \c mobility then the following
    !!      parameters are:
    !!      <ul>
    !!      <li> \b temp (real, unit K): the temperate at which the
    !!           mobility diameters were measured
    !!      <li> \b pressure (real, unit Pa): the pressure at which the
    !!           mobility diameters were measured
    !!      </ul>
    !! <li> \b mode_type (string): the functional form of the mode ---
    !!      must be one of: \c log_normal for a log-normal
    !!      distribution; \c exp for an exponential distribution; \c
    !!      mono for a mono-disperse distribution; or \c sampled for a
    !!      sampled distribution
    !! <li> if \c mode_type is \c log_normal then the mode distribution
    !!      shape is
    !!      \f[ n(\log D) {\rm d}\log D
    !!      = \frac{N_{\rm total}}{\sqrt{2\pi} \log \sigma_{\rm g}}
    !!      \exp\left(\frac{(\log D - \log D_{\rm gn})^2}{2 \log ^2
    !!      \sigma_{\rm g}}\right)
    !!      {\rm d}\log D \f]
    !!      and the following parameters are:
    !!      <ul>
    !!      <li> \b num_conc (real, unit 1/m^3): the total number
    !!           concentration \f$N_{\rm total}\f$ of the mode
    !!      <li> \b geom_mean_diam (real, unit m): the geometric mean
    !!           diameter \f$D_{\rm gn}\f$
    !!      <li> \b log10_geom_std_dev (real, dimensionless):
    !!           \f$\log_{10}\f$ of the geometric standard deviation
    !!           \f$\sigma_{\rm g}\f$ of the diameter
    !!      </ul>
    !! <li> if \c mode_type is \c exp then the mode distribution shape is
    !!      \f[ n(v) {\rm d}v = \frac{N_{\rm total}}{v_{\rm \mu}}
    !!      \exp\left(- \frac{v}{v_{\rm \mu}}\right)
    !!      {\rm d}v \f]
    !!      and the following parameters are:
    !!      <ul>
    !!      <li> \b num_conc (real, unit 1/m^3): the total number
    !!           concentration \f$N_{\rm total}\f$ of the mode
    !!      <li> \b diam_at_mean_vol (real, unit m): the diameter
    !!           \f$D_{\rm \mu}\f$ such that \f$v_{\rm \mu}
    !!           = \frac{\pi}{6} D^3_{\rm \mu}\f$
    !!      </ul>
    !! <li> if \c mode_type is \c mono then the mode distribution shape
    !!      is a delta distribution at diameter \f$D_0\f$ and the
    !!      following parameters are:
    !!      <ul>
    !!      <li> \b num_conc (real, unit 1/m^3): the total number
    !!           concentration \f$N_{\rm total}\f$ of the mode
    !!      <li> \b radius (real, unit m): the radius \f$R_0\f$ of the
    !!           particles, so that \f$D_0 = 2 R_0\f$
    !!      </ul>
    !! <li> if \c mode_type is \c sampled then the mode distribution
    !!      shape is piecewise constant (in log-diameter coordinates)
    !!      and the following parameters are:
    !!      <ul>
    !!      <li> \b size_dist (string): name of file from which to
    !!           read the size distribution --- the file format should
    !!           be \subpage input_format_size_dist
    !!      </ul>
    !! </ul>
    !!
    !! Example:
    !! <pre>
    !! mode_name diesel          # mode name (descriptive only)
    !! mass_frac comp_diesel.dat # mass fractions in each aerosol particle
    !! mode_type log_normal      # type of distribution
    !! num_conc 1.6e8            # particle number density (#/m^3)
    !! geom_mean_diam 2.5e-8     # geometric mean diameter (m)
    !! log10_geom_std_dev 0.24   # log_10 of geometric standard deviation
    !! </pre>
    !!
    !! See also:
    !!   - \ref spec_file_format --- the input file text format
    !!   - \ref input_format_aero_dist --- the format for a complete
    !!     aerosol distribution with several modes
    !!   - \ref input_format_mass_frac --- the format for the mass
    !!     fractions file

    call spec_file_read_line(file, line, eof)
    if (.not. eof) then
       call spec_file_check_line_name(file, line, "mode_name")
       call spec_file_check_line_length(file, line, 1)
       tmp_str = line%data(1) ! hack to avoid gfortran warning
       aero_mode%name = tmp_str(1:AERO_MODE_NAME_LEN)
       aero_mode%source = aero_data_source_by_name(aero_data, aero_mode%name)

       call spec_file_read_string(file, 'mass_frac', mass_frac_filename)
       call spec_file_open(mass_frac_filename, mass_frac_file)
       call spec_file_read_vol_frac(mass_frac_file, aero_data, &
            aero_mode%vol_frac, aero_mode%vol_frac_std)
       call spec_file_close(mass_frac_file)

       call spec_file_read_string(file, 'diam_type', diam_type_str)
       if (trim(diam_type_str) == 'geometric') then
          diam_type = AERO_MODE_DIAM_TYPE_GEOMETRIC
       elseif (trim(diam_type_str) == 'mobility') then
          diam_type = AERO_MODE_DIAM_TYPE_MOBILITY
          call spec_file_read_real(file, 'temp', temp)
          call spec_file_read_real(file, 'pressure', pressure)
       else
          call spec_file_die_msg(804343794, file, &
               "Unknown diam_type: " // trim(diam_type_str))
       end if

       call spec_file_read_string(file, 'mode_type', mode_type)
       aero_mode%sample_radius = [ real(kind=dp) :: ]
       aero_mode%sample_num_conc = [ real(kind=dp) :: ]
       if (trim(mode_type) == 'log_normal') then
          aero_mode%type = AERO_MODE_TYPE_LOG_NORMAL
          call spec_file_read_real(file, 'num_conc', aero_mode%num_conc)
          call spec_file_read_real(file, 'geom_mean_diam', diam)
          call spec_file_read_real(file, 'log10_geom_std_dev', &
               aero_mode%log10_std_dev_radius)
          if (diam_type == AERO_MODE_DIAM_TYPE_GEOMETRIC) then
             aero_mode%char_radius = diam2rad(diam)
          elseif (diam_type == AERO_MODE_DIAM_TYPE_MOBILITY) then
             aero_mode%char_radius &
                  = aero_data_mobility_rad_to_geometric_rad(aero_data, &
                  diam2rad(diam), temp, pressure)

             ! Convert log10_std_dev_radius from mobility to geometric radius.
             ! We do this by finding points +/- one std dev in mobility
             ! radius, converting them to geometric radius, and then using
             ! the distance between them as a measure of 2 * std_dev in
             ! geometric radius.
             log10_r1_mob = log10(diam2rad(diam))
             log10_r0_mob = log10_r1_mob - aero_mode%log10_std_dev_radius
             log10_r2_mob = log10_r1_mob + aero_mode%log10_std_dev_radius
             r0_mob = 10**log10_r0_mob
             r2_mob = 10**log10_r2_mob
             r0_geom = aero_data_mobility_rad_to_geometric_rad(aero_data, &
                  r0_mob, temp, pressure)
             r2_geom = aero_data_mobility_rad_to_geometric_rad(aero_data, &
                  r2_mob, temp, pressure)
             log10_r0_geom = log10(r0_geom)
             log10_r2_geom = log10(r2_geom)
             aero_mode%log10_std_dev_radius &
                  = (log10_r2_geom - log10_r0_geom) / 2d0
          else
             call die_msg(532966100, "Diameter type not handled: " &
                  // integer_to_string(diam_type))
          end if
       elseif (trim(mode_type) == 'exp') then
          aero_mode%type = AERO_MODE_TYPE_EXP
          call spec_file_read_real(file, 'num_conc', aero_mode%num_conc)
          call spec_file_read_real(file, 'diam_at_mean_vol', diam)
          if (diam_type == AERO_MODE_DIAM_TYPE_GEOMETRIC) then
             aero_mode%char_radius = diam2rad(diam)
          elseif (diam_type == AERO_MODE_DIAM_TYPE_MOBILITY) then
             aero_mode%char_radius &
                  = aero_data_mobility_rad_to_geometric_rad(aero_data, &
                  diam2rad(diam), temp, pressure)
          else
             call die_msg(585104460, "Diameter type not handled: " &
                  // integer_to_string(diam_type))
          end if
       elseif (trim(mode_type) == 'mono') then
          aero_mode%type = AERO_MODE_TYPE_MONO
          call spec_file_read_real(file, 'num_conc', aero_mode%num_conc)
          call spec_file_read_real(file, 'diam', diam)
          if (diam_type == AERO_MODE_DIAM_TYPE_GEOMETRIC) then
             aero_mode%char_radius = diam2rad(diam)
          elseif (diam_type == AERO_MODE_DIAM_TYPE_MOBILITY) then
             aero_mode%char_radius &
                  = aero_data_mobility_rad_to_geometric_rad(aero_data, &
                  diam2rad(diam), temp, pressure)
          else
             call die_msg(902864269, "Diameter type not handled: " &
                  // integer_to_string(diam_type))
          end if
       elseif (trim(mode_type) == 'sampled') then
          aero_mode%type = AERO_MODE_TYPE_SAMPLED
          call spec_file_read_string(file, 'size_dist', size_dist_filename)
          call spec_file_open(size_dist_filename, size_dist_file)
          call spec_file_read_size_dist(size_dist_file, &
               aero_mode%sample_radius, aero_mode%sample_num_conc)
          call spec_file_close(size_dist_file)
          if (diam_type == AERO_MODE_DIAM_TYPE_GEOMETRIC) then
             ! do nothing
          elseif (diam_type == AERO_MODE_DIAM_TYPE_MOBILITY) then
             do i_radius = 1,size(aero_mode%sample_radius)
                aero_mode%sample_radius(i_radius) &
                     = aero_data_mobility_rad_to_geometric_rad(aero_data, &
                     aero_mode%sample_radius(i_radius), temp, pressure)
             end do
          else
             call die_msg(239088838, "Diameter type not handled: " &
                  // integer_to_string(diam_type))
          end if
       else
          call spec_file_die_msg(729472928, file, &
               "Unknown distribution mode type: " // trim(mode_type))
       end if
    end if

  end subroutine spec_file_read_aero_mode

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Determines the number of bytes required to pack the given value.
  integer function pmc_mpi_pack_size_aero_mode(val)

    !> Value to pack.
    type(aero_mode_t), intent(in) :: val

    pmc_mpi_pack_size_aero_mode = &
         pmc_mpi_pack_size_string(val%name) &
         + pmc_mpi_pack_size_integer(val%type) &
         + pmc_mpi_pack_size_real(val%char_radius) &
         + pmc_mpi_pack_size_real(val%log10_std_dev_radius) &
         + pmc_mpi_pack_size_real_array(val%sample_radius) &
         + pmc_mpi_pack_size_real_array(val%sample_num_conc) &
         + pmc_mpi_pack_size_real(val%num_conc) &
         + pmc_mpi_pack_size_real_array(val%vol_frac) &
         + pmc_mpi_pack_size_real_array(val%vol_frac_std) &
         + pmc_mpi_pack_size_integer(val%source)

  end function pmc_mpi_pack_size_aero_mode

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Packs the given value into the buffer, advancing position.
  subroutine pmc_mpi_pack_aero_mode(buffer, position, val)

    !> Memory buffer.
    character, intent(inout) :: buffer(:)
    !> Current buffer position.
    integer, intent(inout) :: position
    !> Value to pack.
    type(aero_mode_t), intent(in) :: val

#ifdef PMC_USE_MPI
    integer :: prev_position

    prev_position = position
    call pmc_mpi_pack_string(buffer, position, val%name)
    call pmc_mpi_pack_integer(buffer, position, val%type)
    call pmc_mpi_pack_real(buffer, position, val%char_radius)
    call pmc_mpi_pack_real(buffer, position, val%log10_std_dev_radius)
    call pmc_mpi_pack_real_array(buffer, position, val%sample_radius)
    call pmc_mpi_pack_real_array(buffer, position, val%sample_num_conc)
    call pmc_mpi_pack_real(buffer, position, val%num_conc)
    call pmc_mpi_pack_real_array(buffer, position, val%vol_frac)
    call pmc_mpi_pack_real_array(buffer, position, val%vol_frac_std)
    call pmc_mpi_pack_integer(buffer, position, val%source)
    call assert(497092471, &
         position - prev_position <= pmc_mpi_pack_size_aero_mode(val))
#endif

  end subroutine pmc_mpi_pack_aero_mode

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Unpacks the given value from the buffer, advancing position.
  subroutine pmc_mpi_unpack_aero_mode(buffer, position, val)

    !> Memory buffer.
    character, intent(inout) :: buffer(:)
    !> Current buffer position.
    integer, intent(inout) :: position
    !> Value to pack.
    type(aero_mode_t), intent(inout) :: val

#ifdef PMC_USE_MPI
    integer :: prev_position

    prev_position = position
    call pmc_mpi_unpack_string(buffer, position, val%name)
    call pmc_mpi_unpack_integer(buffer, position, val%type)
    call pmc_mpi_unpack_real(buffer, position, val%char_radius)
    call pmc_mpi_unpack_real(buffer, position, val%log10_std_dev_radius)
    call pmc_mpi_unpack_real_array(buffer, position, val%sample_radius)
    call pmc_mpi_unpack_real_array(buffer, position, val%sample_num_conc)
    call pmc_mpi_unpack_real(buffer, position, val%num_conc)
    call pmc_mpi_unpack_real_array(buffer, position, val%vol_frac)
    call pmc_mpi_unpack_real_array(buffer, position, val%vol_frac_std)
    call pmc_mpi_unpack_integer(buffer, position, val%source)
    call assert(874467577, &
         position - prev_position <= pmc_mpi_pack_size_aero_mode(val))
#endif

  end subroutine pmc_mpi_unpack_aero_mode

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_aero_mode
