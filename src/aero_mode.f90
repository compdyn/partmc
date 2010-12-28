! Copyright (C) 2005-2010 Nicole Riemer and Matthew West
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

  !> An aerosol size distribution mode.
  !!
  !! Each mode is assumed to be fully internally mixed so that every
  !! particle has the same composition. The \c num_conc array then
  !! stores the number concentration distribution.
  type aero_mode_t
     !> Mode name, used to track particle sources.
     character(len=AERO_MODE_NAME_LEN) :: name
     !> Mode type (given by module constants).
     integer :: type
     !> Mean radius of mode (m).
     real(kind=dp) :: mean_radius
     !> Log base 10 of geometric standard deviation of radius, if
     !> necessary (m).
     real(kind=dp) :: log10_std_dev_radius
     !> Total number concentration of mode (#/m^3).
     real(kind=dp) :: num_conc
     !> Species fractions by volume [length \c aero_data%%n_spec] (1).
     real(kind=dp), pointer :: vol_frac(:)
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
    else
       aero_mode_type_to_string = "unknown"
    end if

  end function aero_mode_type_to_string

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Allocates an aero_mode.
  subroutine aero_mode_allocate(aero_mode)

    !> Aerosol mode.
    type(aero_mode_t), intent(out) :: aero_mode

    allocate(aero_mode%vol_frac(0))

  end subroutine aero_mode_allocate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Allocates an aero_mode of the given size.
  subroutine aero_mode_allocate_size(aero_mode, n_spec)

    !> Aerosol mode.
    type(aero_mode_t), intent(out) :: aero_mode
    !> Number of species.
    integer, intent(in) :: n_spec

    allocate(aero_mode%vol_frac(n_spec))

  end subroutine aero_mode_allocate_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Free all storage.
  subroutine aero_mode_deallocate(aero_mode)

    !> Aerosol mode.
    type(aero_mode_t), intent(inout) :: aero_mode

    deallocate(aero_mode%vol_frac)

  end subroutine aero_mode_deallocate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Copy an aero_mode.
  subroutine aero_mode_copy(aero_mode_from, aero_mode_to)

    !> Aerosol mode original.
    type(aero_mode_t), intent(in) :: aero_mode_from
    !> Aerosol mode copy.
    type(aero_mode_t), intent(inout) :: aero_mode_to

    call aero_mode_deallocate(aero_mode_to)
    call aero_mode_allocate_size(aero_mode_to, size(aero_mode_from%vol_frac))
    aero_mode_to%name = aero_mode_from%name
    aero_mode_to%type = aero_mode_from%type
    aero_mode_to%mean_radius = aero_mode_from%mean_radius
    aero_mode_to%log10_std_dev_radius = aero_mode_from%log10_std_dev_radius
    aero_mode_to%num_conc = aero_mode_from%num_conc
    aero_mode_to%vol_frac = aero_mode_from%vol_frac

  end subroutine aero_mode_copy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Compute a log-normal distribution, normalized so that
  !> sum(num_conc(k) * dlnr) = 1
  subroutine num_conc_log_normal(mean_radius, log_sigma, bin_grid, num_conc)
    
    !> Geometric mean radius (m).
    real(kind=dp), intent(in) :: mean_radius
    !> log_10(geom. std dev) (1).
    real(kind=dp), intent(in) :: log_sigma
    !> Bin grid.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Number concentration (#(ln(r))d(ln(r))).
    real(kind=dp), intent(out) :: num_conc(bin_grid%n_bin)
    
    integer :: k
    
    do k = 1,bin_grid%n_bin
       num_conc(k) = 1d0 / (sqrt(2d0 * const%pi) * log_sigma) * &
            dexp(-(dlog10(vol2rad(bin_grid%v(k))) &
            - dlog10(mean_radius))**2d0 &
            / (2d0 * log_sigma**2d0)) / dlog(10d0)
    end do
    
    ! The formula above was originally for a distribution in
    ! log_10(r), while we are using log_e(r) for our bin grid. The
    ! division by dlog(10) at the end corrects for this.

    ! Remember that log_e(r) = log_10(r) * log_e(10).
    
  end subroutine num_conc_log_normal
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Compute a log-normal distribution in volume.
  subroutine vol_conc_log_normal(mean_radius, log_sigma, bin_grid, vol_conc)
    
    !> Geometric mean radius (m).
    real(kind=dp), intent(in) :: mean_radius
    !> log_10(geom. std dev) (1).
    real(kind=dp), intent(in) :: log_sigma
    !> Bin grid.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Volume concentration (V(ln(r))d(ln(r))).
    real(kind=dp), intent(out) :: vol_conc(bin_grid%n_bin)
    
    real(kind=dp) :: num_conc(bin_grid%n_bin)

    call num_conc_log_normal(mean_radius, log_sigma, bin_grid, num_conc)
    vol_conc = num_conc * bin_grid%v
    
  end subroutine vol_conc_log_normal
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Exponential distribution in volume
  !> \f[ n(v) = \frac{1}{\rm mean-vol} \exp(- v / {\rm mean-vol}) \f]
  !> Normalized so that sum(num_conc(k) * dlnr) = 1
  subroutine num_conc_exp(mean_radius, bin_grid, num_conc)
    
    !> Mean radius (m).
    real(kind=dp), intent(in) :: mean_radius
    !> Bin grid.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Number concentration (#(ln(r))d(ln(r))).
    real(kind=dp), intent(out) :: num_conc(bin_grid%n_bin)
    
    integer :: k
    real(kind=dp) :: mean_vol, num_conc_vol
    
    mean_vol = rad2vol(mean_radius)
    do k = 1,bin_grid%n_bin
       num_conc_vol = 1d0 / mean_vol * exp(-(bin_grid%v(k) / mean_vol))
       call vol_to_lnr(vol2rad(bin_grid%v(k)), num_conc_vol, num_conc(k))
    end do
    
  end subroutine num_conc_exp
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Exponential distribution in volume.
  subroutine vol_conc_exp(mean_radius, bin_grid, vol_conc)
    
    !> Mean radius (m).
    real(kind=dp), intent(in) :: mean_radius
    !> Bin grid.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Volume concentration (V(ln(r))d(ln(r))).
    real(kind=dp), intent(out) :: vol_conc(bin_grid%n_bin)
    
    real(kind=dp) :: num_conc(bin_grid%n_bin)

    call num_conc_exp(mean_radius, bin_grid, num_conc)
    vol_conc = num_conc * bin_grid%v
    
  end subroutine vol_conc_exp
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Mono-disperse distribution.
  !> Normalized so that sum(num_conc(k) * dlnr) = 1
  subroutine num_conc_mono(radius, bin_grid, num_conc)
    
    !> Radius of each particle (m^3).
    real(kind=dp), intent(in) :: radius
    !> Bin grid.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Number concentration (#(ln(r))d(ln(r))).
    real(kind=dp), intent(out) :: num_conc(bin_grid%n_bin)
    
    integer :: k

    num_conc = 0d0
    k = bin_grid_particle_in_bin(bin_grid, rad2vol(radius))
    num_conc(k) = 1d0 / bin_grid%dlnr
    
  end subroutine num_conc_mono
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Mono-disperse distribution in volume.
  subroutine vol_conc_mono(radius, bin_grid, vol_conc)
    
    !> Radius of each particle (m^3).
    real(kind=dp), intent(in) :: radius
    !> Bin grid.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Volume concentration (V(ln(r))d(ln(r))).
    real(kind=dp), intent(out) :: vol_conc(bin_grid%n_bin)
    
    integer :: k

    vol_conc = 0d0
    k = bin_grid_particle_in_bin(bin_grid, rad2vol(radius))
    vol_conc(k) = 1d0 / bin_grid%dlnr * rad2vol(radius)
    
  end subroutine vol_conc_mono
  
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
    real(kind=dp), intent(out) :: num_conc(bin_grid%n_bin)

    if (aero_mode%type == AERO_MODE_TYPE_LOG_NORMAL) then
       call num_conc_log_normal(aero_mode%mean_radius, &
            aero_mode%log10_std_dev_radius, bin_grid, num_conc)
    elseif (aero_mode%type == AERO_MODE_TYPE_EXP) then
       call num_conc_exp(aero_mode%mean_radius, bin_grid, num_conc)
    elseif (aero_mode%type == AERO_MODE_TYPE_MONO) then
       call num_conc_mono(aero_mode%mean_radius, bin_grid, num_conc)
    else
       call die_msg(719625922, "unknown aero_mode type: " &
            // integer_to_string(aero_mode%type))
    end if
    num_conc = num_conc * aero_mode%num_conc

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
    real(kind=dp), intent(out) :: vol_conc(bin_grid%n_bin, aero_data%n_spec)

    integer :: i_spec
    real(kind=dp) :: vol_conc_total(bin_grid%n_bin)

    if (aero_mode%type == AERO_MODE_TYPE_LOG_NORMAL) then
       call vol_conc_log_normal(aero_mode%mean_radius, &
            aero_mode%log10_std_dev_radius, bin_grid, vol_conc_total)
    elseif (aero_mode%type == AERO_MODE_TYPE_EXP) then
       call vol_conc_exp(aero_mode%mean_radius, bin_grid, &
            vol_conc_total)
    elseif (aero_mode%type == AERO_MODE_TYPE_MONO) then
       call vol_conc_mono(aero_mode%mean_radius, bin_grid, &
            vol_conc_total)
    else
       call die_msg(314169653, "Unknown aero_mode type: " &
            // integer_to_string(aero_mode%type))
    end if
    vol_conc_total = vol_conc_total * aero_mode%num_conc
    do i_spec = 1,aero_data%n_spec
       vol_conc(:,i_spec) = vol_conc_total &
            * aero_mode%vol_frac(i_spec)
    end do

  end subroutine aero_mode_vol_conc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Return the weighted number concentration for an \c aero_mode.
  real(kind=dp) function aero_mode_weighted_num_conc(aero_mode, aero_weight)

    !> Aero_mode to sample radius from.
    type(aero_mode_t), intent(in) :: aero_mode
    !> Aero weight.
    type(aero_weight_t), intent(in) :: aero_weight

    real(kind=dp) :: x_mean_prime

    if (aero_weight%type == AERO_WEIGHT_TYPE_NONE) then
       aero_mode_weighted_num_conc = aero_mode%num_conc
    elseif ((aero_weight%type == AERO_WEIGHT_TYPE_POWER) &
         .or. (aero_weight%type == AERO_WEIGHT_TYPE_MFA)) then
       if (aero_mode%type == AERO_MODE_TYPE_LOG_NORMAL) then
          x_mean_prime = log10(aero_mode%mean_radius) &
               - aero_weight%exponent * aero_mode%log10_std_dev_radius**2 &
               * log(10d0)
          aero_mode_weighted_num_conc = aero_mode%num_conc &
               * aero_weight%ref_radius**aero_weight%exponent &
               * exp((x_mean_prime**2 - log10(aero_mode%mean_radius)**2) &
               / (2d0 * aero_mode%log10_std_dev_radius**2))
       elseif (aero_mode%type == AERO_MODE_TYPE_EXP) then
          call die_msg(822252601, "exp/power unimplemented")
       elseif (aero_mode%type == AERO_MODE_TYPE_MONO) then
          aero_mode_weighted_num_conc = aero_mode%num_conc &
               / aero_weight_value(aero_weight, aero_mode%mean_radius)
       else
          call die_msg(901140225, "unknown aero_mode type: " &
               // integer_to_string(aero_mode%type))
       end if
    else
       call die_msg(742383510, "unknown aero_weight type: " &
            // integer_to_string(aero_weight%type))
    end if

  end function aero_mode_weighted_num_conc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Return a radius randomly sampled from the mode distribution.
  subroutine aero_mode_sample_radius(aero_mode, aero_weight, radius)

    !> Aero_mode to sample radius from.
    type(aero_mode_t), intent(in) :: aero_mode
    !> Aero weight.
    type(aero_weight_t), intent(in) :: aero_weight
    !> Sampled radius (m).
    real(kind=dp), intent(out) :: radius

    real(kind=dp) :: x_mean_prime

    if (aero_weight%type == AERO_WEIGHT_TYPE_NONE) then
       if (aero_mode%type == AERO_MODE_TYPE_LOG_NORMAL) then
          radius = 10d0**rand_normal(log10(aero_mode%mean_radius), &
               aero_mode%log10_std_dev_radius)
       elseif (aero_mode%type == AERO_MODE_TYPE_EXP) then
          radius = vol2rad(- rad2vol(aero_mode%mean_radius) * log(pmc_random()))
       elseif (aero_mode%type == AERO_MODE_TYPE_MONO) then
          radius = aero_mode%mean_radius
       else
          call die_msg(749122931, "Unknown aero_mode type: " &
               // integer_to_string(aero_mode%type))
       end if
    elseif ((aero_weight%type == AERO_WEIGHT_TYPE_POWER) &
         .or. (aero_weight%type == AERO_WEIGHT_TYPE_MFA)) then
       if (aero_mode%type == AERO_MODE_TYPE_LOG_NORMAL) then
          x_mean_prime = log10(aero_mode%mean_radius) &
               - aero_weight%exponent * aero_mode%log10_std_dev_radius**2 &
               * log(10d0)
          radius = 10d0**rand_normal(x_mean_prime, &
               aero_mode%log10_std_dev_radius)
       elseif (aero_mode%type == AERO_MODE_TYPE_EXP) then
          call die_msg(111024862, "exp/power unimplemented")
       elseif (aero_mode%type == AERO_MODE_TYPE_MONO) then
          radius = aero_mode%mean_radius
       else
          call die_msg(886417976, "unknown aero_mode type: " &
               // integer_to_string(aero_mode%type))
       end if
    else
       call die_msg(863127819, "unknown aero_weight type: " &
            // integer_to_string(aero_weight%type))
    end if

  end subroutine aero_mode_sample_radius

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read volume fractions from a data file.
  subroutine spec_file_read_vol_frac(file, aero_data, vol_frac)

    !> Spec file to read mass fractions from.
    type(spec_file_t), intent(inout) :: file
    !> Aero_data data.
    type(aero_data_t), intent(in) :: aero_data
    !> Aerosol species volume fractions.
    real(kind=dp), intent(inout) :: vol_frac(:)

    integer :: n_species, species, i
    character(len=SPEC_LINE_MAX_VAR_LEN), pointer :: species_name(:)
    real(kind=dp), pointer :: species_data(:,:)
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
    !! indicating that the diesel particles are 30% organic carbon and
    !! 70% black carbon.
    !!
    !! See also:
    !!   - \ref spec_file_format --- the input file text format
    !!   - \ref input_format_aero_dist --- the format for a complete
    !!     aerosol distribution with several modes
    !!   - \ref input_format_aero_mode --- the format for each mode
    !!     of an aerosol distribution

    ! read the aerosol data from the specified file
    allocate(species_name(0))
    allocate(species_data(0,0))
    call spec_file_read_real_named_array(file, 0, species_name, &
         species_data)

    ! check the data size
    n_species = size(species_data, 1)
    if (n_species < 1) then
       call die_msg(427666881, 'file ' // trim(file%name) &
            // ' must contain at least one line of data')
    end if
    if (size(species_data, 2) /= 1) then
       call die_msg(427666881, 'each line in file ' &
            // trim(file%name) // ' must contain exactly one data value')
    end if

    ! copy over the data
    vol_frac = 0d0
    do i = 1,n_species
       species = aero_data_spec_by_name(aero_data, species_name(i))
       if (species == 0) then
          call die_msg(775942501, 'unknown species ' &
               // trim(species_name(i)) // ' in file ' &
               // trim(file%name))
       end if
       vol_frac(species) = species_data(i,1)
    end do
    deallocate(species_name)
    deallocate(species_data)

    ! convert mass fractions to volume fractions
    vol_frac = vol_frac / aero_data%density
    
    ! normalize
    tot_vol_frac = sum(vol_frac)
    if ((minval(vol_frac) < 0d0) .or. (tot_vol_frac <= 0d0)) then
       call die_msg(356648030, 'vol_frac in ' // trim(file%name) &
            // ' is not positive')
    end if
    vol_frac = vol_frac / tot_vol_frac

  end subroutine spec_file_read_vol_frac

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read one mode of an aerosol distribution (number concentration,
  !> volume fractions, and mode shape).
  subroutine spec_file_read_aero_mode(file, aero_data, aero_mode, eof)

    !> Spec file.
    type(spec_file_t), intent(inout) :: file
    !> Aero_data data.
    type(aero_data_t), intent(in) :: aero_data
    !> Aerosol mode.
    type(aero_mode_t), intent(inout) :: aero_mode
    !> If eof instead of reading data.
    logical :: eof

    character(len=SPEC_LINE_MAX_VAR_LEN) :: tmp_str, mode_type
    character(len=SPEC_LINE_MAX_VAR_LEN) :: mass_frac_filename
    type(spec_line_t) :: line
    type(spec_file_t) :: mass_frac_file

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
    !! <li> \b num_conc (real, unit 1/m^3): the total number
    !!      concentration \f$N_{\rm total}\f$ of the mode
    !! <li> \b mode_type (string): the functional form of the mode ---
    !!      must be one of: \c log_normal for a log-normal distribution;
    !!      \c exp for an exponential distribution; or \c mono for a
    !!      mono-disperse distribution
    !! <li> if \c mode_type is \c log_normal then the mode distribution
    !!      shape is
    !!      \f[ n(\log D) {\rm d}\log D
    !!      = \frac{N_{\rm total}}{\sqrt{2\pi} \log \sigma_{\rm g}}
    !!      \exp\left(\frac{(\log D - \log D_{\rm \mu g})^2}{2 \log ^2 \sigma_{\rm g}}\right)
    !!      {\rm d}\log D \f]
    !!      and the following parameters are:
    !!      <ul>
    !!      <li> \b mean_radius (real, unit m): the geometric mean radius
    !!           \f$R_{\rm \mu g}\f$ such that \f$D_{\rm \mu g} = 2 R_{\rm
    !!           \mu g}\f$
    !!      <li> \b log_std_dev (real, dimensionless): \f$\log_{10}\f$ of the geometric
    !!           standard deviation \f$\sigma_{\rm g}\f$
    !!      </ul>
    !! <li> if \c mode_type is \c exp then the mode distribution shape is
    !!      \f[ n(v) {\rm d}v = \frac{N_{\rm total}}{v_{\rm \mu}}
    !!      \exp\left(- \frac{v}{v_{\rm \mu}}\right)
    !!      {\rm d}v \f]
    !!      and the following parameters are:
    !!      <ul>
    !!      <li> \b mean_radius (real, unit m): the mean radius
    !!           \f$R_{\rm \mu}\f$ such that \f$v_{\rm \mu} = \frac{4}{3} \pi R^3_{\rm
    !!           \mu}\f$
    !!      </ul>
    !! <li> if \c mode_type is \c mono then the mode distribution shape
    !!      is a delta distribution at diameter \f$D_0\f$ and the
    !!      following parameters are:
    !!      <ul>
    !!      <li> \b radius (real, unit m): the radius \f$R_0\f$ of the
    !!           particles, so that \f$D_0 = 2 R_0\f$
    !!      </ul>
    !! </ul>
    !!
    !! Example:
    !! <pre>
    !! mode_name diesel          # mode name (descriptive only)
    !! mass_frac comp_diesel.dat # mass fractions in each aerosol particle
    !! num_conc 1.6e8            # particle number density (#/m^3)
    !! mode_type log_normal      # type of distribution
    !! mean_radius 2.5e-08       # mean radius (m)
    !! log_std_dev 0.24          # log_10 of geometric standard deviation
    !! </pre>
    !!
    !! See also:
    !!   - \ref spec_file_format --- the input file text format
    !!   - \ref input_format_aero_dist --- the format for a complete
    !!     aerosol distribution with several modes
    !!   - \ref input_format_mass_frac --- the format for the mass
    !!     fractions file

    call spec_line_allocate(line)
    call spec_file_read_line(file, line, eof)
    if (.not. eof) then
       call spec_file_check_line_name(file, line, "mode_name")
       call spec_file_check_line_length(file, line, 1)
       call aero_mode_deallocate(aero_mode)
       call aero_mode_allocate_size(aero_mode, aero_data%n_spec)
       tmp_str = line%data(1) ! hack to avoid gfortran warning
       aero_mode%name = tmp_str(1:AERO_MODE_NAME_LEN)
       call spec_file_read_string(file, 'mass_frac', mass_frac_filename)
       call spec_file_open(mass_frac_filename, mass_frac_file)
       call spec_file_read_vol_frac(mass_frac_file, aero_data, &
            aero_mode%vol_frac)
       call spec_file_close(mass_frac_file)
       call spec_file_read_real(file, 'num_conc', aero_mode%num_conc)
       call spec_file_read_string(file, 'mode_type', mode_type)
       if (trim(mode_type) == 'log_normal') then
          aero_mode%type = AERO_MODE_TYPE_LOG_NORMAL
          call spec_file_read_real(file, 'mean_radius', aero_mode%mean_radius)
          call spec_file_read_real(file, 'log_std_dev', &
               aero_mode%log10_std_dev_radius)
       elseif (trim(mode_type) == 'exp') then
          aero_mode%type = AERO_MODE_TYPE_EXP
          call spec_file_read_real(file, 'mean_radius', aero_mode%mean_radius)
       elseif (trim(mode_type) == 'mono') then
          aero_mode%type = AERO_MODE_TYPE_MONO
          call spec_file_read_real(file, 'radius', aero_mode%mean_radius)
       else
          call spec_file_die_msg(729472928, file, &
               "Unknown distribution mode type: " // trim(mode_type))
       end if
    end if
    call spec_line_deallocate(line)

  end subroutine spec_file_read_aero_mode

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Determines the number of bytes required to pack the given value.
  integer function pmc_mpi_pack_size_aero_mode(val)

    !> Value to pack.
    type(aero_mode_t), intent(in) :: val

    pmc_mpi_pack_size_aero_mode = &
         pmc_mpi_pack_size_string(val%name) &
         + pmc_mpi_pack_size_integer(val%type) &
         + pmc_mpi_pack_size_real(val%mean_radius) &
         + pmc_mpi_pack_size_real(val%log10_std_dev_radius) &
         + pmc_mpi_pack_size_real(val%num_conc) &
         + pmc_mpi_pack_size_real_array(val%vol_frac)

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
    call pmc_mpi_pack_real(buffer, position, val%mean_radius)
    call pmc_mpi_pack_real(buffer, position, val%log10_std_dev_radius)
    call pmc_mpi_pack_real(buffer, position, val%num_conc)
    call pmc_mpi_pack_real_array(buffer, position, val%vol_frac)
    call assert(579699255, &
         position - prev_position == pmc_mpi_pack_size_aero_mode(val))
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
    call pmc_mpi_unpack_real(buffer, position, val%mean_radius)
    call pmc_mpi_unpack_real(buffer, position, val%log10_std_dev_radius)
    call pmc_mpi_unpack_real(buffer, position, val%num_conc)
    call pmc_mpi_unpack_real_array(buffer, position, val%vol_frac)
    call assert(874467577, &
         position - prev_position == pmc_mpi_pack_size_aero_mode(val))
#endif

  end subroutine pmc_mpi_unpack_aero_mode

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_aero_mode
