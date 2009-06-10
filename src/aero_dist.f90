! Copyright (C) 2005-2009 Nicole Riemer and Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_aero_dist module.

!> The aero_dist_t structure and associated subroutines.
!!
!! The initial size distributions are computed as number densities, so
!! they can be used for both sectional and particle-resolved
!! simulations. The routine dist_to_n() converts a number density
!! distribution to an actual number of particles ready for a
!! particle-resolved simulation.
!!
!! Initial distributions should be normalized so that <tt>sum(n_den) =
!! 1/dlnr</tt>.
module pmc_aero_dist

  use pmc_bin_grid
  use pmc_util
  use pmc_constants
  use pmc_spec_read
  use pmc_aero_data
  use pmc_mpi
  use pmc_rand
#ifdef PMC_USE_MPI
  use mpi
#endif

  !> Maximum length of an aero_dist mode name.
  integer, parameter :: AERO_DIST_NAME_LEN = 300
  !> Maximum length of an aero_dist mode type.
  integer, parameter :: AERO_DIST_TYPE_LEN = 100

  !> An aerosol size distribution mode.
  !!
  !! Each mode is assumed to be fully internally mixed so that every
  !! particle has the same composition. The \c num_den array then
  !! stores the number density distribution.
  type aero_mode_t
     !> Mode name, used to track particle sources.
     character(len=AERO_DIST_NAME_LEN) :: name
     !> Mode type ("log_normal", "exp", or "mono").
     character(len=AERO_DIST_TYPE_LEN) :: type
     !> Mean radius of mode (m).
     real*8 :: mean_radius
     !> Log base 10 of geometric standard deviation of radius, if necessary (m).
     real*8 :: log10_std_dev_radius
     !> Total number density of mode (#/m^3).
     real*8 :: num_den
     !> Species fractions by volume [length \c aero_data%%n_spec] (1).
     real*8, pointer :: vol_frac(:)
  end type aero_mode_t

  !> A complete aerosol distribution con
  type aero_dist_t
     !> Number of modes.
     integer :: n_mode
     !> Internally mixed modes [length \c n_mode].
     type(aero_mode_t), pointer :: mode(:)
  end type aero_dist_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Allocates an aero_mode.
  subroutine aero_mode_alloc(aero_mode, n_spec)

    !> Aerosol mode.
    type(aero_mode_t), intent(out) :: aero_mode
    !> Number of species.
    integer, intent(in) :: n_spec

    allocate(aero_mode%vol_frac(n_spec))

  end subroutine aero_mode_alloc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Free all storage.
  subroutine aero_mode_free(aero_mode)

    !> Aerosol mode.
    type(aero_mode_t), intent(inout) :: aero_mode

    deallocate(aero_mode%vol_frac)

  end subroutine aero_mode_free

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Copy an aero_mode.
  subroutine aero_mode_copy(aero_mode_from, aero_mode_to)

    !> Aerosol mode original.
    type(aero_mode_t), intent(in) :: aero_mode_from
    !> Aerosol mode copy.
    type(aero_mode_t), intent(inout) :: aero_mode_to

    call aero_mode_free(aero_mode_to)
    call aero_mode_alloc(aero_mode_to, size(aero_mode_from%vol_frac))
    aero_mode_to%name = aero_mode_from%name
    aero_mode_to%type = aero_mode_from%type
    aero_mode_to%mean_radius = aero_mode_from%mean_radius
    aero_mode_to%log10_std_dev_radius = aero_mode_from%log10_std_dev_radius
    aero_mode_to%num_den = aero_mode_from%num_den
    aero_mode_to%vol_frac = aero_mode_from%vol_frac

  end subroutine aero_mode_copy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Allocates an aero_dist.
  subroutine aero_dist_alloc(aero_dist, n_mode, n_spec)

    !> Aerosol distribution.
    type(aero_dist_t), intent(out) :: aero_dist
    !> Number of modes.
    integer, intent(in) :: n_mode
    !> Number of species.
    integer, intent(in) :: n_spec

    integer :: i

    aero_dist%n_mode = n_mode
    allocate(aero_dist%mode(n_mode))
    do i = 1,n_mode
       call aero_mode_alloc(aero_dist%mode(i), n_spec)
    end do

  end subroutine aero_dist_alloc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Free all storage.
  subroutine aero_dist_free(aero_dist)

    !> Aerosol distribution.
    type(aero_dist_t), intent(inout) :: aero_dist

    integer :: i

    do i = 1,aero_dist%n_mode
       call aero_mode_free(aero_dist%mode(i))
    end do
    deallocate(aero_dist%mode)

  end subroutine aero_dist_free

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Copy an aero_dist.
  subroutine aero_dist_copy(aero_dist_from, aero_dist_to)

    !> Aero_dist original.
    type(aero_dist_t), intent(in) :: aero_dist_from
    !> Aero_dist copy.
    type(aero_dist_t), intent(inout) :: aero_dist_to

    integer :: n_spec, i

    if (aero_dist_from%n_mode > 0) then
       n_spec = size(aero_dist_from%mode(1)%vol_frac)
    else
       n_spec = 0
    end if
    call aero_dist_free(aero_dist_to)
    call aero_dist_alloc(aero_dist_to, aero_dist_from%n_mode, n_spec)
    do i = 1,aero_dist_from%n_mode
       call aero_mode_copy(aero_dist_from%mode(i), aero_dist_to%mode(i))
    end do

  end subroutine aero_dist_copy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the total number concentration in #/m^3 of a distribution.
  !> (#/m^3)
  real*8 function aero_dist_total_num_den(aero_dist)

    !> Aerosol distribution.
    type(aero_dist_t), intent(in) :: aero_dist

    aero_dist_total_num_den = sum(aero_dist%mode%num_den)

  end function aero_dist_total_num_den

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Compute a log-normal distribution, normalized so that
  !> sum(num_den(k) * dlnr) = 1
  subroutine num_den_log_normal(mean_radius, log_sigma, bin_grid, num_den)
    
    !> Geometric mean radius (m).
    real*8, intent(in) :: mean_radius
    !> log_10(geom. std dev) (1).
    real*8, intent(in) :: log_sigma
    !> Bin grid.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Normalized number density (#(ln(r))d(ln(r))).
    real*8, intent(out) :: num_den(bin_grid%n_bin)
    
    integer :: k
    
    do k = 1,bin_grid%n_bin
       num_den(k) = 1d0 / (sqrt(2d0 * const%pi) * log_sigma) * &
            dexp(-(dlog10(vol2rad(bin_grid%v(k))) &
            - dlog10(mean_radius))**2d0 &
            / (2d0 * log_sigma**2d0)) / dlog(10d0)
    end do
    
    ! The formula above was originally for a distribution in
    ! log_10(r), while we are using log_e(r) for our bin grid. The
    ! division by dlog(10) at the end corrects for this.

    ! Remember that log_e(r) = log_10(r) * log_e(10).
    
  end subroutine num_den_log_normal
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Compute a log-normal distribution in volume.
  subroutine vol_den_log_normal(mean_radius, log_sigma, bin_grid, vol_den)
    
    !> Geometric mean radius (m).
    real*8, intent(in) :: mean_radius
    !> log_10(geom. std dev) (1).
    real*8, intent(in) :: log_sigma
    !> Bin grid.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Normalized volume density (#(ln(r))d(ln(r))).
    real*8, intent(out) :: vol_den(bin_grid%n_bin)
    
    real*8 :: num_den(bin_grid%n_bin)

    call num_den_log_normal(mean_radius, log_sigma, bin_grid, num_den)
    vol_den = num_den * bin_grid%v
    
  end subroutine vol_den_log_normal
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Exponential distribution in volume
  !> n(v) = 1 / mean_vol * exp(- v / mean_vol)
  !> Normalized so that sum(num_den(k) * dlnr) = 1
  subroutine num_den_exp(mean_radius, bin_grid, num_den)
    
    !> Mean radius (m).
    real*8, intent(in) :: mean_radius
    !> Bin grid.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Num den (#(ln(r))d(ln(r))).
    real*8, intent(out) :: num_den(bin_grid%n_bin)
    
    integer :: k
    real*8 :: mean_vol, num_den_vol
    
    mean_vol = rad2vol(mean_radius)
    do k = 1,bin_grid%n_bin
       num_den_vol = 1d0 / mean_vol * exp(-(bin_grid%v(k) / mean_vol))
       call vol_to_lnr(vol2rad(bin_grid%v(k)), num_den_vol, num_den(k))
    end do
    
  end subroutine num_den_exp
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Exponential distribution in volume.
  subroutine vol_den_exp(mean_radius, bin_grid, vol_den)
    
    !> Mean radius (m).
    real*8, intent(in) :: mean_radius
    !> Bin grid.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Num den (#(ln(r))d(ln(r))).
    real*8, intent(out) :: vol_den(bin_grid%n_bin)
    
    real*8 :: num_den(bin_grid%n_bin)

    call num_den_exp(mean_radius, bin_grid, num_den)
    vol_den = num_den * bin_grid%v
    
  end subroutine vol_den_exp
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Mono-disperse distribution.
  !> Normalized so that sum(num_den(k) * dlnr) = 1
  subroutine num_den_mono(radius, bin_grid, num_den)
    
    !> Radius of each particle (m^3).
    real*8, intent(in) :: radius
    !> Bin grid.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Num den (#(ln(r))d(ln(r))).
    real*8, intent(out) :: num_den(bin_grid%n_bin)
    
    integer :: k

    num_den = 0d0
    k = bin_grid_particle_in_bin(bin_grid, rad2vol(radius))
    num_den(k) = 1d0 / bin_grid%dlnr
    
  end subroutine num_den_mono
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Mono-disperse distribution in volume.
  subroutine vol_den_mono(radius, bin_grid, vol_den)
    
    !> Radius of each particle (m^3).
    real*8, intent(in) :: radius
    !> Bin grid.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Num den (#(ln(r))d(ln(r))).
    real*8, intent(out) :: vol_den(bin_grid%n_bin)
    
    integer :: k

    vol_den = 0d0
    k = bin_grid_particle_in_bin(bin_grid, rad2vol(radius))
    vol_den(k) = 1d0 / bin_grid%dlnr * rad2vol(radius)
    
  end subroutine vol_den_mono
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Return a radius randomly sampled from the mode distribution.
  subroutine aero_mode_sample_radius(aero_mode, radius)

    !> Aero_mode to sample radius from.
    type(aero_mode_t), intent(in) :: aero_mode
    !> Sampled radius (m).
    real*8, intent(out) :: radius

    if (aero_mode%type == "log_normal") then
       radius = 10d0**rand_normal(log10(aero_mode%mean_radius), &
            aero_mode%log10_std_dev_radius)
    elseif (aero_mode%type == "exp") then
       radius = vol2rad(- rad2vol(aero_mode%mean_radius) * log(pmc_random()))
    elseif (aero_mode%type == "mono") then
       radius = aero_mode%mean_radius
    else
       call die_msg(749122931, "Unknown aero_mode type")
    end if

  end subroutine aero_mode_sample_radius

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Determine the current aero_dist and rate by interpolating at the
  !> current time with the lists of aero_dists and rates.
  subroutine aero_dist_interp_1d(aero_dist_list, time_list, &
         rate_list, time, aero_dist, rate)

    !> Gas states.
    type(aero_dist_t), intent(in) :: aero_dist_list(:)
    !> Times (s).
    real*8, intent(in) :: time_list(size(aero_dist_list))
    !> Rates (s^{-1}).
    real*8, intent(in) :: rate_list(size(aero_dist_list))
    !> Current time (s).
    real*8, intent(in) :: time
    !> Current gas state.
    type(aero_dist_t), intent(inout) :: aero_dist
    !> Current rate (s^{-1}).
    real*8, intent(out) :: rate

    integer :: n, p, n_bin, n_spec, i, i_new
    real*8 :: y, alpha

    n = size(aero_dist_list)
    p = find_1d(n, time_list, time)
    if (p == 0) then
       ! before the start, just use the first state and rate
       call aero_dist_copy(aero_dist_list(1), aero_dist)
       rate = rate_list(1)
    elseif (p == n) then
       ! after the end, just use the last state and rate
       call aero_dist_copy(aero_dist_list(n), aero_dist)
       rate = rate_list(n)
    else
       ! in the middle, use the previous dist
       call aero_dist_copy(aero_dist_list(p), aero_dist)
       rate = rate_list(p)
    end if

  end subroutine aero_dist_interp_1d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read volume fractions from a data file.
  subroutine spec_read_vol_frac(file, aero_data, vol_frac)

    !> Spec file.
    type(spec_file_t), intent(inout) :: file
    !> Aero_data data.
    type(aero_data_t), intent(in) :: aero_data
    !> Aerosol species volume fractions.
    real*8, intent(out) :: vol_frac(:)

    integer :: n_species, species, i
    character(len=MAX_VAR_LEN) :: read_name
    type(spec_file_t) :: read_file
    character(len=MAX_VAR_LEN), pointer :: species_name(:)
    character(len=MAX_VAR_LEN) :: frac_type
    real*8, pointer :: species_data(:,:)
    real*8 :: tot_vol_frac

    ! read the aerosol data from the specified file
    call spec_read_string(file, 'frac_type', frac_type)
    if ((frac_type /= "volume") &
         .and. (frac_type /= "mass") &
         .and. (frac_type /= "mole")) then
       call spec_read_die_msg(842919812,file, "unknown frac_type " &
            // "(should be 'volume', 'mass', or 'mole'): " // trim(frac_type))
    end if
    call spec_read_string(file, 'frac', read_name)
    call spec_read_open(read_name, read_file)
    call spec_read_real_named_array(read_file, 0, species_name, species_data)
    call spec_read_close(read_file)

    ! check the data size
    n_species = size(species_data, 1)
    if (n_species < 1) then
       call die_msg(427666881, 'file ' // trim(read_name) &
            // ' must contain at least one line of data')
    end if
    if (size(species_data, 2) /= 1) then
       call die_msg(427666881, 'each line in file ' &
            // trim(read_name) // ' must contain exactly one data value')
    end if

    ! copy over the data
    vol_frac = 0d0
    do i = 1,n_species
       species = aero_data_spec_by_name(aero_data, species_name(i))
       if (species == 0) then
          call die_msg(775942501, 'unknown species ' &
               // trim(species_name(i)) // ' in file ' &
               // trim(read_name))
       end if
       vol_frac(species) = species_data(i,1)
    end do
    deallocate(species_name)
    deallocate(species_data)

    ! convert to appropriate fraction type
    if (frac_type == "mass") then
       vol_frac = vol_frac / aero_data%density
    elseif (frac_type == "mole") then
       vol_frac = vol_frac / aero_data%density * aero_data%molec_weight
    elseif (frac_type /= "volume") then
       call die_msg(719123834, "unknown frac_type: " // trim(frac_type))
    end if
    
    ! normalize
    tot_vol_frac = sum(vol_frac)
    if ((minval(vol_frac) < 0d0) .or. (tot_vol_frac <= 0d0)) then
       call die_msg(356648030, 'vol_frac in ' // trim(read_name) &
            // ' is not positive')
    end if
    vol_frac = vol_frac / tot_vol_frac

  end subroutine spec_read_vol_frac

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read the shape (number density profile) of one mode of an aerosol
  !> distribution.
  subroutine spec_read_aero_mode_shape(file, aero_mode)

    !> Spec file.
    type(spec_file_t), intent(inout) :: file
    !> Aerosol mode.
    type(aero_mode_t), intent(inout) :: aero_mode

    character(len=MAX_VAR_LEN) :: mode_type

    call spec_read_real(file, 'num_den', aero_mode%num_den)
    call spec_read_string(file, 'mode_type', mode_type)
    if (len_trim(mode_type) < AERO_DIST_TYPE_LEN) then
       aero_mode%type = mode_type(1:AERO_DIST_TYPE_LEN)
    else
       call spec_read_die_msg(284789262, file, "mode_type string too long")
    end if
    if (trim(mode_type) == 'log_normal') then
       call spec_read_real(file, 'mean_radius', aero_mode%mean_radius)
       call spec_read_real(file, 'log_std_dev', aero_mode%log10_std_dev_radius)
    elseif (trim(mode_type) == 'exp') then
       call spec_read_real(file, 'mean_radius', aero_mode%mean_radius)
    elseif (trim(mode_type) == 'mono') then
       call spec_read_real(file, 'radius', aero_mode%mean_radius)
    else
       call spec_read_die_msg(729472928, file, "Unknown distribution type")
    end if

  end subroutine spec_read_aero_mode_shape

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read one mode of an aerosol distribution (number density and
  !> volume fractions).
  subroutine spec_read_aero_mode(file, aero_data, aero_mode, eof)

    !> Spec file.
    type(spec_file_t), intent(inout) :: file
    !> Aero_data data.
    type(aero_data_t), intent(in) :: aero_data
    !> Aerosol mode (will be allocated).
    type(aero_mode_t), intent(inout) :: aero_mode
    !> If eof instead of reading data.
    logical :: eof

    character(len=MAX_VAR_LEN) :: tmp_str
    type(spec_line_t) :: line

    call spec_read_line(file, line, eof)
    if (.not. eof) then
       call spec_read_check_line_name(file, line, "mode_name")
       call spec_read_check_line_length(file, line, 1)
       call aero_mode_alloc(aero_mode, aero_data%n_spec)
       tmp_str = line%data(1) ! hack to avoid gfortran warning
       aero_mode%name = tmp_str(1:AERO_DIST_NAME_LEN)
       allocate(aero_mode%vol_frac(aero_data%n_spec))
       call spec_read_vol_frac(file, aero_data, aero_mode%vol_frac)
       call spec_read_aero_mode_shape(file, aero_mode)
       call spec_line_free(line)
    end if

  end subroutine spec_read_aero_mode

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read continuous aerosol distribution composed of several modes.
  subroutine spec_read_aero_dist(file, aero_data, aero_dist)

    !> Spec file.
    type(spec_file_t), intent(inout) :: file
    !> Aero_data data.
    type(aero_data_t), intent(in) :: aero_data
    !> Aerosol dist, will be allocated.
    type(aero_dist_t), intent(inout) :: aero_dist

    type(aero_mode_t), pointer :: new_aero_mode_list(:)
    type(aero_mode_t) :: aero_mode
    integer :: i, j
    logical :: eof
    
    aero_dist%n_mode = 0
    allocate(aero_dist%mode(aero_dist%n_mode))
    call spec_read_aero_mode(file, aero_data, aero_mode, eof)
    do while (.not. eof)
       aero_dist%n_mode = aero_dist%n_mode + 1
       allocate(new_aero_mode_list(aero_dist%n_mode))
       do i = 1,aero_dist%n_mode
          call aero_mode_alloc(new_aero_mode_list(i), 0)
       end do
       call aero_mode_copy(aero_mode, &
            new_aero_mode_list(aero_dist%n_mode))
       call aero_mode_free(aero_mode)
       do i = 1,(aero_dist%n_mode - 1)
          call aero_mode_copy(aero_dist%mode(i), &
               new_aero_mode_list(i))
          call aero_mode_free(aero_dist%mode(i))
       end do
       deallocate(aero_dist%mode)
       aero_dist%mode => new_aero_mode_list
       nullify(new_aero_mode_list)
       call spec_read_aero_mode(file, aero_data, aero_mode, eof)
    end do

  end subroutine spec_read_aero_dist

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read aerosol distribution from filename on line in file.
  subroutine spec_read_aero_dist_filename(file, aero_data, bin_grid, &
       name, aero_dist)

    !> Spec file.
    type(spec_file_t), intent(inout) :: file
    !> Aero_data data.
    type(aero_data_t), intent(in) :: aero_data
    !> Bin grid.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Name of data line for filename.
    character(len=*), intent(in) :: name
    !> Aerosol distribution.
    type(aero_dist_t), intent(inout) :: aero_dist

    character(len=MAX_VAR_LEN) :: read_name
    type(spec_file_t) :: read_file

    ! read the aerosol data from the specified file
    call spec_read_string(file, name, read_name)
    call spec_read_open(read_name, read_file)
    call spec_read_aero_dist(read_file, aero_data, aero_dist)
    call spec_read_close(read_file)

  end subroutine spec_read_aero_dist_filename
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read an array of aero_dists with associated times and rates from
  !> the given file.
  subroutine spec_read_aero_dists_times_rates(file, aero_data, &
       bin_grid, name, times, rates, aero_dists)

    !> Spec file.
    type(spec_file_t), intent(inout) :: file
    !> Aero data.
    type(aero_data_t), intent(in) :: aero_data
    !> Bin grid.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Name of data line for filename.
    character(len=*), intent(in) :: name
    !> Times (s).
    real*8, pointer :: times(:)
    !> Rates (s^{-1}).
    real*8, pointer :: rates(:)
    !> Aero dists.
    type(aero_dist_t), pointer :: aero_dists(:)

    character(len=MAX_VAR_LEN) :: read_name
    type(spec_file_t) :: read_file
    type(spec_line_t) :: aero_dist_line
    integer :: n_time, i_time
    character(len=MAX_VAR_LEN), pointer :: names(:)
    real*8, pointer :: data(:,:)

    ! read the filename then read the data from that file
    call spec_read_string(file, name, read_name)
    call spec_read_open(read_name, read_file)
    call spec_read_real_named_array(read_file, 2, names, data)
    call spec_read_line_no_eof(read_file, aero_dist_line)
    call spec_read_check_line_name(read_file, aero_dist_line, "dist")
    call spec_read_check_line_length(read_file, aero_dist_line, size(data, 2))
    call spec_read_close(read_file)

    ! check the data size
    if (trim(names(1)) /= 'time') then
       call die_msg(570205795, 'row 1 in ' // trim(read_name) &
            // ' must start with: time not: ' // trim(names(1)))
    end if
    if (trim(names(2)) /= 'rate') then
       call die_msg(221270915, 'row 2 in ' // trim(read_name) &
            // ' must start with: rate not: ' // trim(names(1)))
    end if
    n_time = size(data, 2)
    if (n_time < 1) then
       call die_msg(457229710, 'each line in ' // trim(read_name) &
            // ' must contain at least one data value')
    end if

    ! copy over the data
    allocate(aero_dists(n_time))
    allocate(times(n_time))
    allocate(rates(n_time))
    do i_time = 1,n_time
       call spec_read_open(aero_dist_line%data(i_time), read_file)
       call spec_read_aero_dist(read_file, aero_data, aero_dists(i_time))
       call spec_read_close(read_file)
       times(i_time) = data(1,i_time)
       rates(i_time) = data(2,i_time)
    end do
    deallocate(names)
    deallocate(data)
    call spec_line_free(aero_dist_line)

  end subroutine spec_read_aero_dists_times_rates

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Determines the number of bytes required to pack the given value.
  integer function pmc_mpi_pack_size_aero_mode(val)

    !> Value to pack.
    type(aero_mode_t), intent(in) :: val

    pmc_mpi_pack_size_aero_mode = &
         pmc_mpi_pack_size_string(val%name) &
         + pmc_mpi_pack_size_string(val%type) &
         + pmc_mpi_pack_size_real(val%mean_radius) &
         + pmc_mpi_pack_size_real(val%log10_std_dev_radius) &
         + pmc_mpi_pack_size_real(val%num_den) &
         + pmc_mpi_pack_size_real_array(val%vol_frac)

  end function pmc_mpi_pack_size_aero_mode

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Determines the number of bytes required to pack the given value.
  integer function pmc_mpi_pack_size_aero_dist(val)

    !> Value to pack.
    type(aero_dist_t), intent(in) :: val

    integer :: i, total_size

    total_size = pmc_mpi_pack_size_integer(val%n_mode)
    do i = 1,size(val%mode)
       total_size = total_size + pmc_mpi_pack_size_aero_mode(val%mode(i))
    end do
    pmc_mpi_pack_size_aero_dist = total_size

  end function pmc_mpi_pack_size_aero_dist

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
    call pmc_mpi_pack_string(buffer, position, val%type)
    call pmc_mpi_pack_real(buffer, position, val%mean_radius)
    call pmc_mpi_pack_real(buffer, position, val%log10_std_dev_radius)
    call pmc_mpi_pack_real(buffer, position, val%num_den)
    call pmc_mpi_pack_real_array(buffer, position, val%vol_frac)
    call assert(579699255, &
         position - prev_position == pmc_mpi_pack_size_aero_mode(val))
#endif

  end subroutine pmc_mpi_pack_aero_mode

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Packs the given value into the buffer, advancing position.
  subroutine pmc_mpi_pack_aero_dist(buffer, position, val)

    !> Memory buffer.
    character, intent(inout) :: buffer(:)
    !> Current buffer position.
    integer, intent(inout) :: position
    !> Value to pack.
    type(aero_dist_t), intent(in) :: val

#ifdef PMC_USE_MPI
    integer :: prev_position, i

    prev_position = position
    call pmc_mpi_pack_integer(buffer, position, val%n_mode)
    do i = 1,size(val%mode)
       call pmc_mpi_pack_aero_mode(buffer, position, val%mode(i))
    end do
    call assert(440557910, &
         position - prev_position == pmc_mpi_pack_size_aero_dist(val))
#endif

  end subroutine pmc_mpi_pack_aero_dist

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Unpacks the given value from the buffer, advancing position.
  subroutine pmc_mpi_unpack_aero_mode(buffer, position, val)

    !> Memory buffer.
    character, intent(inout) :: buffer(:)
    !> Current buffer position.
    integer, intent(inout) :: position
    !> Value to pack.
    type(aero_mode_t), intent(out) :: val

#ifdef PMC_USE_MPI
    integer :: prev_position

    prev_position = position
    call pmc_mpi_unpack_string(buffer, position, val%name)
    call pmc_mpi_unpack_string(buffer, position, val%type)
    call pmc_mpi_unpack_real(buffer, position, val%mean_radius)
    call pmc_mpi_unpack_real(buffer, position, val%log10_std_dev_radius)
    call pmc_mpi_unpack_real(buffer, position, val%num_den)
    call pmc_mpi_unpack_real_array(buffer, position, val%vol_frac)
    call assert(874467577, &
         position - prev_position == pmc_mpi_pack_size_aero_mode(val))
#endif

  end subroutine pmc_mpi_unpack_aero_mode

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Unpacks the given value from the buffer, advancing position.
  subroutine pmc_mpi_unpack_aero_dist(buffer, position, val)

    !> Memory buffer.
    character, intent(inout) :: buffer(:)
    !> Current buffer position.
    integer, intent(inout) :: position
    !> Value to pack.
    type(aero_dist_t), intent(out) :: val

#ifdef PMC_USE_MPI
    integer :: prev_position, i

    prev_position = position
    call pmc_mpi_unpack_integer(buffer, position, val%n_mode)
    allocate(val%mode(val%n_mode))
    do i = 1,size(val%mode)
       call pmc_mpi_unpack_aero_mode(buffer, position, val%mode(i))
    end do
    call assert(742535268, &
         position - prev_position == pmc_mpi_pack_size_aero_dist(val))
#endif

  end subroutine pmc_mpi_unpack_aero_dist

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_aero_dist
