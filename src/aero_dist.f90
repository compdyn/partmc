! Copyright (C) 2005-2008 Nicole Riemer and Matthew West
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
  use pmc_inout
  use pmc_aero_data
  use pmc_mpi
#ifdef PMC_USE_MPI
  use mpi
#endif

  !> An aerosol size distribution mode.
  !!
  !! Each mode is assumed to be fully internally mixed so that every
  !! particle has the same composition. The \c num_den array then
  !! stores the number density distribution.
  type aero_mode_t
     !> Number density [length bin_grid%%n_bin] (#/m^3).
     real*8, pointer :: num_den(:)
     !> Species fractions [length \c aero_data%%n_spec] (1).
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
  subroutine aero_mode_alloc(aero_mode, n_bin, n_spec)

    !> Number of bins.
    integer, intent(in) :: n_bin
    !> Number of species.
    integer, intent(in) :: n_spec
    !> Aerosol mode.
    type(aero_mode_t), intent(out) :: aero_mode

    allocate(aero_mode%num_den(n_bin))
    allocate(aero_mode%vol_frac(n_spec))

  end subroutine aero_mode_alloc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Free all storage.
  subroutine aero_mode_free(aero_mode)

    !> Aerosol mode.
    type(aero_mode_t), intent(inout) :: aero_mode

    deallocate(aero_mode%num_den)
    deallocate(aero_mode%vol_frac)

  end subroutine aero_mode_free

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Copy an aero_mode.
  subroutine aero_mode_copy(aero_mode_from, aero_mode_to)

    !> Aerosol mode original.
    type(aero_mode_t), intent(in) :: aero_mode_from
    !> Aerosol mode copy.
    type(aero_mode_t), intent(inout) :: aero_mode_to

    aero_mode_to%num_den = aero_mode_from%num_den
    aero_mode_to%vol_frac = aero_mode_from%vol_frac

  end subroutine aero_mode_copy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> aero_mode += aero_mode_delta
  subroutine aero_mode_add(aero_mode, aero_mode_delta)

    !> Aerosol mode.
    type(aero_mode_t), intent(inout) :: aero_mode
    !> Increment.
    type(aero_mode_t), intent(in) :: aero_mode_delta

    aero_mode%num_den = aero_mode%num_den + aero_mode_delta%num_den
    aero_mode%vol_frac = aero_mode%vol_frac + aero_mode_delta%vol_frac

  end subroutine aero_mode_add

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Scale an aero_mode.
  subroutine aero_mode_scale(aero_mode, alpha)

    !> Aerosol mode.
    type(aero_mode_t), intent(inout) :: aero_mode
    !> Scale factor.
    real*8, intent(in) :: alpha

    aero_mode%num_den = aero_mode%num_den * alpha
    aero_mode%vol_frac = aero_mode%vol_frac * alpha

  end subroutine aero_mode_scale

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Allocates an aero_dist.
  subroutine aero_dist_alloc(aero_dist, n_mode, n_bin, n_spec)

    !> Number of modes.
    integer, intent(in) :: n_mode
    !> Number of bins.
    integer, intent(in) :: n_bin
    !> Number of species.
    integer, intent(in) :: n_spec
    !> Aerosol distribution.
    type(aero_dist_t), intent(out) :: aero_dist

    integer :: i

    aero_dist%n_mode = n_mode
    allocate(aero_dist%mode(n_mode))
    do i = 1,n_mode
       call aero_mode_alloc(aero_dist%mode(i), n_bin, n_spec)
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

    integer :: n_bin, n_spec, i

    if (aero_dist_from%n_mode > 0) then
       n_bin = size(aero_dist_from%mode(1)%num_den)
       n_spec = size(aero_dist_from%mode(1)%vol_frac)
    else
       n_bin = 0
       n_spec = 0
    end if
    call aero_dist_free(aero_dist_to)
    call aero_dist_alloc(aero_dist_to, aero_dist_from%n_mode, n_bin, n_spec)
    do i = 1,aero_dist_from%n_mode
       call aero_mode_copy(aero_dist_from%mode(i), aero_dist_to%mode(i))
    end do

  end subroutine aero_dist_copy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> aero_dist += aero_dist_delta
  subroutine aero_dist_add(aero_dist, aero_dist_delta)

    !> Aero_dist.
    type(aero_dist_t), intent(inout) :: aero_dist
    !> Increment.
    type(aero_dist_t), intent(in) :: aero_dist_delta

    integer :: i

    do i = 1,aero_dist%n_mode
       call aero_mode_add(aero_dist%mode(i), aero_dist_delta%mode(i))
    end do

  end subroutine aero_dist_add

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> aero_dist *= alpha
  subroutine aero_dist_scale(aero_dist, alpha)

    !> Aero_dist.
    type(aero_dist_t), intent(inout) :: aero_dist
    !> Scale factor.
    real*8, intent(in) :: alpha

    integer :: i

    do i = 1,aero_dist%n_mode
       call aero_mode_scale(aero_dist%mode(i), alpha)
    end do

  end subroutine aero_dist_scale

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the total number concentration in #/m^3 of a distribution.
  !> (#/m^3)
  real*8 function aero_dist_total_num_den(bin_grid, aero_dist)

    !> Bin grid.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Aerosol distribution.
    type(aero_dist_t), intent(in) :: aero_dist

    integer :: i
    
    aero_dist_total_num_den = 0d0
    do i = 1,aero_dist%n_mode
       aero_dist_total_num_den = aero_dist_total_num_den &
            + sum(aero_dist%mode(i)%num_den)
    end do
    aero_dist_total_num_den = aero_dist_total_num_den * bin_grid%dlnr

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
!DEBUG
    write(*,*) 'radius = ', radius
    write(*,*) 'volume = ', rad2vol(radius)
    write(*,*) 'k = ', k
    write(*,*) 'bin_grid%v(k) = ', bin_grid%v(k)
    write(*,'(e50.40)') bin_grid%v(k)
    write(*,*) 'bin_edge(k - 2) = ', bin_edge(bin_grid, k - 2)
    write(*,*) 'bin_edge(k - 1) = ', bin_edge(bin_grid, k - 1)
    write(*,*) 'bin_edge(k) = ', bin_edge(bin_grid, k)
    write(*,*) 'bin_edge(k + 1) = ', bin_edge(bin_grid, k + 1)
    write(*,*) 'bin_edge(k + 2) = ', bin_edge(bin_grid, k + 2)
    write(*,*) 'r bin_edge(k - 2) = ', vol2rad(bin_edge(bin_grid, k - 2))
    write(*,*) 'r bin_edge(k - 1) = ', vol2rad(bin_edge(bin_grid, k - 1))
    write(*,*) 'r bin_edge(k) = ', vol2rad(bin_edge(bin_grid, k))
    write(*,*) 'r bin_edge(k + 1) = ', vol2rad(bin_edge(bin_grid, k + 1))
    write(*,*) 'r bin_edge(k + 2) = ', vol2rad(bin_edge(bin_grid, k + 2))
!DEBUG
    num_den(k) = 1d0 / bin_grid%dlnr
    
  end subroutine num_den_mono
  
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

  !> Write full state.
  subroutine inout_write_aero_mode(file, aero_mode)
    
    !> File to write to.
    type(inout_file_t), intent(inout) :: file
    !> Aero_mode to write.
    type(aero_mode_t), intent(in) :: aero_mode

    call inout_write_comment(file, "begin aero_mode")
    call inout_write_real_array(file, "num_dens(num/m^3)", aero_mode%num_den)
    call inout_write_real_array(file, "volume_frac(1)", aero_mode%vol_frac)
    call inout_write_comment(file, "end aero_mode")

  end subroutine inout_write_aero_mode

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Write full state.
  subroutine inout_write_aero_dist(file, aero_dist)
    
    !> File to write to.
    type(inout_file_t), intent(inout) :: file
    !> Aero_dist to write.
    type(aero_dist_t), intent(in) :: aero_dist

    integer :: i
    
    call inout_write_comment(file, "begin aero_dist")
    call inout_write_integer(file, "n_modes", aero_dist%n_mode)
    do i = 1,aero_dist%n_mode
       call inout_write_integer(file, "mode_number", i)
       call inout_write_aero_mode(file, aero_dist%mode(i))
    end do
    call inout_write_comment(file, "end aero_dist")

  end subroutine inout_write_aero_dist

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read full state.
  subroutine inout_read_aero_mode(file, aero_mode)
    
    !> File to read from.
    type(inout_file_t), intent(inout) :: file
    !> Aero_mode to read.
    type(aero_mode_t), intent(out) :: aero_mode

    call inout_check_comment(file, "begin aero_mode")
    call inout_read_real_array(file, "num_dens(num/m^3)", aero_mode%num_den)
    call inout_read_real_array(file, "volume_frac(1)", aero_mode%vol_frac)
    call inout_check_comment(file, "end aero_mode")

  end subroutine inout_read_aero_mode

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read full state.
  subroutine inout_read_aero_dist(file, aero_dist)
    
    !> File to read from.
    type(inout_file_t), intent(inout) :: file
    !> Aero_dist to read.
    type(aero_dist_t), intent(out) :: aero_dist

    integer :: i, check_i
    
    call inout_check_comment(file, "begin aero_dist")
    call inout_read_integer(file, "n_modes", aero_dist%n_mode)
    allocate(aero_dist%mode(aero_dist%n_mode))
    do i = 1,aero_dist%n_mode
       call inout_read_integer(file, "mode_number", check_i)
       call inout_check_index(file, i, check_i)
       call inout_read_aero_mode(file, aero_dist%mode(i))
    end do
    call inout_check_comment(file, "end aero_dist")

  end subroutine inout_read_aero_dist

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read volume fractions from a data file.
  subroutine spec_read_vol_frac(file, aero_data, vol_frac)

    !> Inout file.
    type(inout_file_t), intent(inout) :: file
    !> Aero_data data.
    type(aero_data_t), intent(in) :: aero_data
    !> Aerosol species volume fractions.
    real*8, intent(out) :: vol_frac(:)

    integer :: n_species, species, i
    character(len=MAX_VAR_LEN) :: read_name
    type(inout_file_t) :: read_file
    character(len=MAX_VAR_LEN), pointer :: species_name(:)
    real*8, pointer :: species_data(:,:)
    real*8 :: tot_vol_frac

    ! read the aerosol data from the specified file
    call inout_read_string(file, 'vol_frac', read_name)
    call inout_open_read(read_name, read_file)
    call inout_read_real_named_array(read_file, 0, species_name, species_data)
    call inout_close(read_file)

    ! check the data size
    n_species = size(species_data, 1)
    if (n_species < 1) then
       write(0,*) 'ERROR: file ', trim(read_name), &
            ' must contain at least one line of data'
       call exit(1)
    end if
    if (size(species_data, 2) /= 1) then
       write(0,*) 'ERROR: each line in ', trim(read_name), &
            ' should contain exactly one data value'
       call exit(1)
    end if

    ! copy over the data
    vol_frac = 0d0
    do i = 1,n_species
       species = aero_data_spec_by_name(aero_data, species_name(i))
       if (species == 0) then
          write(0,*) 'ERROR: unknown species ', trim(species_name(i)), &
               ' in file ', trim(read_name)
          call exit(1)
       end if
       vol_frac(species) = species_data(i,1)
    end do
    deallocate(species_name)
    deallocate(species_data)
    
    ! normalize
    tot_vol_frac = sum(vol_frac)
    if ((minval(vol_frac) < 0d0) .or. (tot_vol_frac <= 0d0)) then
       write(0,*) 'ERROR: vol_frac in ', trim(read_name), &
            ' is not positive'
       call exit(1)
    end if
    vol_frac = vol_frac / tot_vol_frac

  end subroutine spec_read_vol_frac

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read the shape (number density) of one mode of an aerosol
  !> distribution.
  subroutine spec_read_aero_mode_shape(file, aero_data, bin_grid, num_den)

    !> Inout file.
    type(inout_file_t), intent(inout) :: file
    !> Aero_data data.
    type(aero_data_t), intent(in) :: aero_data
    !> Bin grid.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Mode density.
    real*8 :: num_den(bin_grid%n_bin)

    real*8 :: tot_num_den
    character(len=MAX_VAR_LEN) :: mode_type
    real*8 :: mean_radius, log_std_dev, radius

    call inout_read_real(file, 'num_den', tot_num_den)
    call inout_read_string(file, 'mode_type', mode_type)
    if (trim(mode_type) == 'log_normal') then
       call inout_read_real(file, 'mean_radius', mean_radius)
       call inout_read_real(file, 'log_std_dev', log_std_dev)
       call num_den_log_normal(mean_radius, log_std_dev, bin_grid, num_den)
    elseif (trim(mode_type) == 'exp') then
       call inout_read_real(file, 'mean_radius', mean_radius)
       call num_den_exp(mean_radius, bin_grid, num_den)
    elseif (trim(mode_type) == 'mono') then
       call inout_read_real(file, 'radius', radius)
       call num_den_mono(radius, bin_grid, num_den)
    else
       write(0,'(a,a,a,a,a,i3)') 'ERROR: Unknown distribution type ', &
            trim(mode_type), ' in file ', trim(file%name), &
            ' at line ', file%line_num
       call exit(1)
    end if

    num_den = num_den * tot_num_den

  end subroutine spec_read_aero_mode_shape

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read one mode of an aerosol distribution (number density and
  !> volume fractions).
  subroutine spec_read_aero_mode(file, aero_data, bin_grid, aero_mode)

    !> Inout file.
    type(inout_file_t), intent(inout) :: file
    !> Aero_data data.
    type(aero_data_t), intent(in) :: aero_data
    !> Bin grid.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Aerosol mode (will be allocated).
    type(aero_mode_t), intent(inout) :: aero_mode

    allocate(aero_mode%num_den(bin_grid%n_bin))
    allocate(aero_mode%vol_frac(aero_data%n_spec))
    call spec_read_vol_frac(file, aero_data, aero_mode%vol_frac)
    call spec_read_aero_mode_shape(file, aero_data, bin_grid, &
         aero_mode%num_den)

  end subroutine spec_read_aero_mode

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read continuous aerosol distribution composed of several modes.
  subroutine spec_read_aero_dist(file, aero_data, bin_grid, aero_dist)

    !> Inout file.
    type(inout_file_t), intent(inout) :: file
    !> Aero_data data.
    type(aero_data_t), intent(in) :: aero_data
    !> Bin grid.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Aerosol dist,.
    type(aero_dist_t), intent(inout) :: aero_dist
                                                  ! will be allocated

    integer :: i

    call inout_read_integer(file, 'n_modes', aero_dist%n_mode)
    allocate(aero_dist%mode(aero_dist%n_mode))
    do i = 1,aero_dist%n_mode
       call spec_read_aero_mode(file, aero_data, bin_grid, aero_dist%mode(i))
    end do

  end subroutine spec_read_aero_dist

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read aerosol distribution from filename on line in file.
  subroutine spec_read_aero_dist_filename(file, aero_data, bin_grid, &
       name, dist)


    !> Inout file.
    type(inout_file_t), intent(inout) :: file
    !> Aero_data data.
    type(aero_data_t), intent(in) :: aero_data
    !> Bin grid.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Name of data line for filename.
    character(len=*), intent(in) :: name
    !> Aerosol distribution.
    type(aero_dist_t), intent(inout) :: dist

    character(len=MAX_VAR_LEN) :: read_name
    type(inout_file_t) :: read_file

    ! read the aerosol data from the specified file
    call inout_read_string(file, name, read_name)
    call inout_open_read(read_name, read_file)
    call spec_read_aero_dist(read_file, aero_data, bin_grid, dist)
    call inout_close(read_file)

  end subroutine spec_read_aero_dist_filename
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read an array of aero_dists with associated times and rates from
  !> the given file.
  subroutine spec_read_aero_dists_times_rates(file, aero_data, &
       bin_grid, name, times, rates, aero_dists)

    !> Inout file.
    type(inout_file_t), intent(inout) :: file
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
    type(inout_file_t) :: read_file
    type(inout_line_t) :: aero_dist_line
    integer :: n_time, i_time
    character(len=MAX_VAR_LEN), pointer :: names(:)
    real*8, pointer :: data(:,:)

    ! read the filename then read the data from that file
    call inout_read_string(file, name, read_name)
    call inout_open_read(read_name, read_file)
    call inout_read_real_named_array(read_file, 2, names, data)
    call inout_read_line_no_eof(read_file, aero_dist_line)
    call inout_check_line_name(read_file, aero_dist_line, "dist")
    call inout_check_line_length(read_file, aero_dist_line, size(data, 2))
    call inout_close(read_file)

    ! check the data size
    if (trim(names(1)) /= 'time') then
       write(0,*) 'ERROR: row 1 in ', trim(read_name), &
            ' must start with: time'
       call exit(1)
    end if
    if (trim(names(2)) /= 'rate') then
       write(0,*) 'ERROR: row 2 in ', trim(read_name), &
            ' must start with: rate'
       call exit(1)
    end if
    n_time = size(data, 2)
    if (n_time < 1) then
       write(0,*) 'ERROR: each line in ', trim(read_name), &
            ' must contain at least one data value'
       call exit(1)
    end if

    ! copy over the data
    allocate(aero_dists(n_time))
    allocate(times(n_time))
    allocate(rates(n_time))
    do i_time = 1,n_time
       call inout_open_read(aero_dist_line%data(i_time), read_file)
       call spec_read_aero_dist(read_file, aero_data, bin_grid, &
            aero_dists(i_time))
       call inout_close(read_file)
       times(i_time) = data(1,i_time)
       rates(i_time) = data(2,i_time)
    end do
    deallocate(names)
    deallocate(data)
    call inout_line_free(aero_dist_line)

  end subroutine spec_read_aero_dists_times_rates

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Computes the average of an array of aero_mode.
  subroutine aero_mode_average(aero_mode_vec, aero_mode_avg)

    !> Array of aero_mode.
    type(aero_mode_t), intent(in) :: aero_mode_vec(:)
    !> Avg of aero_mode_vec.
    type(aero_mode_t), intent(out) :: aero_mode_avg

    integer :: n_bin, n_spec, i_bin, i_spec, i, n

    n_bin = size(aero_mode_vec(1)%num_den)
    n_spec = size(aero_mode_vec(1)%vol_frac)
    call aero_mode_alloc(aero_mode_avg, n_bin, n_spec)
    n = size(aero_mode_vec)
    do i_bin = 1,n_bin
       call average_real((/(aero_mode_vec(i)%num_den(i_bin),i=1,n)/), &
            aero_mode_avg%num_den(i_bin))
    end do
    do i_spec = 1,n_spec
       call average_real((/(aero_mode_vec(i)%vol_frac(i_spec),i=1,n)/), &
            aero_mode_avg%vol_frac(i_spec))
    end do
    
  end subroutine aero_mode_average
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Computes the average of an array of aero_dist.
  subroutine aero_dist_average(aero_dist_vec, aero_dist_avg)

    !> Array of aero_dist.
    type(aero_dist_t), intent(in) :: aero_dist_vec(:)
    !> Avg of aero_dist_vec.
    type(aero_dist_t), intent(out) :: aero_dist_avg

    integer :: n_modes, i_mode, i, n

    n_modes = aero_dist_vec(1)%n_mode
    call aero_dist_alloc(aero_dist_avg, n_modes, 0, 0)
    n = size(aero_dist_vec)
    do i_mode = 1,n_modes
       call aero_mode_average((/(aero_dist_vec(i)%mode(i_mode),i=1,n)/), &
            aero_dist_avg%mode(i_mode))
    end do
    
  end subroutine aero_dist_average
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Determines the number of bytes required to pack the given value.
  integer function pmc_mpi_pack_size_aero_mode(val)

    !> Value to pack.
    type(aero_mode_t), intent(in) :: val

    pmc_mpi_pack_size_aero_mode = &
         pmc_mpi_pack_size_real_array(val%num_den) &
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
    call pmc_mpi_pack_real_array(buffer, position, val%num_den)
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
    call pmc_mpi_unpack_real_array(buffer, position, val%num_den)
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
