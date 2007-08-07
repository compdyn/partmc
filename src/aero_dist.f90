! Copyright (C) 2005-2007 Nicole Riemer and Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.
!
! Aerosol size distributions.
!
! The initial size distributions are computed as number densities, so
! they can be used for both sectional and particle-resolved
! simulations. The routine dist_to_n() converts a number density
! distribution to an actual number of particles ready for a
! particle-resolved simulation.
!
! Initial distributions should be normalized so that sum(n_den) = 1/dlnr.

module pmc_aero_dist

  type aero_mode_t
     real*8, pointer :: num_den(:)      ! len n_bin, number density (#/m^3)
     real*8, pointer :: vol_frac(:)     ! len n_spec, species fractions (1)
  end type aero_mode_t

  type aero_dist_t
     integer :: n_mode
     type(aero_mode_t), pointer :: mode(:) ! len n_mode, internally mixed modes
  end type aero_dist_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine aero_mode_alloc(aero_mode, n_bin, n_spec)

    ! Allocates an aero_mode.

    integer, intent(in) :: n_bin        ! number of bins
    integer, intent(in) :: n_spec       ! number of species
    type(aero_mode_t), intent(out) :: aero_mode ! aerosol mode

    allocate(aero_mode%num_den(n_bin))
    allocate(aero_mode%vol_frac(n_spec))

  end subroutine aero_mode_alloc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine aero_mode_free(aero_mode)

    ! Free all storage.

    type(aero_mode_t), intent(inout) :: aero_mode ! aerosol mode

    deallocate(aero_mode%num_den)
    deallocate(aero_mode%vol_frac)

  end subroutine aero_mode_free

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine aero_mode_copy(aero_mode_from, aero_mode_to)

    ! Copy an aero_mode.

    type(aero_mode_t), intent(in) :: aero_mode_from ! aerosol mode original
    type(aero_mode_t), intent(inout) :: aero_mode_to ! aerosol mode copy

    aero_mode_to%num_den = aero_mode_from%num_den
    aero_mode_to%vol_frac = aero_mode_from%vol_frac

  end subroutine aero_mode_copy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine aero_mode_scale(aero_mode, alpha)

    ! Scale an aero_mode.

    type(aero_mode_t), intent(inout) :: aero_mode ! aerosol mode
    real*8, intent(in) :: alpha         ! scale factor

    aero_mode%num_den = aero_mode%num_den * alpha

  end subroutine aero_mode_scale

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine aero_dist_alloc(aero_dist, n_mode, n_bin, n_spec)

    ! Allocates an aero_dist.

    integer, intent(in) :: n_mode       ! number of modes
    integer, intent(in) :: n_bin        ! number of bins
    integer, intent(in) :: n_spec       ! number of species
    type(aero_dist_t), intent(out) :: aero_dist ! aerosol distribution

    integer :: i

    aero_dist%n_mode = n_mode
    allocate(aero_dist%mode(n_mode))
    do i = 1,n_mode
       call aero_mode_alloc(aero_dist%mode(i), n_bin, n_spec)
    end do

  end subroutine aero_dist_alloc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine aero_dist_free(aero_dist)

    ! Free all storage.

    type(aero_dist_t), intent(inout) :: aero_dist ! aerosol distribution

    integer :: i

    do i = 1,aero_dist%n_mode
       call aero_mode_free(aero_dist%mode(i))
    end do
    deallocate(aero_dist%mode)

  end subroutine aero_dist_free

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine aero_dist_copy(aero_dist_from, aero_dist_to)

    ! Copy an aero_dist.

    type(aero_dist_t), intent(in) :: aero_dist_from ! aero_dist original
    type(aero_dist_t), intent(inout) :: aero_dist_to ! aero_dist copy

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

  real*8 function aero_dist_total_num_den(bin_grid, aero_dist) ! (#/m^3)

    ! Returns the total number concentration in #/m^3 of a distribution.

    use pmc_bin_grid

    type(bin_grid_t), intent(in) :: bin_grid ! bin grid
    type(aero_dist_t), intent(in) :: aero_dist ! aerosol distribution

    integer :: i
    
    aero_dist_total_num_den = 0d0
    do i = 1,aero_dist%n_mode
       aero_dist_total_num_den = aero_dist_total_num_den &
            + sum(aero_dist%mode(i)%num_den)
    end do
    aero_dist_total_num_den = aero_dist_total_num_den * bin_grid%dlnr

  end function aero_dist_total_num_den

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine num_den_log_normal(d_mean, log_sigma, bin_grid, num_den)

    ! Compute a log-normal distribution.
    
    use pmc_bin_grid
    use pmc_util
    use pmc_constants
    
    real*8, intent(in) :: d_mean        ! geometric mean diameter (m)
    real*8, intent(in) :: log_sigma     ! log_10(geom. std dev) (1)
    type(bin_grid_t), intent(in) :: bin_grid ! bin grid
    real*8, intent(out) :: num_den(bin_grid%n_bin) ! num den (#(ln(r))d(ln(r)))
                                        ! (normalized)
    
    integer k
    
    do k = 1,bin_grid%n_bin
       num_den(k) = 1d0 / (sqrt(2d0 * const%pi) * log_sigma) * &
            dexp(-(dlog10(vol2rad(bin_grid%v(k))) - dlog10(d_mean/2d0))**2d0 &
            / (2d0 * log_sigma**2d0)) / dlog(10d0)
    end do
    
    ! The formula above was originally for a distribution in
    ! log_10(r), while we are using log_e(r). The division by dlog(10)
    ! at the end corrects for this.
    
  end subroutine num_den_log_normal
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine num_den_exp(mean_vol, bin_grid, num_den)
    
    ! Exponential distribution in volume
    ! n(v) = 1 / mean_vol * exp(- v / mean_vol)
    
    use pmc_bin_grid
    use pmc_util
    
    real*8, intent(in) :: mean_vol      ! mean volume (m^3)
    type(bin_grid_t), intent(in) :: bin_grid ! bin grid
    real*8, intent(out) :: num_den(bin_grid%n_bin) ! num den (#(ln(r))d(ln(r)))
    
    integer k
    real*8 num_den_vol
    
    do k = 1,bin_grid%n_bin
       num_den_vol = 1d0 / mean_vol * exp(-(bin_grid%v(k) / mean_vol))
       call vol_to_lnr(vol2rad(bin_grid%v(k)), num_den_vol, num_den(k))
    end do
    
  end subroutine num_den_exp
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine num_den_mono(radius, bin_grid, num_den)
    
    ! Mono-disperse distribution.
    
    use pmc_bin_grid
    use pmc_util
    
    real*8, intent(in) :: radius         ! radius of each particle (m^3)
    type(bin_grid_t), intent(in) :: bin_grid ! bin grid
    real*8, intent(out) :: num_den(bin_grid%n_bin) ! num den (#(ln(r))d(ln(r)))
    
    integer :: k

    num_den = 0d0
    k = bin_grid_particle_in_bin(bin_grid, rad2vol(radius))
    num_den(k) = 1d0 / bin_grid%dlnr
    
  end subroutine num_den_mono
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine aero_dist_interp_1d(aero_dist_list, time_list, &
         rate_list, time, aero_dist, rate)

    ! Determine the current aero_dist and rate by interpolating at the
    ! current time with the lists of aero_dists and rates.

    use pmc_util

    type(aero_dist_t), intent(in) :: aero_dist_list(:) ! gas states
    real*8, intent(in) :: time_list(size(aero_dist_list)) ! times (s)
    real*8, intent(in) :: rate_list(size(aero_dist_list)) ! rates (s^{-1})
    real*8, intent(in) :: time          ! current time (s)
    type(aero_dist_t), intent(inout) :: aero_dist ! current gas state
    real*8, intent(out) :: rate         ! current rate (s^{-1})

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

  subroutine inout_write_aero_mode(file, aero_mode)
    
    ! Write full state.
    
    use pmc_inout
    
    type(inout_file_t), intent(inout) :: file ! file to write to
    type(aero_mode_t), intent(in) :: aero_mode ! aero_mode to write

    call inout_write_real_array(file, "num_dens(num/m^3)", aero_mode%num_den)
    call inout_write_real_array(file, "volume_frac(1)", aero_mode%vol_frac)

  end subroutine inout_write_aero_mode

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine inout_write_aero_dist(file, aero_dist)
    
    ! Write full state.
    
    use pmc_inout
    
    type(inout_file_t), intent(inout) :: file ! file to write to
    type(aero_dist_t), intent(in) :: aero_dist ! aero_dist to write

    integer :: i
    
    call inout_write_integer(file, "n_modes", aero_dist%n_mode)
    do i = 1,aero_dist%n_mode
       call inout_write_integer(file, "mode_number", i)
       call inout_write_aero_mode(file, aero_dist%mode(i))
    end do

  end subroutine inout_write_aero_dist

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine inout_read_aero_mode(file, aero_mode)
    
    ! Read full state.
    
    use pmc_inout
    
    type(inout_file_t), intent(inout) :: file ! file to read from
    type(aero_mode_t), intent(out) :: aero_mode ! aero_mode to read

    call inout_read_real_array(file, "num_dens(num/m^3)", aero_mode%num_den)
    call inout_read_real_array(file, "volume_frac(1)", aero_mode%vol_frac)

  end subroutine inout_read_aero_mode

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine inout_read_aero_dist(file, aero_dist)
    
    ! Read full state.
    
    use pmc_inout
    
    type(inout_file_t), intent(inout) :: file ! file to read from
    type(aero_dist_t), intent(out) :: aero_dist ! aero_dist to read

    integer :: i, check_i
    
    call inout_read_integer(file, "n_modes", aero_dist%n_mode)
    allocate(aero_dist%mode(aero_dist%n_mode))
    do i = 1,aero_dist%n_mode
       call inout_read_integer(file, "mode_number", check_i)
       call inout_check_index(file, i, check_i)
       call inout_read_aero_mode(file, aero_dist%mode(i))
    end do

  end subroutine inout_read_aero_dist

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine spec_read_vol_frac(file, aero_data, vol_frac)

    ! Read volume fractions from a data file.

    use pmc_inout
    use pmc_aero_data

    type(inout_file_t), intent(inout) :: file ! inout file
    type(aero_data_t), intent(in) :: aero_data   ! aero_data data
    real*8, intent(out) :: vol_frac(:)  ! aerosol species volume fractions

    integer :: n_species, species, i
    character(len=MAX_CHAR_LEN) :: read_name
    type(inout_file_t) :: read_file
    character(len=MAX_CHAR_LEN), pointer :: species_name(:)
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

  subroutine spec_read_aero_mode_shape(file, aero_data, bin_grid, num_den)

    ! Read the shape (number density) of one mode of an aerosol
    ! distribution.

    use pmc_inout
    use pmc_bin_grid
    use pmc_aero_data

    type(inout_file_t), intent(inout) :: file ! inout file
    type(aero_data_t), intent(in) :: aero_data ! aero_data data
    type(bin_grid_t), intent(in) :: bin_grid ! bin grid
    real*8 :: num_den(bin_grid%n_bin)   ! mode density

    real*8 :: num_conc
    character(len=MAX_CHAR_LEN) :: mode_type
    real*8 :: mean_vol, std_dev, radius, small_vol, big_vol, big_num

    call inout_read_real(file, 'num_conc', num_conc)
    call inout_read_string(file, 'mode_type', mode_type)
    if (trim(mode_type) == 'log_normal') then
       call inout_read_real(file, 'dist_mean_diam', mean_vol)
       call inout_read_real(file, 'dist_std_dev', std_dev)
       call num_den_log_normal(mean_vol, std_dev, bin_grid, num_den)
    elseif (trim(mode_type) == 'exp') then
       call inout_read_real(file, 'mean_vol', mean_vol)
       call num_den_exp(mean_vol, bin_grid, num_den)
    elseif (trim(mode_type) == 'mono') then
       call inout_read_real(file, 'radius', radius)
       call num_den_mono(radius, bin_grid, num_den)
    else
       write(0,'(a,a,a,a,a,i3)') 'ERROR: Unknown distribution type ', &
            trim(mode_type), ' in file ', trim(file%name), &
            ' at line ', file%line_num
       call exit(1)
    end if

    num_den = num_den * num_conc

  end subroutine spec_read_aero_mode_shape

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine spec_read_aero_mode(file, aero_data, bin_grid, aero_mode)

    ! Read one mode of an aerosol distribution (number density and
    ! volume fractions).

    use pmc_inout
    use pmc_bin_grid
    use pmc_aero_data

    type(inout_file_t), intent(inout) :: file ! inout file
    type(aero_data_t), intent(in) :: aero_data   ! aero_data data
    type(bin_grid_t), intent(in) :: bin_grid ! bin grid
    type(aero_mode_t), intent(inout) :: aero_mode ! aerosol mode,
                                                  ! will be allocated

    allocate(aero_mode%num_den(bin_grid%n_bin))
    allocate(aero_mode%vol_frac(aero_data%n_spec))
    call spec_read_vol_frac(file, aero_data, aero_mode%vol_frac)
    call spec_read_aero_mode_shape(file, aero_data, bin_grid, aero_mode%num_den)

  end subroutine spec_read_aero_mode

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine spec_read_aero_dist(file, aero_data, bin_grid, aero_dist)

    ! Read continuous aerosol distribution composed of several modes.

    use pmc_inout
    use pmc_bin_grid
    use pmc_aero_data

    type(inout_file_t), intent(inout) :: file ! inout file
    type(aero_data_t), intent(in) :: aero_data   ! aero_data data
    type(bin_grid_t), intent(in) :: bin_grid ! bin grid
    type(aero_dist_t), intent(inout) :: aero_dist ! aerosol dist,
                                                  ! will be allocated

    integer :: i

    call inout_read_integer(file, 'n_modes', aero_dist%n_mode)
    allocate(aero_dist%mode(aero_dist%n_mode))
    do i = 1,aero_dist%n_mode
       call spec_read_aero_mode(file, aero_data, bin_grid, aero_dist%mode(i))
    end do

  end subroutine spec_read_aero_dist

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine spec_read_aero_dist_filename(file, aero_data, bin_grid, name, dist)

    ! Read aerosol distribution from filename on line in file.

    use pmc_inout
    use pmc_bin_grid
    use pmc_aero_data

    type(inout_file_t), intent(inout) :: file ! inout file
    type(aero_data_t), intent(in) :: aero_data   ! aero_data data
    type(bin_grid_t), intent(in) :: bin_grid ! bin grid
    character(len=*), intent(in) :: name ! name of data line for filename
    type(aero_dist_t), intent(inout) :: dist ! aerosol distribution

    character(len=MAX_CHAR_LEN) :: read_name
    type(inout_file_t) :: read_file

    ! read the aerosol data from the specified file
    call inout_read_string(file, name, read_name)
    call inout_open_read(read_name, read_file)
    call spec_read_aero_dist(read_file, aero_data, bin_grid, dist)
    call inout_close(read_file)

  end subroutine spec_read_aero_dist_filename
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine spec_read_aero_dists_times_rates(file, aero_data, &
       bin_grid, name, times, rates, aero_dists)

    ! Read an array of aero_dists with associated times and rates from
    ! the given file.

    use pmc_inout
    use pmc_aero_data
    use pmc_bin_grid

    type(inout_file_t), intent(inout) :: file ! inout file
    type(aero_data_t), intent(in) :: aero_data ! aero data
    type(bin_grid_t), intent(in) :: bin_grid ! bin grid
    character(len=*), intent(in) :: name ! name of data line for filename
    real*8, pointer :: times(:)         ! times (s)
    real*8, pointer :: rates(:)         ! rates (s^{-1})
    type(aero_dist_t), pointer :: aero_dists(:) ! aero dists

    character(len=MAX_CHAR_LEN) :: read_name
    type(inout_file_t) :: read_file
    type(inout_line_t) :: aero_dist_line
    integer :: n_time, i_time
    character(len=MAX_CHAR_LEN), pointer :: names(:)
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

  end subroutine spec_read_aero_dists_times_rates

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine aero_mode_average(aero_mode_vec, aero_mode_avg)
    
    ! Computes the average of an array of aero_mode.

    use pmc_util

    type(aero_mode_t), intent(in) :: aero_mode_vec(:) ! array of aero_mode
    type(aero_mode_t), intent(out) :: aero_mode_avg   ! average of aero_mode_vec

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

  subroutine aero_dist_average(aero_dist_vec, aero_dist_avg)
    
    ! Computes the average of an array of aero_dist.

    type(aero_dist_t), intent(in) :: aero_dist_vec(:) ! array of aero_dist
    type(aero_dist_t), intent(out) :: aero_dist_avg   ! average of aero_dist_vec

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

  integer function pmc_mpi_pack_aero_mode_size(val)

    ! Determines the number of bytes required to pack the given value.

    use pmc_mpi

    type(aero_mode_t), intent(in) :: val ! value to pack

    pmc_mpi_pack_aero_mode_size = &
         pmc_mpi_pack_real_array_size(val%num_den) &
         + pmc_mpi_pack_real_array_size(val%vol_frac)

  end function pmc_mpi_pack_aero_mode_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function pmc_mpi_pack_aero_dist_size(val)

    ! Determines the number of bytes required to pack the given value.

    use pmc_mpi

    type(aero_dist_t), intent(in) :: val ! value to pack

    integer :: i, total_size

    total_size = pmc_mpi_pack_integer_size(val%n_mode)
    do i = 1,size(val%mode)
       total_size = total_size + pmc_mpi_pack_aero_mode_size(val%mode(i))
    end do
    pmc_mpi_pack_aero_dist_size = total_size

  end function pmc_mpi_pack_aero_dist_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine pmc_mpi_pack_aero_mode(buffer, position, val)

    ! Packs the given value into the buffer, advancing position.

#ifdef PMC_USE_MPI
    use mpi
    use pmc_mpi
    use pmc_util
#endif

    character, intent(inout) :: buffer(:) ! memory buffer
    integer, intent(inout) :: position  ! current buffer position
    type(aero_mode_t), intent(in) :: val ! value to pack

#ifdef PMC_USE_MPI
    integer :: prev_position

    prev_position = position
    call pmc_mpi_pack_real_array(buffer, position, val%num_den)
    call pmc_mpi_pack_real_array(buffer, position, val%vol_frac)
    call assert(position - prev_position == pmc_mpi_pack_aero_mode_size(val))
#endif

  end subroutine pmc_mpi_pack_aero_mode

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine pmc_mpi_pack_aero_dist(buffer, position, val)

    ! Packs the given value into the buffer, advancing position.

#ifdef PMC_USE_MPI
    use mpi
    use pmc_mpi
    use pmc_util
#endif

    character, intent(inout) :: buffer(:) ! memory buffer
    integer, intent(inout) :: position  ! current buffer position
    type(aero_dist_t), intent(in) :: val ! value to pack

#ifdef PMC_USE_MPI
    integer :: prev_position, i

    prev_position = position
    call pmc_mpi_pack_integer(buffer, position, val%n_mode)
    do i = 1,size(val%mode)
       call pmc_mpi_pack_aero_mode(buffer, position, val%mode(i))
    end do
    call assert(position - prev_position == pmc_mpi_pack_aero_dist_size(val))
#endif

  end subroutine pmc_mpi_pack_aero_dist

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine pmc_mpi_unpack_aero_mode(buffer, position, val)

    ! Unpacks the given value from the buffer, advancing position.

#ifdef PMC_USE_MPI
    use mpi
    use pmc_mpi
    use pmc_util
#endif

    character, intent(inout) :: buffer(:) ! memory buffer
    integer, intent(inout) :: position  ! current buffer position
    type(aero_mode_t), intent(out) :: val ! value to pack

#ifdef PMC_USE_MPI
    integer :: prev_position

    prev_position = position
    call pmc_mpi_unpack_real_array(buffer, position, val%num_den)
    call pmc_mpi_unpack_real_array(buffer, position, val%vol_frac)
    call assert(position - prev_position == pmc_mpi_pack_aero_mode_size(val))
#endif

  end subroutine pmc_mpi_unpack_aero_mode

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine pmc_mpi_unpack_aero_dist(buffer, position, val)

    ! Unpacks the given value from the buffer, advancing position.

#ifdef PMC_USE_MPI
    use mpi
    use pmc_mpi
    use pmc_util
#endif

    character, intent(inout) :: buffer(:) ! memory buffer
    integer, intent(inout) :: position  ! current buffer position
    type(aero_dist_t), intent(out) :: val ! value to pack

#ifdef PMC_USE_MPI
    integer :: prev_position, i

    prev_position = position
    call pmc_mpi_unpack_integer(buffer, position, val%n_mode)
    allocate(val%mode(val%n_mode))
    do i = 1,size(val%mode)
       call pmc_mpi_unpack_aero_mode(buffer, position, val%mode(i))
    end do
    call assert(position - prev_position == pmc_mpi_pack_aero_dist_size(val))
#endif

  end subroutine pmc_mpi_unpack_aero_dist

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_aero_dist
