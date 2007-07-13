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

module mod_aero_dist

  type aero_mode_t
     real*8, pointer :: n_den(:)        ! len n_bin, number density (#/m^3)
     real*8, pointer :: vol_frac(:)     ! len n_spec, species fractions (1)
  end type aero_mode_t

  type aero_dist_t
     integer :: n_modes
     type(aero_mode_t), pointer :: modes(:) ! len n_mode, internally mixed modes
  end type aero_dist_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine alloc_aero_mode(n_bin, n_spec, aero_mode)

    ! Allocates an aero_mode.

    integer, intent(in) :: n_bin        ! number of bins
    integer, intent(in) :: n_spec       ! number of species
    type(aero_mode_t), intent(out) :: aero_mode ! aerosol mode

    allocate(aero_mode%n_den(n_bin))
    allocate(aero_mode%vol_frac(n_spec))

  end subroutine alloc_aero_mode

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine alloc_aero_dist(n_modes, n_bin, n_spec, aero_dist)

    ! Allocates an aero_dist.

    integer, intent(in) :: n_modes      ! number of modes
    integer, intent(in) :: n_bin        ! number of bins
    integer, intent(in) :: n_spec       ! number of species
    type(aero_dist_t), intent(out) :: aero_dist ! aerosol distribution

    integer :: i

    aero_dist%n_modes = n_modes
    allocate(aero_dist%modes(n_modes))
    do i = 1,n_modes
       call alloc_aero_mode(n_bin, n_spec, aero_dist%modes(i))
    end do

  end subroutine alloc_aero_dist

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine aero_dist_total_n_den(aero_dist, n_den)
      
    ! Compute the total number density of an aerosol distribution.
    
    type(aero_dist_t), intent(in) :: aero_dist ! aerosol distribution
    real*8, intent(out) :: n_den(:)     ! total number density (#/m^3)

    integer :: i

    n_den = 0d0
    do i = 1,aero_dist%n_modes
       n_den = n_den + aero_dist%modes(i)%n_den
    end do

  end subroutine aero_dist_total_n_den

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real*8 function aero_dist_total_num_den(bin_grid, aero_dist) ! #/m^3

    ! Returns the total number concentration in #/m^3 of a distribution.

    use mod_bin_grid

    type(bin_grid_t), intent(in) :: bin_grid ! bin grid
    type(aero_dist_t), intent(in) :: aero_dist ! aerosol distribution

    integer :: i
    
    aero_dist_total_num_den = 0d0
    do i = 1,aero_dist%n_modes
       aero_dist_total_num_den = aero_dist_total_num_den &
            + sum(aero_dist%modes(i)%n_den)
    end do
    aero_dist_total_num_den = aero_dist_total_num_den * bin_grid%dlnr

  end function aero_dist_total_num_den

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine init_log_normal(d_mean, log_sigma, bin_grid, n_den)

    ! Compute a log-normal distribution.
    
    use mod_bin_grid
    use mod_util
    use mod_constants
    
    real*8, intent(in) :: d_mean        ! geometric mean diameter of initial number dist (m)
    real*8, intent(in) :: log_sigma     ! log_10(geom. std dev(init dist)) (1)
    type(bin_grid_t), intent(in) :: bin_grid ! bin grid
    real*8,  intent(out) :: n_den(bin_grid%n_bin) ! init number den (#(ln(r))d(ln(r)))
                                         ! (normalized)
    
    integer k
    
    do k = 1,bin_grid%n_bin
       n_den(k) = 1d0 / (sqrt(2d0 * const%pi) * log_sigma) * &
            dexp(-(dlog10(vol2rad(bin_grid%v(k))) - dlog10(d_mean/2d0))**2d0 &
            / (2d0 * log_sigma**2d0)) / dlog(10d0)
    end do
    
    ! The formula above was originally for a distribution in
    ! log_10(r), while we are using log_e(r). The division by dlog(10)
    ! at the end corrects for this.
    
  end subroutine init_log_normal
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine init_exp(mean_vol, bin_grid, n_den)
    
    ! Exponential distribution in volume
    ! n(v) = 1 / mean_vol * exp(- v / mean_vol)
    
    use mod_bin_grid
    use mod_util
    
    real*8, intent(in) :: mean_vol      ! mean volume of init dist (m^3)
    type(bin_grid_t), intent(in) :: bin_grid ! bin grid
    real*8, intent(out) :: n_den(bin_grid%n_bin) ! init number density (#(ln(r))d(ln(r)))
    
    integer k
    real*8 n_den_vol
    
    do k = 1,bin_grid%n_bin
       n_den_vol = 1d0 / mean_vol * exp(-(bin_grid%v(k) / mean_vol))
       call vol_to_lnr(vol2rad(bin_grid%v(k)), n_den_vol, n_den(k))
    end do
    
  end subroutine init_exp
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine init_mono(vol, bin_grid, n_den)
    
    ! Mono-disperse distribution at mean_vol
    
    use mod_bin_grid
    use mod_util
    
    real*8, intent(in) :: vol           ! volume of each particle (m^3)
    type(bin_grid_t), intent(in) :: bin_grid ! bin grid
    real*8, intent(out) :: n_den(bin_grid%n_bin) ! init number density (#(ln(r))d(ln(r)))
    
    integer k

    n_den = 0d0
    call particle_in_bin(vol, bin_grid, k)
    n_den(k) = 1d0 / bin_grid%dlnr
    
  end subroutine init_mono
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine inout_write_aero_mode(file, aero_mode)
    
    ! Write full state.
    
    use mod_inout
    
    type(inout_file_t), intent(inout) :: file ! file to write to
    type(aero_mode_t), intent(in) :: aero_mode ! aero_mode to write

    call inout_write_real_array(file, "num_dens(num/m^3)", aero_mode%n_den)
    call inout_write_real_array(file, "volume_frac(1)", aero_mode%vol_frac)

  end subroutine inout_write_aero_mode

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine inout_write_aero_dist(file, aero_dist)
    
    ! Write full state.
    
    use mod_inout
    
    type(inout_file_t), intent(inout) :: file ! file to write to
    type(aero_dist_t), intent(in) :: aero_dist ! aero_dist to write

    integer :: i
    
    call inout_write_integer(file, "n_modes", aero_dist%n_modes)
    do i = 1,aero_dist%n_modes
       call inout_write_integer(file, "mode_number", i)
       call inout_write_aero_mode(file, aero_dist%modes(i))
    end do

  end subroutine inout_write_aero_dist

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine inout_read_aero_mode(file, aero_mode)
    
    ! Read full state.
    
    use mod_inout
    
    type(inout_file_t), intent(inout) :: file ! file to read from
    type(aero_mode_t), intent(out) :: aero_mode ! aero_mode to read

    call inout_read_real_array(file, "num_dens(num/m^3)", aero_mode%n_den)
    call inout_read_real_array(file, "volume_frac(1)", aero_mode%vol_frac)

  end subroutine inout_read_aero_mode

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine inout_read_aero_dist(file, aero_dist)
    
    ! Read full state.
    
    use mod_inout
    
    type(inout_file_t), intent(inout) :: file ! file to read from
    type(aero_dist_t), intent(out) :: aero_dist ! aero_dist to read

    integer :: i, check_i
    
    call inout_read_integer(file, "n_modes", aero_dist%n_modes)
    allocate(aero_dist%modes(aero_dist%n_modes))
    do i = 1,aero_dist%n_modes
       call inout_read_integer(file, "mode_number", check_i)
       call inout_check_index(file, i, check_i)
       call inout_read_aero_mode(file, aero_dist%modes(i))
    end do

  end subroutine inout_read_aero_dist

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine spec_read_vol_frac(file, aero_data, vol_frac)

    ! Read volume fractions from a data file.

    use mod_inout
    use mod_aero_data

    type(inout_file_t), intent(inout) :: file ! inout file
    type(aero_data_t), intent(in) :: aero_data   ! aero_data data
    real*8, intent(out) :: vol_frac(:)  ! aerosol species volume fractions

    integer :: n_species, species, i
    character(len=MAX_CHAR_LEN) :: read_name
    type(inout_file_t) :: read_file
    character(len=MAX_CHAR_LEN), pointer :: species_name(:)
    real*8, pointer :: species_data(:,:)
    integer :: species_data_shape(2)
    real*8 :: tot_vol_frac

    ! read the aerosol data from the specified file
    call inout_read_string(file, 'vol_frac', read_name)
    call inout_open_read(read_name, read_file)
    call inout_read_real_named_array(read_file, 0, species_name, species_data)
    call inout_close(read_file)

    ! check the data size
    species_data_shape = shape(species_data)
    n_species = species_data_shape(1)
    if (n_species < 1) then
       write(0,*) 'ERROR: file ', trim(read_name), &
            ' must contain at least one line of data'
       call exit(1)
    end if
    if (species_data_shape(2) /= 1) then
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

  subroutine spec_read_aero_mode_shape(file, aero_data, bin_grid, n_den)

    ! Read the shape (number density) of one mode of an aerosol
    ! distribution.

    use mod_inout
    use mod_bin_grid
    use mod_aero_data

    type(inout_file_t), intent(inout) :: file ! inout file
    type(aero_data_t), intent(in) :: aero_data ! aero_data data
    type(bin_grid_t), intent(in) :: bin_grid ! bin grid
    real*8 :: n_den(bin_grid%n_bin)     ! mode density

    real*8 :: num_conc
    character(len=MAX_CHAR_LEN) :: mode_type
    real*8 :: mean_vol, std_dev, vol, small_vol, big_vol, big_num

    call inout_read_real(file, 'num_conc', num_conc)
    call inout_read_string(file, 'mode_type', mode_type)
    if (trim(mode_type) == 'log_normal') then
       call inout_read_real(file, 'dist_mean_diam', mean_vol)
       call inout_read_real(file, 'dist_std_dev', std_dev)
       call init_log_normal(mean_vol, std_dev, bin_grid, n_den)
    elseif (trim(mode_type) == 'exp') then
       call inout_read_real(file, 'mean_vol', mean_vol)
       call init_exp(mean_vol, bin_grid, n_den)
    elseif (trim(mode_type) == 'mono') then
       call inout_read_real(file, 'vol', vol)
       call init_mono(vol, bin_grid, n_den)
    else
       write(0,'(a,a,a,a,a,i3)') 'ERROR: Unknown distribution type ', &
            trim(mode_type), ' in file ', trim(file%name), &
            ' at line ', file%line_num
       call exit(1)
    end if

    n_den = n_den * num_conc

  end subroutine spec_read_aero_mode_shape

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine spec_read_aero_mode(file, aero_data, bin_grid, aero_mode)

    ! Read one mode of an aerosol distribution (number density and
    ! volume fractions).

    use mod_inout
    use mod_bin_grid
    use mod_aero_data

    type(inout_file_t), intent(inout) :: file ! inout file
    type(aero_data_t), intent(in) :: aero_data   ! aero_data data
    type(bin_grid_t), intent(in) :: bin_grid ! bin grid
    type(aero_mode_t), intent(inout) :: aero_mode ! aerosol mode,
                                                  ! will be allocated

    allocate(aero_mode%n_den(bin_grid%n_bin))
    allocate(aero_mode%vol_frac(aero_data%n_spec))
    call spec_read_vol_frac(file, aero_data, aero_mode%vol_frac)
    call spec_read_aero_mode_shape(file, aero_data, bin_grid, aero_mode%n_den)

  end subroutine spec_read_aero_mode

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine spec_read_aero_dist(file, aero_data, bin_grid, aero_dist)

    ! Read continuous aerosol distribution composed of several modes.

    use mod_inout
    use mod_bin_grid
    use mod_aero_data

    type(inout_file_t), intent(inout) :: file ! inout file
    type(aero_data_t), intent(in) :: aero_data   ! aero_data data
    type(bin_grid_t), intent(in) :: bin_grid ! bin grid
    type(aero_dist_t), intent(inout) :: aero_dist ! aerosol dist,
                                                  ! will be allocated

    integer :: i

    call inout_read_integer(file, 'n_modes', aero_dist%n_modes)
    allocate(aero_dist%modes(aero_dist%n_modes))
    do i = 1,aero_dist%n_modes
       call spec_read_aero_mode(file, aero_data, bin_grid, aero_dist%modes(i))
    end do

  end subroutine spec_read_aero_dist

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine spec_read_aero_dist_filename(file, aero_data, bin_grid, name, dist)

    ! Read aerosol distribution from filename on line in file.

    use mod_inout
    use mod_bin_grid
    use mod_aero_data

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

  subroutine average_aero_mode(aero_mode_vec, aero_mode_avg)
    
    ! Computes the average of an array of aero_mode.

    use mod_util

    type(aero_mode_t), intent(in) :: aero_mode_vec(:) ! array of aero_mode
    type(aero_mode_t), intent(out) :: aero_mode_avg   ! average of aero_mode_vec

    integer :: n_bin, n_spec, i_bin, i_spec, i, n

    n_bin = size(aero_mode_vec(1)%n_den)
    n_spec = size(aero_mode_vec(1)%vol_frac)
    call alloc_aero_mode(n_bin, n_spec, aero_mode_avg)
    n = size(aero_mode_vec)
    do i_bin = 1,n_bin
       call average_real((/(aero_mode_vec(i)%n_den(i_bin),i=1,n)/), &
            aero_mode_avg%n_den(i_bin))
    end do
    do i_spec = 1,n_spec
       call average_real((/(aero_mode_vec(i)%vol_frac(i_spec),i=1,n)/), &
            aero_mode_avg%vol_frac(i_spec))
    end do
    
  end subroutine average_aero_mode
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine average_aero_dist(aero_dist_vec, aero_dist_avg)
    
    ! Computes the average of an array of aero_dist.

    type(aero_dist_t), intent(in) :: aero_dist_vec(:) ! array of aero_dist
    type(aero_dist_t), intent(out) :: aero_dist_avg   ! average of aero_dist_vec

    integer :: n_modes, i_mode, i, n

    n_modes = aero_dist_vec(1)%n_modes
    call alloc_aero_dist(n_modes, 0, 0, aero_dist_avg)
    n = size(aero_dist_vec)
    do i_mode = 1,n_modes
       call average_aero_mode((/(aero_dist_vec(i)%modes(i_mode),i=1,n)/), &
            aero_dist_avg%modes(i_mode))
    end do
    
  end subroutine average_aero_dist
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module mod_aero_dist
