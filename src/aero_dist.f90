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
     real*8, pointer :: vol_frac(:)     ! len aero_data%n_spec, species fractions (1)
  end type aero_mode_t

  type aero_dist_t
     integer :: n_modes
     type(aero_mode_t), pointer :: modes(:) ! len n_mode, internally mixed modes
  end type aero_dist_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine dist_total_n_den(bin_grid, aero_data, dist, n_den)
      
    ! Compute the total number density of an aerosol distribution.
    
    use mod_aero_data
    use mod_bin
    use mod_util
    
    type(bin_grid_t), intent(in) :: bin_grid ! bin grid
    type(aero_data_t), intent(in) :: aero_data   ! aero_data data
    type(aero_dist_t), intent(in) :: dist ! aerosol distribution
    real*8, intent(out) :: n_den(bin_grid%n_bin) ! total number density (#/m^3)

    integer :: i

    n_den = 0d0
    do i = 1,dist%n_modes
       n_den = n_den + dist%modes(i)%n_den
    end do

  end subroutine dist_total_n_den

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real*8 function dist_num_conc(bin_grid, dist) ! #/m^3

    ! Returns the total number concentration in #/m^3 of a distribution.

    use mod_bin

    type(bin_grid_t), intent(in) :: bin_grid ! bin grid
    type(aero_dist_t), intent(in) :: dist ! aerosol distribution

    integer :: i
    
    dist_num_conc = 0d0
    do i = 1,dist%n_modes
       dist_num_conc = dist_num_conc + sum(dist%modes(i)%n_den)
    end do
    dist_num_conc = dist_num_conc * bin_grid%dlnr

  end function dist_num_conc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine init_log_normal(d_mean, log_sigma, bin_grid, n_den)

    ! Compute a log-normal distribution.
    
    use mod_bin
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
    
    use mod_bin
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
    
    use mod_bin
    use mod_util
    
    real*8, intent(in) :: vol           ! volume of each particle (m^3)
    type(bin_grid_t), intent(in) :: bin_grid ! bin grid
    real*8, intent(out) :: n_den(bin_grid%n_bin) ! init number density (#(ln(r))d(ln(r)))
    
    integer k

    n_den = 0d0
    call particle_in_bin(vol, bin_grid%n_bin, bin_grid%v, k)
    n_den(k) = 1d0 / bin_grid%dlnr
    
  end subroutine init_mono
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine inout_write_aero_mode(file, aero_mode)
    
    ! Write full state.
    
    use mod_inout
    
    type(inout_file_t), intent(inout) :: file ! file to write to
    type(aero_mode_t), intent(in) :: aero_mode ! aero_mode to write

    call inout_write_real_array(file, "num_density(#/m^3)", aero_mode%n_den)
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

    call inout_read_real_array(file, "num_density(#/m^3)", aero_mode%n_den)
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

end module mod_aero_dist
