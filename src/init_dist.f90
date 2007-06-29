! Copyright (C) 2005-2007 Nicole Riemer and Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.
!
! Initial size distributions.
!
! The initial size distributions are computed as number densities, so
! they can be used for both sectional and particle-resolved
! simulations. The routine dist_to_n() converts a number density
! distribution to an actual number of particles ready for a
! particle-resolved simulation.

module mod_init_dist

  type aero_mode_t
     real*8, pointer :: n_den(:)        ! len n_bin, number density (#/m^3)
     real*8, pointer :: vol_frac(:)     ! len mat%n_spec, species fractions (1)
  end type aero_mode_t

  type aero_dist_t
     integer :: n_modes
     type(aero_mode_t), pointer :: modes(:) ! len n_mode, internally mixed modes
  end type aero_dist_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine dist_total_n_den(bin_grid, mat, dist, n_den)
      
    ! Compute the total number density of an aerosol distribution.
    
    use mod_material
    use mod_bin
    use mod_array
    use mod_util
    
    type(bin_grid_t), intent(in) :: bin_grid ! bin grid
    type(material), intent(in) :: mat   ! material data
    type(aero_dist_t), intent(in) :: dist ! aerosol distribution
    real*8, intent(out) :: n_den(bin_grid%n_bin) ! total number density (#/m^3)

    integer :: i

    n_den = 0d0
    do i = 1,dist%n_modes
       n_den = n_den + dist%modes(i)%n_den
    end do

  end subroutine dist_total_n_den

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine dist_to_part(bin_grid, mat, dist, n_part, aero)

    ! Convert a continuous distribution into particles.
    
    use mod_material
    use mod_bin
    use mod_array
    use mod_util

    type(bin_grid_t), intent(in) :: bin_grid ! bin grid
    type(material), intent(in) :: mat   ! material data
    type(aero_dist_t), intent(in) :: dist ! aerosol distribution
    integer, intent(in) :: n_part       ! total number of particles
    type(aerosol), intent(out) :: aero  ! aerosol distribution, will be alloced

    integer :: i
    real*8 :: total_n_den
    real*8 :: mode_n_dens(dist%n_modes)
    integer :: mode_n_parts(dist%n_modes)
    integer :: num_per_bin(bin_grid%n_bin)

    ! find the total number density of each mode
    total_n_den = 0d0
    do i = 1,dist%n_modes
       mode_n_dens(i) = sum(dist%modes(i)%n_den)
    end do
    total_n_den = sum(mode_n_dens)

    ! allocate particles to modes proportional to their number densities
    call vec_cts_to_disc(dist%n_modes, mode_n_dens, n_part, mode_n_parts)

    ! allocate particles within each mode in proportion to mode shape
    call allocate_aerosol(bin_grid%n_bin, mat%n_spec, aero)
    do i = 1,dist%n_modes
       call vec_cts_to_disc(bin_grid%n_bin, dist%modes(i)%n_den, &
            mode_n_parts(i), num_per_bin)
       call add_particles(bin_grid, mat, dist%modes(i)%vol_frac, &
            num_per_bin, aero)
    end do

  end subroutine dist_to_part
  
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
  
  subroutine init_dist(dist_type, dist_args, n_bin, bin_v, n_den)

    ! Multiplexer to make an initial distribution based on its
    ! name. These should be normalized so that sum(n_den) = 1/dlnr.

    character(len=*), intent(in) :: dist_type ! type of distribution
    real*8, intent(in) :: dist_args(:)  ! distribution parameters
    integer, intent(in) :: n_bin        ! number of bins
    real*8, intent(in) :: bin_v(n_bin)  ! volume of particles in bins (m^3)
    real*8, intent(out) :: n_den(n_bin) ! initial number density 
                                        ! (#(ln(r))d(ln(r))) (normalized)

    if (trim(dist_type) == 'log_normal') then
       call init_log_normal(dist_args(1), dist_args(2), n_bin, &
            bin_v, n_den)
    elseif (trim(dist_type) == 'exp') then
       call init_exp(dist_args(1), n_bin, bin_v, n_den)
    elseif (trim(dist_type) == 'mono') then
       call init_mono(dist_args(1), n_bin, bin_v, n_den)
    else
       write(0,*) 'ERROR: unknown distribution type: ', trim(dist_type)
       call exit(1)
    end if

  end subroutine init_dist
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine init_bidisperse(M, small_vol, big_vol, big_num, &
       n_bin, bin_v, aero)
    
    ! This is not like the other init functions. It does not produce a
    ! number density. Instead it generates MH and VH directly.

    use mod_array
    use mod_bin

    integer, intent(in) :: M            ! total number of particles
    real*8, intent(in) :: small_vol     ! volume of small particles (m^3)
    real*8, intent(in) :: big_vol       ! volume of big particle (m^3)
    real*8, intent(in) :: big_num       ! number of big particle
    integer, intent(in) :: n_bin        ! number of bins
    real*8, intent(in) :: bin_v(n_bin)  ! volume of particles in bins (m^3)
    type(aerosol), intent(inout) :: aero ! aerosol
    
    integer i_small, i_big, n_big

    n_big = nint(big_num)
    call particle_in_bin(small_vol, n_bin, bin_v, i_small)
    call particle_in_bin(big_vol, n_bin, bin_v, i_big)
    aero%n(i_small) = M - n_big
    aero%n(i_big) = n_big
    call enlarge_bin_to(aero%v(i_small), aero%n(i_small))
    call enlarge_bin_to(aero%v(i_big), aero%n(i_big))
    aero%v(i_small)%p = small_vol
    aero%v(i_big)%p = big_vol

  end subroutine init_bidisperse
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine init_log_normal(d_mean, log_sigma, n_bin, &
       bin_v, n_den)

    ! Compute a log-normal distribution.
    
    use mod_util
    use mod_constants
    
    real*8, intent(in) :: d_mean        ! geometric mean diameter of initial number dist (m)
    real*8, intent(in) :: log_sigma     ! log_10(geom. std dev(init dist)) (1)
    integer, intent(in) :: n_bin        ! number of bins
    real*8, intent(in) :: bin_v(n_bin)  ! volume of particles in bins (m^3)
    real*8,  intent(out) :: n_den(n_bin) ! init number den (#(ln(r))d(ln(r)))
                                         ! (normalized)
    
    integer k
    
    do k = 1,n_bin
       n_den(k) = 1d0 / (sqrt(2d0 * const%pi) * log_sigma) * &
            dexp(-(dlog10(vol2rad(bin_v(k))) - dlog10(d_mean/2d0))**2d0 &
            / (2d0 * log_sigma**2d0)) / dlog(10d0)
    end do
    
    ! The formula above was originally for a distribution in
    ! log_10(r), while we are using log_e(r). The division by dlog(10)
    ! at the end corrects for this.
    
  end subroutine init_log_normal
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine init_exp(mean_vol, n_bin, bin_v, n_den)
    
    ! Exponential distribution in volume
    ! n(v) = 1 / mean_vol * exp(- v / mean_vol)
    
    use mod_bin
    use mod_util
    
    real*8, intent(in) :: mean_vol      ! mean volume of init dist (m^3)
    integer, intent(in) :: n_bin        ! number of bins
    real*8, intent(in) :: bin_v(n_bin)  ! volume of particles in bins (m^3)
    real*8, intent(out) :: n_den(n_bin) ! init number density (#(ln(r))d(ln(r)))
    
    integer k
    real*8 n_den_vol
    
    do k = 1,n_bin
       n_den_vol = 1d0 / mean_vol * exp(-(bin_v(k) / mean_vol))
       call vol_to_lnr(vol2rad(bin_v(k)), n_den_vol, n_den(k))
    end do
    
  end subroutine init_exp
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine init_mono(vol, n_bin, bin_v, n_den)
    
    ! Mono-disperse distribution at mean_vol
    
    use mod_bin
    use mod_util
    
    real*8, intent(in) :: vol           ! volume of each particle (m^3)
    integer, intent(in) :: n_bin        ! number of bins
    real*8, intent(in) :: bin_v(n_bin)  ! volume of particles in bins (m^3)
    real*8, intent(out) :: n_den(n_bin) ! init number density (#(ln(r))d(ln(r)))
    
    integer k

    n_den = 0d0
    call particle_in_bin(vol, n_bin, bin_v, k)
    n_den(k) = 1d0
    ! FIXME: should really be 1/dlnr rather than 1
    
  end subroutine init_mono
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
end module mod_init_dist
