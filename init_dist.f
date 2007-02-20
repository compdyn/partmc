! -*- mode: f90; -*-
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
contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine init_dist(dist_type, dist_args, n_bin, bin_v, n_den)

    ! Multiplexer to make an initial distribution based on its name.

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
    else
       write(0,*) 'ERROR: unknown distribution type: ', trim(dist_type)
       call exit(1)
    end if

  end subroutine init_dist
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine init_bidisperse(M, small_vol, big_vol, big_num, &
       n_bin, bin_v, MH, VH)
    
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
    integer, intent(inout) :: MH(n_bin) ! number of particles per bin
    type(bin_p), intent(inout) :: VH(n_bin) ! particle volumes (m^3)
    
    integer i_small, i_big, n_big

    n_big = nint(big_num)
    call particle_in_bin(small_vol, n_bin, bin_v, i_small)
    call particle_in_bin(big_vol, n_bin, bin_v, i_big)
    MH(i_small) = M - n_big
    MH(i_big) = n_big
    call enlarge_bin_to(VH(i_small), MH(i_small))
    call enlarge_bin_to(VH(i_big), MH(i_big))
    VH(i_small)%p = small_vol
    VH(i_big)%p = big_vol

  end subroutine init_bidisperse
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine init_log_normal(d_mean, log_sigma, n_bin, &
       bin_v, n_den)

    ! Compute a log-normal distribution.
    
    use mod_util
    use mod_constants
    
    real*8, intent(in) :: d_mean        ! mean diameter of initial dist (m)
    real*8, intent(in) :: log_sigma     ! log_e(geom. std dev(init dist)) (1)
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
  
  subroutine init_exp(V_0, n_bin, bin_v, n_den)
    
    ! Exponential distribution in volume
    ! n(v) = N_0 / V_0 exp(- v / V_0)
    
    use mod_bin
    use mod_util
    
    real*8, intent(in) :: V_0           ! mean volume of init dist (m^3)
    integer, intent(in) :: n_bin        ! number of bins
    real*8, intent(in) :: bin_v(n_bin)  ! volume of particles in bins (m^3)
    real*8, intent(out) :: n_den(n_bin) ! init number density (#(ln(r))d(ln(r)))
    
    integer k
    real*8 n_den_vol
    
    do k = 1,n_bin
       n_den_vol = 1d0 / V_0 * exp(-(bin_v(k) / V_0))
       call vol_to_lnr(vol2rad(bin_v(k)), n_den_vol, n_den(k))
    end do
    
  end subroutine init_exp
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine dist_to_n(N, dlnr, n_bin, bin_v, n_den, bin_n)
    
    ! Convert a number density (in ln(r)) to actual number of particles
    ! in each bin.
    
    integer, intent(in) :: N            ! total number of particles (approx)
    real*8, intent(in) :: dlnr          ! bin scale factor
    integer, intent(in) :: n_bin        ! number of bins
    real*8, intent(in) :: bin_v(n_bin)  ! volume of particles in bins (m^3)
    real*8, intent(in) :: n_den(n_bin)  ! init number density (#(ln(r))d(ln(r)))
                                        ! n_den(n_bin) has to be normalized
    integer, intent(out) :: bin_n(n_bin) ! number distribution
    
    integer k
    
    do k = 1,n_bin
       bin_n(k) = int(dble(N) * n_den(k) * dlnr)
    end do
    
  end subroutine dist_to_n
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
end module mod_init_dist
