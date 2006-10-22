! -*- mode: f90; -*-
! Copyright (C) 2005,2006 Nicole Riemer and Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.
!
! Initial size distributions.

module mod_init_dist
contains
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine init_bidisperse(MM, n_bin, n_ini)
    
    integer MM           ! INPUT: physical dimension of V
    integer n_bin        ! INPUT: number of bins
    integer n_ini(n_bin) ! OUTPUT: initial number distribution
    
    integer k
    
    if (MM .lt. 126) then
       write(*,*)'ERROR: MM too small for bidisperse'
       call exit(2)
    endif
    
    do k = 1,n_bin
       n_ini(k) = 0
    enddo
    n_ini(97) = MM - 1
    n_ini(126) = 1
    
    return
  end subroutine init_bidisperse
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine init_log_normal(MM, d_mean, log_sigma, dlnr, n_bin, &
       bin_v, bin_r, n_ini)
    
    ! FIXME: make this return a number density
    
    integer MM           ! INPUT: physical dimension of V
    real*8 d_mean        ! INPUT: mean diameter of initial distribution (m)
    real*8 log_sigma     ! INPUT: log_e of the geometric standard
    ! deviation of initial distribution (1)
    real*8 dlnr          ! INPUT: bin scale factor
    integer n_bin        ! INPUT: number of bins
    real*8 bin_v(n_bin)  ! INPUT: volume of particles in bins (m^3)
    real*8 bin_r(n_bin)  ! INPUT: radius of particles in bins (m)
    integer n_ini(n_bin) ! OUTPUT: initial number distribution
    
    real*8 pi
    parameter (pi = 3.14159265358979323846d0)
    
    integer k
    
    do k = 1,n_bin
       n_ini(k) = int(dble(MM) / (sqrt(2d0 * pi) * log_sigma) * &
            dexp(-(dlog10(bin_r(k)) - dlog10(d_mean/2d0))**2d0 &
            / (2d0 * log_sigma**2d0)) * dlnr / dlog(10d0))
    enddo
    
    ! The formula above was originally for a distribution in
    ! log_10(r), while we are using log_e(r). The division by dlog(10)
    ! at the end corrects for this.
    
  end subroutine init_log_normal
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine dist_to_n(N, dlnr, n_bin, bin_v, bin_r, n_den, bin_n)
    
    ! Convert a number density (in ln(r)) to actual number of particles
    ! in each bin.
    
    integer N            ! INPUT: total number of particles (approximate)
    real*8 dlnr          ! INPUT: bin scale factor
    integer n_bin        ! INPUT: number of bins
    real*8 bin_v(n_bin)  ! INPUT: volume of particles in bins (m^3)
    real*8 bin_r(n_bin)  ! INPUT: radius of particles in bins (m)
    real*8 n_den(n_bin)  ! INPUT: initial number density (#(ln(r))d(ln(r)))
    integer bin_n(n_bin) ! OUTPUT: number distribution
    
    integer k
    
    do k = 1,n_bin
       bin_n(k) = int(dble(N) * n_den(k) * dlnr)
    enddo
    
  end subroutine dist_to_n
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine init_exp(V_0, n_bin, bin_v, bin_r, n_den)
    
    ! Exponential distribution in volume (or mass)
    ! n(v) = N_0 / V_0 exp(- v / V_0)
    
    use mod_bin
    
    real*8 V_0           ! INPUT: mean volume of initial distribution (m^3)
    integer n_bin        ! INPUT: number of bins
    real*8 bin_v(n_bin)  ! INPUT: volume of particles in bins (m^3)
    real*8 bin_r(n_bin)  ! INPUT: radius of particles in bins (m)
    real*8 n_den(n_bin)  ! OUTPUT: initial number density (#(ln(r))d(ln(r)))
    
    integer k
    real*8 n_den_vol
    
    do k = 1,n_bin
       n_den_vol = 1d0 / V_0 * exp(-(bin_v(k) / V_0))
       call vol_to_lnr(bin_r(k), n_den_vol, n_den(k))
    enddo
    
  end subroutine init_exp
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
end module mod_init_dist
