! Copyright (C) 2005-2007 Nicole Riemer and Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.
!
! Functions that deal with the bin grid.
!
! The grid of bins is logarithmically spaced, an assumption that is
! quite heavily incorporated into the code. At some point in the
! future it would be nice to relax this assumption.

module mod_bin

  type bin_grid_t
     integer :: n_bin
     real*8, pointer :: v(:)            ! len n_bin, bin center volumes (m^3)
     real*8, pointer :: dlnr            ! bin scale factor (1)
  end type bin_grid_t

  type bin_dist_t
     real*8, pointer :: v(:)            ! len n_bin, volume per bin (m^3)
     real*8, pointer :: vs(:,:)         ! n_bin x n_spec, vol per bin&spec (m^3)
     integer, pointer :: n(:)           ! len n_bin, number per bin (1)
  end type bin_dist_t

contains
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine vol_to_lnr(r, f_vol, f_lnr)
    
    ! Convert a density f(vol)d(vol) to f(ln(r))d(ln(r))
    ! where vol = 4 pi r^3.
    
    use mod_constants
    
    real*8, intent(in) :: r             ! radius (m)
    real*8, intent(in) :: f_vol         ! density as a function of volume
    real*8, intent(out) :: f_lnr        ! density as a function of ln(r)
    
    f_lnr = f_vol * 4d0 * const%pi * r**3
    
  end subroutine vol_to_lnr
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine make_bin_grid(n_bin, scal, v_min, bin_v, dlnr)

    ! Generates the bin grid, given the minimum volume, number of grid
    ! points, and the scale factor.
    
    integer, intent(in) :: n_bin        ! number of bins
    integer, intent(in) :: scal         ! scale factor
    real*8, intent(in) :: v_min         ! minimum volume (m^3)
    real*8, intent(out) :: bin_v(n_bin) ! volume of particles in bins (m^3)
    real*8, intent(out) :: dlnr         ! scale factor
    
    integer i
    real*8 ax
    
    dlnr = dlog(2d0) / (3d0 * dble(scal)) ! ln(r(i) / r(i-1))
    ax = 2d0**(1d0 / dble(scal)) ! ratio bin_v(i)/bin_v(i-1)
    
    do i = 1,n_bin
       bin_v(i) = v_min * 0.5d0 * (ax + 1d0) * ax**(i - 1)  ! (m^3)
    end do
    
  end subroutine make_bin_grid
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine bin_edge(n_bin, bin_v, i, v_edge)

    ! Given a bin grid (which stores the center points of the bins),
    ! find the given edge volume.
    
    integer, intent(in) :: n_bin        ! number of bins
    real*8, intent(in) :: bin_v(n_bin)  ! volume of particles in bins (m^3)
    integer, intent(in) :: i            ! edge number (1 <= i <= n_bin + 1)
    real*8, intent(out) :: v_edge       ! volume at edge
    
    if (i .eq. 1) then
       v_edge = bin_v(1) - (bin_v(2) - bin_v(1)) / 2d0
    elseif (i .eq. (n_bin + 1)) then
       v_edge = bin_v(n_bin) + (bin_v(n_bin) - bin_v(n_bin - 1)) / 2d0
    else
       v_edge = (bin_v(i - 1) + bin_v(i)) / 2d0
    end if
    
  end subroutine bin_edge
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine particle_in_bin(v, n_bin, bin_v, k)
    
    ! Find the bin number that contains a given particle.
    
    real*8, intent(in) :: v             ! volume of particle
    integer, intent(in) :: n_bin        ! number of bins
    real*8, intent(in) :: bin_v(n_bin)  ! volume of particles in bins (m^3)
    integer, intent(out) :: k           ! bin number containing particle
    
    ! FIXME: for log-spaced bins we can do this without search, but we
    ! plan to switch to arbitrary bins at some point, so maybe just
    ! leave it like it is for now.

    k = 0
300 k = k + 1
    if (k .lt. n_bin) then
       if (v .gt. (bin_v(k) + bin_v(k+1)) / 2d0) then
          goto 300
       end if
    end if
    
  end subroutine particle_in_bin
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine bin_n_to_g(n_bin, bin_v, bin_n, bin_g)

    ! Convert a number of particles per bin to the total volume per
    ! bin.
    
    integer, intent(in) :: n_bin        ! number of bins
    real*8, intent(in) :: bin_v(n_bin)  ! volume of particles in bins (m^3)
    integer, intent(in) :: bin_n(n_bin) ! number in bins
    real*8, intent(out) :: bin_g(n_bin) ! volume in bins
    
    integer i
    
    do i = 1,n_bin
       bin_g(i) = dble(bin_n(i)) * bin_v(i)
    end do
    
  end subroutine bin_n_to_g
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
end module mod_bin
