! -*- mode: f90; -*-
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
contains
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine bin_kernel(n_bin, bin_v, kernel, env, k)
   
    use mod_environ
    
    ! Computes an array of kernel values for each bin pair. k(i,j) is
    ! the kernel value at the midpoint of bins i and j.
    
    integer, intent(in) :: n_bin        ! number of bins
    real*8, intent(in) :: bin_v(n_bin)  ! volume of particles in bins (m^3)
    real*8, intent(out) :: k(n_bin,n_bin) ! kernel values
    type(environ), intent(in) :: env    ! environment state

    interface
       subroutine kernel(v1, v2, env, k)
         use mod_environ
         real*8, intent(in) :: v1
         real*8, intent(in) :: v2
         type(environ), intent(in) :: env
         real*8, intent(out) :: k
       end subroutine kernel
    end interface
    
    integer i, j
    
    do i = 1,n_bin
       do j = 1,n_bin
          call kernel(bin_v(i), bin_v(j), env, k(i,j))
       end do
    end do
    
  end subroutine bin_kernel
  
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
    
    dlnr = dlog(2d0) / (3d0 * dble(scal))
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
    if ((k .lt. n_bin) .and. &
         (v .gt. (bin_v(k) + bin_v(k+1)) / 2d0)) goto 300
    
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
  
  subroutine est_k_max_for_bin(n_bin, bin_v, kernel, b1, b2, env, k_max)

    ! Samples within bins b1 and b2 to find the maximum value of the
    ! kernel between particles from the two bins.
   
    use mod_environ
 
    integer, intent(in) :: n_bin        ! number of bins
    real*8, intent(in) :: bin_v(n_bin)  ! volume of particles in bins (m^3)
    integer, intent(in) :: b1           ! first bin
    integer, intent(in) :: b2           ! second bin
    type(environ), intent(in) :: env    ! environment state    
    real*8, intent(out) :: k_max        ! maximum kernel values
    
    interface
       subroutine kernel(v1, v2, env, k)
         use mod_environ
         real*8, intent(in) :: v1
         real*8, intent(in) :: v2
         type(environ), intent(in) :: env
         real*8, intent(out) :: k
       end subroutine kernel
    end interface
    
    real*8 v1, v2, v1_high, v1_low, v2_high, v2_low, k
    integer i, j
    
    integer, parameter :: n_sample = 10  ! number of sample points per bin
    
    ! v1_low < bin_v(b1) < v1_high
    call bin_edge(n_bin, bin_v, b1, v1_low)
    call bin_edge(n_bin, bin_v, b1 + 1, v1_high)
    
    ! v2_low < bin_v(b2) < v2_high
    call bin_edge(n_bin, bin_v, b2, v2_low)
    call bin_edge(n_bin, bin_v, b2 + 1, v2_high)
    
    k_max = 0d0
    do i = 1,n_sample
       do j = 1,n_sample
          v1 = v1_high * dble(n_sample - i) / dble(n_sample - 1) + &
               v1_low * dble(i - 1) / dble(n_sample - 1)
          v2 = v2_high * dble(n_sample - j) / dble(n_sample - 1) + &
               v2_low * dble(j - 1) / dble(n_sample - 1)
          call kernel(v1, v2, env, k)
          if (k .gt. k_max) k_max = k
       end do
    end do
    
  end subroutine est_k_max_for_bin
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine est_k_max_binned(n_bin, bin_v, kernel, env, k_max)

    ! Computes an array of maximum kernel values. Given particles v1
    ! in bin b1 and v2 in bin b2, it is approximately true that
    ! kernel(v1,v2) <= k_max(b1,b2).

    use mod_environ
    
    integer, intent(in) :: n_bin        ! number of bins
    real*8, intent(in) :: bin_v(n_bin)  ! volume of particles in bins (m^3)
    type(environ), intent(in) :: env  ! environment state
    real*8, intent(out) :: k_max(n_bin,n_bin) ! maximum kernel values
    
    interface
       subroutine kernel(v1, v2, env, k)
         use mod_environ
         real*8, intent(in) :: v1
         real*8, intent(in) :: v2
         type(environ), intent(in) :: env
         real*8, intent(out) :: k
       end subroutine kernel
    end interface
    
    integer i, j
    
    do i = 1,n_bin
       do j = 1,n_bin
          call est_k_max_for_bin(n_bin, bin_v, kernel, i, j, &
               env, k_max(i,j))
       end do
    end do
    
  end subroutine est_k_max_binned
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine output_open(output_unit, output_name, n_loop, n_bin, &
       n_spec, n_time)

    ! Open an output file for writing.

    integer, intent(in) :: output_unit  ! unit number to output to
    character(len=300) :: output_name   ! name of output
    integer, intent(in) :: n_loop       ! number of loops
    integer, intent(in) :: n_bin        ! number of bins
    integer, intent(in) :: n_spec       ! number of species
    integer, intent(in) :: n_time       ! number of times

    character(len=300) :: out_file_name

    write(out_file_name, '(a,a,a)') 'out_', trim(output_name), '.d'
    open(unit=output_unit, file=out_file_name)

    write(output_unit,'(a10,i10)') 'n_loop', n_loop
    write(output_unit,'(a10,i10)') 'n_bin', n_bin
    write(output_unit,'(a10,i10)') 'n_time', n_time
    write(output_unit,'(a10,i10)') 'n_spec', n_spec
    
  end subroutine output_open
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine output_info(output_unit, time, n_bin, n_spec, &
       bin_v, bin_g, bin_gs, bin_n, dlnr, env, mat, i_loop)

    ! Write the current binned data to the output file. This version
    ! of the function takes absolute number and absolute volume
    ! per-bin (as produced by a particle-resolved code, for example).
    
    use mod_material
    use mod_environ
    
    integer, intent(in) :: output_unit  ! unit number to output to
    real*8, intent(in) :: time          ! simulation time
    integer, intent(in) :: n_bin        ! number of bins
    integer, intent(in) :: n_spec       ! number of species
    real*8, intent(in) :: bin_v(n_bin)  ! volume of particles in bins (m^3)
    real*8, intent(in) :: bin_g(n_bin)  ! volume in bins (m^3)
    real*8, intent(in) :: bin_gs(n_bin,n_spec) ! species volume in bins
    integer, intent(in) :: bin_n(n_bin) ! number in bins (dimensionless)
    real*8, intent(in) :: dlnr          ! bin scale factor
    type(environ), intent(in) :: env    ! environment state
    type(material), intent(in) :: mat   ! material properties
    integer, intent(in) :: i_loop       ! current loop number
    
    real*8 bin_g_den(n_bin), bin_gs_den(n_bin,n_spec)
    real*8 bin_n_den(n_bin)
    
    bin_g_den = bin_g / env%V_comp / dlnr
    bin_gs_den = bin_gs / env%V_comp / dlnr
    bin_n_den = dble(bin_n) / env%V_comp / dlnr
    call output_info_density(output_unit, time, n_bin, n_spec, bin_v, &
         bin_g_den, bin_gs_den, bin_n_den, env, mat, i_loop)
    
  end subroutine output_info
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine output_info_density(output_unit, time, n_bin, n_spec, bin_v, &
       bin_g_den, bin_gs_den, bin_n_den, env, mat, i_loop)

    ! Write the current binned data to the output file. This version
    ! of the function takes number and volume densities (as produced
    ! by a sectional code, for example).
    
    use mod_material
    use mod_environ
    use mod_util
    
    integer, intent(in) :: output_unit  ! unit number to output to
    real*8, intent(in) :: time          ! simulation time
    integer, intent(in) :: n_bin        ! number of bins
    integer, intent(in) :: n_spec       ! number of species
    real*8, intent(in) :: bin_v(n_bin)  ! volume of particles in bins (m^3)
    real*8, intent(in) :: bin_g_den(n_bin) ! volume density in bins (1)
    real*8, intent(in) :: bin_gs_den(n_bin,n_spec) ! spec vol den in bins (1)
    real*8, intent(in) :: bin_n_den(n_bin) ! number density in bins (1/m^3)
    type(environ), intent(in) :: env    ! environment state
    type(material), intent(in) :: mat   ! material properties
    integer, intent(in) :: i_loop       ! current loop number
    
    integer k

    write(output_unit,'(a10,i20)') 'loop_num', i_loop
    write(output_unit,'(a10,e20.10)') 'time(s)', time
    write(output_unit,'(a10,e20.10)') 'temp(K)', env%T
    write(output_unit,'(a10,e20.10)') 'RH(1)', env%RH
    write(output_unit,'(a1,a9,a20,a20,a20,a30)') '#', 'bin_num', &
         'radius(m)', 'tot num (#/m^3)', 'tot vol (m^3/m^3)', &
         'vol per species (m^3/m^3)'
    do k = 1,n_bin
       write(output_unit, '(i10,20e20.10)') k, vol2rad(bin_v(k)), &
            bin_n_den(k), bin_g_den(k), bin_gs_den(k,:)
    end do
    
  end subroutine output_info_density
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
end module mod_bin
