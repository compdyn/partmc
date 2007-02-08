! -*- mode: f90; -*-
! Copyright (C) 2005-2007 Nicole Riemer and Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.
!
! Functions that deal with the bin grid.

module mod_bin
contains
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine bin_kernel(n_bin, bin_v, kernel, env, k)
   
    use mod_environ
    
    ! Computes the kernel for each bin pair
    
    integer, intent(in) :: n_bin            ! number of bins
    real*8, intent(in) :: bin_v(n_bin)      ! volume of particles in bins (m^3)
    real*8, intent(out) :: k(n_bin,n_bin)   ! kernel values
    type(environ), intent(in) :: env        ! environment state

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
    
    real*8, intent(in) :: r                  ! radius (m)
    real*8, intent(in) :: f_vol              ! density as a function of volume
    real*8, intent(out) :: f_lnr              ! density as a function of ln(r)
    
    f_lnr = f_vol * 4d0 * const%pi * r**3
    
  end subroutine vol_to_lnr
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine make_bin_grid(n_bin, scal, v_min, bin_v, &
       dlnr)
    
    integer, intent(in) :: n_bin        ! number of bins
    integer, intent(in) :: scal         ! scale factor
    real*8, intent(in) :: v_min         ! minimum volume (m^3)
    real*8, intent(out) :: bin_v(n_bin)  ! volume of particles in bins (m^3)
    real*8, intent(out) :: dlnr          ! scale factor
    
    integer i
    real*8 ax
    
    dlnr = dlog(2d0) / (3d0 * dble(scal))
    ax = 2d0**(1d0 / dble(scal)) ! ratio bin_v(i)/bin_v(i-1)
    
    do i = 1,n_bin
       ! volume (m^3)
       bin_v(i) = v_min * 0.5d0 * (ax + 1d0) * ax**(i - 1)
    end do
    
  end subroutine make_bin_grid
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine bin_edge(n_bin, bin_v, i, v_edge)
    
    integer, intent(in) :: n_bin             ! number of bins
    real*8, intent(in) :: bin_v(n_bin)       ! volume of particles in bins (m^3)
    integer, intent(in) :: i                 ! edge number (1 <= i <= n_bin + 1)
    real*8, intent(out) :: v_edge             ! volume at edge
    
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
    ! FIXME: for log-spaced bins we can do this without search
    
    real*8, intent(in) :: v             ! volume of particle
    integer, intent(in) :: n_bin        ! number of bins
    real*8, intent(in) :: bin_v(n_bin)  ! volume of particles in bins (m^3)
    integer, intent(out) :: k            ! bin number containing particle
    
    k = 0
300 k = k + 1
    if ((k .lt. n_bin) .and. &
         (v .gt. (bin_v(k) + bin_v(k+1)) / 2d0)) goto 300
    
  end subroutine particle_in_bin
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine bin_n_to_g(n_bin, bin_v, bin_n, bin_g)
    
    integer, intent(in) :: n_bin             ! number of bins
    real*8, intent(in) :: bin_v(n_bin)       ! volume of particles in bins (m^3)
    integer, intent(in) :: bin_n(n_bin)      ! number in bins
    real*8, intent(out) :: bin_g(n_bin)      ! volume in bins
    
    integer i
    
    do i = 1,n_bin
       bin_g(i) = dble(bin_n(i)) * bin_v(i)
    end do
    
  end subroutine bin_n_to_g
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine est_k_max_for_bin(n_bin, bin_v, kernel, b1, b2, env, k_max)
   
    use mod_environ
 
    integer, intent(in) :: n_bin             ! number of bins
    real*8, intent(in) :: bin_v(n_bin)       ! volume of particles in bins (m^3)
    integer, intent(in) :: b1                ! first bin
    integer, intent(in) :: b2                ! second bin
    real*8, intent(out) :: k_max              ! maximum kernel values

    type(environ), intent(in) :: env  ! environment state    
    
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

    use mod_environ
    
    integer, intent(in) :: n_bin             ! number of bins
    real*8, intent(in) :: bin_v(n_bin)       ! volume of particles in bins (m^3)
    real*8, intent(out) :: k_max(n_bin,n_bin) ! maximum kernel values
    type(environ), intent(in) :: env  ! environment state
    
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

    integer, intent(in) :: output_unit ! unit number to output to
    character(len=300) :: output_name ! name of output
    integer, intent(in) :: n_loop  ! number of loops
    integer, intent(in) :: n_bin   ! number of bins
    integer, intent(in) :: n_time  ! number of times
    integer, intent(in) :: n_spec  ! number of species

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
       bin_v, bin_g, bin_gs, bin_n, dlnr, env, mat)
    
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
    
    real*8 bin_g_den(n_bin), bin_gs_den(n_bin,n_spec)
    real*8 bin_n_den(n_bin)
    
    bin_g_den = bin_g / env%V_comp / dlnr
    bin_gs_den = bin_gs / env%V_comp / dlnr
    bin_n_den = dble(bin_n) / env%V_comp / dlnr
    call output_info_density(output_unit, time, n_bin, n_spec, bin_v, &
         bin_g_den, bin_gs_den, bin_n_den, env, mat)
    
  end subroutine output_info
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine output_info_density(output_unit, time, n_bin, n_spec, bin_v, &
       bin_g_den, bin_gs_den, bin_n_den, env, mat)
    
    use mod_material
    use mod_environ
    use mod_util
    
    integer, intent(in) :: output_unit   ! unit number to output to
    real*8, intent(in) :: time           ! simulation time
    integer, intent(in) :: n_bin         ! number of bins
    integer, intent(in) :: n_spec        ! number of species
    real*8, intent(in) :: bin_v(n_bin)   ! volume of particles in bins (m^3)
    real*8, intent(in) :: bin_g_den(n_bin) ! volume density in bins (dimensionless)
    real*8, intent(in) :: bin_gs_den(n_bin,n_spec) ! species volume density in bins
    ! (dimensionless)
    real*8, intent(in) :: bin_n_den(n_bin)  ! number density in bins (1/m^3)
    type(environ), intent(in) :: env  ! environment state
    type(material), intent(in) :: mat    ! material properties
    
    integer k
    
    write(output_unit,'(a10,e20.10)') 'time(s)', time
    write(output_unit,'(a10,e20.10)') 'temp(K)', env%T
    write(output_unit,'(a10,e20.10)') 'rh(1)', env%RH
    do k = 1,n_bin
       write(output_unit, '(i10,20e20.10)') k, vol2rad(bin_v(k)), &
            bin_n_den(k), bin_g_den(k), bin_gs_den(k,:)
    end do
    
  end subroutine output_info_density
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
end module mod_bin
