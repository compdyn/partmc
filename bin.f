! -*- mode: f90; -*-
! Copyright (C) 2005,2006 Nicole Riemer and Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.
C     Functions that deal with the bin grid.

      module mod_bin
      contains

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine bin_kernel(n_bin, bin_v, kernel, k)

C     Computes the kernel for each bin pair

      integer n_bin             ! INPUT: number of bins
      real*8 bin_v(n_bin)       ! INPUT: volume of particles in bins (m^3)
      real*8 k(n_bin, n_bin)    ! OUTPUT: kernel values

      interface
         subroutine kernel(v1, v2, k)
         real*8, intent(in) :: v1
         real*8, intent(in) :: v2
         real*8, intent(out) :: k
         end subroutine
      end interface

      integer i, j

      do i = 1,n_bin
         do j = 1,n_bin
            call kernel(bin_v(i), bin_v(j), k(i,j))
         enddo
      enddo
      
      end subroutine

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine vol_to_lnr(r, f_vol, f_lnr)

C     Convert a density f(vol)d(vol) to f(ln(r))d(ln(r))
C     where vol = 4 pi r^3.

      real*8 r                  ! INPUT: radius (m)
      real*8 f_vol              ! INPUT: density as a function of volume
      real*8 f_lnr              ! OUTPUT: density as a function of ln(r)
      
      real*8 pi
      parameter (pi = 3.14159265358979323846d0)

      f_lnr = f_vol * 4d0 * pi * r**3

      end subroutine

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine make_bin_grid(n_bin, scal, v_min, bin_v, bin_r,
     $     dlnr)

      integer n_bin        ! INPUT: number of bins
      integer scal         ! INPUT: scale factor
      real*8 v_min         ! INPUT: minimum volume (m^3)
      real*8 bin_v(n_bin)  ! OUTPUT: volume of particles in bins (m^3)
      real*8 bin_r(n_bin)  ! OUTPUT: radius of particles in bins (m)
      real*8 dlnr          ! OUTPUT: scale factor

      integer i
      real*8 ax

      real*8 pi
      parameter (pi = 3.14159265358979323846d0)

      dlnr = dlog(2d0) / (3d0 * dble(scal))
      ax = 2d0**(1d0 / dble(scal)) ! ratio bin_v(i)/bin_v(i-1)

      do i = 1,n_bin
         ! volume (m^3)
         bin_v(i) = v_min * 0.5d0 * (ax + 1d0) * ax**(i - 1)
         ! radius (m)
         bin_r(i) = dexp(dlog(3d0 * bin_v(i) / 
     &             (4d0 * pi)) / 3d0)
      enddo

      end subroutine

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine bin_edge(n_bin, bin_v, i, v_edge)

      integer n_bin             ! INPUT: number of bins
      real*8 bin_v(n_bin)       ! INPUT: volume of particles in bins (m^3)
      integer i                 ! INPUT: edge number (1 <= i <= n_bin + 1)
      real*8 v_edge             ! OUTPUT: volume at edge

      if (i .eq. 1) then
         v_edge = bin_v(1) - (bin_v(2) - bin_v(1)) / 2d0
      elseif (i .eq. (n_bin + 1)) then
         v_edge = bin_v(n_bin) + (bin_v(n_bin) - bin_v(n_bin - 1)) / 2d0
      else
         v_edge = (bin_v(i - 1) + bin_v(i)) / 2d0
      endif

      end subroutine

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine particle_in_bin(v, n_bin, bin_v, k)
      ! FIXME: for log-spaced bins we can do this without search

      real*8 v             ! INPUT: volume of particle
      integer n_bin        ! INPUT: number of bins
      real*8 bin_v(n_bin)  ! INPUT: volume of particles in bins (m^3)
      integer k            ! OUTPUT: bin number containing particle

      k = 0
 300  k = k + 1
      if ((k .lt. n_bin) .and.
     &     (v .gt. (bin_v(k) + bin_v(k+1)) / 2d0)) goto 300

      end subroutine

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine bin_n_to_g(n_bin, bin_v, bin_n, bin_g)

      integer n_bin             ! INPUT: number of bins
      real*8 bin_v(n_bin)       ! INPUT: volume of particles in bins (m^3)
      integer bin_n(n_bin)      ! INPUT: number in bins
      real*8 bin_g(n_bin)       ! OUTPUT: mass in bins

      integer i

      do i = 1,n_bin
         bin_g(i) = dble(bin_n(i)) * bin_v(i)
      enddo

      end subroutine

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine est_k_max_for_bin(n_bin, bin_v, kernel, b1, b2, k_max)

      integer n_bin             ! INPUT: number of bins
      real*8 bin_v(n_bin)       ! INPUT: volume of particles in bins (m^3)
      integer b1                ! INPUT: first bin
      integer b2                ! INPUT: second bin
      real*8 k_max              ! OUTPUT: maximum kernel values

      interface
         subroutine kernel(v1, v2, k)
         real*8, intent(in) :: v1
         real*8, intent(in) :: v2
         real*8, intent(out) :: k
         end subroutine
      end interface

      real*8 v1, v2, v1_high, v1_low, v2_high, v2_low, k
      integer i, j

      integer n_sample
      parameter (n_sample = 10)  ! number of sample points per bin

      ! v1_low < bin_v(b1) < v1_high
      call bin_edge(n_bin, bin_v, b1, v1_low)
      call bin_edge(n_bin, bin_v, b1 + 1, v1_high)

      ! v2_low < bin_v(b2) < v2_high
      call bin_edge(n_bin, bin_v, b2, v2_low)
      call bin_edge(n_bin, bin_v, b2 + 1, v2_high)

      k_max = 0d0
      do i = 1,n_sample
         do j = 1,n_sample
            v1 = v1_high * dble(n_sample - i) / dble(n_sample - 1) +
     $           v1_low * dble(i - 1) / dble(n_sample - 1)
            v2 = v2_high * dble(n_sample - j) / dble(n_sample - 1) +
     $           v2_low * dble(j - 1) / dble(n_sample - 1)
            call kernel(v1, v2, k)
            if (k .gt. k_max) k_max = k
         enddo
      enddo

      end subroutine

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine est_k_max_binned(n_bin, bin_v, kernel, k_max)

      integer n_bin             ! INPUT: number of bins
      real*8 bin_v(n_bin)       ! INPUT: volume of particles in bins (m^3)
      real*8 k_max(n_bin,n_bin) ! OUTPUT: maximum kernel values

      interface
         subroutine kernel(v1, v2, k)
         real*8, intent(in) :: v1
         real*8, intent(in) :: v2
         real*8, intent(out) :: k
         end subroutine
      end interface

      integer i, j

      do i = 1,n_bin
         do j = 1,n_bin
            call est_k_max_for_bin(n_bin, bin_v, kernel, i, j,
     &           k_max(i,j))
         enddo
      enddo

      end subroutine

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      
      subroutine print_header(n_loop, n_bin, n_spec, n_time)

      integer n_loop  ! INPUT: number of loops
      integer n_bin   ! INPUT: number of bins
      integer n_time  ! INPUT: number of times
      integer n_spec  ! INPUT: number of species
 
      write(30,'(a10,i10)') 'n_loop', n_loop
      write(30,'(a10,i10)') 'n_bin', n_bin
      write(30,'(a10,i10)') 'n_time', n_time
      write(30,'(a10,i10)') 'n_spec', n_spec

      end subroutine

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C FIXME: delete
C      
C      subroutine print_info(time, V_comp, n_spec,
C     &     n_bin, bin_v, bin_r, bin_g, bin_gs, bin_n, dlnr, env, mat)
C
C      use mod_material
C      use mod_environ
C
C      real*8 time          ! INPUT: cubin_rent simulation time
C      real*8 V_comp        ! INPUT: computational volume
C      
C      integer n_bin        ! INPUT: number of bins
C      integer n_spec       ! INPUT: number of species
C      real*8 bin_v(n_bin)  ! INPUT: volume of particles in bins (m^3)
C      real*8 bin_r(n_bin)  ! INPUT: radius of particles in bins (m)
C      real*8 bin_g(n_bin)  ! INPUT: mass in bins (???)
C      real*8 bin_gs(n_bin,n_spec) !INPUT: species mass in bins
C      integer bin_n(n_bin) ! INPUT: number in bins
C      real*8 dlnr          ! INPUT: bin scale factor
C      type(environ), intent(inout) :: env  ! environment state
C      type(material), intent(in) :: mat    ! material properties
C
C      integer k, i
C
C      write(30,'(a10,e20.10)') 'time(s)', time
C      write(30,'(a10,e20.10)') 'temp(K)', env%T
C      write(30,'(a10,e20.10)') 'rh(1)', env%RH
C      do k = 1,n_bin
C         write(30, '(i10,20e20.10)')k, bin_r(k),
C     &        dble(bin_n(k)) / V_comp / dlnr,
C     &        bin_g(k) / V_comp / dlnr,
C     &        (bin_gs(k,i) / V_comp / dlnr,i=1,n_spec) 
C      enddo
C
C      end subroutine
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      
      subroutine print_info(time, V_comp, n_spec,
     &     n_bin, bin_v, bin_r, bin_g, bin_gs, bin_n, dlnr, env, mat)

      use mod_material
      use mod_environ

      real*8 time          ! INPUT: simulation time
      real*8 V_comp        ! INPUT: computational volume
      
      integer n_bin        ! INPUT: number of bins
      integer n_spec       ! INPUT: number of species
      real*8 bin_v(n_bin)  ! INPUT: volume of particles in bins (m^3)
      real*8 bin_r(n_bin)  ! INPUT: radius of particles in bins (m)
      real*8 bin_g(n_bin)  ! INPUT: volume in bins (m^3)
      real*8 bin_gs(n_bin,n_spec) !INPUT: species mass in bins
      integer bin_n(n_bin) ! INPUT: number in bins (dimensionless)
      real*8 dlnr          ! INPUT: bin scale factor
      type(environ), intent(inout) :: env  ! environment state
      type(material), intent(in) :: mat    ! material properties

      real*8 bin_g_den(n_bin), bin_gs_den(n_bin,n_spec)
      real*8 bin_n_den(n_bin)

      bin_g_den = bin_g / V_comp / dlnr
      bin_gs_den = bin_gs / V_comp / dlnr
      bin_n_den = dble(bin_n) / V_comp / dlnr
      call print_info_density(time, n_bin, n_spec, bin_v, bin_r,
     $     bin_g_den, bin_gs_den, bin_n_den, env, mat)

      end subroutine

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      
      subroutine print_info_density(time, n_bin, n_spec, bin_v, bin_r,
     &     bin_g_den, bin_gs_den, bin_n_den, env, mat)

      use mod_material
      use mod_environ

      real*8 time              ! INPUT: simulation time
      
      integer n_bin            ! INPUT: number of bins
      integer n_spec       ! INPUT: number of species
      real*8 bin_v(n_bin)      ! INPUT: volume of particles in bins (m^3)
      real*8 bin_r(n_bin)      ! INPUT: radius of particles in bins (m)
      real*8 bin_g_den(n_bin)  ! INPUT: volume density in bins (dimensionless)
      real*8 bin_gs_den(n_bin,n_spec)  ! INPUT: species volume density in bins
                               ! (dimensionless)
      real*8 bin_n_den(n_bin)  ! INPUT: number density in bins (1/m^3)
      type(environ), intent(inout) :: env  ! environment state
      type(material), intent(in) :: mat    ! material properties

      integer k

      write(30,'(a10,e20.10)') 'time(s)', time
      write(30,'(a10,e20.10)') 'temp(K)', env%T
      write(30,'(a10,e20.10)') 'rh(1)', env%RH
      do k = 1,n_bin
         write(30, '(i10,20e20.10)') k, bin_r(k), bin_n_den(k),
     &        bin_g_den(k), bin_gs_den(k,:)
      enddo

      end subroutine

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      end module
