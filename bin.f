C     Functions that deal with the bin grid.

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine bin_kernel(n_bin, bin_v, kernel, k)

C     Computes the kernel for each bin pair

      integer n_bin             ! INPUT: number of bins
      real*8 bin_v(n_bin)       ! INPUT: volume of particles in bins (m^3)
      external kernel           ! INPUT: kernel function
      real*8 k(n_bin, n_bin)    ! OUTPUT: kernel values

      integer i, j

      do i = 1,n_bin
         do j = 1,n_bin
            call kernel(bin_v(i), bin_v(j), k(i,j))
         enddo
      enddo
      
      end

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

      end

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

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine bin_n_to_g(n_bin, bin_v, bin_n, bin_g)

      integer n_bin             ! INPUT: number of bins
      real*8 bin_v(n_bin)       ! INPUT: volume of particles in bins (m^3)
      integer bin_n(n_bin)      ! INPUT: number in bins
      real*8 bin_g(n_bin)       ! OUTPUT: mass in bins

      integer i

      do i = 1,n_bin
         bin_g(i) = bin_n(i) * bin_v(i)
      enddo

      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine est_k_max_for_bin(n_bin, bin_v, kernel, b1, b2, k_max)

      integer n_bin             ! INPUT: number of bins
      real*8 bin_v(n_bin)       ! INPUT: volume of particles in bins (m^3)
      external kernel           ! INPUT: kernel function
      integer b1                ! INPUT: first bin
      integer b2                ! INPUT: second bin
      real*8 k_max              ! OUTPUT: maximum kernel values

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

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine est_k_max_binned(n_bin, bin_v, kernel, k_max)

      integer n_bin             ! INPUT: number of bins
      real*8 bin_v(n_bin)       ! INPUT: volume of particles in bins (m^3)
      external kernel           ! INPUT: kernel function
      real*8 k_max(n_bin,n_bin) ! OUTPUT: maximum kernel values

      integer i, j

      do i = 1,n_bin
         do j = 1,n_bin
            call est_k_max_for_bin(n_bin, bin_v, kernel, i, j, k_max(i,j
     $           ))
         enddo
      enddo

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      
      subroutine print_header(n_loop, n_bin, n_time)

      integer n_loop  ! INPUT: number of loops
      integer n_bin   ! INPUT: number of bins
      integer n_time  ! INPUT: number of times

      write(30,'(a10,i10)') 'n_loop', n_loop
      write(30,'(a10,i10)') 'n_bin', n_bin
      write(30,'(a10,i10)') 'n_time', n_time

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      
      subroutine print_info(time, V_comp,
     &     n_bin, bin_v, bin_r, bin_g, bin_n, dlnr)

      real*8 time          ! INPUT: cubin_rent simulation time
      real*8 V_comp        ! INPUT: computational volume
      
      integer n_bin        ! INPUT: number of bins
      real*8 bin_v(n_bin)  ! INPUT: volume of particles in bins (m^3)
      real*8 bin_r(n_bin)  ! INPUT: radius of particles in bins (m)
      real*8 bin_g(n_bin)  ! INPUT: mass in bins (???)
      integer bin_n(n_bin) ! INPUT: number in bins
      real*8 dlnr          ! INPUT: bin scale factor

      integer k

      write(30,'(a10,e14.5)') 'time', time
      do k = 1,n_bin
         write(30, '(i8,3e14.5)') k, bin_r(k), bin_n(k) / V_comp / dlnr,
     &        bin_g(k) / V_comp / dlnr
      enddo

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
