C     Common utility functions.

      module mod_util
      contains

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      real*8 function vol2rad(v) ! radius (m)

      real*8, intent(in) :: v  ! volume (m^3)

      real*8 pi
      parameter (pi = 3.14159265358979323846d0)

      vol2rad = (v / (4d0 / 3d0 * pi))**(1d0/3d0)

      end function

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      real*8 function vol2diam(v) ! diameter (m)

      real*8, intent(in) :: v  ! volume (m^3)

      real*8 pi
      parameter (pi = 3.14159265358979323846d0)

      vol2diam = 2d0 * (v / (4d0 / 3d0 * pi))**(1d0/3d0)

      end function

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      real*8 function rad2vol(r) ! volume (m^3)

      real*8, intent(in) :: r  ! radius (m)

      real*8 pi
      parameter (pi = 3.14159265358979323846d0)

      rad2vol = 4d0 / 3d0 * pi * r**3d0

      end function

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      real*8 function diam2vol(d) ! volume (m^3)

      real*8, intent(in) :: d  ! diameter (m)

      real*8 pi
      parameter (pi = 3.14159265358979323846d0)

      diam2vol = 4d0 / 3d0 * pi * (d / 2d0)**3d0

      end function

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine sum_int_1d(size1, array, sum)

      integer size1             ! INPUT: size of the array
      integer array(size1)      ! INPUT: array of numbers
      integer sum               ! OUTPUT: sum of array(i)
      
      integer i
      
      sum = 0
      do i = 1,size1
         sum = sum + array(i)
      enddo
      
      end subroutine

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine max_int_2d(size1, size2, array, max)

      integer size1               ! INPUT: leading size of the array
      integer size2               ! INPUT: trailing size of the array
      integer array(size1,size2)  ! INPUT: array of numbers
      integer max                 ! OUTPUT: max of array(i,j)

      integer i, j

      max = 0
      do i = 1,size1
         do j = 1,size2
            if (array(i, j) .gt. max) then
               max = array(i, j)
            endif
         enddo
      enddo

      end subroutine

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      logical function almost_equal(d1, d2)

      ! Tests whether two real numbers are almost equal

      real*8 d1 ! INPUT: first number to compare
      real*8 d2 ! INPUT: second number to compare

      real*8 eps
      parameter (eps = 1d-8) ! relative tolerance

C     handle the 0.0 case
      if (d1 .eq. d2) then
         almost_equal = .true.
      else
         if (abs(d1 - d2) / (abs(d1) + abs(d2)) .lt. eps) then
            almost_equal = .true.
         else
            almost_equal = .false.
         end if
      end if

      end function

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      end module
