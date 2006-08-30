C     Common utility functions.

      module mod_util
      contains

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine vol2rad(v, r)

      real*8, intent(in) :: v  ! volume (m^3)
      real*8, intent(out) :: r ! radius (m)

      real*8 pi
      parameter (pi = 3.14159265358979323846d0)

      r = (v / (4d0 / 3d0 * pi))**(1d0/3d0)

      end subroutine

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine vol2diam(v, d)

      real*8, intent(in) :: v  ! volume (m^3)
      real*8, intent(out) :: d ! diameter (m)

      real*8 pi
      parameter (pi = 3.14159265358979323846d0)

      d = 2d0 * (v / (4d0 / 3d0 * pi))**(1d0/3d0)

      end subroutine

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine rad2vol(r, v)

      real*8, intent(in) :: r  ! radius (m)
      real*8, intent(out) :: v ! volume (m^3)

      real*8 pi
      parameter (pi = 3.14159265358979323846d0)

      v = 4d0 / 3d0 * pi * r**3d0

      end subroutine

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine diam2vol(d, v)

      real*8, intent(in) :: d  ! diameter (m)
      real*8, intent(out) :: v ! volume (m^3)

      real*8 pi
      parameter (pi = 3.14159265358979323846d0)

      v = 4d0 / 3d0 * pi * (d / 2d0)**3d0

      end subroutine

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

      end module
