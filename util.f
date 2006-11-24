! -*- mode: f90; -*-
! Copyright (C) 2005,2006 Nicole Riemer and Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.
!
! Common utility functions.

module mod_util
contains
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real*8 function util_rand()

    ! returns a random number between 0 and 1
#ifdef USE_F95_RAND
    real*8 rnd
    call random_number(rnd)
    util_rand = rnd
#else
    stop
    util_rand = dble(rand())
#endif

  end function util_rand
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  real*8 function vol2rad(v) ! radius (m)
    
    real*8, intent(in) :: v  ! volume (m^3)
    
    real*8 pi
    parameter (pi = 3.14159265358979323846d0)
    
    vol2rad = (v / (4d0 / 3d0 * pi))**(1d0/3d0)
    
  end function vol2rad
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  real*8 function vol2diam(v) ! diameter (m)
    
    real*8, intent(in) :: v  ! volume (m^3)
    
    real*8 pi
    parameter (pi = 3.14159265358979323846d0)
    
    vol2diam = 2d0 * (v / (4d0 / 3d0 * pi))**(1d0/3d0)
    
  end function vol2diam
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  real*8 function rad2vol(r) ! volume (m^3)
    
    real*8, intent(in) :: r  ! radius (m)
    
    real*8 pi
    parameter (pi = 3.14159265358979323846d0)
    
    rad2vol = 4d0 / 3d0 * pi * r**3d0
    
  end function rad2vol
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  real*8 function diam2vol(d) ! volume (m^3)
    
    real*8, intent(in) :: d  ! diameter (m)
    
    real*8 pi
    parameter (pi = 3.14159265358979323846d0)
    
    diam2vol = 4d0 / 3d0 * pi * (d / 2d0)**3d0
    
  end function diam2vol
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine sum_int_1d(size1, array, sum)
    
    integer, intent(in) :: size1             !  size of the array
    integer, intent(in) :: array(size1)      !  array of numbers
    integer, intent(out) :: sum               !  sum of array(i)
    
    integer i
    
    sum = 0
    do i = 1,size1
       sum = sum + array(i)
    enddo
    
  end subroutine sum_int_1d
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine max_int_2d(size1, size2, array, max)
    
    integer, intent(in) :: size1               !  leading size of the array
    integer, intent(in) :: size2               !  trailing size of the array
    integer, intent(in) :: array(size1,size2)  !  array of numbers
    integer, intent(out) :: max                 !  max of array(i,j)
    
    integer i, j
    
    max = 0
    do i = 1,size1
       do j = 1,size2
          if (array(i, j) .gt. max) then
             max = array(i, j)
          endif
       enddo
    enddo
    
  end subroutine max_int_2d
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  logical function almost_equal(d1, d2)
    
    ! Tests whether two real numbers are almost equal
    
    real*8, intent(in) :: d1 !  first number to compare
    real*8, intent(in) :: d2 !  second number to compare
    
    real*8 eps
    parameter (eps = 1d-8) ! relative tolerance
    
    ! handle the 0.0 case
    if (d1 .eq. d2) then
       almost_equal = .true.
    else
       if (abs(d1 - d2) / (abs(d1) + abs(d2)) .lt. eps) then
          almost_equal = .true.
       else
          almost_equal = .false.
       end if
    end if
    
  end function almost_equal
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
end module mod_util
