! -*- mode: f90; -*-
! Copyright (C) 2005-2007 Nicole Riemer and Matthew West
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
    util_rand = dble(rand())
#endif

  end function util_rand
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  real*8 function vol2rad(v) ! radius (m)

    use mod_constants
    
    real*8, intent(in) :: v  ! volume (m^3)
    
    vol2rad = (v / (4d0 / 3d0 * const%pi))**(1d0/3d0)
    
  end function vol2rad
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  real*8 function vol2diam(v) ! diameter (m)
    
    use mod_constants
    
    real*8, intent(in) :: v  ! volume (m^3)
    
    vol2diam = 2d0 * (v / (4d0 / 3d0 * const%pi))**(1d0/3d0)
    
  end function vol2diam
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  real*8 function rad2vol(r) ! volume (m^3)
    
    use mod_constants
    
    real*8, intent(in) :: r  ! radius (m)
    
    rad2vol = 4d0 / 3d0 * const%pi * r**3d0
    
  end function rad2vol
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  real*8 function diam2vol(d) ! volume (m^3)
    
    use mod_constants
    
    real*8, intent(in) :: d  ! diameter (m)
    
    diam2vol = 4d0 / 3d0 * const%pi * (d / 2d0)**3d0
    
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
  
  logical function almost_equal_abs(d1, d2, abs_tol)
    
    ! Tests whether two real numbers are almost equal
    
    real*8, intent(in) :: d1 !  first number to compare
    real*8, intent(in) :: d2 !  second number to compare
    real*8, intent(in) :: abs_tol ! tolerance for when d1 equals d2
    
    real*8 eps
    parameter (eps = 1d-8) ! relative tolerance
    
    ! handle the 0.0 case
    if (d1 .eq. d2) then
       almost_equal_abs = .true.
    else
       if ((abs(d1 - d2) .lt. abs_tol) .or. &
            (abs(d1 - d2) / (abs(d1) + abs(d2)) .lt. eps)) then
          almost_equal_abs = .true.
       else
          almost_equal_abs = .false.
       end if
    end if
    
  end function almost_equal_abs
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine check_event(time, timestep, interval, last_time, &
       do_event)
    
    ! Computes whether an event is scheduled to take place. The events
    ! should occur ideally at times 0, interval, 2*interval, etc. The
    ! events are guaranteed to occur at least interval * (1 -
    ! tolerance) apart, and if at least interval time has passed then
    ! the next call is guaranteed to do the event. Otherwise the
    ! timestep is used to guess whether to do the event.
    
    real*8, intent(in) :: time       !  current time
    real*8, intent(in) :: timestep   !  an estimate of the time to the next call
    real*8, intent(in) :: interval   !  how often the event should be done
    real*8, intent(inout) :: last_time  !  when the event was last done
    logical, intent(out) :: do_event  !  whether the event should be done
    
    real*8, parameter :: tolerance = 1d-6 ! fuzz for event occurance
    
    real*8 closest_interval_time
    
    ! if we are at time 0 then do the event unconditionally
    if (time .eq. 0d0) then
       do_event = .true.
    else
       ! if we are too close to the last time then don't do it
       if ((time - last_time) .lt. interval * (1d0 - tolerance)) then
          do_event = .false.
       else
          ! if it's been too long since the last time then do it
          if ((time - last_time) .ge. interval) then
             do_event = .true.
          else
             ! gray area -- if we are closer than we will be next
             ! time then do it
             closest_interval_time = anint(time / interval) * interval
             if (abs(time - closest_interval_time) &
                  .lt. abs(time + timestep - closest_interval_time)) &
                  then
                do_event = .true.
             else
                do_event = .false.
             end if
          end if
       end if
    end if
    
    if (do_event) then
       last_time = time
    end if
    
  end subroutine check_event
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
end module mod_util
