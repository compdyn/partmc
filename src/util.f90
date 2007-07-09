! Copyright (C) 2005-2007 Nicole Riemer and Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.
!
! Common utility functions.

module mod_util

  integer, parameter :: max_units = 200
  integer, parameter :: unit_offset = 10
  logical, save :: unit_used(max_units) = .false.

contains
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function get_unit()
    
    ! Returns an available unit number. This should be freed by free_unit().

    integer i
    logical found_unit

    found_unit = .false.
    do i = 1,max_units
       if (.not. unit_used(i)) then
          found_unit = .true.
          exit
       end if
    end do
    if (.not. found_unit) then
       write(0,*) 'ERROR: no more units available - need to free_unit()'
       call exit(1)
    end if
    unit_used(i) = .true.
    get_unit = i + unit_offset

  end function get_unit

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine free_unit(unit)

    ! Frees a unit number returned by get_unit().

    integer, intent(in) :: unit

    unit_used(unit - unit_offset) = .false.

  end subroutine free_unit

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real*8 function util_rand()

    ! Returns a random number between 0 and 1. Call this function
    ! rather than rand() or random_number() directly, as some
    ! compilers have trouble with one or the other.

!#ifdef USE_F95_RAND
    real*8 rnd
    call random_number(rnd)
    util_rand = rnd
!#else
!    util_rand = dble(rand())
!#endif

  end function util_rand
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function util_rand_disc(n)

    ! Returns a random integer between 1 and n.

    integer, intent(in) :: n            ! maximum random number to generate

    util_rand_disc = mod(int(util_rand() * dble(n)), n) + 1

  end function util_rand_disc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  real*8 function vol2rad(v)            ! radius (m)

    ! Convert volume to radius.

    use mod_constants
    
    real*8, intent(in) :: v             ! volume (m^3)
    
    vol2rad = (v / (4d0 / 3d0 * const%pi))**(1d0/3d0)
    
  end function vol2rad
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  real*8 function vol2diam(v)           ! diameter (m)
    
    ! Convert volume to diameter.

    use mod_constants
    
    real*8, intent(in) :: v             ! volume (m^3)
    
    vol2diam = 2d0 * (v / (4d0 / 3d0 * const%pi))**(1d0/3d0)
    
  end function vol2diam
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  real*8 function rad2vol(r)            ! volume (m^3)
    
    ! Convert radius to volume.

    use mod_constants
    
    real*8, intent(in) :: r             ! radius (m)
    
    rad2vol = 4d0 / 3d0 * const%pi * r**3d0
    
  end function rad2vol
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  real*8 function diam2vol(d)           ! volume (m^3)
    
    ! Convert diameter to volume.

    use mod_constants
    
    real*8, intent(in) :: d             ! diameter (m)
    
    diam2vol = 4d0 / 3d0 * const%pi * (d / 2d0)**3d0
    
  end function diam2vol
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  logical function almost_equal(d1, d2)
    
    ! Tests whether two real numbers are almost equal using only a
    ! relative tolerance.
    
    real*8, intent(in) :: d1            ! first number to compare
    real*8, intent(in) :: d2            ! second number to compare
    
    real*8, parameter :: eps = 1d-8     ! relative tolerance
    
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
    
    ! Tests whether two real numbers are almost equal using an
    ! absolute and relative tolerance.
    
    real*8, intent(in) :: d1            ! first number to compare
    real*8, intent(in) :: d2            ! second number to compare
    real*8, intent(in) :: abs_tol       ! tolerance for when d1 equals d2
    
    real*8, parameter :: eps = 1d-8     ! relative tolerance
    
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
    
    real*8, intent(in) :: time          ! current time
    real*8, intent(in) :: timestep      ! estimate of the time to the next call
    real*8, intent(in) :: interval      ! how often the event should be done
    real*8, intent(inout) :: last_time  ! when the event was last done
    logical, intent(out) :: do_event    ! whether the event should be done
    
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

  subroutine open_existing(unit, filename)

    ! Opens a file for reading that must already exist, checking for
    ! errors.

    integer, intent(in) :: unit         ! unit to open with
    character(len=*), intent(in) :: filename ! filename of file to open

    integer ios

    open(unit=unit, file=filename, status='old', iostat=ios)
    if (ios /= 0) then
       write(0,*) 'ERROR: unable to open file ', trim(filename), &
            ': IOSTAT = ', ios
       call exit(1)
    end if

  end subroutine open_existing

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function find_1d(n, x_vals, x) ! position of x

    ! Takes an array of x_vals, and a single x value, and returns the
    ! position p such that x_vals(p) <= x < x_vals(p+1). If p == 0
    ! then x < x_vals(1) and if p == n then x_vals(n) <= x. x_vals
    ! must be sorted in increasing order.

    integer, intent(in) :: n            ! number of values
    real*8, intent(in) :: x_vals(n)     ! x value array, must be sorted
    real*8, intent(in) :: x             ! value to interpolate at

    integer p

    if (n == 0) then
       find_1d = 0
       return
    end if
    p = 1
    do while (x >= x_vals(p))
       p = p + 1
       if (p > n) then
          exit
       end if
    end do
    p = p - 1
    find_1d = p

  end function find_1d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real*8 function interp_1d(n, x_vals, y_vals, x) ! y value at x

    ! Takes an array of x and y, and a single x value, and returns the
    ! corresponding y using linear interpolation. x_vals must be
    ! sorted.

    integer, intent(in) :: n            ! number of values
    real*8, intent(in) :: x_vals(n)     ! x value array, must be sorted
    real*8, intent(in) :: y_vals(n)     ! y value array
    real*8, intent(in) :: x             ! value to interpolate at

    integer p
    real*8 y, alpha

    p = find_1d(n, x_vals, x)
    if (p == 0) then
       y = y_vals(1)
    elseif (p == n) then
       y = y_vals(n)
    else
       alpha = (x - x_vals(p)) / (x_vals(p+1) - x_vals(p))
       y = (1d0 - alpha) * y_vals(p) + alpha * y_vals(p+1)
    end if
    interp_1d = y

  end function interp_1d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function string_to_integer(string)

    ! Convert a string to an integer.

    character(len=*), intent(in) :: string ! string to convert
    
    integer :: val
    integer :: ios

    read(string, '(i20)', iostat=ios) val
    if (ios /= 0) then
       write(0,'(a,a,a,i3)') 'Error converting ', string, &
            'to integer: IOSTAT = ', ios
       call exit(1)
    end if
    string_to_integer = val

  end function string_to_integer

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real*8 function string_to_real(string)

    ! Convert a string to an real.

    character(len=*), intent(in) :: string ! string to convert
    
    real*8 :: val
    integer :: ios

    read(string, '(f20.0)', iostat=ios) val
    if (ios /= 0) then
       write(0,'(a,a,a,i3)') 'Error converting ', string, &
            'to real: IOSTAT = ', ios
       call exit(1)
    end if
    string_to_real = val

  end function string_to_real

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  logical function string_to_logical(string)

    ! Convert a string to an logical.

    character(len=*), intent(in) :: string ! string to convert
    
    logical :: val
    integer :: ios

    val = .false.
    if ((trim(string) == 'yes') &
         .or. (trim(string) == 'y') &
         .or. (trim(string) == 'true') &
         .or. (trim(string) == 't') &
         .or. (trim(string) == '1')) then
       val = .true.
    elseif ((trim(string) == 'no') &
         .or. (trim(string) == 'n') &
         .or. (trim(string) == 'false') &
         .or. (trim(string) == 'f') &
         .or. (trim(string) == '0')) then
       val = .false.
    else
       write(0,'(a,a,a)') 'Error converting ', string, 'to logical'
       call exit(1)
    end if
    string_to_logical = val

  end function string_to_logical

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function sample_cts_pdf(n, pdf)

    ! Sample the given probability density function. That is,
    ! return a number k = 1,...,n such that prob(k) = pdf(k) / sum(pdf).

    integer, intent(in) :: n            ! number of entries
    real*8, intent(in) :: pdf(n)        ! probability density function
                                        ! (not normalized)

    real*8 :: pdf_max
    integer :: k
    logical :: found

    ! use accept-reject
    pdf_max = maxval(pdf)
    if (minval(pdf) < 0d0) then
       write(0,*) 'ERROR: pdf contains negative values'
       call exit(1)
    end if
    if (pdf_max <= 0d0) then
       write(*,*) 'ERROR: pdf is not positive'
       call exit(1)
    end if
    found = .false.
    do while (.not. found)
       k = util_rand_disc(n)
       if (util_rand() < pdf(k) / pdf_max) then
          found = .true.
       end if
    end do
    sample_cts_pdf = k

  end function sample_cts_pdf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function sample_disc_pdf(n, pdf)

    ! Sample the given probability density function. That is,
    ! return a number k = 1,...,n such that prob(k) = pdf(k) / sum(pdf).

    integer, intent(in) :: n            ! number of entries
    integer, intent(in) :: pdf(n)       ! probability density function

    integer :: pdf_max, k
    logical :: found

    ! use accept-reject
    pdf_max = maxval(pdf)
    if (minval(pdf) < 0) then
       write(0,*) 'ERROR: pdf contains negative values'
       call exit(1)
    end if
    if (pdf_max <= 0) then
       write(*,*) 'ERROR: pdf is not positive'
       call exit(1)
    end if
    found = .false.
    do while (.not. found)
       k = util_rand_disc(n)
       if (util_rand() < dble(pdf(k)) / dble(pdf_max)) then
          found = .true.
       end if
    end do
    sample_disc_pdf = k

  end function sample_disc_pdf
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine sample_vec_cts_to_disc(n, vec_cts, n_samp, vec_disc)
    
    ! Convert a continuous vector into a discrete vector with n_samp
    ! samples. Each discrete entry is sampled with a PDF given by
    ! vec_cts. This is very slow for large n_samp or large n.
    
    integer, intent(in) :: n          ! number of entries in vector
    real*8, intent(in) :: vec_cts(n)  ! continuous vector
    integer, intent(in) :: n_samp     ! number of discrete samples to use
    integer, intent(out) :: vec_disc(n) ! discretized vector

    integer :: i_samp, k

    vec_disc = 0
    do i_samp = 1,n_samp
       k = sample_cts_pdf(n, vec_cts)
       vec_disc(k) = vec_disc(k) + 1
    end do
    
  end subroutine sample_vec_cts_to_disc
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine vec_cts_to_disc(n, vec_cts, n_samp, vec_disc)
    
    ! Convert a continuous vector into a discrete vector with n_samp
    ! samples. Returns discrete vector whose relative entry sizes
    ! are \ell_1 closest to the continuous vector.
    
    integer, intent(in) :: n          ! number of entries in vector
    real*8, intent(in) :: vec_cts(n)  ! continuous vector
    integer, intent(in) :: n_samp     ! number of discrete samples to use
    integer, intent(out) :: vec_disc(n) ! discretized vector
    
    integer :: k(1)
    real*8 :: vec_tot
    
    vec_tot = sum(vec_cts)
    
    ! assign a best guess for each bin independently
    vec_disc = nint(vec_cts / vec_tot * dble(n_samp))
    
    ! if we have too few particles then add more
    do while (sum(vec_disc) < n_samp)
       k = minloc(abs(dble(vec_disc + 1) - vec_cts) &
            - abs(dble(vec_disc) - vec_cts))
       vec_disc(k) = vec_disc(k) + 1
    end do
    
    ! if we have too many particles then remove some
    do while (sum(vec_disc) > n_samp)
       k = minloc(abs(dble(vec_disc - 1) - vec_cts) &
            - abs(dble(vec_disc) - vec_cts))
       vec_disc(k) = vec_disc(k) - 1
    end do
    
    ! asserts
    if (sum(vec_disc) /= n_samp) then
       write(0,*) 'ERROR: generated incorrect number of samples'
       call exit(1)
    end if
    
  end subroutine vec_cts_to_disc
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine average_integer(int_vec, int_avg)
    
    ! Computes the average of an array of integer numbers.
    
    integer, intent(in) :: int_vec(:)   ! array of integer numbers
    integer, intent(out) :: int_avg     ! average of int_vec
    
    int_avg = sum(int_vec) / size(int_vec)
    
  end subroutine average_integer
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine average_real(real_vec, real_avg)
    
    ! Computes the average of an array of real numbers.
    
    real*8, intent(in) :: real_vec(:) ! array of real numbers
    real*8, intent(out) :: real_avg   ! average of real_vec
    
    real_avg = sum(real_vec) / dble(size(real_vec))
    
  end subroutine average_real
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
end module mod_util
