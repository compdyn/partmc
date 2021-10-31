! Copyright (C) 2005-2016 Nicole Riemer and Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_util module.

!> Common utility subroutines.
module pmc_util

  use pmc_constants
#ifdef PMC_USE_MPI
  use mpi
#endif
#ifdef PMC_USE_C_SORT
  use iso_c_binding
#endif

  !> Maximum number of IO units usable simultaneously.
  integer, parameter :: max_units = 200
  !> Minimum unit number to allocate.
  integer, parameter :: unit_offset = 10
  !> Table of unit numbers storing allocation status.
  logical, save :: unit_used(max_units) = .false.

  !> Length of string for converting numbers.
  integer, parameter :: PMC_UTIL_CONVERT_STRING_LEN = 100
  !> Maximum length of filenames.
  integer, parameter :: PMC_MAX_FILENAME_LEN = 300

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Prints a warning message.
  subroutine warn_msg(code, warning_msg, already_warned)

    !> Status code to use.
    integer, intent(in) :: code
    !> Message to display.
    character(len=*), intent(in) :: warning_msg
    !> Flag to control warning only once (should be a save variable).
    logical, intent(inout), optional :: already_warned

    if (present(already_warned)) then
       if (already_warned) return
       ! set already_warned so next time we will immediately return
       already_warned = .true.
    end if
    write(0,'(a)') 'WARNING (PartMC-' // trim(integer_to_string(code)) &
         // '): ' // trim(warning_msg)

  end subroutine warn_msg

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Prints a warning message if condition_ok is false.
  subroutine warn_assert_msg(code, condition_ok, warning_msg)

    !> Status code to use.
    integer, intent(in) :: code
    !> Whether the assertion is ok.
    logical, intent(in) :: condition_ok
    !> Message to display.
    character(len=*), intent(in) :: warning_msg

    if (.not. condition_ok) then
       call warn_msg(code, warning_msg)
    end if

  end subroutine warn_assert_msg

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Errors unless condition_ok is true.
  subroutine assert_msg(code, condition_ok, error_msg)

    !> Status code to use if assertion fails.
    integer, intent(in) :: code
    !> Whether the assertion is ok.
    logical, intent(in) :: condition_ok
    !> Msg if assertion fails.
    character(len=*), intent(in) :: error_msg

    integer :: ierr

    if (.not. condition_ok) then
       write(0,'(a)') 'ERROR (PartMC-' // trim(integer_to_string(code)) &
            // '): ' // trim(error_msg)
#ifdef PMC_USE_MPI
       call mpi_abort(MPI_COMM_WORLD, code, ierr)
#else
       stop 3
#endif
    end if

  end subroutine assert_msg

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Errors unless condition_ok is true.
  subroutine assert(code, condition_ok)

    !> Status code to use if assertion fails.
    integer, intent(in) :: code
    !> Whether the assertion is ok.
    logical, intent(in) :: condition_ok

    if (.not. condition_ok) then
       ! much faster with gfortran 4.4.5 to do this extra test
       ! FIXME: is it still faster now that assert_msg doesn't
       ! unconditionally construct a code_str?
       call assert_msg(code, condition_ok, 'assertion failed')
    end if

  end subroutine assert

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Error immediately.
  subroutine die(code)

    !> Status code to use if assertion fails.
    integer, intent(in) :: code

    call assert(code, .false.)

  end subroutine die

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Error immediately.
  subroutine die_msg(code, error_msg)

    !> Status code to use if assertion fails.
    integer, intent(in) :: code
    !> Msg if assertion fails.
    character(len=*), intent(in) :: error_msg

    call assert_msg(code, .false., error_msg)

  end subroutine die_msg

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns an available unit number. This should be freed by free_unit().
  integer function get_unit()

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
       call die_msg(690355443, &
            'no more units available - need to free_unit()')
    end if
    unit_used(i) = .true.
    get_unit = i + unit_offset

  end function get_unit

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Frees a unit number returned by get_unit().
  subroutine free_unit(unit)

    integer, intent(in) :: unit

    unit_used(unit - unit_offset) = .false.

  end subroutine free_unit

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Open a file for reading with an automatically assigned unit and
  !> test that it succeeds. The file should be closed with
  !> close_file().
  subroutine open_file_read(filename, unit)

    !> Filename to open.
    character(len=*), intent(in) :: filename
    !> Unit assigned to file.
    integer, intent(out) :: unit

    integer :: ios

    unit = get_unit()
    open(unit=unit, file=filename, status='old', action='read', iostat=ios)
    call assert_msg(544344918, ios == 0, 'unable to open file ' &
         // trim(filename) // ' for reading: ' // integer_to_string(ios))

  end subroutine open_file_read

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Open a file for writing with an automatically assigned unit and
  !> test that it succeeds. The file should be closed with
  !> close_file().
  subroutine open_file_write(filename, unit)

    !> Filename to open.
    character(len=*), intent(in) :: filename
    !> Unit assigned to file.
    integer, intent(out) :: unit

    integer :: ios

    unit = get_unit()
    open(unit=unit, file=filename, status='replace', action='write', &
         iostat=ios)
    call assert_msg(609624199, ios == 0, 'unable to open file ' &
         // trim(filename) // ' for writing: ' // integer_to_string(ios))

  end subroutine open_file_write

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Close a file and de-assign the unit.
  subroutine close_file(unit)

    !> Unit to close.
    integer, intent(in) :: unit

    close(unit)
    call free_unit(unit)

  end subroutine close_file

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Convert mass-equivalent volume \f$V\f$ (m^3) to geometric radius
  !> \f$R_{\rm geo}\f$ (m) for spherical particles.
  real(kind=dp) elemental function sphere_vol2rad(v)

    !> Particle mass-equivalent volume (m^3).
    real(kind=dp), intent(in) :: v

    sphere_vol2rad = (3d0 * v / 4d0 / const%pi)**(1d0 / 3d0)

  end function sphere_vol2rad

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Convert radius (m) to diameter (m).
  real(kind=dp) elemental function rad2diam(r)

    !> Radius (m).
    real(kind=dp), intent(in) :: r

    rad2diam = 2d0 * r

  end function rad2diam

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Convert geometric radius \f$R_{\rm geo}\f$ (m) to mass-equivalent volume
  !> \f$V\f$ (m^3) for spherical particles.
  real(kind=dp) elemental function sphere_rad2vol(r)

    !> Geometric radius (m).
    real(kind=dp), intent(in) :: r

    sphere_rad2vol = 4d0 * const%pi * r**3 / 3d0

  end function sphere_rad2vol

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Convert diameter (m) to radius (m).
  real(kind=dp) elemental function diam2rad(d)

    !> Diameter (m).
    real(kind=dp), intent(in) :: d

    diam2rad = d / 2d0

  end function diam2rad

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Calculate air molecular mean free path \f$l\f$ (m).
  real(kind=dp) function air_mean_free_path(temp, pressure)

    !> Temperature (K).
    real(kind=dp), intent(in) :: temp
    !> Pressure (Pa).
    real(kind=dp), intent(in) :: pressure

    real(kind=dp) :: boltz, avogad, mwair, rgas, rhoair, viscosd, &
         viscosk, gasspeed

    boltz = const%boltzmann
    avogad = const%avagadro
    mwair = const%air_molec_weight
    rgas = const%univ_gas_const

    rhoair = (pressure * mwair) / (rgas * temp)

    viscosd = (1.8325d-5 * (296.16d0 + 120d0) / (temp + 120d0)) &
         * (temp / 296.16d0)**1.5d0
    viscosk = viscosd / rhoair
    gasspeed = sqrt(8d0 * boltz * temp * avogad / (const%pi * mwair))
    air_mean_free_path = 2d0 * viscosk / gasspeed

  end function air_mean_free_path

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Tests whether two real numbers are almost equal using only a
  !> relative tolerance.
  logical function almost_equal(d1, d2)

    !> First number to compare.
    real(kind=dp), intent(in) :: d1
    !> Second number to compare.
    real(kind=dp), intent(in) :: d2

    !> Relative tolerance.
    real(kind=dp), parameter :: eps = 1d-8

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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Tests whether two real numbers are almost equal using an
  !> absolute and relative tolerance.
  logical function almost_equal_abs(d1, d2, abs_tol)

    !> First number to compare.
    real(kind=dp), intent(in) :: d1
    !> Second number to compare.
    real(kind=dp), intent(in) :: d2
    !> Tolerance for when d1 equals d2.
    real(kind=dp), intent(in) :: abs_tol

    !> Relative tolerance.
    real(kind=dp), parameter :: eps = 1d-8

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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Check that the first time interval is close to an integer
  !> multiple of the second, and warn if it is not.
  subroutine check_time_multiple(first_name, first_time, &
       second_name, second_time)

    !> Name of the first time variable.
    character(len=*), intent(in) :: first_name
    !> First time variable (s).
    real(kind=dp), intent(in) :: first_time
    !> Name of the second time variable.
    character(len=*), intent(in) :: second_name
    !> Second time variable (s).
    real(kind=dp), intent(in) :: second_time

    real(kind=dp) :: ratio

    ratio = first_time / second_time
    if (abs(ratio - aint(ratio)) > 1d-6) then
       call warn_msg(952299377, trim(first_name) &
            // " is not an integer multiple of " // trim(second_name))
    end if

  end subroutine check_time_multiple

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Computes whether an event is scheduled to take place.
  !!
  !! The events should occur ideally at times 0, interval, 2*interval,
  !! etc. The events are guaranteed to occur at least interval * (1 -
  !! tolerance) apart, and if at least interval time has passed then
  !! the next call is guaranteed to do the event. Otherwise the
  !! timestep is used to guess whether to do the event.
  subroutine check_event(time, timestep, interval, last_time, &
       do_event)

    !> Current time.
    real(kind=dp), intent(in) :: time
    !> Estimate of the time to the next call.
    real(kind=dp), intent(in) :: timestep
    !> How often the event should be done.
    real(kind=dp), intent(in) :: interval
    !> When the event was last done.
    real(kind=dp), intent(inout) :: last_time
    !> Whether the event should be done.
    logical, intent(out) :: do_event

    !> Fuzz for event occurance.
    real(kind=dp), parameter :: tolerance = 1d-6

    real(kind=dp) closest_interval_time

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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Makes a linearly spaced array from min to max.
  function linspace(min_x, max_x, n)

    !> Minimum array value.
    real(kind=dp), intent(in) :: min_x
    !> Maximum array value.
    real(kind=dp), intent(in) :: max_x
    !> Length of array to create.
    integer, intent(in) :: n

    !> Return value.
    real(kind=dp), allocatable :: linspace(:)

    integer :: i
    real(kind=dp) :: a

    allocate(linspace(n))
    call assert(999299119, n >= 0)
    do i = 2, (n - 1)
       a = real(i - 1, kind=dp) / real(n - 1, kind=dp)
       linspace(i) = (1d0 - a) * min_x + a * max_x
    end do
    if (n > 0) then
       ! make sure these values are exact
       linspace(1) = min_x
       linspace(n) = max_x
    end if

  end function linspace

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Makes a logarithmically spaced array of length n from min to max.
  function logspace(min_x, max_x, n)

    !> Minimum array value.
    real(kind=dp), intent(in) :: min_x
    !> Maximum array value.
    real(kind=dp), intent(in) :: max_x
    !> Length of array to create.
    integer, intent(in) :: n

    !> Return value.
    real(kind=dp), allocatable :: logspace(:)

    real(kind=dp), allocatable :: log_x(:)

    allocate(logspace(n))

    call assert(804623592, n >= 0)
    if (n == 0) return
    call assert(548290438, min_x > 0d0)
    call assert(805259035, max_x > 0d0)
    log_x = linspace(log(min_x), log(max_x), n)
    logspace = exp(log_x)
    ! make sure these values are exact
    logspace(1) = min_x
    logspace(n) = max_x

  end function logspace

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Find the position of a real number in a 1D linear array.
  !!
  !! If xa is the array allocated by linspace(min_x, max_x, xa) then i
  !! = linspace_find(min_x, max_x, n, x) returns the index i
  !! satisfying xa(i) <= x < xa(i+1) for min_x <= x < max_x. If x ==
  !! max_x then i = n - 1.  If x > max_x then i = n. If x < min_x then
  !! i = 0. Thus 0 <= i <= n. Here n is the length of xa.
  !!
  !! This is equivalent to using find_1d() but much faster if the
  !! array is linear.
  integer function linspace_find(min_x, max_x, n, x)

    !> Minimum array value.
    real(kind=dp), intent(in) :: min_x
    !> Maximum array value.
    real(kind=dp), intent(in) :: max_x
    !> Number of entries.
    integer, intent(in) :: n
    !> Value.
    real(kind=dp), intent(in) :: x

    if (x == max_x) then
       linspace_find = n - 1
       return
    end if
    linspace_find = floor((x - min_x) / (max_x - min_x) &
         * real(n - 1, kind=dp)) + 1
    linspace_find = min(linspace_find, n)
    linspace_find = max(linspace_find, 0)

  end function linspace_find

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Find the position of a real number in a 1D logarithmic array.
  !!
  !! If xa is the array allocated by logspace(min_x, max_x, xa) then i
  !! = logspace_find(min_x, max_x, n, x) returns the index i
  !! satisfying xa(i) <= x < xa(i+1) for min_x <= x < max_x. If x >=
  !! max_x then i = n.  If x < min_x then i = 0. Thus 0 <= i <=
  !! n. Here n is the length of xa.
  !!
  !! This is equivalent to using find_1d() but much faster if the
  !! array is logarithmic.
  integer function logspace_find(min_x, max_x, n, x)

    !> Minimum array value.
    real(kind=dp), intent(in) :: min_x
    !> Maximum array value.
    real(kind=dp), intent(in) :: max_x
    !> Number of entries.
    integer, intent(in) :: n
    !> Value.
    real(kind=dp), intent(in) :: x

    logspace_find = linspace_find(log(min_x), log(max_x), n, log(x))

  end function logspace_find

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Find the position of a real number in an arbitrary 1D array.
  !!
  !! Takes an array of x_vals, and a single x value, and returns the
  !! position p such that x_vals(p) <= x < x_vals(p+1). If p == 0 then
  !! x < x_vals(1) and if p == n then x_vals(n) <= x. x_vals must be
  !! sorted in increasing order.
  !!
  !! If the array is known to be linearly or logarithmically spaced
  !! then linspace_find() or logspace_find() will be much faster.
  integer function find_1d(n, x_vals, x)

    !> Number of values.
    integer, intent(in) :: n
    !> X value array, must be sorted.
    real(kind=dp), intent(in) :: x_vals(n)
    !> Value to interpolate at.
    real(kind=dp), intent(in) :: x

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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> 1D linear interpolation.
  !!
  !! Takes an array of x and y, and a single x value, and returns the
  !! corresponding y using linear interpolation. x_vals must be
  !! sorted.
  real(kind=dp) function interp_1d(x_vals, y_vals, x)

    !> X value array, must be sorted.
    real(kind=dp), intent(in) :: x_vals(:)
    !> Y value array.
    real(kind=dp), intent(in) :: y_vals(size(x_vals))
    !> Value to interpolate at.
    real(kind=dp), intent(in) :: x

    integer :: n, p
    real(kind=dp) :: y, alpha

    n = size(x_vals)
    p = find_1d(n, x_vals, x)
    if (p == 0) then
       y = y_vals(1)
    elseif (p == n) then
       y = y_vals(n)
    else
       if (y_vals(p) == y_vals(p+1)) then
          y = y_vals(p)
       else
          alpha = (x - x_vals(p)) / (x_vals(p+1) - x_vals(p))
          y = (1d0 - alpha) * y_vals(p) + alpha * y_vals(p+1)
       end if
    end if
    interp_1d = y

  end function interp_1d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Linear interpolation over discrete indices.
  !!
  !! Takes real values \c x_1 and \c x_n and integers \c n and \c i
  !! and returns the linear interpolation so that \c x_1 is returned
  !! when \c i = 1 and \c x_n is returned when \c i = \c n.
  real(kind=dp) function interp_linear_disc(x_1, x_n, n, i)

    !> Value of \c x when \c i = 1.
    real(kind=dp), intent(in) :: x_1
    !> Value of \c x when \c i = n.
    real(kind=dp), intent(in) :: x_n
    !> Number of points to interpolate over.
    integer, intent(in) :: n
    !> Index to interpolate at.
    integer, intent(in) :: i

    real(kind=dp) :: alpha

    if (x_1 == x_n) then
       interp_linear_disc = x_1
    else
       alpha = real(i - 1, kind=dp) / real(n - 1, kind=dp)
       interp_linear_disc = (1d0 - alpha) * x_1 + alpha * x_n
    end if

  end function interp_linear_disc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Convert a string to an integer.
  integer function string_to_integer(string)

    !> String to convert.
    character(len=*), intent(in) :: string

    integer :: val
    integer :: ios
    character(len=len(string)+300) :: error_msg

    call assert_msg(447772570, len_trim(string) <= 20, &
         'error converting "' // trim(string) &
         // '" to integer: string too long')
    read(string, '(i20)', iostat=ios) val
    call assert_msg(895881873, ios == 0, &
         'error converting "' // trim(string) &
         // '" to integer: IOSTAT = ' // trim(integer_to_string(ios)))
    string_to_integer = val

  end function string_to_integer

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Convert a string to a real.
  real(kind=dp) function string_to_real(string)

    !> String to convert.
    character(len=*), intent(in) :: string

    real(kind=dp) :: val
    integer :: ios
    character(len=len(string)+300) :: error_msg

    call assert_msg(733728030, len_trim(string) <= 30, &
         'error converting "' // trim(string) // '" to real: string too long')
    read(string, '(f30.0)', iostat=ios) val
    call assert_msg(727430097, ios == 0, &
         'error converting "' // trim(string) &
         // '" to real: IOSTAT = ' // trim(integer_to_string(ios)))
    string_to_real = val

  end function string_to_real

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Convert a string to a logical.
  logical function string_to_logical(string)

    !> String to convert.
    character(len=*), intent(in) :: string

    logical :: val
    integer :: ios
    character(len=len(string)+300) :: error_msg

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
       call die_msg(985010153, 'error converting "' // trim(string) &
            // '" to logical')
    end if
    string_to_logical = val

  end function string_to_logical

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Convert an integer to a string format.
  character(len=PMC_UTIL_CONVERT_STRING_LEN) function integer_to_string(val)

    !> Value to convert.
    integer, intent(in) :: val

    character(len=PMC_UTIL_CONVERT_STRING_LEN) :: ret_val

    ret_val = ""
    write(ret_val, '(i30)') val
    integer_to_string = adjustl(ret_val)

  end function integer_to_string

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Convert a real to a string format.
  character(len=PMC_UTIL_CONVERT_STRING_LEN) function real_to_string(val)

    !> Value to convert.
    real(kind=dp), intent(in) :: val

    character(len=PMC_UTIL_CONVERT_STRING_LEN) :: ret_val

    ret_val = ""
    write(ret_val, '(g30.20)') val
    real_to_string = adjustl(ret_val)

  end function real_to_string

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Convert a logical to a string format.
  character(len=PMC_UTIL_CONVERT_STRING_LEN) function logical_to_string(val)

    !> Value to convert.
    logical, intent(in) :: val

    character(len=PMC_UTIL_CONVERT_STRING_LEN) :: ret_val

    ret_val = ""
    if (val) then
       ret_val = "TRUE"
    else
       ret_val = "FALSE"
    end if
    logical_to_string = ret_val

  end function logical_to_string

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Convert a complex to a string format.
  character(len=PMC_UTIL_CONVERT_STRING_LEN) function complex_to_string(val)

    !> Value to convert.
    complex(kind=dc), intent(in) :: val

    character(len=PMC_UTIL_CONVERT_STRING_LEN) :: ret_val

    ret_val = ""
    ret_val = "(" // trim(real_to_string(real(val))) &
         // ", " // trim(real_to_string(aimag(val))) // ")"
    complex_to_string = adjustl(ret_val)

  end function complex_to_string

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Convert an integer to a string format of maximum length.
  character(len=PMC_UTIL_CONVERT_STRING_LEN) &
       function integer_to_string_max_len(val, max_len)

    !> Value to convert.
    integer, intent(in) :: val
    !> Maximum length of resulting string.
    integer, intent(in) :: max_len

    character(len=PMC_UTIL_CONVERT_STRING_LEN) :: ret_val

    ret_val = integer_to_string(val)
    if (len_trim(ret_val) > max_len) then
       ret_val = real_to_string_max_len(real(val, kind=dp), max_len)
    end if
    integer_to_string_max_len = adjustl(ret_val)

  end function integer_to_string_max_len

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Convert a real to a string format of maximum length.
  character(len=PMC_UTIL_CONVERT_STRING_LEN) &
       function real_to_string_max_len(val, max_len)

    !> Value to convert.
    real(kind=dp), intent(in) :: val
    !> Maximum length of resulting string.
    integer, intent(in) :: max_len

    character(len=PMC_UTIL_CONVERT_STRING_LEN) :: ret_val, exp_str, frac_str
    integer :: exp_val, exp_len, frac_len, use_frac_len, min_frac_len, i
    real(kind=dp) :: frac_val

    ret_val = ""
    if (val == 0d0) then
       if (max_len >= 3) then
          ret_val = "0e0"
       else
          do i = 1,max_len
             ret_val(i:i) = "*"
          end do
       end if
       real_to_string_max_len = adjustl(ret_val)
       return
    end if

    exp_val = floor(log10(abs(val)))
    frac_val = val / 10d0**exp_val
    exp_str = integer_to_string(exp_val)
    frac_str = real_to_string(frac_val)

    exp_len = len_trim(exp_str)
    frac_len = len_trim(frac_str)
    use_frac_len = max_len - 1 - exp_len
    if (use_frac_len > frac_len) then
       use_frac_len = frac_len
    end if
    if (val < 0d0) then
       min_frac_len = 2
    else
       min_frac_len = 1
    end if
    if (use_frac_len < min_frac_len) then
       do i = 1,max_len
          ret_val(i:i) = "*"
       end do
    else
       ret_val = frac_str(1:use_frac_len) // "e" // trim(exp_str)
    end if
    real_to_string_max_len = adjustl(ret_val)

  end function real_to_string_max_len

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Convert a time to a string format of maximum length.
  character(len=PMC_UTIL_CONVERT_STRING_LEN) &
       function time_to_string_max_len(time, max_len)

    !> Time to convert (s).
    real(kind=dp), intent(in) :: time
    !> Maximum length of resulting string.
    integer, intent(in) :: max_len

    integer, dimension(4), parameter :: scale  = (/   1,  60,  60,  24 /)
    character, dimension(4), parameter :: unit = (/ "s", "m", "h", "d" /)

    character(len=PMC_UTIL_CONVERT_STRING_LEN) :: ret_val
    integer :: i
    logical :: len_ok
    real(kind=dp) :: scaled_time

    scaled_time = time
    len_ok = .false.
    do i = 1,4
       scaled_time = scaled_time / real(scale(i), kind=dp)
       ret_val = trim(integer_to_string(nint(scaled_time))) // unit(i)
       if (len_trim(ret_val) <= max_len) then
          len_ok = .true.
          exit
       end if
    end do
    if (.not. len_ok) then
       ret_val = ""
       do i = 1,max_len
          ret_val(i:i) = "*"
       end do
    end if
    time_to_string_max_len = adjustl(ret_val)

  end function time_to_string_max_len

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Convert a real-valued vector into an integer-valued vector.
  !!
  !! Use n_samp samples. Returns discrete vector whose relative entry
  !! sizes are \f$ \ell_1 \f$ closest to the continuous vector.
  subroutine vec_cts_to_disc(n, vec_cts, n_samp, vec_disc)

    !> Number of entries in vectors.
    integer, intent(in) :: n
    !> Continuous vector.
    real(kind=dp), intent(in) :: vec_cts(n)
    !> Number of discrete samples to use.
    integer, intent(in) :: n_samp
    !> Discretized vector.
    integer, intent(out) :: vec_disc(n)

    integer :: k(1)
    real(kind=dp) :: vec_tot

    vec_tot = sum(vec_cts)

    ! assign a best guess for each bin independently
    vec_disc = nint(vec_cts / vec_tot * real(n_samp, kind=dp))

    ! if we have too few particles then add more
    do while (sum(vec_disc) < n_samp)
       k = minloc(abs(real(vec_disc + 1, kind=dp) - vec_cts) &
            - abs(real(vec_disc, kind=dp) - vec_cts))
       vec_disc(k) = vec_disc(k) + 1
    end do

    ! if we have too many particles then remove some
    do while (sum(vec_disc) > n_samp)
       k = minloc(abs(real(vec_disc - 1, kind=dp) - vec_cts) &
            - abs(real(vec_disc, kind=dp) - vec_cts))
       vec_disc(k) = vec_disc(k) - 1
    end do

    call assert_msg(323412496, sum(vec_disc) == n_samp, &
         'generated incorrect number of samples')

  end subroutine vec_cts_to_disc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Computes the average of an array of integer numbers.
  subroutine average_integer(int_vec, int_avg)

    !> Array of integer numbers.
    integer, intent(in) :: int_vec(:)
    !> Average of int_vec.
    integer, intent(out) :: int_avg

    int_avg = sum(int_vec) / size(int_vec)

  end subroutine average_integer

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Computes the average of an array of real numbers.
  subroutine average_real(real_vec, real_avg)

    !> Array of real numbers.
    real(kind=dp), intent(in) :: real_vec(:)
    !> Average of real_vec.
    real(kind=dp), intent(out) :: real_avg

    real_avg = sum(real_vec) / real(size(real_vec), kind=dp)

  end subroutine average_real

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Return the index of the first occurance of the given value in the
  !> array, or 0 if it is not present.
  integer function string_array_find(array, val)

    !> Array of values.
    character(len=*), intent(in) :: array(:)
    !> Value to find.
    character(len=*), intent(in) :: val

    do string_array_find = 1,size(array)
       if (trim(array(string_array_find)) == trim(val)) return
    end do
    string_array_find = 0

  end function string_array_find

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Allocate or reallocate the given array to ensure it is of the
  !> given size, preserving any data and/or initializing to 0.
  subroutine ensure_real_array_size(x, n, only_grow)

    !> Array of real numbers.
    real(kind=dp), intent(inout), allocatable :: x(:)
    !> Desired size of array.
    integer, intent(in) :: n
    !> Whether to only increase the array size (default .true.).
    logical, intent(in), optional :: only_grow

    integer :: new_n
    real(kind=dp), allocatable :: tmp_x(:)

    if (allocated(x)) then
       new_n = n
       if (present(only_grow)) then
          new_n = max(new_n, size(x))
       end if
       if (size(x) /= new_n) then
          allocate(tmp_x(new_n))
          tmp_x = 0d0
          tmp_x(1:min(new_n, size(x))) = x(1:min(new_n, size(x)))
          call move_alloc(tmp_x, x)
       end if
    else
       allocate(x(n))
       x = 0d0
    end if

  end subroutine ensure_real_array_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Allocate or reallocate the given array to ensure it is of the
  !> given size, preserving any data and/or initializing to 0.
  subroutine ensure_real_array_2d_size(x, n1, n2, only_grow)

    !> Array of real numbers.
    real(kind=dp), intent(inout), allocatable :: x(:, :)
    !> Desired first size of array.
    integer, intent(in) :: n1
    !> Desired second size of array.
    integer, intent(in) :: n2
    !> Whether to only increase the array size (default .true.).
    logical, intent(in), optional :: only_grow

    integer :: new_n1, new_n2, n1_min, n2_min
    real(kind=dp), allocatable :: tmp_x(:, :)

    if (allocated(x)) then
       new_n1 = n1
       new_n2 = n2
       if (present(only_grow)) then
          new_n1 = max(new_n1, size(x, 1))
          new_n2 = max(new_n2, size(x, 2))
       end if
       if ((size(x, 1) /= new_n1) .or. (size(x, 2) /= new_n2)) then
          allocate(tmp_x(new_n1, new_n2))
          n1_min = min(new_n1, size(x, 1))
          n2_min = min(new_n2, size(x, 2))
          tmp_x = 0d0
          tmp_x(1:n1_min, 1:n2_min) = x(1:n1_min, 1:n2_min)
          call move_alloc(tmp_x, x)
       end if
    else
       allocate(x(n1, n2))
       x = 0d0
    end if

  end subroutine ensure_real_array_2d_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Allocate or reallocate the given array to ensure it is of the
  !> given size, preserving any data and/or initializing to 0.
  subroutine ensure_integer_array_size(x, n, only_grow)

    !> Array of integer numbers.
    integer, intent(inout), allocatable :: x(:)
    !> Desired size of array.
    integer, intent(in) :: n
    !> Whether to only increase the array size (default .true.).
    logical, intent(in), optional :: only_grow

    integer :: new_n
    integer, allocatable :: tmp_x(:)

    if (allocated(x)) then
       new_n = n
       if (present(only_grow)) then
          new_n = max(new_n, size(x))
       end if
       if (size(x) /= new_n) then
          allocate(tmp_x(new_n))
          tmp_x = 0
          tmp_x(1:min(new_n, size(x))) = x(1:min(new_n, size(x)))
          call move_alloc(tmp_x, x)
       end if
    else
       allocate(x(n))
       x = 0
    end if

  end subroutine ensure_integer_array_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Allocate or reallocate the given array to ensure it is of the
  !> given size, preserving any data and/or initializing to 0.
  subroutine ensure_integer_array_2d_size(x, n1, n2, only_grow)

    !> Array of integer numbers.
    integer, intent(inout), allocatable :: x(:, :)
    !> Desired first size of array.
    integer, intent(in) :: n1
    !> Desired second size of array.
    integer, intent(in) :: n2
    !> Whether to only increase the array size (default .true.).
    logical, intent(in), optional :: only_grow

    integer :: new_n1, new_n2, n1_min, n2_min
    integer, allocatable :: tmp_x(:, :)

    if (allocated(x)) then
       new_n1 = n1
       new_n2 = n2
       if (present(only_grow)) then
          new_n1 = max(new_n1, size(x, 1))
          new_n2 = max(new_n2, size(x, 2))
       end if
       if ((size(x, 1) /= new_n1) .or. (size(x, 2) /= new_n2)) then
          allocate(tmp_x(new_n1, new_n2))
          n1_min = min(new_n1, size(x, 1))
          n2_min = min(new_n2, size(x, 2))
          tmp_x = 0
          tmp_x(1:n1_min, 1:n2_min) = x(1:n1_min, 1:n2_min)
          call move_alloc(tmp_x, x)
       end if
    else
       allocate(x(n1, n2))
       x = 0
    end if

  end subroutine ensure_integer_array_2d_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Allocate or reallocate the given array to ensure it is of the
  !> given size, without preserving data.
  subroutine ensure_string_array_size(x, n)

    !> Array of strings numbers.
    character(len=*), intent(inout), allocatable :: x(:)
    !> Desired size of array.
    integer, intent(in) :: n

    if (allocated(x)) then
       if (size(x) /= n) then
          deallocate(x)
          allocate(x(n))
       end if
    else
       allocate(x(n))
    end if

  end subroutine ensure_string_array_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Strip the extension to find the basename.
  !!
  !! The filename is assumed to be of the form basename.extension,
  !! where extension contains no periods. If there is no period in
  !! filename then basename is the whole filename.
  subroutine get_basename(filename, basename)

    !> Full filename.
    character(len=*), intent(in) :: filename
    !> Basename.
    character(len=*), intent(out) :: basename

    integer :: i
    logical :: found_period

    basename = filename
    i = len_trim(basename)
    found_period = .false.
    do while ((i > 0) .and. (.not. found_period))
       ! Fortran .and. does not short-circuit, so we can't do the
       ! obvious do while ((i > 0) .and. (basename(i:i) /= ".")),
       ! instead we have to use this hack with the found_period
       ! logical variable.
       if (basename(i:i) == ".") then
          found_period = .true.
       else
          i = i - 1
       end if
    end do
    if (i > 0) then
       basename(i:) = ""
    end if

  end subroutine get_basename

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Current date and time in ISO 8601 format.
  subroutine iso8601_date_and_time(date_time)

    !> Result string.
    character(len=*), intent(out) :: date_time

    character(len=10) :: date, time, zone

    call assert_msg(893219839, len(date_time) >= 29, &
         "date_time string must have length at least 29")
    call date_and_time(date, time, zone)
    date_time = ""
    write(date_time, '(14a)') date(1:4), "-", date(5:6), "-", &
         date(7:8), "T", time(1:2), ":", time(3:4), ":", &
         time(5:10), zone(1:3), ":", zone(4:5)

  end subroutine iso8601_date_and_time

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Convert degrees to radians.
  real(kind=dp) function deg2rad(deg)

    !> Input degrees.
    real(kind=dp), intent(in) :: deg

    deg2rad = deg / 180d0 * const%pi

  end function deg2rad

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Convert radians to degrees.
  real(kind=dp) function rad2deg(rad)

    !> Input radians.
    real(kind=dp), intent(in) :: rad

    rad2deg = rad / const%pi * 180d0

  end function rad2deg

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Sort the given data array and return the permutation defining the sort.
  !!
  !! If \c data is the original data array, \c new_data is the sorted
  !! value of data, and \c perm is the returned permutation, then
  !! <tt>new_data(i) = data(perm(i))</tt>.
  subroutine integer_sort(data, perm)

    !> Data array to sort, sorted on exit.
    integer, intent(inout) :: data(:)
    !> Permutation defining the sort: <tt>new_data(i) = data(perm(i))</tt>.
    integer, intent(out) :: perm(size(data))

#ifdef PMC_USE_C_SORT
    integer(kind=c_int) :: n_c
    integer(kind=c_int), target :: data_c(size(data))
    integer(kind=c_int), target :: perm_c(size(data))
    type(c_ptr) :: data_ptr, perm_ptr

#ifndef DOXYGEN_SKIP_DOC
    interface
       subroutine integer_sort_c(n_c, data_ptr, perm_ptr) bind(c)
         use iso_c_binding
         integer(kind=c_int), value :: n_c
         type(c_ptr), value :: data_ptr, perm_ptr
       end subroutine integer_sort_c
    end interface
#endif

    data_c = int(data, kind=c_int)
    perm_c = 0_c_int
    n_c = int(size(data), kind=c_int)
    data_ptr = c_loc(data_c)
    perm_ptr = c_loc(perm_c)
    call integer_sort_c(n_c, data_ptr, perm_ptr)
    data = int(data_c)
    perm = int(perm_c)
#endif

  end subroutine integer_sort

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Load a real array from a text file.
  subroutine loadtxt(filename, data)

    !> Filename to read from.
    character(len=*), intent(in) :: filename
    !> Array to store data in.
    real(kind=dp), intent(inout), allocatable :: data(:,:)

    integer :: unit, row, col
    logical :: eol, eof
    character(len=1000) :: word

    deallocate(data)
    allocate(data(1,0))
    call open_file_read(filename, unit)

    eof = .false.
    row = 1
    col = 1
    do while (.not. eof)
       call read_word_raw(unit, word, eol, eof)
       if (len_trim(word) > 0) then
          if (row == 1) then
             if (col > size(data, 2)) then
                call reallocate_real_array2d(data, 1, 2 * col)
             end if
          else
             if (col > size(data, 2)) then
                call assert_msg(516120334, col <= size(data, 2), &
                     trim(filename) // ": line " &
                     // trim(integer_to_string(row)) // " too long")
             end if
          end if
          if (row > size(data, 1)) then
             call reallocate_real_array2d(data, 2 * row, size(data, 2))
          end if
          data(row, col) = string_to_real(word)
          col = col + 1
       end if
       if (eol .or. eof) then
          if (row == 1) then
             call reallocate_real_array2d(data, 1, col - 1)
          else
             call assert_msg(266924956, &
                  (col == 1) .or. (col == size(data, 2) + 1), &
                  trim(filename) // ": line " &
                  // trim(integer_to_string(row)) // " too short")
          end if
       end if
       if (eol) then
          row = row + 1
          col = 1
       end if
    end do

    if (col == 1) then
       call reallocate_real_array2d(data, row - 1, size(data, 2))
    end if
    call close_file(unit)

  end subroutine loadtxt

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Write a real 1D array to a text file.
  subroutine savetxt_1d(filename, data)

    !> Filename to write to.
    character(len=*), intent(in) :: filename
    !> Array of data to write.
    real(kind=dp), intent(in) :: data(:)

    integer :: unit, i

    call open_file_write(filename, unit)
    do i = 1,size(data)
       write(unit, '(e30.15e3)') data(i)
    end do
    call close_file(unit)

  end subroutine savetxt_1d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Write a real 2D array to a text file.
  subroutine savetxt_2d(filename, data)

    !> Filename to write to.
    character(len=*), intent(in) :: filename
    !> Array of data to write.
    real(kind=dp), intent(in) :: data(:,:)

    integer :: unit, i, j

    call open_file_write(filename, unit)
    do i = 1,size(data, 1)
       do j = 1,size(data, 2)
          write(unit, '(e30.15e3)', advance='no') data(i, j)
       end do
       write(unit, *) ''
    end do
    call close_file(unit)

  end subroutine savetxt_2d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Reallocate a 2D real array to the given size, preserving the contents.
  subroutine reallocate_real_array2d(data, rows, cols)

    !> Array to reallocate.
    real(kind=dp), allocatable, intent(inout) :: data(:,:)
    !> New leading dimension.
    integer, intent(in) :: rows
    !> New trailing dimension.
    integer, intent(in) :: cols

    real(kind=dp) :: tmp_data(rows, cols)
    integer :: data_rows, data_cols

    data_rows = min(rows, size(data, 1))
    data_cols = min(cols, size(data, 2))
    tmp_data(1:data_rows, 1:data_cols) = data(1:data_rows, 1:data_cols)
    deallocate(data)
    allocate(data(rows, cols))
    data(1:data_rows, 1:data_cols) = tmp_data(1:data_rows, 1:data_cols)

  end subroutine reallocate_real_array2d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read a single character from a file, signaling if we have hit
  !> end-of-line (EOL) or end-of-file (EOF). If EOL or EOF are true
  !> then the character value should be ignored.
  !!
  !! Testing with gfortran 4.5.3 and different length files shows:
  !!
  !! Empty file (total file length 0):
  !! * eol = .false., eof = .true.
  !! * subsequent calls error
  !!
  !! File containing a single 'A' character (total file length 1):
  !! * char = 'A', eol = .false., eof = .false.
  !! * eol = .false., eof = .true.
  !! * subsequent calls error
  !!
  !! File containing a single newline '\n' (total file length 1):
  !! * eol = .true., eof = .false.
  !! * eol = .false., eof = .true.
  !! * subsequent calls error
  !!
  !! File containing a character and newline 'A\n' (total file length 2):
  !! * char = 'A', eol = .false., eof = .false.
  !! * eol = .true., eof = .false.
  !! * eol = .false., eof = .true.
  !! * subsequent calls error
  subroutine read_char_raw(unit, char, eol, eof)

    !> Unit number to read from.
    integer, intent(in) :: unit
    !> Character read.
    character, intent(out) :: char
    !> True if at EOL (end of line).
    logical, intent(out) :: eol
    !> True if at EOF (end of file).
    logical, intent(out) :: eof

    integer :: ios, n_read
    character(len=1) :: read_char

    eol = .false.
    eof = .false.
    char = " " ! shut up uninitialized variable warnings
    read_char = "" ! needed for pgf95 for reading blank lines
    read(unit=unit, fmt='(a)', advance='no', end=100, eor=110, &
         iostat=ios) read_char
    if (ios /= 0) then
       write(0,*) 'ERROR: reading file: IOSTAT = ', ios
       stop 2
    end if
    ! only reach here if we didn't hit end-of-record (end-of-line) in
    ! the above read
    char = read_char
    goto 120

100 eof = .true. ! goto here if end-of-file was encountered immediately
    goto 120

110 eol = .true. ! goto here if end-of-record, meaning end-of-line

120 return

  end subroutine read_char_raw

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read a white-space delimited word from a file, signaling if we
  !> have EOL or EOF. If EOL or EOF are true then the word will still
  !> be meaningful data. If there was no data to be read then
  !> len(word) will be 0.
  subroutine read_word_raw(unit, word, eol, eof)

    !> Unit number to read from.
    integer, intent(in) :: unit
    !> Word read.
    character(len=*), intent(out) :: word
    !> True if at EOL (end of line).
    logical, intent(out) :: eol
    !> True if at EOF (end of file).
    logical, intent(out) :: eof

    integer :: i
    character :: char

    word = ""

    ! skip over spaces
    call read_char_raw(unit, char, eol, eof)
    do while (((ichar(char) == 9) .or. (ichar(char) == 32)) &
         .and. (.not. eol) .and. (.not. eof))
       call read_char_raw(unit, char, eol, eof)
    end do
    if (eol .or. eof) return

    ! char is now the first word character
    i = 1
    word(i:i) = char
    call read_char_raw(unit, char, eol, eof)
    do while ((ichar(char) /= 9) .and. (ichar(char) /= 32) &
         .and. (.not. eol) .and. (.not. eof) .and. (i < len(word)))
       i = i + 1
       word(i:i) = char
       if (i < len(word)) then
          call read_char_raw(unit, char, eol, eof)
       end if
    end do

  end subroutine read_word_raw

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Checks whether a string starts with a given other string.
  !!
  !! <tt>starts_with(A, B)</tt> returns \c true if string \c A starts
  !! with string \c B.
  logical function starts_with(string, start_string)

    !> String to test.
    character(len=*), intent(in) :: string
    !> Starting string.
    character(len=*), intent(in) :: start_string

    if (len(string) < len(start_string)) then
       starts_with = .false.
       return
    end if
    if (string(1:len(start_string)) == start_string) then
       starts_with = .true.
    else
       starts_with = .false.
    end if

  end function starts_with

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Compute the harmonic mean of two numbers.
  elemental real(kind=dp) function harmonic_mean(x1, x2)

    !> First number to average.
    real(kind=dp), intent(in) :: x1
    !> Second number to average.
    real(kind=dp), intent(in) :: x2

    harmonic_mean = 1d0 / (1d0 / x1 + 1d0 / x2)

  end function harmonic_mean

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Compute \f$ - p \ln p\f$ for computing entropy.
  elemental real(kind=dp) function nplogp(p)

    !> Probability \f$p\f$.
    real(kind=dp), intent(in) :: p

    if (p <= 0d0) then
       nplogp = 0d0
    else
       nplogp = - p * log(p)
    end if

  end function nplogp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Compute the entropy of a probability mass function (non
  !> necessarily normalized).
  real(kind=dp) function entropy(p)

    !> Probability mass function, will be normalized before use.
    real(kind=dp), intent(in) :: p(:)

    entropy = sum(nplogp(p / sum(p)))

  end function entropy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Return the least power-of-2 that is at least equal to n.
  integer function pow2_above(n)

    !> Lower bound on the power of 2.
    integer, intent(in) :: n

    if (n <= 0) then
       pow2_above = 0
       return
    end if

    ! LEADZ is in Fortran 2008
    pow2_above = ibset(0, bit_size(n) - leadz(n - 1))

  end function pow2_above

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_util
