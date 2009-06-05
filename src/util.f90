! Copyright (C) 2005-2009 Nicole Riemer and Matthew West
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

  !> Maximum number of IO units usable simultaneously.
  integer, parameter :: max_units = 200
  !> Minimum unit number to allocate.
  integer, parameter :: unit_offset = 10
  !> Table of unit numbers storing allocation status.
  logical, save :: unit_used(max_units) = .false.

contains
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Errors unless condition_ok is true.
  subroutine assert_msg(code, condition_ok, error_msg)

    !> Status code to use if assertion fails.
    integer, intent(in) :: code
    !> Whether the assertion is ok.
    logical, intent(in) :: condition_ok
    !> Msg if assertion fails.
    character(len=*), intent(in) :: error_msg

    integer :: ierr
    character(len=100) :: code_str

    write(code_str,*) code
    code_str = adjustl(code_str)
    if (.not. condition_ok) then
       write(0,*) 'ERROR (PartMC-', trim(code_str), '): ', trim(error_msg)
#ifdef PMC_USE_MPI
       call mpi_abort(MPI_COMM_WORLD, code, ierr)
#else
       call exit(3)
#endif
    end if

  end subroutine assert_msg

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Errors unless condition_ok is true.
  subroutine assert(code, condition_ok)

    !> Status code to use if assertion fails.
    integer, intent(in) :: code
    !> Whether the assertion is ok.
    logical, intent(in) :: condition_ok

    call assert_msg(code, condition_ok, 'assertion failed')

  end subroutine assert

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Error immediately.
  subroutine die(code)

    !> Status code to use if assertion fails.
    integer, intent(in) :: code

    call assert(code, .false.)

  end subroutine die

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Error immediately.
  subroutine die_msg(code, error_msg)

    !> Status code to use if assertion fails.
    integer, intent(in) :: code
    !> Msg if assertion fails.
    character(len=*), intent(in) :: error_msg

    call assert_msg(code, .false., error_msg)

  end subroutine die_msg

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
       write(0,*) 'ERROR: no more units available - need to free_unit()'
       call exit(1)
    end if
    unit_used(i) = .true.
    get_unit = i + unit_offset

  end function get_unit

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Frees a unit number returned by get_unit().
  subroutine free_unit(unit)

    integer, intent(in) :: unit

    unit_used(unit - unit_offset) = .false.

  end subroutine free_unit

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Allocate a new unit and open it with a filename given by
  !> basename + suffix.
  subroutine open_output(basename, suffix, out_unit)

    !> Basename of the output file.
    character(len=*), intent(in) :: basename
    !> Suffix of the output file.
    character(len=*), intent(in) :: suffix
    !> Unit for the file.
    integer, intent(out) :: out_unit

    character(len=len(basename)+len(suffix)) :: filename
    
    filename = basename
    filename((len_trim(filename)+1):) = suffix
    out_unit = get_unit()
    open(out_unit, file=filename)
    write(*,'(a,a)') 'Writing ', trim(filename)
    
  end subroutine open_output

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Initializes the random number generator to the state defined by
  !> the given seed. If the seed is 0 then a seed is auto-generated
  !> from the current time.
  subroutine pmc_srand(seed)

    !> Random number generator seed.
    integer, intent(in) :: seed

    integer :: i, n, clock
    integer, allocatable :: seed_vec(:)
    ! FIXME: HACK
    real*8 :: r
    ! FIXME: end HACK

    call random_seed(size = n)
    allocate(seed_vec(n))
    if (seed == 0) then
       call system_clock(count = clock)
    else
       clock = seed
    end if
    i = 0 ! HACK to shut up gfortran warning
    seed_vec = clock + 37 * (/ (i - 1, i = 1, n) /)
    call random_seed(put = seed_vec)
    deallocate(seed_vec)
    ! FIXME: HACK for bad rng behavior from pgf90
    do i = 1,1000
       r = pmc_random()
    end do
    ! FIXME: end HACK

  end subroutine pmc_srand

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns a random number between 0 and 1.
  real*8 function pmc_random()

    real*8 rnd

    call random_number(rnd)
    pmc_random = rnd

  end function pmc_random
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns a random integer between 1 and n.
  integer function pmc_rand_int(n)

    !> Maximum random number to generate.
    integer, intent(in) :: n

    pmc_rand_int = mod(int(pmc_random() * dble(n)), n) + 1
    call assert(515838689, pmc_rand_int >= 1)
    call assert(802560153, pmc_rand_int <= n)

  end function pmc_rand_int

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Round val to floor(val) or ceil(val) with probability
  !> proportional to the relative distance from val. That is,
  !> Prob(prob_round(val) == floor(val)) = ceil(val) - val.
  integer function prob_round(val)

    !> Value to round.
    real*8, intent(in) :: val
    
    prob_round = int(val)
    if (pmc_random() .lt. mod(val, 1d0)) then
       prob_round = prob_round + 1
    endif

  end function prob_round

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Convert volume (m^3) to radius (m).
  real*8 function vol2rad(v)

    !> Volume (m^3).
    real*8, intent(in) :: v
    
    vol2rad = (v / (4d0 / 3d0 * const%pi))**(1d0/3d0)
    
  end function vol2rad
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Convert volume (m^3) to diameter (m).
  real*8 function vol2diam(v)

    !> Volume (m^3).
    real*8, intent(in) :: v
    
    vol2diam = 2d0 * (v / (4d0 / 3d0 * const%pi))**(1d0/3d0)
    
  end function vol2diam
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Convert radius (m) to volume (m^3).
  real*8 function rad2vol(r)

    !> Radius (m).
    real*8, intent(in) :: r
    
    rad2vol = 4d0 / 3d0 * const%pi * r**3d0
    
  end function rad2vol
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Convert diameter (m) to volume (m^3).
  real*8 function diam2vol(d)

    !> Diameter (m).
    real*8, intent(in) :: d
    
    diam2vol = 4d0 / 3d0 * const%pi * (d / 2d0)**3d0
    
  end function diam2vol
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Tests whether two real numbers are almost equal using only a
  !> relative tolerance.
  logical function almost_equal(d1, d2)
    
    !> First number to compare.
    real*8, intent(in) :: d1
    !> Second number to compare.
    real*8, intent(in) :: d2
    
    !> Relative tolerance.
    real*8, parameter :: eps = 1d-8
    
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

  !> Tests whether two real numbers are almost equal using an
  !> absolute and relative tolerance.
  logical function almost_equal_abs(d1, d2, abs_tol)
    
    !> First number to compare.
    real*8, intent(in) :: d1
    !> Second number to compare.
    real*8, intent(in) :: d2
    !> Tolerance for when d1 equals d2.
    real*8, intent(in) :: abs_tol
    
    !> Relative tolerance.
    real*8, parameter :: eps = 1d-8
    
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
    real*8, intent(in) :: time
    !> Estimate of the time to the next call.
    real*8, intent(in) :: timestep
    !> How often the event should be done.
    real*8, intent(in) :: interval
    !> When the event was last done.
    real*8, intent(inout) :: last_time
    !> Whether the event should be done.
    logical, intent(out) :: do_event
    
    !> Fuzz for event occurance.
    real*8, parameter :: tolerance = 1d-6
    
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

  !> Opens a file for reading that must already exist, checking for
  !> errors.
  subroutine open_existing(unit, filename)

    !> Unit to open with.
    integer, intent(in) :: unit
    !> Filename of file to open.
    character(len=*), intent(in) :: filename

    integer ios

    open(unit=unit, file=filename, status='old', iostat=ios)
    if (ios /= 0) then
       write(0,*) 'ERROR: unable to open file ', trim(filename), &
            ': IOSTAT = ', ios
       call exit(1)
    end if

  end subroutine open_existing

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Makes a linearly spaced array of length n from min to max.
  subroutine linspace(min_x, max_x, n, x)

    !> Minimum array value.
    real*8, intent(in) :: min_x
    !> Maximum array value.
    real*8, intent(in) :: max_x
    !> Number of entries.
    integer, intent(in) :: n
    !> Array.
    real*8, intent(out) :: x(n)

    integer :: i
    real*8 :: a

    do i = 2, (n - 1)
       a = dble(i - 1) / dble(n - 1)
       x(i) = (1d0 - a) * min_x + a * max_x
    end do
    if (n > 0) then
       ! make sure these values are exact
       x(1) = min_x
       x(n) = max_x
    end if
    
  end subroutine linspace

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Makes a logarithmically spaced array of length n from min to max.
  subroutine logspace(min_x, max_x, n, x)

    !> Minimum array value.
    real*8, intent(in) :: min_x
    !> Maximum array value.
    real*8, intent(in) :: max_x
    !> Number of entries.
    integer, intent(in) :: n
    !> Array.
    real*8, intent(out) :: x(n)

    real*8 :: log_x(n)

    call assert(548290438, min_x > 0d0)
    call assert(805259035, max_x > 0d0)
    call linspace(log(min_x), log(max_x), n, log_x)
    x = exp(log_x)
    if (n > 0) then
       ! make sure these values are exact
       x(1) = min_x
       x(n) = max_x
    end if
    
  end subroutine logspace

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Find the position of a real number in a 1D linear array.
  !!
  !! If xa is the array allocated by linspace(min_x, max_x, n, xa)
  !! then i = linspace_find(min_x, max_x, n, x) returns the index i
  !! satisfying xa(i) <= x < xa(i+1) for min_x <= x < max_x. If
  !! x >= max_x then i = n - 1.  If x < min_x then i = 1. Thus
  !! 1 <= i <= n - 1.
  !!
  !! This is equivalent to using find_1d() but much faster if the
  !! array is linear.
  integer function linspace_find(min_x, max_x, n, x)

    !> Minimum array value.
    real*8, intent(in) :: min_x
    !> Maximum array value.
    real*8, intent(in) :: max_x
    !> Number of entries.
    integer, intent(in) :: n
    !> Value.
    real*8, intent(in) :: x

    linspace_find = floor((x - min_x) / (max_x - min_x) * dble(n - 1)) + 1
    linspace_find = min(linspace_find, n - 1)
    linspace_find = max(linspace_find, 1)
    
  end function linspace_find

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Find the position of a real number in a 1D logarithmic array.
  !!
  !! If xa is the array allocated by logspace(min_x, max_x, n, xa) then
  !! i = logspace_find(min_x, max_x, n, x) returns the index i satisfying
  !! xa(i) <= x < xa(i+1) for min_x <= x < max_x. If x >= max_x then i = n.
  !! If x < min_x then i = 1. Thus 1 <= i <= n.
  !!
  !! This is equivalent to using find_1d() but much faster if the
  !! array is known to be logarithmic.
  integer function logspace_find(min_x, max_x, n, x)

    !> Minimum array value.
    real*8, intent(in) :: min_x
    !> Maximum array value.
    real*8, intent(in) :: max_x
    !> Number of entries.
    integer, intent(in) :: n
    !> Value.
    real*8, intent(in) :: x

    logspace_find = linspace_find(log(min_x), log(max_x), n, log(x))
    
  end function logspace_find

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
    real*8, intent(in) :: x_vals(n)
    !> Value to interpolate at.
    real*8, intent(in) :: x

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

  !> 1D linear interpolation.
  !!
  !! Takes an array of x and y, and a single x value, and returns the
  !! corresponding y using linear interpolation. x_vals must be
  !! sorted.
  real*8 function interp_1d(x_vals, y_vals, x)

    !> X value array, must be sorted.
    real*8, intent(in) :: x_vals(:)
    !> Y value array.
    real*8, intent(in) :: y_vals(size(x_vals))
    !> Value to interpolate at.
    real*8, intent(in) :: x

    integer :: n, p
    real*8 :: y, alpha

    n = size(x_vals)
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

  !> Convert a string to an integer.
  integer function string_to_integer(string)

    !> String to convert.
    character(len=*), intent(in) :: string
    
    integer :: val
    integer :: ios

    read(string, '(i20)', iostat=ios) val
    if (ios /= 0) then
       write(0,'(a,a,a,i3)') 'Error converting ', trim(string), &
            ' to integer: IOSTAT = ', ios
       call exit(1)
    end if
    string_to_integer = val

  end function string_to_integer

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Convert a string to a real.
  real*8 function string_to_real(string)

    !> String to convert.
    character(len=*), intent(in) :: string
    
    real*8 :: val
    integer :: ios

    read(string, '(f20.0)', iostat=ios) val
    if (ios /= 0) then
       write(0,'(a,a,a,i3)') 'Error converting ', trim(string), &
            ' to real: IOSTAT = ', ios
       call exit(1)
    end if
    string_to_real = val

  end function string_to_real

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Convert a string to a logical.
  logical function string_to_logical(string)

    !> String to convert.
    character(len=*), intent(in) :: string
    
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
       write(0,'(a,a,a)') 'Error converting ', trim(string), ' to logical'
       call exit(1)
    end if
    string_to_logical = val

  end function string_to_logical

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Sample the given continuous probability density function.
  !!
  !! That is, return a number k = 1,...,n such that prob(k) = pdf(k) /
  !! sum(pdf).
  integer function sample_cts_pdf(n, pdf)

    !> Number of entries.
    integer, intent(in) :: n
    !> Probability density function (not normalized).
    real*8, intent(in) :: pdf(n)

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
       k = pmc_rand_int(n)
       if (pmc_random() < pdf(k) / pdf_max) then
          found = .true.
       end if
    end do
    sample_cts_pdf = k

  end function sample_cts_pdf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Sample the given discrete probability density function.
  !!
  !! That is, return a number k = 1,...,n such that prob(k) = pdf(k) /
  !! sum(pdf).
  integer function sample_disc_pdf(n, pdf)

    !> Number of entries.
    integer, intent(in) :: n
    !> Probability density function.
    integer, intent(in) :: pdf(n)

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
       k = pmc_rand_int(n)
       if (pmc_random() < dble(pdf(k)) / dble(pdf_max)) then
          found = .true.
       end if
    end do
    sample_disc_pdf = k

  end function sample_disc_pdf
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Convert a real-valued vector into an integer-valued vector by
  !> sampling.
  !!
  !! Use n_samp samples. Each discrete entry is sampled with a PDF
  !! given by vec_cts. This is very slow for large n_samp or large n.
  subroutine sample_vec_cts_to_disc(n, vec_cts, n_samp, vec_disc)
    
    !> Number of entries in vector.
    integer, intent(in) :: n
    !> Continuous vector.
    real*8, intent(in) :: vec_cts(n)
    !> Number of discrete samples to use.
    integer, intent(in) :: n_samp
    !> Discretized vector.
    integer, intent(out) :: vec_disc(n)

    integer :: i_samp, k

    vec_disc = 0
    do i_samp = 1,n_samp
       k = sample_cts_pdf(n, vec_cts)
       vec_disc(k) = vec_disc(k) + 1
    end do
    
  end subroutine sample_vec_cts_to_disc
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Convert a real-valued vector into an integer-valued vector.
  !!
  !! Use n_samp samples. Returns discrete vector whose relative entry
  !! sizes are \f$ \ell_1 \f$ closest to the continuous vector.
  subroutine vec_cts_to_disc(n, vec_cts, n_samp, vec_disc)
    
    !> Number of entries in vectors.
    integer, intent(in) :: n
    !> Continuous vector.
    real*8, intent(in) :: vec_cts(n)
    !> Number of discrete samples to use.
    integer, intent(in) :: n_samp
    !> Discretized vector.
    integer, intent(out) :: vec_disc(n)
    
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

  !> Computes the average of an array of integer numbers.
  subroutine average_integer(int_vec, int_avg)
    
    !> Array of integer numbers.
    integer, intent(in) :: int_vec(:)
    !> Average of int_vec.
    integer, intent(out) :: int_avg
    
    int_avg = sum(int_vec) / size(int_vec)
    
  end subroutine average_integer
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Computes the average of an array of real numbers.
  subroutine average_real(real_vec, real_avg)
    
    !> Array of real numbers.
    real*8, intent(in) :: real_vec(:)
    !> Average of real_vec.
    real*8, intent(out) :: real_avg
    
    real_avg = sum(real_vec) / dble(size(real_vec))
    
  end subroutine average_real
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
end module pmc_util
