! Copyright (C) 2009-2012, 2016 Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The numeric_diff program.

!> Compare two files containing numerical arrays and check whether they
!> are the same as each other, to within the specified tolerance.
!>
!> If the arrays in the two files are of different sizes then they are
!> automatically different. Otherwise they are the same if
!>       \f[ | A_1 - A_2 |_2 < \verb+abs_tol+ \f]
!> and
!>       \f[ \frac{| A_1 - A_2 |_2}{| A_1 |_2 + | A_2 |_2}
!>           < \verb+rel_tol+ \f]
!> and are otherwise different. Setting \c abs_tol or \c rel_tol to zero
!> skips the corresponding test.
!>
!> If the files are the same then "<tt>files match within the given
!> tolerances</tt>" is printed on stdout, otherwise "<tt>files
!> differ</tt>" is printed, followed by the absolute and relative
!> differences, as above, or a message describing the difference. The
!> files will be reported as different if they have a different
!> pattern of end-of-lines and end-of-files, or if they have
!> whitespace in different places (amount of whitespace is
!> irrelevant).
!>
!> The exit status of the program is:
!>   \li 0 if the files are the same
!>   \li 1 if the files are different
!>   \li 3 if an error occurred
program numeric_diff

  use pmc_util
  use pmc_mpi
  use getopt_m

  integer, parameter :: BY_ARRAY = 1
  integer, parameter :: BY_ROW = 2
  integer, parameter :: BY_COL = 3
  integer, parameter :: BY_ELEM = 4

  integer, parameter :: NORM_ONE = 1
  integer, parameter :: NORM_TWO = 2
  integer, parameter :: NORM_SUP = 3

  character(len=PMC_MAX_FILENAME_LEN) :: filename1, filename2
  integer :: by, norm, min_row, max_row, min_col, max_col, n_row, n_col
  real(kind=dp) :: abs_tol, rel_tol
  real(kind=dp), allocatable, target, dimension(:,:) :: data1, data2
  real(kind=dp), allocatable, dimension(:,:) :: diff, norm1, abs_err, rel_err
  real(kind=dp), pointer, dimension(:,:)  :: use_data1, use_data2
  type(option_s) :: opts(9)

  call pmc_mpi_init()

  opts(1) = option_s("help", .false., 'h')
  opts(2) = option_s("abs-tol", .true., 't')
  opts(3) = option_s("rel-tol", .true., 'T')
  opts(4) = option_s("min-row", .true., 'r')
  opts(5) = option_s("max-row", .true., 'R')
  opts(6) = option_s("min-col", .true., 'c')
  opts(7) = option_s("max-col", .true., 'C')
  opts(8) = option_s("by", .true., 'b')
  opts(9) = option_s("norm", .true., 'n')

  abs_tol = 0d0
  rel_tol = 0d0
  min_row = 0
  max_row = 0
  min_col = 0
  max_col = 0
  by = BY_ARRAY
  norm = NORM_TWO

  do
     select case(getopt("ht:T:r:R:c:C:pP", opts))
     case(char(0))
        exit
     case('h')
        call print_help()
        stop
     case('t')
        abs_tol = string_to_real(optarg)
     case('T')
        rel_tol = string_to_real(optarg)
     case('r')
        min_row = string_to_integer(optarg)
     case('R')
        max_row = string_to_integer(optarg)
     case('c')
        min_col = string_to_integer(optarg)
     case('C')
        max_col = string_to_integer(optarg)
     case('b')
        select case(optarg)
        case('array')
           by = BY_ARRAY
        case('row')
           by = BY_ROW
        case('col')
           by = BY_COL
        case('elem')
           by = BY_ELEM
        case default
           call die_msg(526174645, "unknown --by argument: " // trim(optarg))
        end select
     case('n')
        select case(optarg)
        case('one')
           norm = NORM_ONE
        case('two')
           norm = NORM_TWO
        case('sup')
           norm = NORM_SUP
        case default
           call die_msg(568020730, "unknown --norm argument: " // trim(optarg))
        end select
     case( '?' )
        call die_msg(141541134, 'unknown option')
     case default
        call die_msg(816884701, 'unhandled option: ' // trim(optopt))
     end select
  end do

  if (optind /= command_argument_count() - 1) then
     call print_help()
     call die_msg(142676480, &
          'expected exactly two non-option prefix arguments')
  end if

  call get_command_argument(optind, filename1)
  call get_command_argument(optind + 1, filename2)

  allocate(data1(0,0))
  allocate(data2(0,0))
  call loadtxt(filename1, data1)
  call loadtxt(filename2, data2)

  if (min_row <= 0) then
     min_row = 1
  end if
  if (max_row <= 0) then
     call assert_msg(266216891, size(data1, 1) == size(data2, 1), &
          "number of rows differs between input files")
     max_row = size(data1, 1)
  else
     call assert_msg(136425118, max_row <= size(data1, 1), &
          "max-row exceeds the number of rows in " // trim(filename1))
     call assert_msg(279083405, max_row <= size(data2, 1), &
          "max-row exceeds the number of rows in " // trim(filename2))
  end if

  if (min_col <= 0) then
     min_col = 1
  end if
  if (max_col <= 0) then
     call assert_msg(148743161, size(data1, 2) == size(data2, 2), &
          "number of columns differs between input files")
     max_col = size(data1, 2)
  else
     call assert_msg(884008689, max_col <= size(data1, 2), &
          "max-col exceeds the number of columns in " // trim(filename1))
     call assert_msg(553561214, max_col <= size(data2, 2), &
          "max-col exceeds the number of columns in " // trim(filename2))
  end if

  use_data1 => data1(min_row:max_row, min_col:max_col)
  use_data2 => data2(min_row:max_row, min_col:max_col)

  n_row = max_row - min_row + 1
  n_col = max_col - min_col + 1
  diff = use_data1 - use_data2

  select case(by)
  case(BY_ARRAY)
     allocate(norm1(1, 1))
     allocate(abs_err(1, 1))
     select case(norm)
     case(NORM_ONE)
        norm1(1, 1) = sum(abs(use_data1))
        abs_err(1, 1) = sum(abs(diff))
     case(NORM_TWO)
        norm1(1, 1) = sqrt(sum(use_data1**2))
        abs_err(1, 1) = sqrt(sum(diff**2))
     case(NORM_SUP)
        norm1(1, 1) = maxval(abs(use_data1))
        abs_err(1, 1) = maxval(abs(diff))
     case default
        call die(644692127)
     end select
  case(BY_ROW)
     allocate(norm1(size(diff, 1), 1))
     allocate(abs_err(size(diff, 1), 1))
     select case(norm)
     case(NORM_ONE)
        norm1(:, 1) = sum(abs(use_data1), 2)
        abs_err(:, 1) = sum(abs(diff), 2)
     case(NORM_TWO)
        norm1(:, 1) = sqrt(sum(use_data1**2, 2))
        abs_err(:, 1) = sqrt(sum(diff**2, 2))
     case(NORM_SUP)
        norm1(:, 1) = maxval(abs(use_data1), 2)
        abs_err(:, 1) = maxval(abs(diff), 2)
     case default
        call die(698913943)
     end select
  case(BY_COL)
     allocate(norm1(1, size(diff, 2)))
     allocate(abs_err(1, size(diff, 2)))
     select case(norm)
     case(NORM_ONE)
        norm1(1, :) = sum(abs(use_data1), 1)
        abs_err(1, :) = sum(abs(diff), 1)
     case(NORM_TWO)
        norm1(1, :) = sqrt(sum(use_data1**2, 1))
        abs_err(1, :) = sqrt(sum(diff**2, 1))
     case(NORM_SUP)
        norm1(1, :) = maxval(abs(use_data1), 1)
        abs_err(1, :) = maxval(abs(diff), 1)
     case default
        call die(351454435)
     end select
  case(BY_ELEM)
     allocate(norm1(size(diff, 1), size(diff, 2)))
     allocate(abs_err(size(diff, 1), size(diff, 2)))
     norm1(:, :) = abs(use_data1)
     abs_err(:, :) = abs(diff)
  case default
     call die(681575403)
  end select

  if (.not. allocated(norm1)) then
     allocate(norm1(0,0)) ! prevent uninitialized warnings
  end if

  allocate(rel_err(size(abs_err, 1), size(abs_err, 2)))
  where (norm1 > 0d0)
     rel_err = abs_err / norm1
  elsewhere
     rel_err = 0d0
  end where

  if ((abs_tol <= 0d0) .and. (rel_tol <= 0d0)) then
     call print_errors(abs_err, rel_err)
  elseif (((abs_tol <= 0d0) .or. all(abs_err < abs_tol)) &
       .and. ((rel_tol <= 0d0) .or. all(rel_err < rel_tol))) then
     write(*,'(a)') 'files match within the given relative tolerance'
  else
     write(*,'(a)') 'files differ'
     call print_errors(abs_err, rel_err)
     stop 1
  end if

  call pmc_mpi_finalize()

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine print_help()

    write(*,'(a)') 'Usage: numeric_diff [options] <reference_file> <test_file>'
    write(*,'(a)') ''
    write(*,'(a)') 'options are:'
    write(*,'(a)') '  -h, --help          Print this help message.'
    write(*,'(a)') '  -t, --abs-tol <N>   Absolute error tolerance.'
    write(*,'(a)') '  -T, --rel-tol <N>   Relative error tolerance.'
    write(*,'(a)') '  -r, --min-row <N>   Minimum row number of data to use.'
    write(*,'(a)') '  -R, --max-row <N>   Maximum row number of data to use.'
    write(*,'(a)') '  -c, --min-col <N>   Minimum column number of data to ' &
         // 'use.'
    write(*,'(a)') '  -C, --max-col <N>   Maximum column number of data to ' &
         // 'use.'
    write(*,'(a)') '  -b, --by <S>        Compute error by <S>. <S> is one ' &
         // 'of "array", "row",'
    write(*,'(a)') '                      "col", or "elem". Default: "array".'
    write(*,'(a)') '  -n, --norm <S>      Compute error with norm <S>. <S> ' &
         // 'is one of "one",'
    write(*,'(a)') '                      "two", or "sup". Default: "two".'
    write(*,'(a)') ''
    write(*,'(a)') 'Examples:'
    write(*,'(a)') '  numeric_diff --rel-tol 1e-3 ref_data.txt test_data.txt'
    write(*,'(a)') '  numeric_diff --by col --rel-tol 1e-6 --min-col 2 ' &
         // 'ref_data.txt test_data.txt'
    write(*,'(a)') ''

  end subroutine print_help

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine print_errors(abs_err, rel_err)

    !> Absolute errors.
    real(kind=dp) :: abs_err(:,:)
    !> Relative errors.
    real(kind=dp) :: rel_err(:,:)

    integer :: i_row, i_col
    character(len=3) :: advance

    call assert(434301862, (size(abs_err, 1) == size(rel_err, 1)) &
         .and. (size(abs_err, 2) == size(rel_err, 2)))

    if ((size(abs_err, 1) == 1) .and. (size(abs_err, 2) <= 5)) then
       advance = 'no'
    else
       advance = 'yes'
    end if

    write(*,'(a)', advance=advance) 'absolute error: '
    do i_row = 1,size(abs_err, 1)
       do i_col = 1,size(abs_err, 2)
          write(*, '(e12.3)', advance='no') abs_err(i_row, i_col)
       end do
       write(*,'(a)') ''
    end do
    write(*,'(a)', advance=advance) 'relative error: '
    do i_row = 1,size(abs_err, 1)
       do i_col = 1,size(abs_err, 2)
          write(*, '(e12.3)', advance='no') rel_err(i_row, i_col)
       end do
       write(*,'(a)') ''
    end do

    if ((size(abs_err, 1) > 1) .or. (size(abs_err, 2) > 1)) then
       write(*, '(a,e12.3)') 'maximum absolute error: ', maxval(abs_err)
       write(*, '(a,e12.3)') 'maximum relative error: ', maxval(rel_err)
    end if

  end subroutine print_errors

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program numeric_diff
