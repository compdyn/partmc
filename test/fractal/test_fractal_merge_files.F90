! Copyright (C) 2009-2012 Matthew West
! Copyright (C) 2012 Jian Tian
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The merge_dimless_time_files program.

!> Merge too text files from base run and restart run into one
!> file for dimensionless time series.
program merge_dimless_time_files

  use pmc_util
  use pmc_mpi
  use getopt_m

  character(len=PMC_MAX_FILENAME_LEN) :: in_prefix, out_filename
  integer :: row, col, out_unit, old_row
  real(kind=dp), allocatable :: data1(:,:), data2(:,:)
  character(len=PMC_MAX_FILENAME_LEN) :: filename1, filename2
  type(option_s) :: opts(2)

  call pmc_mpi_init()

  opts(1) = option_s("help", .false., 'h')
  opts(2) = option_s("output", .true., 'o')

  out_filename = ""

  do
     select case(getopt("ho:", opts))
     case(char(0))
        exit
     case('h')
        call print_help()
        stop
     case('o')
        out_filename = optarg
     end select
  end do

  if (optind /= command_argument_count()) then
     call print_help()
     call die_msg(967032898, 'expected exactly one non-option prefix argument')
  end if

  call get_command_argument(optind, in_prefix)

  if (out_filename == "") then
     out_filename = "out_dimless_t/" // trim(in_prefix) // "_dimless_t_series.txt"
  end if

  filename1 = "out_dimless_t/" // trim(in_prefix) // "_dimless_time.txt"
  filename2 = "out_dimless_t/restart/" // trim(in_prefix) // "_dimless_time.txt"

  allocate (data1(0,0))
  allocate (data2(0,0))
  call loadtxt(filename1, data1)
  call loadtxt(filename2, data2)
  old_row = size(data1, 1)
  call reallocate_real_array2d(data1, size(data1, 1) + size(data2, 1) - 1, size(data1, 2))
  do row = old_row + 1, size(data1, 1)
     do col = 1, size(data1, 2)
        if (col == 1) then
           data1(row, col) = data2(row - old_row + 1, col) + data1(old_row, 1)
        else
           data1(row, col) = data2(row - old_row + 1, col)
        end if
     end do
  end do

  call open_file_write(out_filename, out_unit)
  do row = 1, size(data1, 1)
     write(out_unit, '(e30.15e3)', advance='no') data1(row, 1)
     write(out_unit, '(e30.15e3)', advance='no') data1(row, 2)
     write(out_unit, '(a)') ''
  end do

  call close_file(out_unit)

  call pmc_mpi_finalize()

contains

  subroutine print_help()

    write(*,'(a)') 'Usage: test_fractal_merge_files [options] <netcdf_prefix>'
    write(*,'(a)') ''
    write(*,'(a)') 'options are:'
    write(*,'(a)') '  -h, --help        Print this help message.'
    write(*,'(a)') '  -o, --out <file>  Output filename.'
    write(*,'(a)') ''
    write(*,'(a)') 'Examples:'
    write(*,'(a)') '  test_fractal_merge_files data_0001'
    write(*,'(a)') ''

  end subroutine print_help

end program merge_dimless_time_files
