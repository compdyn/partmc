! Copyright (C) 2009-2012 Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The extract_sectional_aero_size program.

!> Read NetCDF sectional output files and write out the aerosol size
!> distribution in text format.
program extract_sectional_aero_size

  use pmc_aero_binned
  use pmc_output
  use pmc_mpi
  use getopt_m

  integer, parameter :: DIST_TYPE_NONE = 0
  integer, parameter :: DIST_TYPE_NUM = 1
  integer, parameter :: DIST_TYPE_MASS = 2

  character(len=PMC_MAX_FILENAME_LEN) :: in_prefix, out_filename
  character(len=PMC_MAX_FILENAME_LEN), allocatable :: filename_list(:)
  type(bin_grid_t) :: bin_grid
  type(aero_data_t) :: aero_data
  type(aero_binned_t) :: aero_binned
  integer :: index, out_unit
  integer :: i_file, n_file, i_bin, dist_type
  real(kind=dp) :: time, del_t
  character(len=PMC_UUID_LEN) :: uuid, run_uuid
  real(kind=dp), allocatable :: aero_dist(:,:)
  type(option_s) :: opts(4)

  call pmc_mpi_init()

  opts(1) = option_s("help", .false., 'h')
  opts(2) = option_s("num", .false., 'n')
  opts(3) = option_s("mass", .false., 'm')
  opts(4) = option_s("output", .true., 'o')

  dist_type = DIST_TYPE_NONE
  out_filename = ""

  do
     select case(getopt("hnmo:", opts))
     case(char(0))
        exit
     case('h')
        call print_help()
        stop
     case('n')
        if (dist_type /= DIST_TYPE_NONE) then
           call print_help()
           call die_msg(413086480, 'multiple distribution types selected')
        end if
        dist_type = DIST_TYPE_NUM
     case('m')
        if (dist_type /= DIST_TYPE_NONE) then
           call print_help()
           call die_msg(528866290, 'multiple distribution types selected')
        end if
        dist_type = DIST_TYPE_MASS
     case('o')
        out_filename = optarg
     case( '?' )
        call print_help()
        call die_msg(546118086, 'unknown option: ' // trim(optopt))
     case default
        call print_help()
        call die_msg(720731240, 'unhandled option: ' // trim(optopt))
     end select
  end do

  if (optind /= command_argument_count()) then
     call print_help()
     call die_msg(699147496, 'expected exactly one non-option prefix argument')
  end if

  call get_command_argument(optind, in_prefix)

  if (dist_type == DIST_TYPE_NONE) then
     call print_help()
     call die_msg(576941805, 'must select distribution type (--num or --mass)')
  end if

  if (out_filename == "") then
     if (dist_type == DIST_TYPE_NUM) then
        out_filename = trim(in_prefix) // "_aero_size_num.txt"
     elseif (dist_type == DIST_TYPE_MASS) then
        out_filename = trim(in_prefix) // "_aero_size_mass.txt"
     else
        call die(767619107)
     end if
  end if

  allocate(filename_list(0))
  call input_filename_list(in_prefix, filename_list)
  n_file = size(filename_list)
  call assert_msg(792400289, n_file > 0, &
       "no NetCDF files found with prefix: " // trim(in_prefix))

  call input_sectional(filename_list(1), index, time, del_t, uuid, &
       bin_grid=bin_grid, aero_data=aero_data, aero_binned=aero_binned)
  run_uuid = uuid

  allocate(aero_dist(bin_grid_size(bin_grid), n_file))

  do i_file = 1,n_file
     call input_sectional(filename_list(i_file), index, time, del_t, uuid, &
          bin_grid=bin_grid, aero_data=aero_data, aero_binned=aero_binned)

     call assert_msg(838088000, uuid == run_uuid, &
          "UUID mismatch between " // trim(filename_list(1)) // " and " &
          // trim(filename_list(i_file)))

     if (dist_type == DIST_TYPE_NUM) then
        aero_dist(:, i_file) = aero_binned%num_conc
     elseif (dist_type == DIST_TYPE_MASS) then
        do i_bin = 1,bin_grid_size(bin_grid)
           aero_dist(i_bin, i_file) = sum(aero_binned%vol_conc(i_bin, :) &
                * aero_data%density)
        end do
     else
        call die(141087960)
     end if
  end do

  write(*,'(a)') "Output file: " // trim(out_filename)
  write(*,'(a)') "  Each row of output is one size bin."
  write(*,'(a)') "  The columns of output are:"
  write(*,'(a)') "    column   1: bin center diameter (m)"
  if (dist_type == DIST_TYPE_NUM) then
     write(*,'(a)') "    column j+1: number concentration at time(j) (#/m^3)"
  elseif (dist_type == DIST_TYPE_MASS) then
     write(*,'(a)') "    column j+1: mass concentration at time(j) (kg/m^3)"
  end if
  write(*,'(a)') "  Diameter bins have logarithmic width:"
  write(*,'(a,e20.10)') "    log_width = ln(diam(i+1)) - ln(diam(i)) =", &
       bin_grid%widths(1)

  call open_file_write(out_filename, out_unit)
  do i_bin = 1,bin_grid_size(bin_grid)
     write(out_unit, '(e30.15e3)', advance='no') &
          rad2diam(bin_grid%centers(i_bin))
     do i_file = 1,n_file
        write(out_unit, '(e30.15e3)', advance='no') aero_dist(i_bin, i_file)
     end do
     write(out_unit, '(a)') ''
  end do
  call close_file(out_unit)

  deallocate(filename_list)
  deallocate(aero_dist)

  call pmc_mpi_finalize()

contains

  subroutine print_help()

    write(*,'(a)') 'Usage: extract_sectional_aero_size [options] ' &
         // '<netcdf_prefix>'
    write(*,'(a)') ''
    write(*,'(a)') 'options are:'
    write(*,'(a)') '  -h, --help        Print this help message.'
    write(*,'(a)') '  -n, --num         Output number distribution.'
    write(*,'(a)') '  -m, --mass        Output mass distribution.'
    write(*,'(a)') '  -o, --out <file>  Output filename.'
    write(*,'(a)') ''
    write(*,'(a)') 'Examples:'
    write(*,'(a)') '  extract_sectional_aero_size --num data_0001'
    write(*,'(a)') ''

  end subroutine print_help

end program extract_sectional_aero_size
