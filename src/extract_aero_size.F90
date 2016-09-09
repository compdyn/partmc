! Copyright (C) 2009-2012 Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The extract_aero_size program.

!> Read NetCDF output files and write out the aerosol number or mass
!> size distributions in text format.
program extract_aero_size

  use pmc_aero_state
  use pmc_aero_particle
  use pmc_output
  use pmc_mpi
  use getopt_m

  integer, parameter :: DIST_TYPE_NONE = 0
  integer, parameter :: DIST_TYPE_NUM = 1
  integer, parameter :: DIST_TYPE_MASS = 2

  character(len=PMC_MAX_FILENAME_LEN) :: in_prefix, out_filename
  character(len=PMC_MAX_FILENAME_LEN), allocatable :: filename_list(:)
  character(len=1000) :: tmp_str
  type(bin_grid_t) :: diam_grid
  type(aero_data_t) :: aero_data
  type(aero_state_t) :: aero_state
  integer :: index, i_repeat, i_part, i_spec, out_unit
  integer :: i_file, n_file, i_bin, n_bin, dist_type
  real(kind=dp) :: time, del_t
  real(kind=dp) :: d_min, d_max, diam, volume
  character(len=PMC_UUID_LEN) :: uuid, run_uuid
  real(kind=dp), allocatable :: diameters(:), num_concs(:), masses(:), hist(:)
  real(kind=dp), allocatable :: aero_dist(:,:)
  type(option_s) :: opts(7)

  call pmc_mpi_init()

  opts(1) = option_s("help", .false., 'h')
  opts(2) = option_s("num", .false., 'n')
  opts(3) = option_s("mass", .false., 'm')
  opts(4) = option_s("dmin", .true., 'N')
  opts(5) = option_s("dmax", .true., 'X')
  opts(6) = option_s("nbin", .true., 'b')
  opts(7) = option_s("output", .true., 'o')

  dist_type = DIST_TYPE_NONE
  d_min = 1d-10
  d_max = 1d-3
  n_bin = 100
  out_filename = ""

  do
     select case(getopt("hnmN:X:b:o:", opts))
     case(char(0))
        exit
     case('h')
        call print_help()
        stop
     case('n')
        if (dist_type /= DIST_TYPE_NONE) then
           call print_help()
           call die_msg(525113814, 'multiple distribution types selected')
        end if
        dist_type = DIST_TYPE_NUM
     case('m')
        if (dist_type /= DIST_TYPE_NONE) then
           call print_help()
           call die_msg(155494931, 'multiple distribution types selected')
        end if
        dist_type = DIST_TYPE_MASS
     case('N')
        d_min = string_to_real(optarg)
     case('X')
        d_max = string_to_real(optarg)
     case('b')
        n_bin = string_to_integer(optarg)
     case('o')
        out_filename = optarg
     case( '?' )
        call print_help()
        call die_msg(956456220, 'unknown option: ' // trim(optopt))
     case default
        call print_help()
        call die_msg(203991511, 'unhandled option: ' // trim(optopt))
     end select
  end do

  if (optind /= command_argument_count()) then
     call print_help()
     call die_msg(533171694, 'expected exactly one non-option prefix argument')
  end if

  call get_command_argument(optind, in_prefix)

  if (dist_type == DIST_TYPE_NONE) then
     call print_help()
     call die_msg(540839314, 'must select distribution type (--num or --mass)')
  end if

  if (out_filename == "") then
     if (dist_type == DIST_TYPE_NUM) then
        out_filename = trim(in_prefix) // "_aero_size_num.txt"
     elseif (dist_type == DIST_TYPE_MASS) then
        out_filename = trim(in_prefix) // "_aero_size_mass.txt"
     else
        call die(545030852)
     end if
  end if

  allocate(filename_list(0))
  call input_filename_list(in_prefix, filename_list)
  n_file = size(filename_list)
  call assert_msg(554271458, n_file > 0, &
       "no NetCDF files found with prefix: " // trim(in_prefix))

  call bin_grid_make(diam_grid, BIN_GRID_TYPE_LOG, n_bin, d_min, d_max)
  allocate(aero_dist(n_bin, n_file))

  do i_file = 1,n_file
     call input_state(filename_list(i_file), index, time, del_t, i_repeat, &
          uuid, aero_data=aero_data, aero_state=aero_state)

     if (i_file == 1) then
        run_uuid = uuid
     else
        call assert_msg(657993562, uuid == run_uuid, &
             "UUID mismatch between " // trim(filename_list(1)) // " and " &
             // trim(filename_list(i_file)))
     end if

     diameters = aero_state_diameters(aero_state, aero_data)
     num_concs = aero_state_num_concs(aero_state, aero_data)
     if (dist_type == DIST_TYPE_NUM) then
        hist = bin_grid_histogram_1d(diam_grid, diameters, num_concs)
     elseif (dist_type == DIST_TYPE_MASS) then
        masses = aero_state_masses(aero_state, aero_data)
        hist = bin_grid_histogram_1d(diam_grid, diameters, num_concs * masses)
     else
        call die(123323238)
     end if
     aero_dist(:, i_file) = hist
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
       diam_grid%widths(1)

  call open_file_write(out_filename, out_unit)
  do i_bin = 1,n_bin
     write(out_unit, '(e30.15e3)', advance='no') diam_grid%centers(i_bin)
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

    write(*,'(a)') 'Usage: extract_aero_size [options] <netcdf_prefix>'
    write(*,'(a)') ''
    write(*,'(a)') 'options are:'
    write(*,'(a)') '  -h, --help        Print this help message.'
    write(*,'(a)') '  -n, --num         Output number distribution.'
    write(*,'(a)') '  -m, --mass        Output mass distribution.'
    write(*,'(a)') '  -N, --dmin <D>    Minimum diameter (m).'
    write(*,'(a)') '  -X, --dmax <D>    Maximum diameter (m).'
    write(*,'(a)') '  -b, --nbin <N>    Number of size bins.'
    write(*,'(a)') '  -o, --out <file>  Output filename.'
    write(*,'(a)') ''
    write(*,'(a)') 'Examples:'
    write(*,'(a)') '  extract_aero_size --num data_0001'
    write(*,'(a)') &
         '  extract_aero_size --mass --dmin 1e-8 --dmax 1e-4 data_0001'
    write(*,'(a)') ''

  end subroutine print_help

end program extract_aero_size
