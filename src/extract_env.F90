! Copyright (C) 2009-2012 Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The extract_env program.

!> Read NetCDF output files and write out the environment variables in
!> text format.
program extract_env

  use pmc_env_state
  use pmc_output
  use pmc_mpi
  use getopt_m

  character(len=PMC_MAX_FILENAME_LEN) :: in_prefix, out_filename
  character(len=PMC_MAX_FILENAME_LEN), allocatable :: filename_list(:)
  character(len=1000) :: tmp_str
  type(env_state_t) :: env_state
  integer :: index, i_repeat, i_spec, out_unit
  integer :: i_file, n_file
  real(kind=dp) :: time, del_t
  character(len=PMC_UUID_LEN) :: uuid, run_uuid
  real(kind=dp), allocatable :: times(:), temps(:), rel_humids(:)
  real(kind=dp), allocatable :: pressures(:), mix_heights(:)
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
     case( '?' )
        call print_help()
        call die_msg(909107230, 'unknown option: ' // trim(optopt))
     case default
        call print_help()
        call die_msg(368158543, 'unhandled option: ' // trim(optopt))
     end select
  end do

  if (optind /= command_argument_count()) then
     call print_help()
     call die_msg(410427558, 'expected exactly one non-option prefix argument')
  end if

  call get_command_argument(optind, in_prefix)

  if (out_filename == "") then
     out_filename = trim(in_prefix) // "_env.txt"
  end if

  allocate(filename_list(0))
  call input_filename_list(in_prefix, filename_list)
  n_file = size(filename_list)
  call assert_msg(399220907, n_file > 0, &
       "no NetCDF files found with prefix: " // trim(in_prefix))

  call input_state(filename_list(1), index, time, del_t, i_repeat, uuid, &
       env_state=env_state)
  run_uuid = uuid

  allocate(times(n_file))
  allocate(temps(n_file))
  allocate(rel_humids(n_file))
  allocate(pressures(n_file))
  allocate(mix_heights(n_file))

  do i_file = 1,n_file
     call input_state(filename_list(i_file), index, time, del_t, i_repeat, &
          uuid, env_state=env_state)

     call assert_msg(276800431, uuid == run_uuid, &
          "UUID mismatch between " // trim(filename_list(1)) // " and " &
          // trim(filename_list(i_file)))

     times(i_file) = time
     temps(i_file) = env_state%temp
     rel_humids(i_file) = env_state%rel_humid
     pressures(i_file) = env_state%pressure
     mix_heights(i_file) = env_state%height
  end do

  write(*,'(a,a)') "Output file: ", trim(out_filename)
  write(*,'(a)') "  Each row of output is one time."
  write(*,'(a)') "  The columns of output are:"
  write(*,'(a)') "    column 1: time (s)"
  write(*,'(a)') "    column 2: temperature (K)"
  write(*,'(a)') "    column 3: relative_humidity (1)"
  write(*,'(a)') "    column 4: pressure (Pa)"
  write(*,'(a)') "    column 5: mixing height (m)"

  call open_file_write(out_filename, out_unit)
  do i_file = 1,n_file
     write(out_unit, '(e30.15e3)', advance='no') times(i_file)
     write(out_unit, '(e30.15e3)', advance='no') temps(i_file)
     write(out_unit, '(e30.15e3)', advance='no') rel_humids(i_file)
     write(out_unit, '(e30.15e3)', advance='no') pressures(i_file)
     write(out_unit, '(e30.15e3)', advance='no') mix_heights(i_file)
     write(out_unit, '(a)') ''
  end do
  call close_file(out_unit)

  deallocate(times)
  deallocate(temps)
  deallocate(rel_humids)
  deallocate(pressures)
  deallocate(mix_heights)

  call pmc_mpi_finalize()

contains

  subroutine print_help()

    write(*,'(a)') 'Usage: extract_env [options] <netcdf_prefix>'
    write(*,'(a)') ''
    write(*,'(a)') 'options are:'
    write(*,'(a)') '  -h, --help        Print this help message.'
    write(*,'(a)') '  -o, --out <file>  Output filename.'
    write(*,'(a)') ''
    write(*,'(a)') 'Examples:'
    write(*,'(a)') '  extract_env data_0001'
    write(*,'(a)') ''

  end subroutine print_help

end program extract_env
