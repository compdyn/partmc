! Copyright (C) 2009-2012 Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The extract_gas program.

!> Read NetCDF output files and write out the gas mixing ratios in text
!> format.
program extract_gas

  use pmc_gas_state
  use pmc_output
  use pmc_mpi
  use getopt_m

  character(len=PMC_MAX_FILENAME_LEN) :: in_prefix, out_filename
  character(len=PMC_MAX_FILENAME_LEN), allocatable :: filename_list(:)
  character(len=1000) :: tmp_str
  type(gas_data_t) :: gas_data
  type(gas_state_t) :: gas_state
  integer :: index, i_repeat, i_spec, out_unit
  integer :: i_file, n_file
  real(kind=dp) :: time, del_t
  character(len=PMC_UUID_LEN) :: uuid, run_uuid
  real(kind=dp), allocatable :: times(:), gas_mixing_ratios(:,:)
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
        call die_msg(478715112, 'unknown option: ' // trim(optopt))
     case default
        call print_help()
        call die_msg(935521190, 'unhandled option: ' // trim(optopt))
     end select
  end do

  if (optind /= command_argument_count()) then
     call print_help()
     call die_msg(744333329, 'expected exactly one non-option prefix argument')
  end if

  call get_command_argument(optind, in_prefix)

  if (out_filename == "") then
     out_filename = trim(in_prefix) // "_gas.txt"
  end if

  allocate(filename_list(0))
  call input_filename_list(in_prefix, filename_list)
  n_file = size(filename_list)
  call assert_msg(579059629, n_file > 0, &
       "no NetCDF files found with prefix: " // trim(in_prefix))

  call input_state(filename_list(1), index, time, del_t, i_repeat, uuid, &
       gas_data=gas_data, gas_state=gas_state)
  run_uuid = uuid

  allocate(times(n_file))
  allocate(gas_mixing_ratios(n_file, gas_data_n_spec(gas_data)))

  do i_file = 1,n_file
     call input_state(filename_list(i_file), index, time, del_t, i_repeat, &
          uuid, gas_data=gas_data, gas_state=gas_state)

     call assert_msg(390171757, uuid == run_uuid, &
          "UUID mismatch between " // trim(filename_list(1)) // " and " &
          // trim(filename_list(i_file)))

     times(i_file) = time
     gas_mixing_ratios(i_file, :) = gas_state%mix_rat
  end do

  write(*,'(a,a)') "Output file: ", trim(out_filename)
  write(*,'(a)') "  Each row of output is one time."
  write(*,'(a)') "  The columns of output are:"
  write(*,'(a)') "    column  1: time (s)"
  do i_spec = 1,gas_data_n_spec(gas_data)
     write(*,'(a,i2,a,a,a)') "    column ", i_spec + 1, ": gas ", &
          trim(gas_data%name(i_spec)), " mixing ratio (ppb)"
  end do

  call open_file_write(out_filename, out_unit)
  do i_file = 1,n_file
     write(out_unit, '(e30.15e3)', advance='no') times(i_file)
     do i_spec = 1,gas_data_n_spec(gas_data)
        write(out_unit, '(e30.15e3)', advance='no') &
             gas_mixing_ratios(i_file, i_spec)
     end do
     write(out_unit, '(a)') ''
  end do
  call close_file(out_unit)

  deallocate(times)
  deallocate(gas_mixing_ratios)

  call pmc_mpi_finalize()

contains

  subroutine print_help()

    write(*,'(a)') 'Usage: extract_gas [options] <netcdf_prefix>'
    write(*,'(a)') ''
    write(*,'(a)') 'options are:'
    write(*,'(a)') '  -h, --help        Print this help message.'
    write(*,'(a)') '  -o, --out <file>  Output filename.'
    write(*,'(a)') ''
    write(*,'(a)') 'Examples:'
    write(*,'(a)') '  extract_gas data_0001'
    write(*,'(a)') ''

  end subroutine print_help

end program extract_gas
