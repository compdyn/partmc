!! Copyright (C) 2009-2012 Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The freezing_process program.

!> Read NetCDF output files and write out the time evolution of
!> temperature and frozen fraction in text format.
program freezing_process

  use pmc_aero_state
  use pmc_aero_particle
  use pmc_output
  use pmc_mpi
  use getopt_m

  character(len=PMC_MAX_FILENAME_LEN) :: in_prefix, out_filename
  character(len=PMC_MAX_FILENAME_LEN), allocatable :: filename_list(:)
  character(len=1000) :: tmp_str
  type(aero_data_t) :: aero_data
  type(aero_state_t) :: aero_state
  type(env_state_t) :: env_state
  integer :: index, i_repeat, i_spec, out_unit
  integer :: i_file, n_file
  real(kind=dp) :: time, del_t
  character(len=PMC_UUID_LEN) :: uuid, run_uuid
  real(kind=dp), allocatable :: times(:), temperatures(:), frozen_fractions(:)
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
        call die_msg(514364550, 'unknown option: ' // trim(optopt))
     case default
        call print_help()
        call die_msg(603100341, 'unhandled option: ' // trim(optopt))
     end select
  end do

  if (optind /= command_argument_count()) then
     call print_help()
     call die_msg(967032896, 'expected exactly one non-option prefix argument')
  end if

  call get_command_argument(optind, in_prefix)

  if (out_filename == "") then
     out_filename = trim(in_prefix) // "_aero_time.txt"
  end if

  allocate(filename_list(0))
  call input_filename_list(in_prefix, filename_list)
  n_file = size(filename_list)
  call assert_msg(323514871, n_file > 0, &
       "no NetCDF files found with prefix: " // trim(in_prefix))

  call input_state(filename_list(1), index, time, del_t, i_repeat, uuid, &
       aero_data=aero_data, aero_state=aero_state, env_state=env_state)
  run_uuid = uuid

  allocate(times(n_file))
  allocate(temperatures(n_file))
  allocate(frozen_fractions(n_file))

  do i_file = 1,n_file
     call input_state(filename_list(i_file), index, time, del_t, i_repeat, &
          uuid, aero_data=aero_data, aero_state=aero_state, env_state=env_state)

     call assert_msg(397906326, uuid == run_uuid, &
          "UUID mismatch between " // trim(filename_list(1)) // " and " &
          // trim(filename_list(i_file)))

     times(i_file) = time
     temperatures(i_file) = env_state%temp
     frozen_fractions(i_file) = aero_state_frozen_fraction(aero_state, aero_data)
     
  end do

  write(*,'(a,a)') "Output file: ", trim(out_filename)
  write(*,'(a)') "  Each row of output is one time."
  write(*,'(a)') "  The columns of output are:"
  write(*,'(a)') "    column  1: time (s)"
  write(*,'(a)') "    column  2: temperature (K)"
  write(*,'(a)') "    column  3: frozen fraction (unitless)"

  call open_file_write(out_filename, out_unit)
  do i_file = 1,n_file
     write(out_unit, '(e30.15e3)', advance='no') times(i_file)
     write(out_unit, '(e30.15e3)', advance='no') temperatures(i_file)
     write(out_unit, '(e30.15e3)', advance='no') frozen_fractions(i_file)
     write(out_unit, '(a)') ''
  end do
  call close_file(out_unit)

  deallocate(times)
  deallocate(temperatures)
  deallocate(frozen_fractions)
  deallocate(filename_list)

  call pmc_mpi_finalize()

contains

  subroutine print_help()

    write(*,'(a)') 'Usage: freezing_process [options] <netcdf_prefix>'
    write(*,'(a)') ''
    write(*,'(a)') 'options are:'
    write(*,'(a)') '  -h, --help        Print this help message.'
    write(*,'(a)') '  -o, --out <file>  Output filename.'
    write(*,'(a)') ''
    write(*,'(a)') 'Examples:'
    write(*,'(a)') '  freezing_process freezing_part_0001'
    write(*,'(a)') ''

  end subroutine print_help

end program freezing_process


