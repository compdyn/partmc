! Copyright (C) 2009-2012, 2016 Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The extract_sectional_aero_time program.

!> Read NetCDF sectional output files and write out the time evolution of
!> aerosol number and mass concentrations in text format.
program extract_sectional_aero_time

  use pmc_aero_binned
  use pmc_output
  use pmc_mpi
  use getopt_m

  character(len=PMC_MAX_FILENAME_LEN) :: in_prefix, out_filename
  character(len=PMC_MAX_FILENAME_LEN), allocatable :: filename_list(:)
  type(bin_grid_t) :: bin_grid
  type(aero_data_t) :: aero_data
  type(aero_binned_t) :: aero_binned
  integer :: index, out_unit
  integer :: i_file, n_file, i_spec
  real(kind=dp) :: time, del_t
  character(len=PMC_UUID_LEN) :: uuid, run_uuid
  real(kind=dp), allocatable :: times(:), time_num_concs(:), time_mass_concs(:)
  real(kind=dp), allocatable :: time_species_concs(:,:)
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
        call die_msg(559277549, 'unknown option: ' // trim(optopt))
     case default
        call print_help()
        call die_msg(247146202, 'unhandled option: ' // trim(optopt))
     end select
  end do

  if (optind /= command_argument_count()) then
     call print_help()
     call die_msg(959981535, 'expected exactly one non-option prefix argument')
  end if

  call get_command_argument(optind, in_prefix)

  if (out_filename == "") then
     out_filename = trim(in_prefix) // "_aero_time.txt"
  end if

  allocate(filename_list(0))
  call input_filename_list(in_prefix, filename_list)
  n_file = size(filename_list)
  call assert_msg(875939143, n_file > 0, &
       "no NetCDF files found with prefix: " // trim(in_prefix))

  call input_sectional(filename_list(1), index, time, del_t, uuid, &
       bin_grid=bin_grid, aero_data=aero_data, aero_binned=aero_binned)
  run_uuid = uuid

  allocate(times(n_file))
  allocate(time_num_concs(n_file))
  allocate(time_mass_concs(n_file))
  allocate(time_species_concs(n_file, aero_data_n_spec(aero_data)))

  do i_file = 1,n_file
     call input_sectional(filename_list(i_file), index, time, del_t, uuid, &
          bin_grid=bin_grid, aero_data=aero_data, aero_binned=aero_binned)

     call assert_msg(865522513, uuid == run_uuid, &
          "UUID mismatch between " // trim(filename_list(1)) // " and " &
          // trim(filename_list(i_file)))

     times(i_file) = time
     time_num_concs(i_file) = sum(aero_binned%num_conc * bin_grid%widths)
     ! FIXME: line below assumes bin_grid%widths is constant
     time_species_concs(i_file, :) = sum(aero_binned%vol_conc &
          * bin_grid%widths(1), 1) * aero_data%density
     time_mass_concs(i_file) = sum(time_species_concs(i_file, :))
  end do

  write(*,'(a,a)') "Output file: ", trim(out_filename)
  write(*,'(a)') "  Each row of output is one time."
  write(*,'(a)') "  The columns of output are:"
  write(*,'(a)') "    column  1: time (s)"
  write(*,'(a)') "    column  2: aerosol number concentration (#/m^3)"
  write(*,'(a)') "    column  3: aerosol mass concentration (kg/m^3)"
  do i_spec = 1,aero_data_n_spec(aero_data)
     write(*,'(a,i2,a,a,a)') "    column ", i_spec + 3, ": aerosol ", &
          trim(aero_data%name(i_spec)), " concentration (kg/m^3)"
  end do

  call open_file_write(out_filename, out_unit)
  do i_file = 1,n_file
     write(out_unit, '(e30.15e3)', advance='no') times(i_file)
     write(out_unit, '(e30.15e3)', advance='no') time_num_concs(i_file)
     write(out_unit, '(e30.15e3)', advance='no') time_mass_concs(i_file)
     do i_spec = 1,aero_data_n_spec(aero_data)
        write(out_unit, '(e30.15e3)', advance='no') &
             time_species_concs(i_file, i_spec)
     end do
     write(out_unit, '(a)') ''
  end do
  call close_file(out_unit)

  deallocate(times)
  deallocate(time_num_concs)
  deallocate(time_mass_concs)
  deallocate(time_species_concs)
  deallocate(filename_list)

  call pmc_mpi_finalize()

contains

  subroutine print_help()

    write(*,'(a)') 'Usage: extract_sectional_aero_time [options] ' &
         // '<netcdf_prefix>'
    write(*,'(a)') ''
    write(*,'(a)') 'options are:'
    write(*,'(a)') '  -h, --help        Print this help message.'
    write(*,'(a)') '  -o, --out <file>  Output filename.'
    write(*,'(a)') ''
    write(*,'(a)') 'Examples:'
    write(*,'(a)') '  extract_sectional_aero_time data_0001'
    write(*,'(a)') ''

  end subroutine print_help

end program extract_sectional_aero_time
