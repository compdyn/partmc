! Copyright (C) 2009-2012 Matthew West
! Copyright (C) 2012 Jian Tian
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The extract_sectional_dimless_time program.

!> Read sectional NetCDF output files and write out the time evolution of
!> aerosol dimensionless number concentrations in text format.
program extract_sectional_dimless_time

  use pmc_aero_binned
  use pmc_output
  use pmc_mpi
  use pmc_util
  use pmc_constants
  use getopt_m

  integer, parameter :: REGIME_FREE = 0
  integer, parameter :: REGIME_CONT = 1

  character(len=PMC_MAX_FILENAME_LEN) :: in_prefix, out_filename
  character(len=PMC_MAX_FILENAME_LEN), allocatable :: filename_list(:)
  character(len=1000) :: tmp_str
  type(bin_grid_t) :: bin_grid
  type(aero_data_t) :: aero_data
  type(aero_binned_t) :: aero_binned
  type(env_state_t) :: env_state
  integer :: index, i_spec, out_unit
  integer :: i_file, n_file, regime
  real(kind=dp) :: time, del_t
  character(len=PMC_UUID_LEN) :: uuid, run_uuid
  real(kind=dp), allocatable :: times(:), time_num_concs(:)
  real(kind=dp), allocatable :: dimless_times(:), dimless_time_num_concs(:)
  real(kind=dp), parameter :: N_INIT = 1d14
  real(kind=dp), parameter :: DENSITY = 4200d0
  type(option_s) :: opts(4)

  call pmc_mpi_init()

  opts(1) = option_s("help", .false., 'h')
  opts(2) = option_s("free", .false., 'f')
  opts(3) = option_s("cont", .false., 'c')
  opts(4) = option_s("output", .true., 'o')

  out_filename = ""

  do
     select case(getopt("hfco:", opts))
     case(char(0))
        exit
     case('h')
        call print_help()
        stop
     case('f')
        regime = REGIME_FREE
     case('c')
        regime = REGIME_CONT
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
     out_filename = trim(in_prefix) // "_dimless_time.txt"
  end if

  call bin_grid_allocate(bin_grid)
  call aero_data_allocate(aero_data)
  call env_state_allocate(env_state)
  call aero_binned_allocate(aero_binned)

  allocate(filename_list(0))
  call input_filename_list(in_prefix, filename_list)
  n_file = size(filename_list)
  call assert_msg(323514871, n_file > 0, &
       "no NetCDF files found with prefix: " // trim(in_prefix))

  call input_sectional(filename_list(1), index, time, del_t, uuid, &
       bin_grid=bin_grid, aero_data=aero_data, aero_binned=aero_binned, env_state=env_state)
  run_uuid = uuid

  allocate(times(n_file))
  allocate(dimless_times(n_file))
  allocate(time_num_concs(n_file))
  allocate(dimless_time_num_concs(n_file))

  do i_file = 1,n_file
     call input_sectional(filename_list(i_file), index, time, del_t, uuid, &
       bin_grid=bin_grid, aero_data=aero_data, aero_binned=aero_binned, env_state=env_state)

     call assert_msg(397906326, uuid == run_uuid, &
          "UUID mismatch between " // trim(filename_list(1)) // " and " &
          // trim(filename_list(i_file)))

     times(i_file) = time
     if (regime == REGIME_FREE) then
        dimless_times(i_file) = (6d0 * const%boltzmann * env_state%temp &
             * aero_data%fractal%prime_radius / DENSITY)**(1d0 / 2d0)   &
             * (3d0 / 4d0 / const%pi)**(1d0 / 6d0) * N_INIT * times(i_file)
     elseif (regime == REGIME_CONT) then
        dimless_times(i_file) = (2d0 * const%boltzmann * env_state%temp &
             / 3d0 / const%air_dyn_visc) * N_INIT * times(i_file)
     else
        call die(123323239)
     end if
     time_num_concs(i_file) = sum(aero_binned%num_conc * bin_grid%log_width)
     dimless_time_num_concs(i_file) = time_num_concs(i_file) / N_INIT
  end do

  write(*,'(a,a)') "Output file: ", trim(out_filename)
  write(*,'(a)') "  Each row of output is one time."
  write(*,'(a)') "  The columns of output are:"
  write(*,'(a)') "    column  1: dimensionless time"
  write(*,'(a)') "    column  2: dimensionless aerosol number concentration"

  call open_file_write(out_filename, out_unit)
  do i_file = 1,n_file
     write(out_unit, '(e30.15e3)', advance='no') dimless_times(i_file)
     write(out_unit, '(e30.15e3)', advance='no') dimless_time_num_concs(i_file)
     write(out_unit, '(a)') ''
  end do
  call close_file(out_unit)

  deallocate(times)
  deallocate(dimless_times)
  deallocate(time_num_concs)
  deallocate(dimless_time_num_concs)
  deallocate(filename_list)
  call bin_grid_allocate(bin_grid)
  call aero_data_deallocate(aero_data)
  call aero_binned_deallocate(aero_binned)
  call env_state_deallocate(env_state)

  call pmc_mpi_finalize()

contains

  subroutine print_help()

    write(*,'(a)') 'Usage: extract_sectional_dimless_time [options] <netcdf_prefix>'
    write(*,'(a)') ''
    write(*,'(a)') 'options are:'
    write(*,'(a)') '  -h, --help        Print this help message.'
    write(*,'(a)') '  -o, --out <file>  Output filename.'
    write(*,'(a)') ''
    write(*,'(a)') 'Examples:'
    write(*,'(a)') '  extract_sectional_dimless_time data_0001'
    write(*,'(a)') ''

  end subroutine print_help

end program extract_sectional_dimless_time
