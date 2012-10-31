! Copyright (C) 2009-2012 Matthew West
! Copyright (C) 2012 Jian Tian
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The test_fractal_sectional_dimless_time program.

!> Read sectional NetCDF output files and write out the time evolution of
!> aerosol dimensionless number concentrations in text format.
program test_fractal_sectional_dimless_time

  use pmc_aero_binned
  use pmc_output
  use pmc_mpi
  use pmc_util
  use pmc_constants
  use getopt_m

  integer, parameter :: REGIME_FREE = 1
  integer, parameter :: REGIME_CONT = 2

  character(len=PMC_MAX_FILENAME_LEN) :: in_prefix, out_filename
  character(len=PMC_MAX_FILENAME_LEN), allocatable :: filename_list(:)
  character(len=1000) :: tmp_str
  type(bin_grid_t) :: bin_grid
  type(aero_data_t) :: aero_data
  type(aero_binned_t) :: aero_binned
  type(env_state_t) :: env_state
  integer :: index, i_spec, out_unit
  integer :: i_file, n_file, regime
  real(kind=dp) :: time, del_t, density, n_init
  character(len=PMC_UUID_LEN) :: uuid, run_uuid
  real(kind=dp), allocatable :: times(:), time_num_concs(:)
  real(kind=dp), allocatable :: time_mass_concs(:)
  real(kind=dp), allocatable :: time_vol_concs(:)
  real(kind=dp), allocatable :: dimless_times(:)
  real(kind=dp), allocatable :: dimless_time_num_concs(:)
  real(kind=dp), allocatable :: time_species_concs(:,:)
  real(kind=dp), allocatable :: time_species_vol_concs(:,:)
  type(option_s) :: opts(3)

  call pmc_mpi_init()

  opts(1) = option_s("free", .false., 'f')
  opts(2) = option_s("cont", .false., 'c')
  opts(3) = option_s("output", .true., 'o')

  regime = 0
  n_init = 1d14
  out_filename = ""

  do
     select case(getopt("fco:", opts))
     case(char(0))
        exit
     case('f')
        regime = REGIME_FREE
     case('c')
        regime = REGIME_CONT
     case('o')
        out_filename = optarg
     end select
  end do

  if (regime == 0) then
     call die_msg(377136887, 'missing aerosol regime: ' &
          // 'please specify --free or --cont')
  end if

  if (optind /= command_argument_count()) then
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
       bin_grid=bin_grid, aero_data=aero_data, &
       aero_binned=aero_binned, env_state=env_state)
  run_uuid = uuid

  allocate(times(n_file))
  allocate(dimless_times(n_file))
  allocate(time_num_concs(n_file))
  allocate(time_mass_concs(n_file))
  allocate(time_vol_concs(n_file))
  allocate(dimless_time_num_concs(n_file))
  allocate(time_species_concs(n_file, aero_data%n_spec))
  allocate(time_species_vol_concs(n_file, aero_data%n_spec))

  do i_file = 1,n_file
     call input_sectional(filename_list(i_file), index, time, del_t, uuid, &
          bin_grid=bin_grid, aero_data=aero_data, aero_binned=aero_binned, &
          env_state=env_state)

     call assert_msg(397906326, uuid == run_uuid, &
          "UUID mismatch between " // trim(filename_list(1)) // " and " &
          // trim(filename_list(i_file)))

     times(i_file) = time
     time_num_concs(i_file) = sum(aero_binned%num_conc * bin_grid%log_width)
     dimless_time_num_concs(i_file) = time_num_concs(i_file) / n_init
     time_species_concs(i_file, :) = sum(aero_binned%vol_conc &
          * bin_grid%log_width, 1) * aero_data%density
     time_mass_concs(i_file) = sum(time_species_concs(i_file, :))
     time_species_vol_concs(i_file, :) = sum(aero_binned%vol_conc &
          * bin_grid%log_width, 1)
     time_vol_concs(i_file) = sum(time_species_vol_concs(i_file, :))
     density = time_mass_concs(i_file) / time_vol_concs(i_file)
     if (regime == REGIME_FREE) then
        dimless_times(i_file) = (6d0 * const%boltzmann * env_state%temp &
             * aero_data%fractal%prime_radius / density)**(1d0 / 2d0)   &
             * n_init * times(i_file)
     elseif (regime == REGIME_CONT) then
        dimless_times(i_file) = (2d0 * const%boltzmann * env_state%temp &
             / 3d0 / const%air_dyn_visc) * n_init * times(i_file)
     else
        call die(123323239)
     end if
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
  deallocate(time_mass_concs)
  deallocate(time_vol_concs)
  deallocate(dimless_time_num_concs)
  deallocate(time_species_concs)
  deallocate(time_species_vol_concs)
  deallocate(filename_list)
  call bin_grid_allocate(bin_grid)
  call aero_data_deallocate(aero_data)
  call aero_binned_deallocate(aero_binned)
  call env_state_deallocate(env_state)

  call pmc_mpi_finalize()

end program test_fractal_sectional_dimless_time
