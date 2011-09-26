! Copyright (C) 2009-2011 Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The extract_aero_size_num program.

!> Read NetCDF output files and write out the aerosol number size
!> distributions in text format.
program extract_aero_size_num

  use pmc_aero_state
  use pmc_aero_particle
  use pmc_output

  character(len=PMC_MAX_FILENAME_LEN) :: in_prefix, out_filename
  character(len=PMC_MAX_FILENAME_LEN), allocatable :: filename_list(:)
  character(len=1000) :: tmp_str
  type(bin_grid_t) :: diam_grid
  type(aero_data_t) :: aero_data
  type(aero_state_t) :: aero_state
  type(gas_data_t) :: gas_data
  type(gas_state_t) :: gas_state
  type(env_state_t) :: env_state
  integer :: index, i_repeat, i_part, i_spec, out_unit
  integer :: i_file, n_file, i_bin, n_bin
  real(kind=dp) :: time, del_t
  type(aero_particle_t), pointer :: aero_particle
  real(kind=dp) :: d_min, d_max, diam, volume, log_width
  character(len=PMC_UUID_LEN) :: uuid, run_uuid
  real(kind=dp), allocatable :: diameters(:), num_concs(:), hist(:)
  real(kind=dp), allocatable :: aero_dist(:,:)

  if (command_argument_count() .ne. 5) then
     write(6,*) 'Usage: extract_aero_size_num <d_min> <d_max> <n_bin> ' &
          // '<netcdf_state_prefix> <output_filename>'
     stop 2
  endif
  call get_command_argument(1, tmp_str)
  d_min = string_to_real(tmp_str)
  call get_command_argument(2, tmp_str)
  d_max = string_to_real(tmp_str)
  call get_command_argument(3, tmp_str)
  n_bin = string_to_integer(tmp_str)
  call get_command_argument(4, in_prefix)
  call get_command_argument(5, out_filename)

  call bin_grid_allocate(diam_grid)
  call aero_data_allocate(aero_data)
  call aero_state_allocate(aero_state)
  call gas_data_allocate(gas_data)
  call gas_state_allocate(gas_state)
  call env_state_allocate(env_state)

  allocate(filename_list(0))
  call input_filename_list(in_prefix, filename_list)
  n_file = size(filename_list)
  call assert_msg(875939143, n_file > 0, &
       "no NetCDF files found with prefix: " // trim(in_prefix))

  call bin_grid_make(diam_grid, n_bin, d_min, d_max)
  allocate(aero_dist(n_bin, n_file))
  allocate(hist(n_bin))
  allocate(diameters(0))
  allocate(num_concs(0))

  do i_file = 1,n_file
     call input_state(filename_list(i_file), aero_data, aero_state, gas_data, &
          gas_state, env_state, index, time, del_t, i_repeat, uuid)

     if (i_file == 1) then
        run_uuid = uuid
     else
        call assert_msg(657993562, uuid == run_uuid, &
             "UUID mismatch between " // trim(filename_list(1)) // " and " &
             // trim(filename_list(i_file)))
     end if

     call aero_state_diameters(aero_state, diameters)
     call aero_state_num_concs(aero_state, num_concs)
     call bin_grid_histogram_1d(diam_grid, diameters, num_concs, hist)
     aero_dist(:, i_file) = hist
  end do

  write(*,'(a)') "Output file: " // trim(out_filename)
  write(*,'(a)') "  Each row of output is one size bin."
  write(*,'(a)') "  The columns of output are:"
  write(*,'(a)') "    column   1: bin center diameter (m)"
  write(*,'(a)') "    column j+1: number concentration at time(j) (#/m^3)"
  write(*,'(a)') "  Diameter bins have logarithmic width:"
  write(*,'(a,e20.10)') "    log_width = ln(diam(i+1)) - ln(diam(i)) =", &
       diam_grid%log_width

  call open_file_write(out_filename, out_unit)
  do i_bin = 1,n_bin
     write(out_unit, '(e30.15e3)', advance='no') &
          diam_grid%center_radius(i_bin)
     do i_file = 1,n_file
        write(out_unit, '(e30.15e3)', advance='no') aero_dist(i_bin, i_file)
     end do
     write(out_unit, '(a)') ''
  end do
  call close_file(out_unit)

  deallocate(filename_list)
  deallocate(aero_dist)
  deallocate(hist)
  deallocate(diameters)
  deallocate(num_concs)
  call bin_grid_allocate(diam_grid)
  call aero_data_deallocate(aero_data)
  call aero_state_deallocate(aero_state)
  call gas_data_deallocate(gas_data)
  call gas_state_deallocate(gas_state)
  call env_state_deallocate(env_state)

end program extract_aero_size_num
