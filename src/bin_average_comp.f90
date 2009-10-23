! Copyright (C) 2009 Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The bin_average_comp program.

!> Read a NetCDF file, average the composition of all particles within
!> each bin, and write the data out as another NetCDF file.
program bin_average_comp

  use pmc_aero_state
  use pmc_gas_data
  use pmc_gas_state
  use pmc_env_state
  use pmc_aero_data
  use pmc_output
  use netcdf

  character(len=1000) :: in_filename, out_prefix
  type(bin_grid_t) :: bin_grid
  type(aero_data_t) :: aero_data
  type(aero_state_t) :: aero_state
  type(gas_data_t) :: gas_data
  type(gas_state_t) :: gas_state
  type(env_state_t) :: env_state
  integer :: n_bin, index, i_loop
  real(kind=dp) :: r_min, r_max, time, del_t
  character(len=1000) :: output_type, tmp_str
  logical :: record_removals

  ! process commandline arguments
  if (command_argument_count() .ne. 5) then
     write(6,*) 'Usage: bin_average_comp <r_min> <r_max> <n_bin> ' &
          // '<input_filename> <output_prefix>'
     stop 2
  endif
  call get_command_argument(1, tmp_str)
  r_min = string_to_real(tmp_str)
  call get_command_argument(2, tmp_str)
  r_max = string_to_real(tmp_str)
  call get_command_argument(3, tmp_str)
  n_bin = string_to_integer(tmp_str)
  call get_command_argument(4, in_filename)
  call get_command_argument(5, out_prefix)

  call bin_grid_allocate(bin_grid)
  call aero_data_allocate(aero_data)
  call aero_state_allocate(aero_state)
  call gas_data_allocate(gas_data)
  call gas_state_allocate(gas_state)
  call env_state_allocate(env_state)

  call bin_grid_make(bin_grid, n_bin, rad2vol(r_min), rad2vol(r_max))
  
  call input_state_netcdf(in_filename, bin_grid, aero_data, &
       aero_state, gas_data, gas_state, env_state, index, time, &
       del_t, i_loop)

  call aero_state_bin_average_comp(aero_state, bin_grid, aero_data)

  output_type = "central"
  record_removals = .false.
  call output_state(out_prefix, output_type, bin_grid, aero_data, &
       aero_state, gas_data, gas_state, env_state, index, time, &
       del_t, i_loop, record_removals)

  call bin_grid_deallocate(bin_grid)
  call aero_data_deallocate(aero_data)
  call aero_state_deallocate(aero_state)
  call gas_data_deallocate(gas_data)
  call gas_state_deallocate(gas_state)
  call env_state_deallocate(env_state)

end program bin_average_comp
