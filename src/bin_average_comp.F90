! Copyright (C) 2009-2010 Matthew West
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
  use pmc_aero_weight
  use pmc_output
  use pmc_rand
  use netcdf

  character(len=1000) :: in_filename, out_prefix
  type(bin_grid_t) :: bin_grid
  type(aero_data_t) :: aero_data
  type(aero_weight_t) :: aero_weight
  type(aero_state_t) :: aero_state
  type(gas_data_t) :: gas_data
  type(gas_state_t) :: gas_state
  type(env_state_t) :: env_state
  integer :: n_bin, index, i_repeat, output_type
  real(kind=dp) :: d_min, d_max, time, del_t
  character(len=1000) :: tmp_str
  logical :: record_removals, dry_volume, record_optical
  character(len=PMC_UUID_LEN) :: uuid

  ! process commandline arguments
  if (command_argument_count() .ne. 6) then
     write(6,*) 'Usage: bin_average_comp <d_min> <d_max> <n_bin> ' &
          // '<"wet" or "dry"> <input_filename> <output_prefix>'
     write(6,*) ''
     write(6,*) '  d_min: minimum bin diameter (m)'
     write(6,*) '  d_max: maximum bin diameter (m)'
     write(6,*) '  n_bin: number of bins'
     write(6,*) '  wet/dry: average wet or dry sizes'
     write(6,*) '  input_filename: like scenario_0001_00000001.nc'
     write(6,*) '  output_prefix: like scenario_comp_average'
     stop 2
  endif
  call get_command_argument(1, tmp_str)
  d_min = string_to_real(tmp_str)
  call get_command_argument(2, tmp_str)
  d_max = string_to_real(tmp_str)
  call get_command_argument(3, tmp_str)
  n_bin = string_to_integer(tmp_str)
  call get_command_argument(4, tmp_str)
  if (trim(tmp_str) == "wet") then
     dry_volume = .false.
  elseif (trim(tmp_str) == "dry") then
     dry_volume = .true.
  else
     write(6,*) 'Argument 4 must be "wet" or "dry", not ' &
          // trim(tmp_str)
     stop 1
  end if
  call get_command_argument(5, in_filename)
  call get_command_argument(6, out_prefix)

  call pmc_mpi_init()

  call bin_grid_allocate(bin_grid)
  call aero_data_allocate(aero_data)
  call aero_weight_allocate(aero_weight)
  call aero_state_allocate(aero_state)
  call gas_data_allocate(gas_data)
  call gas_state_allocate(gas_state)
  call env_state_allocate(env_state)

  call bin_grid_make(bin_grid, n_bin, diam2rad(d_min), diam2rad(d_max))

  call input_state(in_filename, bin_grid, aero_data, &
       aero_weight, aero_state, gas_data, gas_state, env_state, &
       index, time, del_t, i_repeat, uuid)

  if (dry_volume) then
     call aero_state_make_dry(aero_state, bin_grid, aero_data)
  end if

  call aero_state_bin_average_comp(aero_state, bin_grid, aero_data, &
       aero_weight, dry_volume)

  output_type = OUTPUT_TYPE_SINGLE
  record_removals = .false.
  record_optical = .true.
  call output_state(out_prefix, output_type, bin_grid, aero_data, &
       aero_weight, aero_state, gas_data, gas_state, env_state, &
       index, time, del_t, i_repeat, record_removals, record_optical, uuid)

  call bin_grid_deallocate(bin_grid)
  call aero_data_deallocate(aero_data)
  call aero_weight_deallocate(aero_weight)
  call aero_state_deallocate(aero_state)
  call gas_data_deallocate(gas_data)
  call gas_state_deallocate(gas_state)
  call env_state_deallocate(env_state)
  
  call pmc_mpi_finalize()

end program bin_average_comp