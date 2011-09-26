! Copyright (C) 2009-2011 Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The extract_aero_particle_mass program.

!> Read NetCDF output files and write out the individual particle
!> masses.
program extract_aero_particle_mass

  use pmc_aero_state
  use pmc_aero_particle
  use pmc_output

  character(len=1000) :: in_filename, out_filename
  type(aero_data_t) :: aero_data
  type(aero_state_t) :: aero_state
  type(gas_data_t) :: gas_data
  type(gas_state_t) :: gas_state
  type(env_state_t) :: env_state
  integer :: index, i_repeat, i_part, i_spec, out_unit
  real(kind=dp) :: time, del_t
  character(len=PMC_UUID_LEN) :: uuid
  type(aero_particle_t), pointer :: aero_particle

  if (command_argument_count() .ne. 2) then
     write(6,*) 'Usage: extract_aero_particle_mass ' &
          // '<netcdf_state_file> <output_filename>'
     stop 2
  endif
  call get_command_argument(1, in_filename)
  call get_command_argument(2, out_filename)

  call aero_data_allocate(aero_data)
  call aero_state_allocate(aero_state)
  call gas_data_allocate(gas_data)
  call gas_state_allocate(gas_state)
  call env_state_allocate(env_state)

  call input_state(in_filename, aero_data, aero_state, gas_data, gas_state, &
       env_state, index, time, del_t, i_repeat, uuid)

  write(*,'(a)') "Output file: " // trim(out_filename)
  write(*,'(a)') "  Output data is for time = " &
       // trim(real_to_string(time)) // " (s)"
  write(*,'(a)') "  Each row of output is one particle."
  write(*,'(a)') "  The columns of output are:"
  write(*,'(a)') "    column  1: particle ID number"
  write(*,'(a)') "    column  2: number concentration (m^{-3})"
  write(*,'(a)') "    column  3: particle diameter (m)"
  write(*,'(a)') "    column  4: particle total mass (kg)"
  do i_spec = 1,aero_data%n_spec
     write(*,'(a,i2,a,a,a,e10.4,a)') '    column ', i_spec + 4, ': ', &
          trim(aero_data%name(i_spec)), ' mass (kg) - density = ', &
          aero_data%density(i_spec), ' (kg/m^3)'
  end do

  call open_file_write(out_filename, out_unit)
  do i_part = 1,aero_state%apa%n_part
     aero_particle => aero_state%apa%particle(i_part)
     write(out_unit, '(i15)', advance='no') aero_particle%id
     write(out_unit, '(e30.15e3)', advance='no') &
          aero_state_particle_num_conc(aero_state, aero_particle)
     write(out_unit, '(e30.15e3)', advance='no') &
          aero_particle_diameter(aero_particle)
     write(out_unit, '(e30.15e3)', advance='no') &
          aero_particle_mass(aero_particle, aero_data)
     do i_spec = 1,aero_data%n_spec
        write(out_unit, '(e30.15e3)', advance='no') &
             aero_particle_species_mass(aero_particle, i_spec, aero_data)
     end do
     write(out_unit, *) ''
  end do
  call close_file(out_unit)

  call aero_data_deallocate(aero_data)
  call aero_state_deallocate(aero_state)
  call gas_data_deallocate(gas_data)
  call gas_state_deallocate(gas_state)
  call env_state_deallocate(env_state)

end program extract_aero_particle_mass
