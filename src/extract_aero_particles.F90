! Copyright (C) 2009-2012 Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The extract_aero_particles program.

!> Read NetCDF output files and write out the individual particle
!> masses.
program extract_aero_particles

  use pmc_aero_state
  use pmc_aero_particle
  use pmc_output
  use pmc_mpi
  use getopt_m

  character(len=PMC_MAX_FILENAME_LEN) :: in_filename, out_filename
  type(aero_data_t) :: aero_data
  type(aero_state_t) :: aero_state
  integer :: index, i_repeat, i_part, i_spec, out_unit, i_char
  real(kind=dp) :: time, del_t
  character(len=PMC_UUID_LEN) :: uuid
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
        call die_msg(885067714, 'unknown option: ' // trim(optopt))
     case default
        call print_help()
        call die_msg(516704561, 'unhandled option: ' // trim(optopt))
     end select
  end do

  if (optind /= command_argument_count()) then
     call print_help()
     call die_msg(605963222, &
          'expected exactly one non-option filename argument')
  end if

  call get_command_argument(optind, in_filename)

  if (out_filename == "") then
     i_char = scan(in_filename, '.', back=.true.)
     if (i_char == 0) then
        out_filename = trim(in_filename) // "_aero_particles.txt"
     else
        out_filename = trim(in_filename(1:(i_char - 1))) &
             // "_aero_particles.txt"
     end if
  end if

  call input_state(in_filename, index, time, del_t, i_repeat, uuid, &
       aero_data=aero_data, aero_state=aero_state)

  write(*,'(a)') "Output file: " // trim(out_filename)
  write(*,'(a)') "  Output data is for time = " &
       // trim(real_to_string(time)) // " (s)"
  write(*,'(a)') "  Each row of output is one particle."
  write(*,'(a)') "  The columns of output are:"
  write(*,'(a)') "    column  1: particle ID number"
  write(*,'(a)') "    column  2: number concentration (m^{-3})"
  write(*,'(a)') "    column  3: particle diameter (m)"
  write(*,'(a)') "    column  4: particle total mass (kg)"
  do i_spec = 1,aero_data_n_spec(aero_data)
     write(*,'(a,i2,a,a,a,e10.4,a)') '    column ', i_spec + 4, ': ', &
          trim(aero_data%name(i_spec)), ' mass (kg) - density = ', &
          aero_data%density(i_spec), ' (kg/m^3)'
  end do

  call open_file_write(out_filename, out_unit)
  do i_part = 1,aero_state_n_part(aero_state)
     write(out_unit, '(i15,e30.15e3,e30.15e3,e30.15e3)', advance='no') &
          aero_state%apa%particle(i_part)%id, &
          aero_state_particle_num_conc(aero_state, &
          aero_state%apa%particle(i_part), aero_data), &
          aero_particle_diameter(aero_state%apa%particle(i_part), aero_data), &
          aero_particle_mass(aero_state%apa%particle(i_part), aero_data)
     do i_spec = 1,aero_data_n_spec(aero_data)
        write(out_unit, '(e30.15e3)', advance='no') &
             aero_particle_species_mass(aero_state%apa%particle(i_part), &
             i_spec, aero_data)
     end do
     write(out_unit, *) ''
  end do
  call close_file(out_unit)

  call pmc_mpi_finalize()

contains

  subroutine print_help()

    write(*,'(a)') 'Usage: extract_aero_particles [options] <netcdf_prefix>'
    write(*,'(a)') ''
    write(*,'(a)') 'options are:'
    write(*,'(a)') '  -h, --help        Print this help message.'
    write(*,'(a)') '  -o, --out <file>  Output filename.'
    write(*,'(a)') ''
    write(*,'(a)') 'Examples:'
    write(*,'(a)') '  extract_aero_particles data_0001_00000001.nc'
    write(*,'(a)') ''

  end subroutine print_help

end program extract_aero_particles
