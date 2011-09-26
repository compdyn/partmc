! Copyright (C) 2009-2011 Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The extract_aero_size program.

!> Read NetCDF output files and write out the aerosol number or mass
!> size distributions in text format.
program extract_aero_size

  use pmc_aero_state
  use pmc_aero_particle
  use pmc_output
  use getopt_m

  character(len=PMC_MAX_FILENAME_LEN) :: in_prefix, out_filename
  character(len=PMC_MAX_FILENAME_LEN), allocatable :: filename_list(:)
  character(len=1000) :: tmp_str
  type(aero_data_t) :: aero_data
  type(aero_state_t) :: aero_state
  type(gas_data_t) :: gas_data
  type(gas_state_t) :: gas_state
  type(env_state_t) :: env_state
  integer :: index, i_repeat, i_spec, out_unit
  integer :: i_file, n_file
  real(kind=dp) :: time, del_t
  character(len=PMC_UUID_LEN) :: uuid, run_uuid
  real(kind=dp), allocatable :: particle_num_concs(:), particle_masses(:)
  real(kind=dp), allocatable :: times(:), time_num_concs(:), time_mass_concs(:)
  real(kind=dp), allocatable :: time_species_concs(:,:)
  type(option_s) :: opts(2)

  opts(1) = option_s("help", .false., 'h')
  opts(2) = option_s("output", .true., 'o')

  out_filename = ""

  do
     select case(getopt("hnmN:X:b:o:", opts))
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
     out_filename = trim(in_prefix) // "_time.txt"
  end if

  call aero_data_allocate(aero_data)
  call aero_state_allocate(aero_state)
  call gas_data_allocate(gas_data)
  call gas_state_allocate(gas_state)
  call env_state_allocate(env_state)

  allocate(filename_list(0))
  call input_filename_list(in_prefix, filename_list)
  n_file = size(filename_list)
  call assert_msg(323514871, n_file > 0, &
       "no NetCDF files found with prefix: " // trim(in_prefix))

  call input_state(filename_list(1), aero_data, aero_state, gas_data, &
       gas_state, env_state, index, time, del_t, i_repeat, uuid)
  run_uuid = uuid

  allocate(times(n_file))
  allocate(time_num_concs(n_file))
  allocate(time_mass_concs(n_file))
  allocate(time_species_concs(n_file, aero_data%n_spec))

  allocate(particle_num_concs(0))
  allocate(particle_masses(0))

  do i_file = 1,n_file
     call input_state(filename_list(i_file), aero_data, aero_state, gas_data, &
          gas_state, env_state, index, time, del_t, i_repeat, uuid)

     call assert_msg(397906326, uuid == run_uuid, &
          "UUID mismatch between " // trim(filename_list(1)) // " and " &
          // trim(filename_list(i_file)))

     times(i_file) = time
     call aero_state_num_concs(aero_state, particle_num_concs)
     time_num_concs(i_file) = sum(particle_num_concs)
     call aero_state_masses(aero_state, aero_data, particle_masses)
     time_mass_concs(i_file) = sum(particle_masses * particle_num_concs)
     do i_spec = 1,aero_data%n_spec
        call aero_state_masses(aero_state, aero_data, particle_masses, &
             include=(/aero_data%name(i_spec)/))
        time_species_concs(i_file, i_spec) &
             = sum(particle_masses * particle_num_concs)
     end do
  end do

  write(*,'(a,a)') "Output file: ", trim(out_filename)
  write(*,'(a)') "  Each row of output is one time."
  write(*,'(a)') "  The columns of output are:"
  write(*,'(a)') "    column  1: time (s)"
  write(*,'(a)') "    column  2: aerosol number concentration (#/m^3)"
  write(*,'(a)') "    column  3: aerosol mass concentration (kg/m^3)"
  do i_spec = 1,aero_data%n_spec
     write(*,'(a,i2,a,a,a)') "    column ", i_spec + 3, ": aerosol ", &
          trim(aero_data%name(i_spec)), " concentration (kg/m^3)"
  end do

  call open_file_write(out_filename, out_unit)
  do i_file = 1,n_file
     write(out_unit, '(e30.15e3)', advance='no') times(i_file)
     write(out_unit, '(e30.15e3)', advance='no') time_num_concs(i_file)
     write(out_unit, '(e30.15e3)', advance='no') time_mass_concs(i_file)
     do i_spec = 1,aero_data%n_spec
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
  deallocate(particle_num_concs)
  deallocate(particle_masses)
  call aero_data_deallocate(aero_data)
  call aero_state_deallocate(aero_state)
  call gas_data_deallocate(gas_data)
  call gas_state_deallocate(gas_state)
  call env_state_deallocate(env_state)

contains

  subroutine print_help()

    write(*,'(a)') 'Usage: extract_aero_time [options] <netcdf_prefix>'
    write(*,'(a)') ''
    write(*,'(a)') 'options are:'
    write(*,'(a)') '  -h, --help        Print this help message.'
    write(*,'(a)') '  -o, --out <file>  Output filename.'
    write(*,'(a)') ''
    write(*,'(a)') 'Examples:'
    write(*,'(a)') '  extract_aero_time data_0001'
    write(*,'(a)') ''

  end subroutine print_help

end program extract_aero_size