!! Copyright (C) 2009-2012 Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The extract_freezing program.

!> Read NetCDF output files and write out the time evolution of
!> frozen fraction average among ensembles in text format.
program extract_freezing

  use pmc_aero_state
  use pmc_aero_particle
  use pmc_output
  use pmc_mpi
  use getopt_m

  integer :: n_ensemble
  character(len=PMC_MAX_FILENAME_LEN) :: in_prefix, in_prefix_ensemble, out_filename
  character(len=PMC_MAX_FILENAME_LEN), allocatable :: filename_list(:)
  character(len=1000) :: tmp_str
  type(aero_data_t) :: aero_data
  type(aero_state_t) :: aero_state
  type(env_state_t) :: env_state
  integer :: index, i_repeat, i_spec, out_unit
  integer :: i_file, n_file, i_ensemble
  real(kind=dp) :: time, del_t
  character(len=PMC_UUID_LEN) :: uuid, run_uuid
  real(kind=dp), allocatable ::  frozen_fractions(:),&
       frozen_fractions_mean(:)
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
     out_filename = trim(in_prefix) // "_frozen_fraction_ensemble_mean.txt"
  end if

 
  call input_n_files(in_prefix, n_ensemble, n_file)
 
  allocate(filename_list(0))

  write(in_prefix_ensemble, "(a,a,i4.4)") trim(in_prefix), "_", 1
  call input_filename_list(in_prefix_ensemble, filename_list)

  allocate(frozen_fractions(n_file))
  allocate(frozen_fractions_mean(n_file))
  frozen_fractions_mean = 0d0

  do i_ensemble = 1,N_ensemble
    write(in_prefix_ensemble, "(a,a,i4.4)") trim(in_prefix), '_', i_ensemble
    call input_filename_list(in_prefix_ensemble, filename_list)
    n_file = size(filename_list)
    call assert_msg(323514871, n_file > 0, &
         "no NetCDF files found with prefix: " // trim(in_prefix_ensemble))

    call input_state(filename_list(1), index, time, del_t, i_repeat, uuid, &
         aero_data=aero_data, aero_state=aero_state, env_state=env_state)
    run_uuid = uuid

  

    do i_file = 1,n_file
       call input_state(filename_list(i_file), index, time, del_t, i_repeat, &
            uuid, aero_data=aero_data, aero_state=aero_state, env_state=env_state)

       call assert_msg(397906326, uuid == run_uuid, &
            "UUID mismatch between " // trim(filename_list(1)) // " and " &
            // trim(filename_list(i_file)))

       frozen_fractions(i_file) = aero_state_frozen_fraction(aero_state, aero_data)
    end do
    frozen_fractions_mean = frozen_fractions_mean + frozen_fractions

  end do
  frozen_fractions_mean = frozen_fractions_mean / N_ensemble

  write(*,'(a,a)') "Output file: ", trim(out_filename)
  write(*,'(a)') "  Each row of output is one time."
  write(*,'(a)') "  The columns of output are:"
  write(*,'(a)') "    column  1: frozen fraction ensemble mean (unitless)"
  call open_file_write(out_filename, out_unit)
  do i_file = 1,n_file
     write(out_unit, '(e30.15e3)', advance='no') frozen_fractions_mean(i_file)
     write(out_unit, '(a)') ''
  end do
  call close_file(out_unit)

  deallocate(frozen_fractions)
  deallocate(frozen_fractions_mean)
  deallocate(filename_list)

  call pmc_mpi_finalize()

contains

  subroutine print_help()

    write(*,'(a)') 'Usage: extract_freezing [options] <netcdf_prefix>'
    write(*,'(a)') ''
    write(*,'(a)') 'options are:'
    write(*,'(a)') '  -h, --help        Print this help message.'
    write(*,'(a)') '  -o, --out <file>  Output filename.'
    write(*,'(a)') ''
    write(*,'(a)') 'Examples:'
    write(*,'(a)') '  extract_freezing freezing_part'
    write(*,'(a)') ''

  end subroutine print_help

end program extract_freezing


