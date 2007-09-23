! Copyright (C) 2005-2007 Nicole Riemer and Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.
!
! Process the saved state files to obtain summary data.

program process_state_new

  use pmc_inout
  use pmc_bin_grid
  use pmc_aero_data
  use pmc_aero_state
  use pmc_gas_data
  use pmc_gas_state
  use pmc_env
  use pmc_output_state
  use pmc_process_state_hist
  use pmc_mosaic
  use pmc_mpi
  use pmc_process_spec
  use pmc_util

  character(len=300) :: filename
  type(process_spec_t), pointer :: process_spec_list(:)
  integer :: i
  type(inout_file_t) :: file

  call pmc_mpi_init()

  if (iargc() < 2) then
     call print_usage()
     call exit(1)
  end if

  call getarg(1, filename)
  call inout_open_read(filename, file)
  call inout_read_process_spec_list(file, process_spec_list)
  call inout_close(file)

  do i = 2,iargc()
     call getarg(i, filename)
     call process_state_file(filename, process_spec_list)
  end do

  do i = 1,size(process_spec_list)
     call process_spec_free(process_spec_list(i))
  end do
  deallocate(process_spec_list)

  call pmc_mpi_finalize()

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine print_usage()

    write(*,*) 'Usage: process_state <process.spec> <state_filenames...>'
    
  end subroutine print_usage

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine process_state_file(filename, process_spec_list)

    character(len=*), intent(in) :: filename ! state filename
    type(process_spec_t), intent(in) :: process_spec_list(:) ! process specs

    character(len=len(filename)) :: basename
    type(bin_grid_t) :: bin_grid
    type(aero_data_t) :: aero_data
    type(aero_state_t) :: aero_state
    type(gas_data_t) :: gas_data
    type(gas_state_t) :: gas_state
    type(env_t) :: env
    real*8 :: time
    integer :: index, i

    call get_basename(filename, basename)

    call inout_read_state(filename, bin_grid, aero_data, aero_state, &
         gas_data, gas_state, env, time, index)
    write(*,'(a,e20.10)') 'time (s) = ', time

    call process_state_spec_list(basename, process_spec_list, &
         bin_grid, aero_data, aero_state, gas_data, gas_state, &
         env, time, index)

  end subroutine process_state_file

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
end program process_state_new
