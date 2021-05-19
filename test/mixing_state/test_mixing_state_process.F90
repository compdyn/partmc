! Copyright (C) 2021 Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The mixing state process program.

!> Read NetCDF output files and process them.
program test_mixing_state_process

  use pmc_output

  character(len=PMC_MAX_FILENAME_LEN), parameter :: prefix &
       = "out/mixing_state"

  character(len=PMC_MAX_FILENAME_LEN) :: in_filename, out_filename
  type(aero_data_t) :: aero_data
  type(aero_state_t) :: aero_state
  type(env_state_t) :: env_state
  integer :: ncid, index, repeat, i_index, i_repeat, unit
  real(kind=dp) :: time, del_t
  real(kind=dp) :: d_alpha, d_gamma, chi
  character(len=PMC_UUID_LEN) :: uuid
  character(len=AERO_NAME_LEN), allocatable :: mixing_state_groups(:,:)

  call pmc_mpi_init()

  allocate(mixing_state_groups(2, 2))
  mixing_state_groups(1,:) = ["A", " "]
  mixing_state_groups(2,:) = ["C", "D"]

  i_index = 1
  i_repeat = 1
  call make_filename(in_filename, prefix, ".nc", i_index, i_repeat)
  call input_state(in_filename, index, time, del_t, repeat, &
       uuid, aero_data=aero_data, aero_state=aero_state, &
       env_state=env_state)

  write(*,'(a20,a24,a24,a24)') 'name', 'd_alpha    ', 'd_gamma    ', 'chi    '

  call aero_state_mixing_state_metrics(aero_state, aero_data, &
       d_alpha, d_gamma, chi)
  write(*,'(a20,g24.16,g24.16,g24.16)') 'all_species', &
       d_alpha, d_gamma, chi
  call open_file_write('out/test_all_species.txt', unit)
  write(unit,'(g30.20,g30.20,g30.20)') d_alpha, d_gamma, chi
  call close_file(unit)

  call aero_state_mixing_state_metrics(aero_state, aero_data, &
       d_alpha, d_gamma, chi, groups=mixing_state_groups)
  write(*,'(a20,g24.16,g24.16,g24.16)') 'groups', &
       d_alpha, d_gamma, chi
  call open_file_write('out/test_groups.txt', unit)
  write(unit,'(g30.20,g30.20,g30.20)') d_alpha, d_gamma, chi
  call close_file(unit)

  call aero_state_mixing_state_metrics(aero_state, aero_data, &
       d_alpha, d_gamma, chi, group=["A", "C"])
  write(*,'(a20,g24.16,g24.16,g24.16)') 'groupAC', &
       d_alpha, d_gamma, chi
  call open_file_write('out/test_groupAC.txt', unit)
  write(unit,'(g30.20,g30.20,g30.20)') d_alpha, d_gamma, chi
  call close_file(unit)

  call aero_state_mixing_state_metrics(aero_state, aero_data, &
       d_alpha, d_gamma, chi, group=["B", "C"], exclude=["D"])
  write(*,'(a20,g24.16,g24.16,g24.16)') 'group_exclude', &
       d_alpha, d_gamma, chi
  call open_file_write('out/test_group_exclude.txt', unit)
  write(unit,'(g30.20,g30.20,g30.20)') d_alpha, d_gamma, chi
  call close_file(unit)

  call aero_state_mixing_state_metrics(aero_state, aero_data, &
       d_alpha, d_gamma, chi, group=["A", "D"], include=["B", "C", "D"])
  write(*,'(a20,g24.16,g24.16,g24.16)') 'group_include', &
       d_alpha, d_gamma, chi
  call open_file_write('out/test_group_include.txt', unit)
  write(unit,'(g30.20,g30.20,g30.20)') d_alpha, d_gamma, chi
  call close_file(unit)

  call pmc_mpi_finalize()

end program test_mixing_state_process
