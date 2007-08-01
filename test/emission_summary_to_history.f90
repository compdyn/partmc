! Copyright (C) 2007 Nicole Riemer and Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.
!
! Process summary to make a time-history of the three non-empty bins.

program emission_summary_to_history

  use pmc_util
  use pmc_inout
  use pmc_bin_grid
  use pmc_aero_binned
  use pmc_gas_state
  use pmc_gas_data
  use pmc_aero_data
  use pmc_env
  
  type(inout_file_t) :: in_file
  integer :: f_out, bin_1, bin_2, bin_3
  integer :: n_loop, n_time, i, i_loop, i_time, i_bin, i_loop_check
  character(len=1000) :: name_in, name_out
  character(len=30) :: num_str
  type(bin_grid_t) :: bin_grid
  type(aero_data_t) :: aero_data
  type(gas_data_t) :: gas_data
 
  real*8, allocatable :: time(:,:), time_avg(:)
  type(env_t), allocatable :: env(:,:), env_avg(:)
  type(aero_binned_t), allocatable :: aero_binned(:,:), aero_binned_avg(:)
  type(gas_state_t), allocatable :: gas_state(:,:), gas_state_avg(:)
 
  ! check there is exactly one commandline argument
  if (iargc() .ne. 1) then
     write(6,*) 'Usage: emission_summary_to_history <filename.d>'
     call exit(2)
  endif
  
  ! get and check first commandline argument (must be "filename.d")
  call getarg(1, name_in)
  i = len_trim(name_in)
  if (i .gt. 40) then
     write(6,*) 'ERROR: filename too long'
     call exit(2)
  endif
  if ((name_in(i:i) .ne. 'd') .or. &
       (name_in((i-1):(i-1)) .ne. '.')) then
     write(6,*) 'ERROR: Filename must end in .d'
     call exit(2)
  endif
  
  ! compute names of output files
  name_out = name_in
  name_out((i-1):) = '_history.d'
  
  write(6,*) 'name_in = ', trim(name_in)
  write(6,*) 'name_out = ', trim(name_out)
  
  ! allocate unit numbers
  f_out = get_unit()

  ! open files
  call inout_open_read(name_in, in_file)
  open(f_out, file=name_out)
  
  ! read header
  call inout_read_integer(in_file, 'n_loop', n_loop)
  call inout_read_integer(in_file, 'n_time', n_time)
  call inout_read_bin_grid(in_file, bin_grid)
  call inout_read_gas_data(in_file, gas_data)
  call inout_read_aero_data(in_file, aero_data)

  allocate(time(n_loop,n_time), time_avg(n_time))
  allocate(env(n_loop,n_time), env_avg(n_time))
  allocate(aero_binned(n_loop,n_time), aero_binned_avg(n_time))
  allocate(gas_state(n_loop,n_time), gas_state_avg(n_time))

  write(*,*) 'n_loop = ', n_loop
  write(*,*) 'n_time = ', n_time
  
  ! read all data
  do i_loop = 1,n_loop
     do i_time = 1,n_time
        call inout_read_integer(in_file, 'loop_num', i_loop_check)
        call inout_check_index(in_file, i_loop, i_loop_check)
        call inout_read_real(in_file, 'time(s)', time(i_loop,i_time))
        call inout_read_env(in_file, env(i_loop,i_time))
        call inout_read_aero_binned(in_file, aero_binned(i_loop,i_time))
        call inout_read_gas_state(in_file, gas_state(i_loop,i_time))
     enddo
  enddo
  
  ! compute simple loop averages
  do i_time = 1,n_time
     call average_real(time(:,i_time), time_avg(i_time))
     call env_average(env(:,i_time), env_avg(i_time))
     call aero_binned_average(aero_binned(:,i_time), aero_binned_avg(i_time))
     call gas_state_average(gas_state(:,i_time), gas_state_avg(i_time))
  enddo

  ! find three non-empty bins
  bin_1 = 1
  do while (aero_binned_avg(2)%num_den(bin_1) == 0d0)
     bin_1 = bin_1 + 1
  end do
  bin_2 = bin_1 + 1
  do while (aero_binned_avg(2)%num_den(bin_2) == 0d0)
     bin_2 = bin_2 + 1
  end do
  bin_3 = bin_2 + 1
  do while (aero_binned_avg(2)%num_den(bin_3) == 0d0)
     bin_3 = bin_3 + 1
  end do

  ! output bin histories
  write(f_out, '(a1,7a20)') '#', 'time(s)', 'num_den_1(#/m^3)', &
       'num_den_2(#/m^3)', 'num_den_3(#/m^3)', &
       'vol_den_1(m^3/m^3)', 'vol_den_2(m^3/m^3)', 'vol_den_3(m^3/m^3)'
  do i_time = 1,n_time
     write(f_out, '(a1,7e20.10)') ' ', time_avg(i_time), &
          aero_binned_avg(i_time)%num_den(bin_1), &
          aero_binned_avg(i_time)%num_den(bin_2), &
          aero_binned_avg(i_time)%num_den(bin_3), &
          aero_binned_avg(i_time)%vol_den(bin_1,1), &
          aero_binned_avg(i_time)%vol_den(bin_2,1), &
          aero_binned_avg(i_time)%vol_den(bin_3,1)
  end do

  call inout_close(in_file)
  close(f_out)

end program emission_summary_to_history
