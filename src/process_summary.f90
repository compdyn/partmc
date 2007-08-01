! Copyright (C) 2005-2007 Nicole Riemer and Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.
!
! Process output data files.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program process_summary

  use pmc_util
  use pmc_inout
  use pmc_bin_grid
  use pmc_aero_binned
  use pmc_gas_state
  use pmc_gas_data
  use pmc_aero_data
  use pmc_env
  
  type(inout_file_t) :: in_file
  integer :: f_out_env, f_out_gas, f_out_aero_binned, f_out_aero_total
  integer :: n_loop, n_time, i, i_loop, i_time, i_bin, i_loop_check
  character(len=1000) :: name_in, name_out_env, name_out_gas
  character(len=1000) :: name_out_aero_binned, name_out_aero_total
  character(len=30) :: num_str
  type(bin_grid_t) :: bin_grid
  type(aero_data_t) :: aero_data
  type(gas_data_t) :: gas_data
 
  real*8, allocatable :: time(:,:), time_avg(:)
  type(env_t), allocatable :: env(:,:), env_avg(:)
  type(aero_binned_t), allocatable :: aero_binned(:,:), aero_binned_avg(:)
  type(gas_state_t), allocatable :: gas_state(:,:), gas_state_avg(:)
  real*8, allocatable :: tot_num_den(:), tot_vol_den(:,:)
 
  ! check there is exactly one commandline argument
  if (iargc() .ne. 1) then
     write(6,*) 'Usage: process_summary <filename.d>'
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
  name_out_env = name_in
  name_out_gas = name_in
  name_out_aero_binned = name_in
  name_out_aero_total = name_in
  name_out_env((i-1):) = '_env.d'
  name_out_gas((i-1):) = '_gas.d'
  name_out_aero_binned((i-1):) = '_aero_binned.d'
  name_out_aero_total((i-1):) = '_aero_total.d'
  
  write(6,*) 'name_in = ', trim(name_in)
  write(6,*) 'name_out_env = ', trim(name_out_env)
  write(6,*) 'name_out_gas = ', trim(name_out_gas)
  write(6,*) 'name_out_aero_binned = ', trim(name_out_aero_binned)
  write(6,*) 'name_out_aero_total = ', trim(name_out_aero_total)
  
  ! allocate unit numbers
  f_out_env = get_unit()
  f_out_gas = get_unit()
  f_out_aero_binned = get_unit()
  f_out_aero_total = get_unit()

  ! open files
  call inout_open_read(name_in, in_file)
  open(f_out_env, file=name_out_env)
  open(f_out_gas, file=name_out_gas)
  open(f_out_aero_binned, file=name_out_aero_binned)
  open(f_out_aero_total, file=name_out_aero_total)
  
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
  allocate(tot_num_den(n_time), tot_vol_den(n_time,aero_data%n_spec))

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
  
  ! compute loop averages
  do i_time = 1,n_time
     call average_real(time(:,i_time), time_avg(i_time))
     call env_average(env(:,i_time), env_avg(i_time))
     call aero_binned_average(aero_binned(:,i_time), aero_binned_avg(i_time))
     call average_gas_state(gas_state(:,i_time), gas_state_avg(i_time))
     tot_num_den(i_time) = sum(aero_binned_avg(i_time)%num_den)
     tot_vol_den(i_time,:) = sum(aero_binned_avg(i_time)%vol_den, 1)
  enddo

  ! output environment
  write(f_out_env, '(a1,5a20)') '#', 'time(s)', &
       'temp(K)', 'rel_hum(1)', 'press(Pa)', &
       'dens(kg/m^3)'
  do i_time = 1,n_time
     write(f_out_env,'(a1,5e20.10)') ' ', time_avg(i_time), &
          env_avg(i_time)%temp, env_avg(i_time)%rel_humid, &
          env_avg(i_time)%pressure, env_avg(i_time)%air_den
  end do
  
  ! output gas
  if (gas_data%n_spec > 0) then
     write(num_str, '(i10)') gas_data%n_spec
     write(f_out_gas, '(a1,a20,'//num_str//'a20)') '#', &
          'time(s)', gas_data%name
     ! FIXME: can we add (ppb) to the end of each name?
     do i_time = 1,n_time
        write(f_out_gas, '(a1,e20.10,'//num_str//'e20.10)') ' ', &
             time_avg(i_time), gas_state_avg(i_time)%conc
     end do
  end if

  ! output aerosol binned
  write(num_str, '(i10)') aero_data%n_spec
  do i_time = 1,n_time
     write(f_out_aero_binned, '(a2,a10,e20.10)') '# ', 'time(s) = ', &
          time_avg(i_time)
     write(f_out_aero_binned, '(a2,a)') '# ', &
          'species are volume density (m^3/m^3)'
     write(f_out_aero_binned, '(a1,a20,a20,'//num_str//'a20)') '#', &
          'radius(m)', 'num_dens(#/m^3)', aero_data%name
     do i_bin = 1,bin_grid%n_bin
        write(f_out_aero_binned, '(a1,e20.10,e20.10,'//num_str//'e20.10)') &
             ' ', vol2rad(bin_grid%v(i_bin)), &
             aero_binned_avg(i_time)%num_den(i_bin), &
             aero_binned_avg(i_time)%vol_den(i_bin,:)
     end do
     write(f_out_aero_binned, '(a)') ''
     write(f_out_aero_binned, '(a)') ''
  end do

  ! output aerosol totals
  write(num_str, '(i10)') aero_data%n_spec
  write(f_out_aero_total, '(a1,a20,a20,'//num_str//'a20)') '#', &
       'time(s)', 'num_dens(#/m^3)', aero_data%name
  ! FIXME: can we add (m^3/m^3) to the end of each name?
  do i_time = 1,n_time
     write(f_out_aero_total, '(a1,e20.10,e20.10,'//num_str//'e20.10)') ' ', &
          time_avg(i_time), tot_num_den(i_time), tot_vol_den(i_time,:)
  end do

  call inout_close(in_file)
  close(f_out_env)
  close(f_out_gas)
  close(f_out_aero_binned)
  close(f_out_aero_total)
  
end program process_summary

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
