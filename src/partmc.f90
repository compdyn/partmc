! Copyright (C) 2007 Nicole Riemer and Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.
!
! Top level driver that reads the .spec file and calls simulation
! routines.

program partmc

  use pmc_inout

  type(inout_file_t) :: file
  character(len=300) :: in_name
  character(len=100) :: run_type
  integer :: i
  
  ! check there is exactly one commandline argument
  if (iargc() .ne. 1) then
     write(6,*) 'Usage: partmc <filename.d>'
     call exit(2)
  end if
  
  ! get and check first commandline argument (must be "filename.spec")
  call getarg(1, in_name)
  i = len_trim(in_name)
  if (in_name((i-4):i) /= '.spec') then
     write(6,*) 'ERROR: input filename must end in .spec'
     call exit(2)
  end if

  call inout_open_read(in_name, file)

  call inout_read_string(file, 'run_type', run_type)

  if (trim(run_type) == 'mc') then
     call partmc_mc(file)
  elseif (trim(run_type) == 'exact') then
     call partmc_exact(file)
  elseif (trim(run_type) == 'sect') then
     call partmc_sect(file)
  else
     write(0,*) 'ERROR: unknown run_type: ', trim(run_type)
     call exit(1)
  end if

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine partmc_mc(file)

    use pmc_bin_grid
    use pmc_aero_state
    use pmc_aero_dist
    use pmc_aero_binned
    use pmc_condensation
    use pmc_kernel_sedi
    use pmc_kernel_golovin
    use pmc_kernel_constant
    use pmc_kernel_brown
    use pmc_kernel_zero
    use pmc_aero_data
    use pmc_env
    use pmc_run_mc
    use pmc_inout
    use pmc_gas_data
    use pmc_gas_state
    use pmc_output_summary
    use pmc_util

    type(inout_file_t), intent(out) :: file ! spec file

    character(len=300) :: summary_name  ! name of output files
    character(len=100) :: kernel_name   ! coagulation kernel name
    type(gas_data_t) :: gas_data        ! gas data
    type(gas_state_t) :: gas_init       ! gas initial condition
    type(gas_state_t) :: gas_state      ! gas current state for run
    type(aero_data_t) :: aero_data      ! aero_data data
    type(aero_dist_t) :: aero_dist_init ! aerosol initial distribution
    type(aero_state_t) :: aero_init     ! aerosol initial condition
    type(aero_state_t) :: aero_state    ! aerosol current state for run
    type(env_t) :: env                ! environment data
    type(bin_grid_t) :: bin_grid        ! bin grid
    type(aero_binned_t) :: aero_binned  ! binned distributions
    type(run_mc_opt_t) :: mc_opt        ! Monte Carlo options
    type(inout_file_t) :: summary_file  ! summary output file

    integer :: i_loop, rand_init
    
    call inout_read_string(file, 'output_file', summary_name)
    call inout_read_string(file, 'state_prefix', mc_opt%state_prefix)
    call inout_read_integer(file, 'n_loop', mc_opt%n_loop)
    call inout_read_integer(file, 'n_part', mc_opt%n_part_max)
    call inout_read_string(file, 'kernel', kernel_name)
    
    call inout_read_real(file, 't_max', mc_opt%t_max)
    call inout_read_real(file, 'del_t', mc_opt%del_t)
    call inout_read_real(file, 't_output', mc_opt%t_output)
    call inout_read_real(file, 't_state', mc_opt%t_state)
    call inout_read_real(file, 't_progress', mc_opt%t_progress)

    call spec_read_bin_grid(file, bin_grid)
    
    call spec_read_gas_data(file, gas_data)
    call spec_read_gas_state(file, gas_data, 'gas_init', gas_init)

    call spec_read_aero_data_filename(file, aero_data)
    call spec_read_aero_dist_filename(file, aero_data, bin_grid, &
         'aerosol_init', aero_dist_init)
    call aero_dist_to_state(bin_grid, aero_data, aero_dist_init, &
         mc_opt%n_part_max, aero_init)
    call aero_dist_free(aero_dist_init)

    call spec_read_env(file, bin_grid, gas_data, aero_data, env)

    call inout_read_integer(file, 'rand_init', rand_init)
    call inout_read_logical(file, 'do_coagulation', mc_opt%do_coagulation)
    call inout_read_logical(file, 'allow_double', mc_opt%allow_double)
    call inout_read_logical(file, 'do_condensation', mc_opt%do_condensation)
    call inout_read_logical(file, 'do_mosaic', mc_opt%do_mosaic)
    call inout_read_logical(file, 'do_restart', mc_opt%do_restart)
    call inout_read_string(file, 'restart_name', mc_opt%restart_name)
    
    call inout_close(file)

    ! finished reading .spec data, now do the run

    call inout_open_write(summary_name, summary_file)
    call output_summary_header(summary_file, bin_grid, gas_data, &
         aero_data, mc_opt%n_loop, nint(mc_opt%t_max / mc_opt%t_output) + 1)
    
    call util_srand(rand_init)

    call aero_binned_alloc(aero_binned, bin_grid%n_bin, aero_data%n_spec)
    call gas_state_alloc(gas_state, gas_data%n_spec)
    call aero_state_alloc(bin_grid%n_bin, aero_data%n_spec, aero_state)
    call cpu_time(mc_opt%t_wall_start)

    do i_loop = 1,mc_opt%n_loop
       mc_opt%i_loop = i_loop
       
       call gas_state_copy(gas_init, gas_state)
       call aero_state_copy(aero_init, aero_state)

       if (mc_opt%do_condensation) then
          call aero_state_equilibriate(bin_grid, env, aero_data, aero_state)
       end if
       
       if (trim(kernel_name) == 'sedi') then
          call run_mc(kernel_sedi, bin_grid, aero_binned, env, aero_data, &
               aero_state, gas_data, gas_state, mc_opt, summary_file)
       elseif (trim(kernel_name) == 'golovin') then
          call run_mc(kernel_golovin, bin_grid, aero_binned, env, aero_data, &
               aero_state, gas_data, gas_state, mc_opt, summary_file)
       elseif (trim(kernel_name) == 'constant') then
          call run_mc(kernel_constant, bin_grid, aero_binned, env, aero_data, &
               aero_state, gas_data, gas_state, mc_opt, summary_file)
       elseif (trim(kernel_name) == 'brown') then
          call run_mc(kernel_brown, bin_grid, aero_binned, env, aero_data, &
               aero_state, gas_data, gas_state, mc_opt, summary_file)
       elseif (trim(kernel_name) == 'zero') then
          call run_mc(kernel_zero, bin_grid, aero_binned, env, aero_data, &
               aero_state, gas_data, gas_state, mc_opt, summary_file)
       else
          write(0,*) 'ERROR: Unknown kernel type; ', trim(kernel_name)
          call exit(1)
       end if
       
    end do

    call inout_close(summary_file)

    call gas_data_free(gas_data)
    call gas_state_free(gas_init)
    call gas_state_free(gas_state)
    call aero_data_free(aero_data)
    call aero_state_free(aero_init)
    call aero_state_free(aero_state)
    call env_free(env)
    call bin_grid_free(bin_grid)
    call aero_binned_free(aero_binned)

  end subroutine partmc_mc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine partmc_exact(file)

    use pmc_util
    use pmc_bin_grid
    use pmc_aero_state
    use pmc_aero_dist
    use pmc_condensation
    use pmc_kernel_golovin
    use pmc_kernel_constant
    use pmc_kernel_zero
    use pmc_aero_data
    use pmc_env
    use pmc_run_exact
    use pmc_inout
    use pmc_gas_data
    use pmc_gas_state
    use pmc_output_summary
    use pmc_kernel_zero

    type(inout_file_t), intent(out) :: file ! spec file

    character(len=300) :: summary_name  ! name of output files
    character(len=100) :: soln_name     ! exact solution name
    type(aero_data_t) :: aero_data      ! aero_data data
    type(env_t) :: env                ! environment data
    type(run_exact_opt_t) :: exact_opt  ! exact solution options
    type(bin_grid_t) :: bin_grid        ! bin grid
    type(inout_file_t) :: summary_file  ! summary output file
    type(gas_data_t) :: gas_data        ! dummy gas data
    
    call inout_read_string(file, 'output_file', summary_name)
    call inout_read_real(file, 'num_conc', exact_opt%num_conc)

    call inout_read_real(file, 't_max', exact_opt%t_max)
    call inout_read_real(file, 't_output', exact_opt%t_output)

    call spec_read_bin_grid(file, bin_grid)
    call spec_read_gas_data(file, gas_data)
    call spec_read_aero_data_filename(file, aero_data)
    call spec_read_env(file, bin_grid, gas_data, aero_data, env)

    call inout_read_string(file, 'soln', soln_name)

    call aero_dist_alloc(exact_opt%aero_dist_init, 0, 0, 0)
    if (trim(soln_name) == 'golovin_exp') then
       call inout_read_real(file, 'mean_vol', exact_opt%mean_vol)
    elseif (trim(soln_name) == 'constant_exp_cond') then
       call inout_read_real(file, 'mean_vol', exact_opt%mean_vol)
    elseif (trim(soln_name) == 'zero') then
       call aero_dist_free(exact_opt%aero_dist_init)
       call spec_read_aero_dist_filename(file, aero_data, bin_grid, &
            'aerosol_init', exact_opt%aero_dist_init)
    else
       write(0,*) 'ERROR: unknown solution type: ', trim(soln_name)
       call exit(1)
    end if
    
    call inout_close(file)

    ! finished reading .spec data, now do the run

    call inout_open_write(summary_name, summary_file)
    call gas_data_alloc(gas_data, 0)
    call output_summary_header(summary_file, bin_grid, gas_data, &
         aero_data, 1, nint(exact_opt%t_max / exact_opt%t_output) + 1)

    if (trim(soln_name) == 'golovin_exp') then
       call run_exact(bin_grid, env, aero_data, exact_opt, &
            soln_golovin_exp, summary_file)
    elseif (trim(soln_name) == 'golovin_exp') then
       call run_exact(bin_grid, env, aero_data, exact_opt, &
            soln_constant_exp_cond, summary_file)
    elseif (trim(soln_name) == 'zero') then
       call run_exact(bin_grid, env, aero_data, exact_opt, &
            soln_zero, summary_file)
    else
       write(0,*) 'ERROR: unknown solution type: ', trim(soln_name)
       call exit(1)
    end if

    call inout_close(summary_file)

    call aero_data_free(aero_data)
    call env_free(env)
    call bin_grid_free(bin_grid)
    call gas_data_free(gas_data)
    call aero_dist_free(exact_opt%aero_dist_init)
    
  end subroutine partmc_exact

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine partmc_sect(file)

    use pmc_util
    use pmc_aero_data
    use pmc_env
    use pmc_run_sect
    use pmc_kernel_sedi
    use pmc_kernel_golovin
    use pmc_kernel_constant
    use pmc_kernel_brown
    use pmc_kernel_zero
    use pmc_bin_grid
    use pmc_aero_state
    use pmc_aero_dist
    use pmc_gas_data
    use pmc_gas_state
    use pmc_output_summary

    type(inout_file_t), intent(out) :: file ! spec file

    character(len=300) :: summary_name  ! name of output files
    character(len=100) :: kernel_name   ! coagulation kernel name
    type(run_sect_opt_t) :: sect_opt    ! sectional code options
    type(aero_data_t) :: aero_data      ! aero_data data
    type(aero_dist_t) :: aero_dist_init ! aerosol initial distribution
    type(aero_state_t) :: aero_init     ! aerosol initial condition
    type(env_t) :: env                ! environment data
    type(bin_grid_t) :: bin_grid        ! bin grid
    type(inout_file_t) :: summary_file  ! summary output file
    type(gas_data_t) :: gas_data        ! dummy gas data

    call inout_read_string(file, 'output_file', summary_name)
    call inout_read_string(file, 'kernel', kernel_name)

    call inout_read_real(file, 't_max', sect_opt%t_max)
    call inout_read_real(file, 'del_t', sect_opt%del_t)
    call inout_read_real(file, 't_output', sect_opt%t_output)
    call inout_read_real(file, 't_progress', sect_opt%t_progress)

    call spec_read_bin_grid(file, bin_grid)

    call spec_read_gas_data(file, gas_data)

    call spec_read_aero_data_filename(file, aero_data)
    call spec_read_aero_dist_filename(file, aero_data, bin_grid, &
         'aerosol_init', aero_dist_init)

    call spec_read_env(file, bin_grid, gas_data, aero_data, env)

    call inout_read_logical(file, 'do_coagulation', sect_opt%do_coagulation)
    
    call inout_close(file)

    ! finished reading .spec data, now do the run

    call inout_open_write(summary_name, summary_file)
    call output_summary_header(summary_file, bin_grid, gas_data, &
         aero_data, 1, nint(sect_opt%t_max / sect_opt%t_output) + 1)

    if (trim(kernel_name) == 'sedi') then
       call run_sect(bin_grid, gas_data, aero_data, aero_dist_init, env, &
            kernel_sedi, sect_opt, summary_file)
    elseif (trim(kernel_name) == 'golovin') then
       call run_sect(bin_grid, gas_data, aero_data, aero_dist_init, env, &
            kernel_golovin, sect_opt, summary_file)
    elseif (trim(kernel_name) == 'constant') then
       call run_sect(bin_grid, gas_data, aero_data, aero_dist_init, env, &
            kernel_constant, sect_opt, summary_file)
    elseif (trim(kernel_name) == 'brown') then
       call run_sect(bin_grid, gas_data, aero_data, aero_dist_init, env, &
            kernel_brown, sect_opt, summary_file)
    elseif (trim(kernel_name) == 'zero') then
       call run_sect(bin_grid, gas_data, aero_data, aero_dist_init, env, &
            kernel_zero, sect_opt, summary_file)
    else
       write(0,*) 'ERROR: Unknown kernel type; ', trim(kernel_name)
       call exit(1)
    end if

    call inout_close(summary_file)

    call aero_data_free(aero_data)
    call aero_dist_free(aero_dist_init)
    call env_free(env)
    call bin_grid_free(bin_grid)
    call gas_data_free(gas_data)
    
  end subroutine partmc_sect

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program partmc
