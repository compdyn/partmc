! Copyright (C) 2007 Nicole Riemer and Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.
!
! Top level driver that reads the .spec file and calls simulation
! routines.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program partmc

  use mod_inout

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

    use mod_bin
    use mod_aero_state
    use mod_aero_dist
    use mod_condensation
    use mod_kernel_sedi
    use mod_kernel_golovin
    use mod_kernel_constant
    use mod_kernel_brown
    use mod_aero_data
    use mod_environ
    use mod_run_mc
    use mod_inout
    use mod_gas_data
    use mod_gas_state
    use mod_output_summary
    use mod_util

    type(inout_file_t), intent(out) :: file ! spec file

    character(len=300) :: output_file   ! name of output files
    character(len=300) :: state_prefix  ! prefix for state files
    character(len=100) :: kernel_name   ! coagulation kernel name
    
    type(gas_data_t) :: gas_data        ! gas data
    type(gas_state_t) :: gas_init       ! gas initial condition
    type(gas_state_t) :: gas_state      ! gas current state for run

    type(aero_data_t) :: aero_data      ! aero_data data
    type(aero_dist_t) :: aero_init_dist ! aerosol initial distribution
    type(aero_state_t) :: aero_init     ! aerosol initial condition
    type(aero_state_t) :: aero_state    ! aerosol current state for run

    type(environ) :: env                ! environment data
    type(bin_grid_t) :: bin_grid        ! bin grid
    type(bin_dist_t) :: bin_dist        ! binned distributions
    type(mc_opt_t) :: mc_opt            ! Monte Carlo options

    integer :: i_loop, rand_init
    
    call inout_read_string(file, 'output_file', output_file)
    call inout_read_string(file, 'state_prefix', mc_opt%state_prefix)
    call inout_read_integer(file, 'n_loop', mc_opt%n_loop)
    call inout_read_integer(file, 'n_part', mc_opt%n_part_max)
    call inout_read_string(file, 'kernel', kernel_name)
    
    call inout_read_real(file, 't_max', mc_opt%t_max)
    call inout_read_real(file, 'del_t', mc_opt%del_t)
    call inout_read_real(file, 't_output', mc_opt%t_output)
    call inout_read_real(file, 't_state', mc_opt%t_state)
    call inout_read_real(file, 't_progress', mc_opt%t_progress)

    call spec_read_environ(file, env)

    call spec_read_bin_grid(file, bin_grid)
    
    call spec_read_gas_data(file, gas_data)
    call spec_read_gas_state(file, gas_data, 'gas_init', gas_init)
    call spec_read_gas_state(file, gas_data, 'gas_emissions', env%gas_emissions)
    call inout_read_real(file, 'gas_emission_rate', env%gas_emission_rate)
    call spec_read_gas_state(file, gas_data, 'gas_background', &
         env%gas_background)
    call inout_read_real(file, 'gas_dilution_rate', env%gas_dilution_rate)

    call spec_read_aero_data_filename(file, aero_data)
    call spec_read_aero_dist_filename(file, aero_data, bin_grid, &
         'aerosol_init', aero_init_dist)
    call aero_dist_to_part(bin_grid, aero_data, aero_init_dist, &
         mc_opt%n_part_max, aero_init)
    call spec_read_aero_dist_filename(file, aero_data, bin_grid, &
         'aerosol_emissions', env%aero_emissions)
    call inout_read_real(file, 'aerosol_emission_rate', env%aero_emission_rate)
    call spec_read_aero_dist_filename(file, aero_data, bin_grid, &
         'aerosol_background', env%aero_background)
    call inout_read_real(file, 'aerosol_dilution_rate', env%aero_dilution_rate)

    call inout_read_integer(file, 'rand_init', rand_init)
    call inout_read_logical(file, 'do_coagulation', mc_opt%do_coagulation)
    call inout_read_logical(file, 'allow_double', mc_opt%allow_double)
    call inout_read_logical(file, 'do_condensation', mc_opt%do_condensation)
    call inout_read_logical(file, 'do_mosaic', mc_opt%do_mosaic)
    call inout_read_logical(file, 'do_restart', mc_opt%do_restart)
    call inout_read_string(file, 'restart_name', mc_opt%restart_name)
    
    call inout_close(file)

    ! finished reading .spec data, now do the run
    
    mc_opt%output_unit = get_unit()
    mc_opt%state_unit = get_unit()

    call output_summary_open(mc_opt%output_unit, output_file, &
         mc_opt%n_loop, bin_grid%n_bin, aero_data%n_spec, &
         nint(mc_opt%t_max / mc_opt%t_output) + 1)
    
    if (rand_init /= 0) then
       call srand(rand_init)
    else
       call srand(time())
    end if

    call allocate_bin_dist(bin_grid%n_bin, aero_data%n_spec, bin_dist)
    call allocate_gas_state(gas_data, gas_state)
    call allocate_aero_state(bin_grid%n_bin, aero_data%n_spec, aero_state)
    call cpu_time(mc_opt%t_wall_start)

    do i_loop = 1,mc_opt%n_loop
       mc_opt%i_loop = i_loop
       
       call copy_gas_state(gas_init, gas_state)
       call copy_aero_state(aero_init, aero_state)

       if (mc_opt%do_condensation) then
          call equilibriate_aero(bin_grid, env, aero_data, aero_state)
       end if
       
       if (trim(kernel_name) == 'sedi') then
          call run_mc(kernel_sedi, bin_grid, bin_dist, env, aero_data, &
               aero_state, gas_data, gas_state, mc_opt)
       elseif (trim(kernel_name) == 'golovin') then
          call run_mc(kernel_golovin, bin_grid, bin_dist, env, aero_data, &
               aero_state, gas_data, gas_state, mc_opt)
       elseif (trim(kernel_name) == 'constant') then
          call run_mc(kernel_constant, bin_grid, bin_dist, env, aero_data, &
               aero_state, gas_data, gas_state, mc_opt)
       elseif (trim(kernel_name) == 'brown') then
          call run_mc(kernel_brown, bin_grid, bin_dist, env, aero_data, &
               aero_state, gas_data, gas_state, mc_opt)
       else
          write(0,*) 'ERROR: Unknown kernel type; ', trim(kernel_name)
          call exit(1)
       end if
       
    end do

  end subroutine partmc_mc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine partmc_exact(file)

    use mod_bin
    use mod_aero_state
    use mod_aero_dist
    use mod_condensation
    use mod_kernel_golovin
    use mod_kernel_constant
    use mod_aero_data
    use mod_environ
    use mod_run_exact
    use mod_inout
    use mod_gas_data
    use mod_gas_state
    use mod_output_summary

    type(inout_file_t), intent(out) :: file     ! spec file

    integer, parameter :: output_unit = 32

    real*8, allocatable :: bin_v(:)
    real*8, allocatable :: bin_g_den(:), bin_gs_den(:,:), bin_n_den(:)
    real*8 :: dlnr
    
    character(len=300) :: output_file   ! name of output files
    real*8 :: num_conc                  ! particle concentration (#/m^3)

    character(len=100) :: soln_name     ! exact solution name
    real*8 :: mean_vol                  ! mean volume of initial distribution
    
    real*8 :: t_max                     ! total simulation time (s)
    real*8 :: t_output                  ! output interval (0 disables) (s)
    
    type(aero_data_t) :: aero_data               ! aero_data data
    type(environ) :: env                ! environment data

    integer :: n_bin                    ! number of bins
    real*8 :: v_min                     ! volume of smallest bin (m^3)
    integer :: scal                     ! scale factor (integer)
    type(bin_grid_t) :: bin_grid        ! bin grid
    
    call inout_read_string(file, 'output_file', output_file)
    call inout_read_real(file, 'num_conc', num_conc)

    call inout_read_string(file, 'soln', soln_name)

    if (trim(soln_name) == 'golovin_exp') then
       call inout_read_real(file, 'mean_vol', mean_vol)
    elseif (trim(soln_name) == 'constant_exp_cond') then
       call inout_read_real(file, 'mean_vol', mean_vol)
    else
       write(0,*) 'ERROR: unknown solution type: ', trim(soln_name)
       call exit(1)
    end if
    
    call inout_read_real(file, 't_max', t_max)
    call inout_read_real(file, 't_output', t_output)
    
    call spec_read_environ(file, env)
    call spec_read_aero_data(file, aero_data)

    call inout_read_integer(file, 'n_bin', n_bin)
    call inout_read_real(file, 'v_min', v_min)
    call inout_read_integer(file, 'scal', scal)

    call inout_close(file)

    ! finished reading .spec data, now do the run
    
    call output_summary_open(output_unit, output_file, 1, n_bin, &
         aero_data%n_spec, nint(t_max / t_output) + 1)

    allocate(bin_v(n_bin), bin_g_den(n_bin), bin_n_den(n_bin))
    allocate(bin_gs_den(n_bin,aero_data%n_spec))
    
    call make_bin_grid(n_bin, scal, v_min, bin_grid)
    ! FIXME: delete following
    bin_v = bin_grid%v
    dlnr = bin_grid%dlnr
    
    if (trim(soln_name) == 'golovin_exp') then
       call run_exact(n_bin, aero_data%n_spec, bin_v, bin_g_den, bin_gs_den, &
            bin_n_den, num_conc, mean_vol, aero_data%rho(1), soln_golovin_exp, &
            t_max, t_output, output_unit, env, aero_data)
    elseif (trim(soln_name) == 'golovin_exp') then
       call run_exact(n_bin, aero_data%n_spec, bin_v, bin_g_den, bin_gs_den, &
            bin_n_den, num_conc, mean_vol, aero_data%rho(1), soln_constant_exp_cond, &
            t_max, t_output, output_unit, env, aero_data)
    else
       write(0,*) 'ERROR: unknown solution type: ', trim(soln_name)
       call exit(1)
    end if
    
  end subroutine partmc_exact

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine partmc_sect(file)

    use mod_aero_data
    use mod_environ
    use mod_run_sect
    use mod_kernel_sedi
    use mod_kernel_golovin
    use mod_kernel_constant
    use mod_kernel_brown
    use mod_bin
    use mod_aero_state
    use mod_aero_dist
    use mod_gas_data
    use mod_gas_state
    use mod_output_summary

    type(inout_file_t), intent(out) :: file     ! spec file

    integer, parameter :: max_dist_args = 10
    integer, parameter :: output_unit = 32

    character(len=300) :: out_file_name
    integer, allocatable :: bin_n(:)
    real*8, allocatable :: bin_v(:), n_den(:), bin_g(:), bin_gs(:,:)
    real*8 :: dlnr
    
    character(len=300) :: output_file   ! name of output files
    real*8 :: num_conc                  ! particle concentration (#/m^3)
    character(len=100) :: kernel_name   ! coagulation kernel name
    
    real*8 :: t_max                     ! total simulation time (s)
    real*8 :: del_t                     ! timestep (s)
    real*8 :: t_output                  ! output interval (0 disables) (s)
    real*8 :: t_progress                ! progress interval (0 disables) (s)
    
    type(aero_data_t) :: aero_data               ! aero_data data
    type(aero_dist_t) :: aero_init_dist ! aerosol initial distribution
    type(aero_state_t) :: aero_init          ! aerosol initial condition
    type(environ) :: env                ! environment data
    
    character(len=100) :: dist_type     ! initial distribution
    real*8 :: dist_args(max_dist_args)  ! distribution arguments
    
    integer :: n_bin                    ! number of bins
    real*8 :: v_min                     ! volume of smallest bin (m^3)
    integer :: scal                     ! scale factor (integer)
    type(bin_grid_t) :: bin_grid        ! bin grid

    call inout_read_string(file, 'output_file', output_file)
    call inout_read_real(file, 'num_conc', num_conc)

    call inout_read_string(file, 'kernel', kernel_name)

    call inout_read_real(file, 't_max', t_max)
    call inout_read_real(file, 'del_t', del_t)
    call inout_read_real(file, 't_output', t_output)
    call inout_read_real(file, 't_progress', t_progress)

    call inout_read_integer(file, 'n_bin', n_bin)
    call inout_read_real(file, 'v_min', v_min)
    call inout_read_integer(file, 'scal', scal)
    allocate(bin_v(n_bin), n_den(n_bin))
    call make_bin_grid(n_bin, scal, v_min, bin_grid)
    ! FIXME: eventually delete following
    bin_v = bin_grid%v
    dlnr = bin_grid%dlnr
    
    call spec_read_environ(file, env)
    call spec_read_aero_data(file, aero_data)

    allocate(bin_g(n_bin), bin_gs(n_bin,aero_data%n_spec), bin_n(n_bin))

    call spec_read_aero_dist_filename(file, aero_data, bin_grid, &
         'aerosol_init', aero_init_dist)
    call dist_total_n_den(bin_grid, aero_data, aero_init_dist, n_den)
    
    call inout_close(file)

    ! finished reading .spec data, now do the run

    call output_summary_open(output_unit, output_file, 1, n_bin, &
         aero_data%n_spec, nint(t_max / t_output) + 1)
    
    if (trim(kernel_name) == 'sedi') then
       call run_sect(n_bin, bin_v, dlnr, n_den, kernel_sedi, &
            t_max, del_t, t_output, t_progress, output_unit, aero_data, env)
    elseif (trim(kernel_name) == 'golovin') then
       call run_sect(n_bin, bin_v, dlnr, n_den, kernel_golovin, &
            t_max, del_t, t_output, t_progress, output_unit, aero_data, env)
    elseif (trim(kernel_name) == 'constant') then
       call run_sect(n_bin, bin_v, dlnr, n_den, kernel_constant, &
            t_max, del_t, t_output, t_progress, output_unit, aero_data, env)
    elseif (trim(kernel_name) == 'brown') then
       call run_sect(n_bin, bin_v, dlnr, n_den, kernel_brown, &
            t_max, del_t, t_output, t_progress, output_unit, aero_data, env)
    else
       write(0,*) 'ERROR: Unknown kernel type; ', trim(kernel_name)
       call exit(1)
    end if
    
  end subroutine partmc_sect

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program partmc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
