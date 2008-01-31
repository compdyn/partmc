! Copyright (C) 2007, 2008 Nicole Riemer and Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.
!
! Top level driver that reads the .spec file and calls simulation
! routines.

program partmc

  use pmc_mpi
  use pmc_process_spec
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
  use pmc_env_data
  use pmc_env_state
  use pmc_run_mc
  use pmc_run_exact
  use pmc_run_sect
  use pmc_inout
  use pmc_gas_data
  use pmc_gas_state
  use pmc_util

  character(len=300) :: spec_name, tmp, process_name, nc_name, action
  
  call pmc_mpi_init()
  if (pmc_mpi_rank() == 0) then
     ! only the root process accesses the commandline

     if (iargc() == 1) then
        call getarg(1, spec_name)
        action = 'run'
     elseif (iargc() >= 4) then
        call getarg(1, tmp)
        if (trim(tmp) /= '-p') then
           write(0,*) 'ERROR: first argument must be -p'
           call print_usage()
           call exit(2)
        end if
        call getarg(2, process_name)
        call getarg(3, nc_name)
        action = 'process'
     elseif (iargc() > 1) then
        write(0,*) 'ERROR: invalid number of arguments'
        call print_usage()
        action = 'die'
     else
        call print_usage()
        action = 'die'
     end if
  end if
  
  call pmc_mpi_bcast_string(action)
  
  if (trim(action) == 'run') then
     call pmc_mpi_bcast_string(spec_name)
     call partmc_run(spec_name)
  elseif (trim(action) == 'process') then
     call pmc_mpi_bcast_string(process_name)
     call pmc_mpi_bcast_string(nc_name)
     call partmc_process(process_name, nc_name)
  elseif (trim(action) == 'die') then
     call pmc_mpi_abort(2)
  else
     if (pmc_mpi_rank() == 0) then
        write(0,*) 'ERROR: unknown action: ', trim(action)
     end if
     call pmc_mpi_abort(1)
  end if
     
  call pmc_mpi_finalize()

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine print_usage()

    write(0,*) 'Usage: partmc <spec-file>'
    write(0,*) 'or:    partmc -p <process.dat> <nc-file> <state-files...>'

  end subroutine print_usage

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine partmc_process(process_name, nc_name)

    character(len=*), intent(in) :: process_name ! process.dat filename
    character(len=*), intent(in) :: nc_name ! output NetCDF filename

    type(inout_file_t) :: file
    type(process_spec_t), pointer :: process_spec_list(:)
    integer :: i
    character(len=300) :: state_name
    integer :: ncid

    if (pmc_mpi_rank() /= 0) then
       return
    end if

    call inout_open_read(process_name, file)
    call inout_read_process_spec_list(file, process_spec_list)
    call inout_close(file)
    
    call output_processed_open(nc_name, -1, ncid)
    do i = 4,iargc()
       call getarg(i, state_name)
       call partmc_process_state_file(ncid, state_name, process_spec_list)
    end do
    call output_processed_close(ncid)
    
    call process_spec_list_free(process_spec_list)

  end subroutine partmc_process

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine partmc_process_state_file(ncid, state_name, process_spec_list)

    integer, intent(in) :: ncid         ! NetCDF file ID, must be open
    character(len=*), intent(in) :: state_name ! state filename
    type(process_spec_t), pointer :: process_spec_list(:) ! process spec

    type(bin_grid_t) :: bin_grid
    type(aero_data_t) :: aero_data
    type(aero_state_t) :: aero_state
    type(gas_data_t) :: gas_data
    type(gas_state_t) :: gas_state
    type(env_state_t) :: env_state
    real*8 :: time, del_t
    integer :: index, i_loop
    
    call inout_read_state(state_name, bin_grid, aero_data, aero_state, &
         gas_data, gas_state, env_state, time, index, del_t, i_loop)
    
    call output_processed(ncid, process_spec_list, &
         bin_grid, aero_data, aero_state, gas_data, gas_state, &
         env_state, index, time, del_t, i_loop)

    call bin_grid_free(bin_grid)
    call aero_data_free(aero_data)
    call aero_state_free(aero_state)
    call gas_data_free(gas_data)
    call gas_state_free(gas_state)
    call env_state_free(env_state)
    
  end subroutine partmc_process_state_file

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine partmc_run(spec_name)
    
    character(len=*), intent(in) :: spec_name ! spec filename

    type(inout_file_t) :: file
    character(len=100) :: run_type
    integer :: i

    ! check filename (must be "filename.spec")
    i = len_trim(spec_name)
    if (spec_name((i-4):i) /= '.spec') then
       write(6,*) 'ERROR: input filename must end in .spec'
       call pmc_mpi_abort(3)
    end if
    
    if (pmc_mpi_rank() == 0) then
       ! only the root process does I/O
       call inout_open_read(spec_name, file)
       call inout_read_string(file, 'run_type', run_type)
    end if
    
    call pmc_mpi_bcast_string(run_type)
    if (trim(run_type) == 'mc') then
       call partmc_mc(file)
    elseif (trim(run_type) == 'exact') then
       call partmc_exact(file)
    elseif (trim(run_type) == 'sect') then
       call partmc_sect(file)
    else
       if (pmc_mpi_rank() == 0) then
          write(0,*) 'ERROR: unknown run_type: ', trim(run_type)
       end if
       call pmc_mpi_abort(1)
    end if

  end subroutine partmc_run

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine partmc_mc(file)

    type(inout_file_t), intent(out) :: file ! spec file

    character(len=300) :: summary_name  ! name of output files
    character(len=100) :: kernel_name   ! coagulation kernel name
    type(gas_data_t) :: gas_data        ! gas data
    type(gas_state_t) :: gas_init       ! gas initial condition
    type(gas_state_t) :: gas_state      ! gas current state for run
    type(aero_data_t) :: aero_data      ! aero_data data
    type(aero_dist_t) :: aero_dist_init ! aerosol initial distribution
    type(aero_state_t) :: aero_state    ! aerosol current state for run
    type(env_data_t) :: env_data        ! environment data
    type(env_state_t) :: env_state                  ! environment state
    type(bin_grid_t) :: bin_grid        ! bin grid
    type(aero_binned_t) :: aero_binned  ! binned distributions
    type(run_mc_opt_t) :: mc_opt        ! Monte Carlo options
    type(process_spec_t), pointer :: process_spec_list(:) ! process specs
    integer :: i_loop                   ! current loop number
    integer :: rand_init                ! random number generator init
    character, allocatable :: buffer(:) ! buffer for MPI
    integer :: buffer_size              ! length of buffer
    integer :: position                 ! current position in buffer
    
    if (pmc_mpi_rank() == 0) then
       ! only the root process does I/O

       call inout_read_string(file, 'output_file', summary_name)
       call inout_read_string(file, 'state_prefix', mc_opt%state_prefix)
       call spec_read_process_spec_list_filename(file, 'process_spec', &
            process_spec_list)
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
       
       call spec_read_env_data(file, bin_grid, gas_data, aero_data, env_data)
       call spec_read_env_state(file, env_state)
       
       call inout_read_integer(file, 'rand_init', rand_init)
       call inout_read_real(file, 'mix_rate', mc_opt%mix_rate)
       call inout_read_logical(file, 'do_coagulation', mc_opt%do_coagulation)
       call inout_read_logical(file, 'allow_double', mc_opt%allow_double)
       call inout_read_logical(file, 'do_condensation', mc_opt%do_condensation)
       call inout_read_logical(file, 'do_mosaic', mc_opt%do_mosaic)
       if (mc_opt%do_mosaic .and. (.not. mosaic_support())) then
          write(0,'(a,i3,a,a,a)') 'ERROR: line ', file%line_num, &
            ' of input file ', trim(file%name), &
            ': cannot use MOSAIC, support is not compiled in'
          call exit(1)
       end if

       call inout_read_logical(file, 'do_restart', mc_opt%do_restart)
       call inout_read_string(file, 'restart_name', mc_opt%restart_name)
       
       call inout_close(file)
    end if

    ! finished reading .spec data, now broadcast data

#ifdef PMC_USE_MPI
    if (pmc_mpi_rank() == 0) then
       ! root process determines size
       buffer_size = 0
       buffer_size = buffer_size + pmc_mpi_pack_size_mc_opt(mc_opt)
       buffer_size = buffer_size + pmc_mpi_pack_size_string(kernel_name)
       buffer_size = buffer_size + pmc_mpi_pack_size_bin_grid(bin_grid)
       buffer_size = buffer_size + pmc_mpi_pack_size_gas_data(gas_data)
       buffer_size = buffer_size + pmc_mpi_pack_size_gas_state(gas_init)
       buffer_size = buffer_size + pmc_mpi_pack_size_aero_data(aero_data)
       buffer_size = buffer_size &
            + pmc_mpi_pack_size_aero_dist(aero_dist_init)
       buffer_size = buffer_size + pmc_mpi_pack_size_env_data(env_data)
       buffer_size = buffer_size + pmc_mpi_pack_size_env_state(env_state)
       buffer_size = buffer_size + pmc_mpi_pack_size_integer(rand_init)
       buffer_size = buffer_size &
            + pmc_mpi_pack_size_process_spec_list(process_spec_list)
    end if

    ! tell everyone the size and allocate buffer space
    call pmc_mpi_bcast_integer(buffer_size)
    allocate(buffer(buffer_size))

    if (pmc_mpi_rank() == 0) then
       ! root process packs data
       position = 0
       call pmc_mpi_pack_mc_opt(buffer, position, mc_opt)
       call pmc_mpi_pack_string(buffer, position, kernel_name)
       call pmc_mpi_pack_bin_grid(buffer, position, bin_grid)
       call pmc_mpi_pack_gas_data(buffer, position, gas_data)
       call pmc_mpi_pack_gas_state(buffer, position, gas_init)
       call pmc_mpi_pack_aero_data(buffer, position, aero_data)
       call pmc_mpi_pack_aero_dist(buffer, position, aero_dist_init)
       call pmc_mpi_pack_env_data(buffer, position, env_data)
       call pmc_mpi_pack_env_state(buffer, position, env_state)
       call pmc_mpi_pack_integer(buffer, position, rand_init)
       call pmc_mpi_pack_process_spec_list(buffer, position, process_spec_list)
    end if

    ! broadcast data to everyone
    call pmc_mpi_bcast_packed(buffer)

    if (pmc_mpi_rank() /= 0) then
       ! non-root processes unpack data
       position = 0
       call pmc_mpi_unpack_mc_opt(buffer, position, mc_opt)
       call pmc_mpi_unpack_string(buffer, position, kernel_name)
       call pmc_mpi_unpack_bin_grid(buffer, position, bin_grid)
       call pmc_mpi_unpack_gas_data(buffer, position, gas_data)
       call pmc_mpi_unpack_gas_state(buffer, position, gas_init)
       call pmc_mpi_unpack_aero_data(buffer, position, aero_data)
       call pmc_mpi_unpack_aero_dist(buffer, position, aero_dist_init)
       call pmc_mpi_unpack_env_data(buffer, position, env_data)
       call pmc_mpi_unpack_env_state(buffer, position, env_state)
       call pmc_mpi_unpack_integer(buffer, position, rand_init)
       call pmc_mpi_unpack_process_spec_list(buffer, position, &
            process_spec_list)
    end if

    ! free the buffer
    deallocate(buffer)
#endif

    call pmc_srand(rand_init + pmc_mpi_rank())

    call aero_binned_alloc(aero_binned, bin_grid%n_bin, aero_data%n_spec)
    call gas_state_alloc(gas_state, gas_data%n_spec)
    call cpu_time(mc_opt%t_wall_start)

    call aero_state_alloc(0, 0, aero_state)
    do i_loop = 1,mc_opt%n_loop
       mc_opt%i_loop = i_loop
       
       call gas_state_copy(gas_init, gas_state)
       call aero_state_free(aero_state)
       call aero_dist_to_state(bin_grid, aero_data, aero_dist_init, &
            mc_opt%n_part_max, aero_state)
       call env_data_init_state(env_data, env_state, 0d0)

       if (mc_opt%do_condensation) then
          call aero_state_equilibriate(bin_grid, env_state, aero_data, aero_state)
       end if
       
       if (trim(kernel_name) == 'sedi') then
          call run_mc(kernel_sedi, bin_grid, aero_binned, env_data, &
               env_state, aero_data, aero_state, gas_data, gas_state, &
               mc_opt, process_spec_list)
       elseif (trim(kernel_name) == 'golovin') then
          call run_mc(kernel_golovin, bin_grid, aero_binned, env_data, &
               env_state, aero_data, aero_state, gas_data, gas_state, &
               mc_opt, process_spec_list)
       elseif (trim(kernel_name) == 'constant') then
          call run_mc(kernel_constant, bin_grid, aero_binned, env_data, &
               env_state, aero_data, aero_state, gas_data, gas_state, &
               mc_opt, process_spec_list)
       elseif (trim(kernel_name) == 'brown') then
          call run_mc(kernel_brown, bin_grid, aero_binned, env_data, &
               env_state, aero_data, aero_state, gas_data, gas_state, &
               mc_opt, process_spec_list)
       elseif (trim(kernel_name) == 'zero') then
          call run_mc(kernel_zero, bin_grid, aero_binned, env_data, &
               env_state, aero_data, aero_state, gas_data, gas_state, &
               mc_opt, process_spec_list)
       else
          if (pmc_mpi_rank() == 0) then
             write(0,*) 'ERROR: Unknown kernel type; ', trim(kernel_name)
          end if
          call pmc_mpi_abort(1)
       end if

    end do

    call gas_data_free(gas_data)
    call gas_state_free(gas_init)
    call gas_state_free(gas_state)
    call aero_data_free(aero_data)
    call aero_dist_free(aero_dist_init)
    call aero_state_free(aero_state)
    call env_data_free(env_data)
    call env_state_free(env_state)
    call bin_grid_free(bin_grid)
    call aero_binned_free(aero_binned)
    call process_spec_list_free(process_spec_list)

  end subroutine partmc_mc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine partmc_exact(file)

    type(inout_file_t), intent(out) :: file ! spec file

    character(len=300) :: summary_name  ! name of output files
    character(len=100) :: soln_name     ! exact solution name
    type(aero_data_t) :: aero_data      ! aero_data data
    type(env_data_t) :: env_data        ! environment data
    type(env_state_t) :: env_state                  ! environment state
    type(run_exact_opt_t) :: exact_opt  ! exact solution options
    type(bin_grid_t) :: bin_grid        ! bin grid
    type(gas_data_t) :: gas_data        ! dummy gas data
    type(process_spec_t), pointer :: process_spec_list(:) ! process specs

    ! only serial code here
    if (pmc_mpi_rank() /= 0) then
       return
    end if
    
    call inout_read_string(file, 'output_file', summary_name)
    call inout_read_string(file, 'output_prefix', exact_opt%prefix)
    call spec_read_process_spec_list_filename(file, 'process_spec', &
         process_spec_list)
    call inout_read_real(file, 'num_den', exact_opt%num_den)

    call inout_read_real(file, 't_max', exact_opt%t_max)
    call inout_read_real(file, 't_output', exact_opt%t_output)

    call spec_read_bin_grid(file, bin_grid)
    call spec_read_gas_data(file, gas_data)
    call spec_read_aero_data_filename(file, aero_data)
    call spec_read_env_data(file, bin_grid, gas_data, aero_data, env_data)
    call spec_read_env_state(file, env_state)

    call inout_read_string(file, 'soln', soln_name)

    call aero_dist_alloc(exact_opt%aero_dist_init, 0, 0, 0)
    if (trim(soln_name) == 'golovin_exp') then
       call inout_read_real(file, 'mean_radius', exact_opt%mean_radius)
    elseif (trim(soln_name) == 'constant_exp_cond') then
       call inout_read_real(file, 'mean_radius', exact_opt%mean_radius)
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

    if (trim(soln_name) == 'golovin_exp') then
       call run_exact(bin_grid, env_data, env_state, aero_data, exact_opt, &
            soln_golovin_exp, process_spec_list)
    elseif (trim(soln_name) == 'golovin_exp') then
       call run_exact(bin_grid, env_data, env_state, aero_data, exact_opt, &
            soln_constant_exp_cond, process_spec_list)
    elseif (trim(soln_name) == 'zero') then
       call run_exact(bin_grid, env_data, env_state, aero_data, exact_opt, &
            soln_zero, process_spec_list)
    else
       write(0,*) 'ERROR: unknown solution type: ', trim(soln_name)
       call exit(1)
    end if

    call aero_data_free(aero_data)
    call env_data_free(env_data)
    call env_state_free(env_state)
    call bin_grid_free(bin_grid)
    call gas_data_free(gas_data)
    call aero_dist_free(exact_opt%aero_dist_init)
    call process_spec_list_free(process_spec_list)
    
  end subroutine partmc_exact

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine partmc_sect(file)

    type(inout_file_t), intent(out) :: file ! spec file

    character(len=300) :: summary_name  ! name of output files
    character(len=100) :: kernel_name   ! coagulation kernel name
    type(run_sect_opt_t) :: sect_opt    ! sectional code options
    type(aero_data_t) :: aero_data      ! aero_data data
    type(aero_dist_t) :: aero_dist_init ! aerosol initial distribution
    type(aero_state_t) :: aero_init     ! aerosol initial condition
    type(env_data_t) :: env_data        ! environment data
    type(env_state_t) :: env_state                  ! environment state
    type(bin_grid_t) :: bin_grid        ! bin grid
    type(gas_data_t) :: gas_data        ! dummy gas data
    type(process_spec_t), pointer :: process_spec_list(:) ! process specs

    ! only serial code here
    if (pmc_mpi_rank() /= 0) then
       return
    end if
    
    call inout_read_string(file, 'output_file', summary_name)
    call inout_read_string(file, 'output_prefix', sect_opt%prefix)
    call spec_read_process_spec_list_filename(file, 'process_spec', &
         process_spec_list)
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

    call spec_read_env_data(file, bin_grid, gas_data, aero_data, env_data)
    call spec_read_env_state(file, env_state)

    call inout_read_logical(file, 'do_coagulation', sect_opt%do_coagulation)
    
    call inout_close(file)

    ! finished reading .spec data, now do the run

    if (trim(kernel_name) == 'sedi') then
       call run_sect(bin_grid, gas_data, aero_data, aero_dist_init, &
            env_data, env_state, kernel_sedi, sect_opt, process_spec_list)
    elseif (trim(kernel_name) == 'golovin') then
       call run_sect(bin_grid, gas_data, aero_data, aero_dist_init, &
            env_data, env_state, kernel_golovin, sect_opt, process_spec_list)
    elseif (trim(kernel_name) == 'constant') then
       call run_sect(bin_grid, gas_data, aero_data, aero_dist_init, &
            env_data, env_state, kernel_constant, sect_opt, process_spec_list)
    elseif (trim(kernel_name) == 'brown') then
       call run_sect(bin_grid, gas_data, aero_data, aero_dist_init, &
            env_data, env_state, kernel_brown, sect_opt, process_spec_list)
    elseif (trim(kernel_name) == 'zero') then
       call run_sect(bin_grid, gas_data, aero_data, aero_dist_init, &
            env_data, env_state, kernel_zero, sect_opt, process_spec_list)
    else
       write(0,*) 'ERROR: Unknown kernel type; ', trim(kernel_name)
       call exit(1)
    end if

    call aero_data_free(aero_data)
    call aero_dist_free(aero_dist_init)
    call env_state_free(env_state)
    call env_data_free(env_data)
    call bin_grid_free(bin_grid)
    call gas_data_free(gas_data)
    call process_spec_list_free(process_spec_list)
    
  end subroutine partmc_sect

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program partmc
