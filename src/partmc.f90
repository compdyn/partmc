! Copyright (C) 2007-2009 Nicole Riemer and Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The partmc program.

!> \mainpage
!>
!> \dotfile partmc_modules.gv

!> Top level driver.
program partmc

  use pmc_mpi
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
  use pmc_spec_file
  use pmc_gas_data
  use pmc_gas_state
  use pmc_util

  character(len=300) :: spec_name
  
  call pmc_mpi_init()

  if (pmc_mpi_rank() == 0) then
     ! only the root process accesses the commandline

     if (iargc() /= 1) then
        call print_usage()
        call die_msg(739173192, "invalid commandline arguments")
     end if

     call getarg(1, spec_name)
  end if

  call pmc_mpi_bcast_string(spec_name)
  call partmc_run(spec_name)

  call pmc_mpi_finalize()

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Print the usage text to stderr.
  subroutine print_usage()

    write(*,*) 'Usage: partmc <spec-file>'

  end subroutine print_usage

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Do a PartMC run.
  subroutine partmc_run(spec_name)
    
    !> Spec filename.
    character(len=*), intent(in) :: spec_name

    type(spec_file_t) :: file
    character(len=100) :: run_type
    integer :: i

    ! check filename (must be "filename.spec")
    i = len_trim(spec_name)
    if (spec_name((i-4):i) /= '.spec') then
       call die_msg(710381938, "input filename must end in .spec")
    end if
    
    if (pmc_mpi_rank() == 0) then
       ! only the root process does I/O
       call spec_file_open(spec_name, file)
       call spec_file_read_string(file, 'run_type', run_type)
    end if
    
    call pmc_mpi_bcast_string(run_type)
    if (trim(run_type) == 'mc') then
       call partmc_mc(file)
    elseif (trim(run_type) == 'exact') then
       call partmc_exact(file)
    elseif (trim(run_type) == 'sect') then
       call partmc_sect(file)
    else
       call die_msg(719261940, "unknown run_type: " // trim(run_type))
    end if

  end subroutine partmc_run

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Run a Monte Carlo simulation.
  subroutine partmc_mc(file)

    !> Spec file.
    type(spec_file_t), intent(out) :: file

    character(len=100) :: kernel_name
    type(gas_data_t) :: gas_data
    type(gas_state_t) :: gas_init
    type(gas_state_t) :: gas_state
    type(aero_data_t) :: aero_data
    type(aero_dist_t) :: aero_dist_init
    type(aero_state_t) :: aero_state
    type(env_data_t) :: env_data
    type(env_state_t) :: env_state
    type(bin_grid_t) :: bin_grid
    type(run_mc_opt_t) :: mc_opt
    integer :: i_loop
    integer :: rand_init
    character, allocatable :: buffer(:)
    integer :: buffer_size
    integer :: position
    
    if (pmc_mpi_rank() == 0) then
       ! only the root process does I/O

       call spec_file_read_string(file, 'output_prefix', mc_opt%output_prefix)
       call spec_file_read_integer(file, 'n_loop', mc_opt%n_loop)
       call spec_file_read_integer(file, 'n_part', mc_opt%n_part_max)
       call spec_file_read_string(file, 'kernel', kernel_name)
       
       call spec_file_read_real(file, 't_max', mc_opt%t_max)
       call spec_file_read_real(file, 'del_t', mc_opt%del_t)
       call spec_file_read_real(file, 't_output', mc_opt%t_output)
       call spec_file_read_real(file, 't_progress', mc_opt%t_progress)

       call bin_grid_allocate(bin_grid)
       call spec_file_read_bin_grid(file, bin_grid)

       call gas_data_allocate(gas_data)
       call spec_file_read_gas_data(file, gas_data)
       call gas_state_allocate(gas_init)
       call spec_file_read_gas_state(file, gas_data, 'gas_init', gas_init)

       call aero_data_allocate(aero_data)
       call spec_file_read_aero_data_filename(file, aero_data)
       call aero_dist_allocate(aero_dist_init)
       call spec_file_read_aero_dist_filename(file, aero_data, bin_grid, &
            'aerosol_init', aero_dist_init)

       call env_data_allocate(env_data)
       call spec_file_read_env_data(file, bin_grid, gas_data, aero_data, env_data)
       call env_state_allocate(env_state)
       call spec_file_read_env_state(file, env_state)
       
       call spec_file_read_integer(file, 'rand_init', rand_init)
       call spec_file_read_real(file, 'mix_rate', mc_opt%mix_rate)
       call spec_file_read_logical(file, 'do_coagulation', mc_opt%do_coagulation)
       call spec_file_read_logical(file, 'allow_doubling', mc_opt%allow_doubling)
       call spec_file_read_logical(file, 'allow_halving', mc_opt%allow_halving)
       call spec_file_read_logical(file, 'do_condensation', mc_opt%do_condensation)
       call spec_file_read_logical(file, 'do_mosaic', mc_opt%do_mosaic)
       if (mc_opt%do_mosaic .and. (.not. mosaic_support())) then
          call spec_file_die_msg(230495365, file, &
               'cannot use MOSAIC, support is not compiled in')
       end if

       call spec_file_read_logical(file, 'record_removals', mc_opt%record_removals)
       
       call spec_file_close(file)
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
    end if

    ! free the buffer
    deallocate(buffer)
#endif

    call pmc_srand(rand_init + pmc_mpi_rank())

    call gas_state_allocate_size(gas_state, gas_data%n_spec)
    call cpu_time(mc_opt%t_wall_start)

    call aero_state_allocate_size(aero_state, bin_grid%n_bin, aero_data%n_spec)
    do i_loop = 1,mc_opt%n_loop
       mc_opt%i_loop = i_loop
       
       call gas_state_copy(gas_init, gas_state)
       call aero_state_deallocate(aero_state)
       call aero_state_allocate_size(aero_state, bin_grid%n_bin, aero_data%n_spec)
       aero_state%comp_vol = dble(mc_opt%n_part_max) / &
            aero_dist_total_num_conc(aero_dist_init)
       call aero_state_add_aero_dist_sample(aero_state, bin_grid, &
            aero_data, aero_dist_init, 1d0, 0d0)
       call env_data_init_state(env_data, env_state, 0d0)

       if (mc_opt%do_condensation) then
          call aero_state_equilibriate(bin_grid, env_state, aero_data, &
               aero_state)
       end if
       
       if (trim(kernel_name) == 'sedi') then
          call run_mc(kernel_sedi, kernel_sedi_max, bin_grid, &
               env_data, env_state, aero_data, &
               aero_state, gas_data, gas_state, mc_opt)
       elseif (trim(kernel_name) == 'golovin') then
          call run_mc(kernel_golovin, kernel_golovin_max, bin_grid, &
               env_data, env_state, aero_data, &
               aero_state, gas_data, gas_state, mc_opt)
       elseif (trim(kernel_name) == 'constant') then
          call run_mc(kernel_constant, kernel_constant_max, bin_grid, &
               env_data, env_state, aero_data, &
               aero_state, gas_data, gas_state, mc_opt)
       elseif (trim(kernel_name) == 'brown') then
          call run_mc(kernel_brown, kernel_brown_max, bin_grid, &
               env_data, env_state, aero_data, &
               aero_state, gas_data, gas_state, mc_opt)
       elseif (trim(kernel_name) == 'zero') then
          call run_mc(kernel_zero, kernel_zero_max, bin_grid, &
               env_data, env_state, aero_data, &
               aero_state, gas_data, gas_state, mc_opt)
       else
          if (pmc_mpi_rank() == 0) then
             call die_msg(727498351, 'unknown kernel type: ' &
                  // trim(kernel_name))
          end if
          call pmc_mpi_abort(1)
       end if

    end do

    call gas_data_deallocate(gas_data)
    call gas_state_deallocate(gas_init)
    call gas_state_deallocate(gas_state)
    call aero_data_deallocate(aero_data)
    call aero_dist_deallocate(aero_dist_init)
    call aero_state_deallocate(aero_state)
    call env_data_deallocate(env_data)
    call env_state_deallocate(env_state)
    call bin_grid_deallocate(bin_grid)

  end subroutine partmc_mc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Run an exact solution simulation.
  subroutine partmc_exact(file)

    !> Spec file.
    type(spec_file_t), intent(out) :: file

    character(len=100) :: soln_name
    type(aero_data_t) :: aero_data
    type(env_data_t) :: env_data
    type(env_state_t) :: env_state
    type(run_exact_opt_t) :: exact_opt
    type(bin_grid_t) :: bin_grid
    type(gas_data_t) :: gas_data

    ! only serial code here
    if (pmc_mpi_rank() /= 0) then
       return
    end if
    
    call bin_grid_allocate(bin_grid)
    call gas_data_allocate(gas_data)
    call aero_data_allocate(aero_data)
    call env_data_allocate(env_data)
    call env_state_allocate(env_state)
    call aero_dist_allocate(exact_opt%aero_dist_init)

    call spec_file_read_string(file, 'output_prefix', exact_opt%prefix)
    call spec_file_read_real(file, 'num_conc', exact_opt%num_conc)

    call spec_file_read_real(file, 't_max', exact_opt%t_max)
    call spec_file_read_real(file, 't_output', exact_opt%t_output)

    call spec_file_read_bin_grid(file, bin_grid)
    call spec_file_read_gas_data(file, gas_data)
    call spec_file_read_aero_data_filename(file, aero_data)
    call spec_file_read_env_data(file, bin_grid, gas_data, aero_data, env_data)
    call spec_file_read_env_state(file, env_state)

    call spec_file_read_string(file, 'soln', soln_name)

    if (trim(soln_name) == 'golovin_exp') then
       call spec_file_read_real(file, 'mean_radius', exact_opt%mean_radius)
    elseif (trim(soln_name) == 'constant_exp_cond') then
       call spec_file_read_real(file, 'mean_radius', exact_opt%mean_radius)
    elseif (trim(soln_name) == 'zero') then
       call aero_dist_deallocate(exact_opt%aero_dist_init)
       call spec_file_read_aero_dist_filename(file, aero_data, bin_grid, &
            'aerosol_init', exact_opt%aero_dist_init)
    else
       call die_msg(955390033, 'unknown solution type: ' &
            // trim(soln_name))
    end if
    
    call spec_file_close(file)

    ! finished reading .spec data, now do the run

    call env_data_init_state(env_data, env_state, 0d0)

    if (trim(soln_name) == 'golovin_exp') then
       call run_exact(bin_grid, env_data, env_state, aero_data, exact_opt, &
            soln_golovin_exp)
    elseif (trim(soln_name) == 'golovin_exp') then
       call run_exact(bin_grid, env_data, env_state, aero_data, exact_opt, &
            soln_constant_exp_cond)
    elseif (trim(soln_name) == 'zero') then
       call run_exact(bin_grid, env_data, env_state, aero_data, exact_opt, &
            soln_zero)
    else
       call die_msg(859292825, 'unknown solution type: ' &
            // trim(soln_name))
    end if

    call aero_data_deallocate(aero_data)
    call env_data_deallocate(env_data)
    call env_state_deallocate(env_state)
    call bin_grid_deallocate(bin_grid)
    call gas_data_deallocate(gas_data)
    call aero_dist_deallocate(exact_opt%aero_dist_init)
    
  end subroutine partmc_exact

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Run a sectional code simulation.
  subroutine partmc_sect(file)

    !> Spec file.
    type(spec_file_t), intent(out) :: file

    character(len=100) :: kernel_name
    type(run_sect_opt_t) :: sect_opt
    type(aero_data_t) :: aero_data
    type(aero_dist_t) :: aero_dist_init
    type(aero_state_t) :: aero_init
    type(env_data_t) :: env_data
    type(env_state_t) :: env_state
    type(bin_grid_t) :: bin_grid
    type(gas_data_t) :: gas_data

    ! only serial code here
    if (pmc_mpi_rank() /= 0) then
       return
    end if
    
    call aero_data_allocate(aero_data)
    call aero_dist_allocate(aero_dist_init)
    call env_state_allocate(env_state)
    call env_data_allocate(env_data)
    call bin_grid_allocate(bin_grid)
    call gas_data_allocate(gas_data)

    call spec_file_read_string(file, 'output_prefix', sect_opt%prefix)
    call spec_file_read_string(file, 'kernel', kernel_name)

    call spec_file_read_real(file, 't_max', sect_opt%t_max)
    call spec_file_read_real(file, 'del_t', sect_opt%del_t)
    call spec_file_read_real(file, 't_output', sect_opt%t_output)
    call spec_file_read_real(file, 't_progress', sect_opt%t_progress)

    call spec_file_read_bin_grid(file, bin_grid)

    call spec_file_read_gas_data(file, gas_data)

    call spec_file_read_aero_data_filename(file, aero_data)
    call spec_file_read_aero_dist_filename(file, aero_data, bin_grid, &
         'aerosol_init', aero_dist_init)

    call spec_file_read_env_data(file, bin_grid, gas_data, aero_data, env_data)
    call spec_file_read_env_state(file, env_state)

    call spec_file_read_logical(file, 'do_coagulation', sect_opt%do_coagulation)
    
    call spec_file_close(file)

    ! finished reading .spec data, now do the run

    call env_data_init_state(env_data, env_state, 0d0)

    if (trim(kernel_name) == 'sedi') then
       call run_sect(bin_grid, gas_data, aero_data, aero_dist_init, &
            env_data, env_state, kernel_sedi, sect_opt)
    elseif (trim(kernel_name) == 'golovin') then
       call run_sect(bin_grid, gas_data, aero_data, aero_dist_init, &
            env_data, env_state, kernel_golovin, sect_opt)
    elseif (trim(kernel_name) == 'constant') then
       call run_sect(bin_grid, gas_data, aero_data, aero_dist_init, &
            env_data, env_state, kernel_constant, sect_opt)
    elseif (trim(kernel_name) == 'brown') then
       call run_sect(bin_grid, gas_data, aero_data, aero_dist_init, &
            env_data, env_state, kernel_brown, sect_opt)
    elseif (trim(kernel_name) == 'zero') then
       call run_sect(bin_grid, gas_data, aero_data, aero_dist_init, &
            env_data, env_state, kernel_zero, sect_opt)
    else
       call die_msg(859292825, 'unknown kernel type: ' &
            // trim(kernel_name))
    end if

    call aero_data_deallocate(aero_data)
    call aero_dist_deallocate(aero_dist_init)
    call env_state_deallocate(env_state)
    call env_data_deallocate(env_data)
    call bin_grid_deallocate(bin_grid)
    call gas_data_deallocate(gas_data)
    
  end subroutine partmc_sect

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program partmc
