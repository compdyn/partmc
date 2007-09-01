! Copyright (C) 2005-2007 Nicole Riemer and Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.
!
! Process the saved state files to obtain summary data.

program process_state

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

  character(len=100) :: filename        ! input filename
  character(len=100) :: basename        ! basename of the input filename
  type(bin_grid_t) :: bin_grid          ! bin_grid structure
  type(aero_data_t) :: aero_data        ! aero_data structure
  type(aero_state_t) :: aero_state      ! aero_state structure
  type(gas_data_t) :: gas_data          ! gas_data structure
  type(gas_state_t) :: gas_state        ! gas_state structure
  type(env_t) :: env                    ! env structure
  real*8 :: time                        ! current time
  character(len=100) :: command         ! process command

  call pmc_mpi_init()

  if (iargc() < 1) then
     call print_usage()
     call exit(1)
  end if
  
  call get_filename(filename, basename)

  call inout_read_state(filename, bin_grid, aero_data, aero_state, &
       gas_data, gas_state, env, time)
  write(*,'(a,e20.10)') 'time (s) = ', time

  if (iargc() == 1) then
     call process_env(env)
     call process_info(bin_grid, aero_data, aero_state)
     call process_moments(basename, bin_grid, aero_data, aero_state, time)
     call process_hist(basename, "_n_orig_part", bin_grid, env, aero_data, &
          aero_state, orig_part_step_comp_grid, orig_part_step_comp, &
          orig_part_particle_func, .false.)
  else
     call getarg(2, command)

     if (command == "comp") then
        call process_hist(basename, "_comp", bin_grid, env, aero_data, &
             aero_state, comp_step_comp_grid, comp_step_comp, &
             comp_particle_func, .false.)

     elseif (command == "kappa") then
        call process_hist(basename, "_kappa", bin_grid, env, aero_data, &
             aero_state, kappa_step_comp_grid, kappa_step_comp, &
             kappa_particle_func, .false.)

     elseif (command == "absorb") then
        call mosaic_init(bin_grid, env, 0d0)
        call mosaic_aero_optical(bin_grid, env, aero_data, &
             aero_state, gas_data, gas_state, time)
        call process_hist(basename, "_absorb", bin_grid, env, aero_data, &
             aero_state, absorb_step_comp_grid, absorb_step_comp, &
             absorb_particle_func, .true.)

     elseif (command == "scatter") then
        call mosaic_init(bin_grid, env, 0d0)
        call mosaic_aero_optical(bin_grid, env, aero_data, &
             aero_state, gas_data, gas_state, time)
        call process_hist(basename, "_scatter", bin_grid, env, aero_data, &
             aero_state, scatter_step_comp_grid, scatter_step_comp, &
             scatter_particle_func, .true.)

     elseif (command == "extinct") then
        call mosaic_init(bin_grid, env, 0d0)
        call mosaic_aero_optical(bin_grid, env, aero_data, &
             aero_state, gas_data, gas_state, time)
        call process_hist(basename, "_extinct", bin_grid, env, aero_data, &
             aero_state, extinct_step_comp_grid, extinct_step_comp, &
             extinct_particle_func, .true.)

     else
        write(0,*) 'ERROR: unknown command'
        call exit(1)
     end if
  end if

  call pmc_mpi_finalize()

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine print_usage()

    write(0,*) 'Usage: process_state <filename.d>'
    write(0,*) '       process_state <filename.d> comp' &
         // ' <n_steps> -a <A species> -b <B species>'
    write(0,*) '       process_state <filename.d> kappa <n_steps>'
    
  end subroutine print_usage

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine get_filename(filename, basename)
    
    character(len=*), intent(out) :: filename ! input filename
    character(len=*), intent(out) :: basename ! basename of the input filename

    integer :: i
    
    ! get and check first commandline argument (must be "filename.d")
    call getarg(1, filename)
    i = len_trim(filename)
    if (i > len(filename)) then
       write(0,*) 'ERROR: filename too long'
       call exit(1)
    end if
    if ((filename(i:i) /= 'd') .or. &
         (filename((i-1):(i-1)) /= '.')) then
       write(0,*) 'ERROR: Filename must end in .d'
       call exit(1)
    end if
    
    ! chop .d off the end of the filename to get the basename
    basename = filename
    basename((i-1):i) = '  '
    
  end subroutine get_filename
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine process_env(env)

    type(env_t), intent(in) :: env      ! environment state

    write(*,'(a,e20.10)') 'temp (K) = ', env%temp
    write(*,'(a,e20.10)') 'rel_humid (1) = ', env%rel_humid
    write(*,'(a,e20.10)') 'pressure (Pa) = ', env%pressure
    write(*,'(a,e20.10)') 'height (m) = ', env%height

  end subroutine process_env

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine process_info(bin_grid, aero_data, aero_state)

    type(bin_grid_t), intent(in) :: bin_grid ! bin_grid structure
    type(aero_data_t), intent(in) :: aero_data ! aero_data structure
    type(aero_state_t), intent(in) :: aero_state ! aero_state structure

    integer :: n_part

    n_part = total_particles(aero_state)
    write(*,'(a,i20)') 'total particles = ', n_part
    write(*,'(a,e20.10)') 'comp_vol (m^3) = ', aero_state%comp_vol
    write(*,'(a,e20.10)') 'num_dens (#/m^3) = ', &
         dble(n_part) / aero_state%comp_vol

  end subroutine process_info

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine process_moments(basename, bin_grid, aero_data, aero_state, &
       time)

    use pmc_util
    use pmc_aero_binned

    character(len=*), intent(in) :: basename ! basename of the input filename
    type(bin_grid_t), intent(in) :: bin_grid ! bin_grid structure
    type(aero_data_t), intent(in) :: aero_data ! aero_data structure
    type(aero_state_t), intent(in) :: aero_state ! aero_state structure
    real*8, intent(in) :: time          ! current time (s)

    type(aero_binned_t) :: aero_binned
    integer :: f_out, i_bin, i_spec

    call aero_binned_alloc(aero_binned, bin_grid%n_bin, aero_data%n_spec)
    call aero_state_to_binned(bin_grid, aero_data, aero_state, aero_binned)
    call open_output(basename, "_aero_binned.d", f_out)
    call aero_binned_write_summary(aero_binned, aero_data, bin_grid, &
         time, 0, f_out)
    close(unit=f_out)
    call aero_binned_free(aero_binned)

  end subroutine process_moments

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
end program process_state
