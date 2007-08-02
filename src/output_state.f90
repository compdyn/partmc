! Copyright (C) 2005-2007 Nicole Riemer and Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.
!
! Save and restore the state of particle-resolved Monte Carlo
! runs. The state file should contain enough data to restart the
! simulation at the point it was written.
!
! Because it contains the full state of every particle, this is also
! the best way to gain complete access to all statistics of the
! simulation.

module pmc_output_state
contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine inout_write_state(state_prefix, bin_grid, aero_data, &
       aero_state, gas_data, gas_state, env, index, time, i_loop)

    ! Write the current state.

    use pmc_bin_grid
    use pmc_aero_data
    use pmc_aero_state
    use pmc_env
    use pmc_util
    use pmc_inout
    use pmc_gas_data
    use pmc_mpi
    
    character(len=*), intent(in) :: state_prefix ! prefix of state file
    type(bin_grid_t), intent(in) :: bin_grid ! bin grid
    type(aero_data_t), intent(in) :: aero_data ! aerosol data
    type(aero_state_t), intent(in) :: aero_state ! aerosol state
    type(gas_data_t), intent(in) :: gas_data ! gas data
    type(gas_state_t), intent(in) :: gas_state ! gas state
    type(env_t), intent(in) :: env      ! environment state
    integer, intent(in) :: index        ! filename index
    real*8, intent(in) :: time          ! current time (s)
    integer, intent(in) :: i_loop       ! current loop number
    
    character*300 :: filename
    type(inout_file_t) :: file

    if (pmc_mpi_rank() == 0) then
       ! FIXME: for now only the root process writes state

       write(filename, '(a,a,i4.4,a,i8.8,a)') trim(state_prefix), &
            '_', i_loop, '_', index, '.d'
       call inout_open_write(filename, file)
       
       call inout_write_real(file, 'time(s)', time)
       call inout_write_integer(file, 'loop', i_loop)
       call inout_write_integer(file, 'index', index)
       
       call inout_write_env(file, env)
       call inout_write_bin_grid(file, bin_grid)
       call inout_write_gas_data(file, gas_data)
       call inout_write_gas_state(file, gas_state)
       call inout_write_aero_data(file, aero_data)
       call inout_write_aero_state(file, aero_state)
       
       call inout_close(file)
    end if
    
  end subroutine inout_write_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine inout_read_state(state_name, bin_grid, aero_data, &
       aero_state, gas_data, gas_state, env, time)

    ! Read the current state.

    use pmc_bin_grid
    use pmc_aero_data
    use pmc_aero_state
    use pmc_env
    use pmc_util
    use pmc_inout
    use pmc_gas_data
    use pmc_mpi
    
    character(len=*), intent(in) :: state_name ! name of state file
    type(bin_grid_t), intent(out) :: bin_grid ! bin grid
    type(aero_data_t), intent(out) :: aero_data ! aerosol data
    type(aero_state_t), intent(out) :: aero_state ! aerosol state
    type(gas_data_t), intent(out) :: gas_data ! gas data
    type(gas_state_t), intent(out) :: gas_state ! gas state
    type(env_t), intent(out) :: env     ! environment state
    real*8, intent(out) :: time         ! current time (s)
    
    type(inout_file_t) :: file
    integer :: dummy_integer

    if (pmc_mpi_rank() /= 0) then
       call pmc_mpi_abort(52115)
    end if

    call inout_open_read(state_name, file)
    
    call inout_read_real(file, 'time(s)', time)
    call inout_read_integer(file, 'loop', dummy_integer)
    call inout_read_integer(file, 'index', dummy_integer)

    call inout_read_env(file, env)
    call inout_read_bin_grid(file, bin_grid)
    call inout_read_gas_data(file, gas_data)
    call inout_read_gas_state(file, gas_state)
    call inout_read_aero_data(file, aero_data)
    call inout_read_aero_state(file, aero_state)

    call inout_close(file)
    
  end subroutine inout_read_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_output_state
