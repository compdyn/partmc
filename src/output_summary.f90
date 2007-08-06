! Copyright (C) 2005-2007 Nicole Riemer and Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

module pmc_output_summary
contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine output_summary_header(file, bin_grid, gas_data, &
       aero_data, n_loop, n_time)

    ! Print summary header.

    use pmc_bin_grid
    use pmc_inout
    use pmc_aero_data
    use pmc_gas_data
    use pmc_mpi

    type(inout_file_t), intent(inout) :: file ! file to output to
    type(bin_grid_t), intent(in) :: bin_grid ! bin grid
    type(gas_data_t), intent(in) :: gas_data ! gas data
    type(aero_data_t), intent(in) :: aero_data ! aerosol data
    integer, intent(in) :: n_loop       ! number of loops
    integer, intent(in) :: n_time       ! number of times

    if (pmc_mpi_rank() == 0) then
       ! only the root process does I/O
       call inout_write_integer(file, 'n_loop', n_loop)
       call inout_write_integer(file, 'n_time', n_time)
       call inout_write_bin_grid(file, bin_grid)
       call inout_write_gas_data(file, gas_data)
       call inout_write_aero_data(file, aero_data)
    end if

  end subroutine output_summary_header
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine output_summary(file, time, bin_grid, aero_data, &
       aero_binned, gas_data, gas_state, env, i_loop)

    ! Write the current binned data to the output file. This version
    ! of the function takes absolute number and absolute volume
    ! per-bin (as produced by a particle-resolved code, for example).
    
    use pmc_bin_grid
    use pmc_aero_data
    use pmc_aero_binned
    use pmc_env
    use pmc_inout
    use pmc_gas_data
    use pmc_gas_state
    use pmc_mpi

    type(inout_file_t), intent(inout) :: file ! file to output to
    real*8, intent(in) :: time          ! simulation time
    type(bin_grid_t), intent(in) :: bin_grid ! bin grid
    type(aero_data_t), intent(in) :: aero_data ! aerosol data
    type(aero_binned_t), intent(in) :: aero_binned ! binned aerosol data
    type(gas_data_t), intent(in) :: gas_data ! gas data
    type(gas_state_t), intent(in) :: gas_state ! gas state
    type(env_t), intent(in) :: env      ! environment state
    integer, intent(in) :: i_loop       ! current loop number

#ifdef PMC_USE_MPI
    type(aero_binned_t) :: aero_binned_avg
    type(gas_state_t) :: gas_state_avg
    type(env_t) :: env_avg

    call aero_binned_alloc(aero_binned_avg, bin_grid%n_bin, aero_data%n_spec)
    call gas_state_alloc(gas_state_avg, gas_data%n_spec)
    call env_alloc(env_avg)

    call pmc_mpi_reduce_avg_aero_binned(aero_binned, aero_binned_avg)
    call pmc_mpi_reduce_avg_gas_state(gas_state, gas_state_avg)
    call pmc_mpi_reduce_avg_env(env, env_avg)
#endif

    if (pmc_mpi_rank() == 0) then
       ! only the root process does I/O
       call inout_write_integer(file, 'loop_num', i_loop)
       call inout_write_real(file, 'time(s)', time)
#ifdef PMC_USE_MPI
       call inout_write_env(file, env_avg)
       call inout_write_aero_binned(file, aero_binned_avg)
       call inout_write_gas_state(file, gas_state_avg)
#else
       call inout_write_env(file, env)
       call inout_write_aero_binned(file, aero_binned)
       call inout_write_gas_state(file, gas_state)
#endif
    end if

  end subroutine output_summary
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_output_summary
