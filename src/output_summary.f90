! Copyright (C) 2005-2007 Nicole Riemer and Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

module pmc_output_summary

  use pmc_process

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
    type(env_t) :: env_avg
    type(aero_binned_t) :: aero_binned_avg
    type(gas_state_t) :: gas_state_avg

    call aero_binned_alloc(aero_binned_avg, bin_grid%n_bin, aero_data%n_spec)
    call gas_state_alloc(gas_state_avg, gas_data%n_spec)

    call pmc_mpi_reduce_avg_env(env, env_avg)
    call pmc_mpi_reduce_avg_aero_binned(aero_binned, aero_binned_avg)
    call pmc_mpi_reduce_avg_gas_state(gas_state, gas_state_avg)
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

#ifdef PMC_USE_MPI
    call env_free(env_avg)
    call aero_binned_free(aero_binned_avg)
    call gas_state_free(gas_state_avg)
#endif

  end subroutine output_summary
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine output_processed_open(prefix, i_loop, ncid)

    ! Open the processed state output file.

    use pmc_bin_grid
    use pmc_aero_data
    use pmc_aero_state
    use pmc_env
    use pmc_util
    use pmc_inout
    use pmc_gas_data
    use pmc_mpi
    use pmc_process_spec
    
    character(len=*), intent(in) :: prefix ! prefix of files to write
    integer, intent(in) :: i_loop       ! current loop number
    integer, intent(out) :: ncid        ! new NetCDF file ID, in data mode

    character(len=len(prefix)+20) :: filename
    character(len=20) :: date_str, time_str, zone_str
    character(len=500) :: history

    write(filename, '(a,a,i4.4,a)') trim(prefix), '_', i_loop, '.nc'
    call pmc_nc_check(nf90_create(filename, NF90_CLOBBER, ncid))

    call pmc_nc_check(nf90_put_att(ncid, NF90_GLOBAL, "title", &
         "PartMC output file"))
    call date_and_time(date_str, time_str, zone_str)
    write(history, '(15a)') date_str(1:4), "-", date_str(5:6), "-", &
         date_str(7:8), "T", time_str(1:2), ":", time_str(3:4), ":", &
         time_str(5:10), zone_str(1:3), ":", zone_str(4:5), &
         " created by PartMC"
    call pmc_nc_check(nf90_put_att(ncid, NF90_GLOBAL, "history", history))

    call pmc_nc_check(nf90_enddef(ncid))

  end subroutine output_processed_open

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine output_processed_close(ncid)

    ! Close the processed state output file.

    integer, intent(out) :: ncid        ! new NetCDF file ID, in data mode

    call pmc_nc_check(nf90_close(ncid))

  end subroutine output_processed_close

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine output_processed(ncid, prefix, process_spec_list, bin_grid, &
       aero_data, aero_state, gas_data, gas_state, env, index, time, &
       del_t, i_loop)

    ! Write the current processed state.

    use pmc_bin_grid
    use pmc_aero_data
    use pmc_aero_state
    use pmc_env
    use pmc_util
    use pmc_inout
    use pmc_gas_data
    use pmc_mpi
    use pmc_process_spec
    
    integer, intent(in) :: ncid         ! NetCDF file ID, in data mode
    character(len=*), intent(in) :: prefix ! prefix of files to write
    type(process_spec_t), intent(in) :: process_spec_list(:) ! processings specs
    type(bin_grid_t), intent(in) :: bin_grid ! bin grid
    type(aero_data_t), intent(in) :: aero_data ! aerosol data
    type(aero_state_t), intent(in) :: aero_state ! aerosol state
    type(gas_data_t), intent(in) :: gas_data ! gas data
    type(gas_state_t), intent(in) :: gas_state ! gas state
    type(env_t), intent(in) :: env      ! environment state
    integer, intent(in) :: index        ! filename index
    real*8, intent(in) :: time          ! current time (s)
    real*8, intent(in) :: del_t         ! current output time-step (s)
    integer, intent(in) :: i_loop       ! current loop number

    character(len=len(prefix)+20) :: basename

    write(basename, '(a,a,i4.4,a,i8.8)') trim(prefix), '_', i_loop, '_', index
    call process_time(ncid, time, index, del_t)
    call process_state_spec_list(ncid, basename, process_spec_list, &
         bin_grid, aero_data, aero_state, gas_data, gas_state, &
         env, time, index, del_t)

  end subroutine output_processed

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine output_binned(ncid, prefix, process_spec_list, bin_grid, &
       aero_data, aero_binned, gas_data, gas_state, env, index, time, del_t)

    ! Write the current binned data.
    
    use pmc_bin_grid
    use pmc_aero_data
    use pmc_aero_binned
    use pmc_env
    use pmc_inout
    use pmc_gas_data
    use pmc_gas_state
    use pmc_mpi
    use pmc_process_spec

    integer, intent(in) :: ncid         ! NetCDF file ID, in data mode
    character(len=*), intent(in) :: prefix ! prefix of files to write
    type(process_spec_t), intent(in) :: process_spec_list(:) ! processings specs
    type(bin_grid_t), intent(in) :: bin_grid ! bin grid
    type(aero_data_t), intent(in) :: aero_data ! aerosol data
    type(aero_binned_t), intent(in) :: aero_binned ! binned aerosol data
    type(gas_data_t), intent(in) :: gas_data ! gas data
    type(gas_state_t), intent(in) :: gas_state ! gas state
    type(env_t), intent(in) :: env      ! environment state
    integer, intent(in) :: index        ! filename index
    real*8, intent(in) :: time          ! current time (s)
    real*8, intent(in) :: del_t         ! current output time-step (s)

    character(len=(len(prefix)+30)) :: basename
    integer :: i

    call process_time(ncid, time, index, del_t)

    write(basename, '(a,a,i8.8)') trim(prefix), '_', index
    do i = 1,size(process_spec_list)
       if (process_spec_list(i)%type == "aero") then
          call output_aero(ncid, basename, process_spec_list(i)%suffix, &
               time, index, bin_grid, aero_data, aero_binned)
       end if
    end do

  end subroutine output_binned
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_output_summary
