! Copyright (C) 2005-2008 Nicole Riemer and Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_output_state module.

!> Save and restore the exact internal state (a checkpoint).
!!
!! The state file will contain enough data to restart the simulation
!! at the point it was written.
!!
!! Because it contains the full state of every particle, this is also
!! the best way to gain complete access to all statistics of the
!! simulation.
module pmc_output_state

  use pmc_bin_grid
  use pmc_aero_data
  use pmc_aero_state
  use pmc_env_state
  use pmc_util
  use pmc_inout
  use pmc_gas_data
  use pmc_mpi
  
contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Write the current state.
  subroutine inout_write_state(state_prefix, bin_grid, aero_data, &
       aero_state, gas_data, gas_state, env_state, index, time, del_t, i_loop)

    !> Prefix of state file.
    character(len=*), intent(in) :: state_prefix
    !> Bin grid.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Aerosol state.
    type(aero_state_t), intent(in) :: aero_state
    !> Gas data.
    type(gas_data_t), intent(in) :: gas_data
    !> Gas state.
    type(gas_state_t), intent(in) :: gas_state
    !> Environment state.
    type(env_state_t), intent(in) :: env_state
    !> Filename index.
    integer, intent(in) :: index
    !> Current time (s).
    real*8, intent(in) :: time
    !> Current timestep (s).
    real*8, intent(in) :: del_t
    !> Current loop number.
    integer, intent(in) :: i_loop
    
    character*300 :: filename
    type(inout_file_t) :: file
    type(env_state_t) :: env_write
    type(gas_state_t) :: gas_state_write
    type(aero_state_t) :: aero_state_write
    integer :: ierr, status, buffer_size, i_proc, position
    character, allocatable :: buffer(:)

    ! only root node actually writes to the file
    if (pmc_mpi_rank() == 0) then
       ! write all the common data
       write(filename, '(a,a,i4.4,a,i8.8,a)') trim(state_prefix), &
            '_', i_loop, '_', index, '.d'
       call inout_open_write(filename, file)
       
       call inout_write_real(file, 'time(s)', time)
       call inout_write_real(file, 'timestep(s)', del_t)
       call inout_write_integer(file, 'loop', i_loop)
       call inout_write_integer(file, 'index', index)

       call inout_write_bin_grid(file, bin_grid)
       call inout_write_gas_data(file, gas_data)
       call inout_write_aero_data(file, aero_data)

       ! write root node's state
       call inout_write_integer(file, 'n_processor', pmc_mpi_size())
       call inout_write_integer(file, 'processor', 0)
       call inout_write_env_state(file, env_state)
       call inout_write_gas_state(file, gas_state)
       call inout_write_aero_state(file, aero_state)
    end if

#ifdef PMC_USE_MPI
    ! write everyone else's state
    do i_proc = 1,(pmc_mpi_size() - 1)
       call pmc_mpi_barrier()
       ! compute and send buffer_size from remote node
       if (pmc_mpi_rank() == i_proc) then
          buffer_size = 0
          buffer_size = buffer_size + pmc_mpi_pack_size_env_state(env_state)
          buffer_size = buffer_size + pmc_mpi_pack_size_gas_state(gas_state)
          buffer_size = buffer_size + pmc_mpi_pack_size_aero_state(aero_state)
          call mpi_send(buffer_size, 1, MPI_INTEGER, 0, 57, &
               MPI_COMM_WORLD, ierr)
          call pmc_mpi_check_ierr(ierr)
       end if
       ! get buffer_size at root node
       if (pmc_mpi_rank() == 0) then
          call mpi_recv(buffer_size, 1, MPI_INTEGER, i_proc, 57, &
               MPI_COMM_WORLD, status, ierr)
          call pmc_mpi_check_ierr(ierr)
       end if
       ! send buffer from remote node
       if (pmc_mpi_rank() == i_proc) then
          allocate(buffer(buffer_size))
          position = 0
          call pmc_mpi_pack_env_state(buffer, position, env_state)
          call pmc_mpi_pack_gas_state(buffer, position, gas_state)
          call pmc_mpi_pack_aero_state(buffer, position, aero_state)
          call assert(542772170, position == buffer_size)
          call mpi_send(buffer, buffer_size, MPI_CHARACTER, 0, 58, &
               MPI_COMM_WORLD, ierr)
          call pmc_mpi_check_ierr(ierr)
          deallocate(buffer)
       end if
       ! get buffer at root node
       if (pmc_mpi_rank() == 0) then
          allocate(buffer(buffer_size))
          call mpi_recv(buffer, buffer_size, MPI_CHARACTER, i_proc, 58, &
               MPI_COMM_WORLD, status, ierr)
          position = 0
          call pmc_mpi_unpack_env_state(buffer, position, env_write)
          call pmc_mpi_unpack_gas_state(buffer, position, gas_state_write)
          call pmc_mpi_unpack_aero_state(buffer, position, aero_state_write)
          call assert(518174881, position == buffer_size)
          deallocate(buffer)
       end if
       ! process state at root node
       if (pmc_mpi_rank() == 0) then
          call inout_write_integer(file, 'processor', i_proc)
          call inout_write_env_state(file, env_write)
          call inout_write_gas_state(file, gas_state_write)
          call inout_write_aero_state(file, aero_state_write)
          
          call env_state_free(env_write)
          call gas_state_free(gas_state_write)
          call aero_state_free(aero_state_write)
       end if
    end do
#endif
       
    if (pmc_mpi_rank() == 0) then
       call inout_close(file)
    end if
    
  end subroutine inout_write_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read the current state.
  subroutine inout_read_state(state_name, bin_grid, aero_data, &
       aero_state, gas_data, gas_state, env_state, time, index, del_t, i_loop)

    !> Name of state file.
    character(len=*), intent(in) :: state_name
    !> Bin grid.
    type(bin_grid_t), intent(out) :: bin_grid
    !> Aerosol data.
    type(aero_data_t), intent(out) :: aero_data
    !> Aerosol state.
    type(aero_state_t), intent(out) :: aero_state
    !> Gas data.
    type(gas_data_t), intent(out) :: gas_data
    !> Gas state.
    type(gas_state_t), intent(out) :: gas_state
    !> Environment state.
    type(env_state_t), intent(out) :: env_state
    !> Current time (s).
    real*8, intent(out) :: time
    !> Current index.
    integer, intent(out) :: index
    !> Current time-step (s).
    real*8, intent(out) :: del_t
    !> Current loop number.
    integer, intent(out) :: i_loop
    
    type(inout_file_t) :: file
    integer :: n_proc, i_proc, check_i_proc
    type(env_state_t) :: env_read
    type(gas_state_t) :: gas_state_read
    type(aero_state_t) :: aero_state_read

    if (pmc_mpi_rank() /= 0) then
       call pmc_mpi_abort(52115)
    end if

    call inout_open_read(state_name, file)
    
    call inout_read_real(file, 'time(s)', time)
!DEBUG
    write(*,*) 'time = ', time
!DEBUG
    call inout_read_real(file, 'timestep(s)', del_t)
!DEBUG
    write(*,*) 'timestep = ', del_t
!DEBUG
    call inout_read_integer(file, 'loop', i_loop)
!DEBUG
    write(*,*) 'loop = ', i_loop
!DEBUG
    call inout_read_integer(file, 'index', index)
!DEBUG
    write(*,*) 'index = ', index
!DEBUG

!DEBUG
    write(*,*) 'about to read bin_grid'
!DEBUG
    call inout_read_bin_grid(file, bin_grid)
!DEBUG
    write(*,*) 'read bin_grid'
!DEBUG
    call inout_read_gas_data(file, gas_data)
    call inout_read_aero_data(file, aero_data)

    ! read root node's state
    call inout_read_integer(file, 'n_processor', n_proc)
    call inout_read_integer(file, 'processor', check_i_proc)
    call inout_check_index(file, 0, check_i_proc)
    call inout_read_env_state(file, env_state)
    call inout_read_gas_state(file, gas_state)
    call inout_read_aero_state(file, aero_state)

    ! read other nodes' states
    do i_proc = 1,(n_proc - 1)
       call inout_read_integer(file, 'processor', check_i_proc)
       call inout_check_index(file, i_proc, check_i_proc)
       call inout_read_env_state(file, env_read)
       call inout_read_gas_state(file, gas_state_read)
       call inout_read_aero_state(file, aero_state_read)

       call env_state_add(env_state, env_read)
       call gas_state_add(gas_state, gas_state_read)
       call aero_state_add(aero_state, aero_state_read)
       
       call env_state_free(env_read)
       call gas_state_free(gas_state_read)
       call aero_state_free(aero_state_read)
    end do

    ! average data
    call env_state_scale(env_state, dble(n_proc))
    call gas_state_scale(gas_state, dble(n_proc))
    
    call inout_close(file)
    
  end subroutine inout_read_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_output_state
