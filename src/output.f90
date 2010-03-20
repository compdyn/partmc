! Copyright (C) 2005-2010 Nicole Riemer and Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_output module.

!> \page output_format Output NetCDF File Format
!!
!! PartMC output files are in the <a
!! href="http://www.unidata.ucar.edu/software/netcdf/">NetCDF Classic
!! Format</a> (also known as NetCDF-3 format). The dimensions and
!! variables in the files will depend on the type of run (particle,
!! analytical solution, etc), and options in the spec file (e.g. \c
!! record_removals). The output data comes in several different
!! groups, as follows:
!!
!! \subpage output_format_general "General Variables"
!!
!! \subpage output_format_env_state "Environment State"
!!
!! \subpage output_format_gas_data "Gas Material Data"
!!
!! \subpage output_format_gas_state "Gas State"
!!
!! \subpage output_format_aero_data "Aerosol Material Data"
!!
!! \subpage output_format_aero_state "Aerosol Particle State"
!!
!! \subpage output_format_aero_removed "Aerosol Particle Removal Information"
!!
!! \subpage output_format_bin_grid "Bin Grid Data"
!!
!! \subpage output_format_aero_binned "Aerosol Binned Sectional State"

!> Write data in NetCDF format.
module pmc_output

  use pmc_bin_grid
  use pmc_aero_data
  use pmc_aero_state
  use pmc_aero_binned
  use pmc_netcdf
  use pmc_env_state
  use pmc_util
  use pmc_gas_data
  use pmc_mpi
#ifdef PMC_USE_MPI
  use mpi
#endif

  integer, parameter :: TAG_OUTPUT_STATE_CENTRAL = 4341
  integer, parameter :: TAG_OUTPUT_STATE_SINGLE  = 4342
  
contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Write the current state.
  subroutine output_state(prefix, output_type, bin_grid, aero_data, &
       aero_weight, aero_state, gas_data, gas_state, env_state, index, &
       time, del_t, i_loop, record_removals)

    !> Prefix of state file.
    character(len=*), intent(in) :: prefix
    !> Output type for parallel runs (central/dist/single).
    character(len=*), intent(in) :: output_type
    !> Bin grid.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Aerosol weight.
    type(aero_weight_t), intent(in) :: aero_weight
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
    real(kind=dp), intent(in) :: time
    !> Current timestep (s).
    real(kind=dp), intent(in) :: del_t
    !> Current loop number.
    integer, intent(in) :: i_loop
    !> Whether to output particle removal info.
    logical, intent(in) :: record_removals
    
    integer :: rank, n_proc
#ifdef PMC_USE_MPI
    type(env_state_t) :: env_state_write
    type(gas_state_t) :: gas_state_write
    type(aero_state_t) :: aero_state_write
    integer :: ierr, status(MPI_STATUS_SIZE), buffer_size, i_proc, position
    character, allocatable :: buffer(:)
#endif

    rank = pmc_mpi_rank()
    n_proc = pmc_mpi_size()
    if (output_type == "central") then
       ! write per-processor data to separate files, but do it by
       ! transferring data to processor 0 and having it do the writes
       if (rank == 0) then
          call output_state_to_file(prefix, bin_grid, aero_data, &
               aero_weight, aero_state, gas_data, gas_state, env_state, &
               index, time, del_t, i_loop, record_removals, rank, n_proc)
#ifdef PMC_USE_MPI
          do i_proc = 1,(n_proc - 1)
             call recv_output_state_central(prefix, bin_grid, &
                  aero_data, aero_weight, gas_data, index, time, del_t, &
                  i_loop, record_removals, i_proc)
          end do
#endif
       else ! rank /= 0
          call send_output_state_central(aero_state, gas_state, env_state)
       end if
    elseif (output_type == "dist") then
       ! have each processor write its own data directly
       call output_state_to_file(prefix, bin_grid, aero_data, &
            aero_weight, aero_state, gas_data, gas_state, env_state, &
            index, time, del_t, i_loop, record_removals, rank, n_proc)
    elseif (output_type == "single") then
       if (n_proc == 1) then
          call output_state_to_file(prefix, bin_grid, aero_data, &
               aero_weight, aero_state, gas_data, gas_state, &
               env_state, index, time, del_t, i_loop, &
               record_removals, rank, n_proc)
       else
#ifdef PMC_USE_MPI
          ! collect all data onto processor 0 and then write it to a single file
          call env_state_allocate(env_state_write)
          call gas_state_allocate(gas_state_write)
          call env_state_copy(env_state, env_state_write)
          call gas_state_copy(gas_state, gas_state_write)
          call env_state_reduce_avg(env_state_write)
          call gas_state_reduce_avg(gas_state_write)
          !FIXME: use aero_state_mpi_gather() here.
          if (rank == 0) then
             call aero_state_allocate(aero_state_write)
             call aero_state_copy(aero_state, aero_state_write)
             do i_proc = 1,(n_proc - 1)
                call recv_output_state_single(aero_state_write, i_proc)
             end do
             call output_state_to_file(prefix, bin_grid, aero_data, &
                  aero_weight, aero_state_write, gas_data, gas_state_write, &
                  env_state_write, index, time, del_t, i_loop, &
                  record_removals, rank, 1)
             call aero_state_deallocate(aero_state_write)
          else ! rank /= 0
             call send_output_state_single(aero_state)
          end if
          call env_state_deallocate(env_state_write)
          call gas_state_deallocate(gas_state_write)
#endif
       end if
    else
       call die_msg(626743323, "Unknown output_type: " // trim(output_type))
    end if

  end subroutine output_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Helper routine to write various global attributes. Do not call
  !> directly.
  subroutine write_header_and_time(ncid, time, del_t, index)

    !> NetCDF file ID, in data mode.
    integer, intent(in) :: ncid
    !> Current time (s).
    real(kind=dp), intent(in) :: time
    !> Current timestep (s).
    real(kind=dp), intent(in) :: del_t
    !> Filename index.
    integer, intent(in) :: index

    character(len=500) :: history

    call pmc_nc_check(nf90_redef(ncid))

    call pmc_nc_check(nf90_put_att(ncid, NF90_GLOBAL, "source", &
         "PartMC version 2.0.0"))
    call iso8601_date_and_time(history)
    history((len_trim(history)+1):) = " created by PartMC"
    call pmc_nc_check(nf90_put_att(ncid, NF90_GLOBAL, "history", history))
    call pmc_nc_check(nf90_put_att(ncid, NF90_GLOBAL, "Conventions", "CF-1.4"))
    
    call pmc_nc_check(nf90_enddef(ncid))
    
    call pmc_nc_write_real(ncid, time, "time", unit="s", &
         description="time elapsed since simulation start")
    call pmc_nc_write_real(ncid, del_t, "timestep", unit="s", &
         description="current timestep size")
    call pmc_nc_write_integer(ncid, index, "timestep_index", &
         description="an integer that is 1 on the first timestep, " &
         // "2 on the second timestep, etc.")

  end subroutine write_header_and_time

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Write the current state for a single processor. Do not call this
  !> subroutine directly, but rather call output_state().
  subroutine output_state_to_file(prefix, bin_grid, aero_data, &
       aero_weight, aero_state, gas_data, gas_state, env_state, index, &
       time, del_t, i_loop, record_removals, write_rank, write_n_proc)

    !> Prefix of state file.
    character(len=*), intent(in) :: prefix
    !> Bin grid.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Aerosol weight.
    type(aero_weight_t), intent(in) :: aero_weight
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
    real(kind=dp), intent(in) :: time
    !> Current timestep (s).
    real(kind=dp), intent(in) :: del_t
    !> Current loop number.
    integer, intent(in) :: i_loop
    !> Whether to output particle removal info.
    logical, intent(in) :: record_removals
    !> Rank to write into file.
    integer, intent(in) :: write_rank
    !> Number of processors to write into file.
    integer, intent(in) :: write_n_proc
    
    character(len=len(prefix)+100) :: filename
    integer :: ncid

#ifdef PMC_USE_MPI
    if (write_n_proc > 1) then
       write(filename, '(a,a,i4.4,a,i4.4,a,i8.8,a)') trim(prefix), &
            '_', i_loop, '_', (write_rank + 1), '_', index, '.nc'
    else
    write(filename, '(a,a,i4.4,a,i8.8,a)') trim(prefix), &
            '_', i_loop, '_', index, '.nc'
    end if
#else
    write(filename, '(a,a,i4.4,a,i8.8,a)') trim(prefix), &
         '_', i_loop, '_', index, '.nc'
#endif
    call pmc_nc_check_msg(nf90_create(filename, NF90_CLOBBER, ncid), &
         "opening " // trim(filename))
    
    !> \page output_format_general Output NetCDF File Format: General Variables
    !!
    !! The general global attributes are:
    !!   - \b title: always set to the string "PartMC output file"
    !!   - \b source: set to the string "PartMC version V.V.V" where V.V.V
    !!     is the PartMC version that created the file
    !!   - \b Conventions: set to the string "CF-1.4", indicating
    !!     compliance with the <a
    !!     href="http://cf-pcmdi.llnl.gov/documents/cf-conventions/1.4">CF
    !!     convention format</a>
    !!   - \b history: set to the string
    !!     "YYYY-MM-DDThh:mm:ss[+-]ZZ:zz created by PartMC" where the first
    !!     term is the file creation time in the
    !!     <a href="http://en.wikipedia.org/wiki/ISO_8601">ISO 8601
    !!     format</a>. For example, noon Pacific Standard Time (PST) on
    !!     February 1st, 2000 would be written 2000-02-01T12:00:00-08:00.
    !!     The date and time variables are:
    !!     - YYYY: four-digit year
    !!     - MM: two-digit month number
    !!     - DD: two-digit day within month
    !!     - T: literal "T" character
    !!     - hh: two-digit hour in 24-hour format
    !!     - mm: two-digit minute
    !!     - ss: two-digit second
    !!     - [+-]: a literal "+" or "-" character giving the time zone
    !!       offset sign
    !!     - ZZ: two-digit hours of the time zone offset from UTC
    !!     - zz: two-digit minutes of the time zone offset from UTC
    !!
    !! The general variables are:
    !!   - \b time (s): time elapsed since the simulation start time, as
    !!     specified in the \ref output_format_env_state section
    !!   - \b timestep (s): the current timestep size
    !!   - \b loop: the loop number of this simulation (starting from 1)
    !!   - \b timestep_index: an integer that is 1 on the first timestep, 2
    !!     on the second timestep, etc.
    !!   - \b processor (MPI only): the processor number (starting from 1)
    !!     that output this data file
    !!   - \b total_processors (MPI only): the total number of processors

    call pmc_nc_check(nf90_put_att(ncid, NF90_GLOBAL, "title", &
         "PartMC output file"))
    call pmc_nc_check(nf90_enddef(ncid))

    call write_header_and_time(ncid, time, del_t, index)
    call pmc_nc_write_integer(ncid, i_loop, "loop", &
         description="loop repeat number of this simulation " &
         // "(starting from 1)")
#ifdef PMC_USE_MPI
    call pmc_nc_write_integer(ncid, write_rank + 1, "processor", &
         description="the processor number (starting from 1) " &
         "that output this data file")
    call pmc_nc_write_integer(ncid, write_n_proc, "total_processors", &
         description="total number of processors")
#endif

    call env_state_output_netcdf(env_state, ncid)
    call gas_data_output_netcdf(gas_data, ncid)
    call gas_state_output_netcdf(gas_state, ncid, gas_data)
    call aero_data_output_netcdf(aero_data, ncid)
    call aero_state_output_netcdf(aero_state, ncid, bin_grid, &
         aero_data, aero_weight, record_removals)

    call pmc_nc_check(nf90_close(ncid))
    
  end subroutine output_state_to_file

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Send the state for the "central" output method to the root processor.
  subroutine send_output_state_central(aero_state, gas_state, env_state)

    !> Aerosol state.
    type(aero_state_t), intent(in) :: aero_state
    !> Gas state.
    type(gas_state_t), intent(in) :: gas_state
    !> Environment state.
    type(env_state_t), intent(in) :: env_state

#ifdef PMC_USE_MPI
    integer :: buffer_size, position, ierr
    character, allocatable :: buffer(:)

    call assert(645797304, pmc_mpi_rank() /= 0)

    buffer_size = 0
    buffer_size = buffer_size + pmc_mpi_pack_size_env_state(env_state)
    buffer_size = buffer_size + pmc_mpi_pack_size_gas_state(gas_state)
    buffer_size = buffer_size + pmc_mpi_pack_size_aero_state(aero_state)
    allocate(buffer(buffer_size))
    position = 0
    call pmc_mpi_pack_env_state(buffer, position, env_state)
    call pmc_mpi_pack_gas_state(buffer, position, gas_state)
    call pmc_mpi_pack_aero_state(buffer, position, aero_state)
    call assert(839343839, position == buffer_size)
    call mpi_send(buffer, buffer_size, MPI_CHARACTER, 0, &
         TAG_OUTPUT_STATE_CENTRAL, MPI_COMM_WORLD, ierr)
    call pmc_mpi_check_ierr(ierr)
    deallocate(buffer)
#endif

  end subroutine send_output_state_central

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Receive the state for the "central" output method on the root
  !> processor.
  subroutine recv_output_state_central(prefix, bin_grid, aero_data, &
       aero_weight, gas_data, index, time, del_t, i_loop, record_removals, &
       remote_proc)

    !> Prefix of state file.
    character(len=*), intent(in) :: prefix
    !> Bin grid.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Aerosol weight.
    type(aero_weight_t), intent(in) :: aero_weight
    !> Gas data.
    type(gas_data_t), intent(in) :: gas_data
    !> Filename index.
    integer, intent(in) :: index
    !> Current time (s).
    real(kind=dp), intent(in) :: time
    !> Current timestep (s).
    real(kind=dp), intent(in) :: del_t
    !> Current loop number.
    integer, intent(in) :: i_loop
    !> Whether to output particle removal info.
    logical, intent(in) :: record_removals
    !> Processor number to receive from.
    integer, intent(in) :: remote_proc

#ifdef PMC_USE_MPI
    type(env_state_t) :: env_state
    type(gas_state_t) :: gas_state
    type(aero_state_t) :: aero_state
    integer :: buffer_size, position, status(MPI_STATUS_SIZE)
    integer :: n_proc, ierr
    character, allocatable :: buffer(:)

    call assert(206980035, pmc_mpi_rank() == 0)
    call assert(291452117, remote_proc /= 0)
    n_proc = pmc_mpi_size()

    ! get buffer size
    call mpi_probe(remote_proc, TAG_OUTPUT_STATE_CENTRAL, MPI_COMM_WORLD, &
         status, ierr)
    call pmc_mpi_check_ierr(ierr)
    call mpi_get_count(status, MPI_CHARACTER, buffer_size, ierr)

    ! get message
    allocate(buffer(buffer_size))
    call mpi_recv(buffer, buffer_size, MPI_CHARACTER, remote_proc, &
         TAG_OUTPUT_STATE_CENTRAL, MPI_COMM_WORLD, status, ierr)
    call pmc_mpi_check_ierr(ierr)

    ! unpack message
    position = 0
    call env_state_allocate(env_state)
    call gas_state_allocate(gas_state)
    call aero_state_allocate(aero_state)
    call pmc_mpi_unpack_env_state(buffer, position, env_state)
    call pmc_mpi_unpack_gas_state(buffer, position, gas_state)
    call pmc_mpi_unpack_aero_state(buffer, position, aero_state)
    call assert(279581330, position == buffer_size)
    deallocate(buffer)
    
    call output_state_to_file(prefix, bin_grid, aero_data, &
         aero_weight, aero_state, gas_data, gas_state, env_state, &
         index, time, del_t, i_loop, record_removals, remote_proc, n_proc)
    
    call env_state_deallocate(env_state)
    call gas_state_deallocate(gas_state)
    call aero_state_deallocate(aero_state)
#endif
    
  end subroutine recv_output_state_central

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Send the state for the "single" output method to the root processor.
  subroutine send_output_state_single(aero_state)

    !> Aerosol state.
    type(aero_state_t), intent(in) :: aero_state

#ifdef PMC_USE_MPI
    integer :: buffer_size, position, ierr
    character, allocatable :: buffer(:)

    call assert(699329210, pmc_mpi_rank() /= 0)

    buffer_size = 0
    buffer_size = buffer_size + pmc_mpi_pack_size_aero_state(aero_state)
    allocate(buffer(buffer_size))
    position = 0
    call pmc_mpi_pack_aero_state(buffer, position, aero_state)
    call assert(269866580, position == buffer_size)
    !>DEBUG
    !write(*,*) 'send_output_state_single: proc/n_part/bs = ', pmc_mpi_rank(), aero_state%n_part, buffer_size
    !if (pmc_mpi_rank() == 1) then
    !   write(*,*) 'proc 1, send buffer(01) = ', ichar(buffer(01))
    !   write(*,*) 'proc 1, send buffer(02) = ', ichar(buffer(02))
    !   write(*,*) 'proc 1, send buffer(03) = ', ichar(buffer(03))
    !   write(*,*) 'proc 1, send buffer(04) = ', ichar(buffer(04))
    !   write(*,*) 'proc 1, send buffer(05) = ', ichar(buffer(05))
    !   write(*,*) 'proc 1, send buffer(06) = ', ichar(buffer(06))
    !   write(*,*) 'proc 1, send buffer(07) = ', ichar(buffer(07))
    !   write(*,*) 'proc 1, send buffer(08) = ', ichar(buffer(08))
    !   write(*,*) 'proc 1, send buffer(09) = ', ichar(buffer(09))
    !   write(*,*) 'proc 1, send buffer(10) = ', ichar(buffer(10))
    !   write(*,*) 'proc 1, send buffer(11) = ', ichar(buffer(11))
    !   write(*,*) 'proc 1, send buffer(12) = ', ichar(buffer(12))
    !end if
    !<DEBUG
    call mpi_send(buffer, buffer_size, MPI_CHARACTER, 0, &
         TAG_OUTPUT_STATE_SINGLE, MPI_COMM_WORLD, ierr)
    call pmc_mpi_check_ierr(ierr)
    deallocate(buffer)
#endif

  end subroutine send_output_state_single

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Receive the state for the "central" output method on the root
  !> processor and add it to the given \c aero_state.
  subroutine recv_output_state_single(aero_state, remote_proc)

    !> Aerosol state to add sent state to.
    type(aero_state_t), intent(inout) :: aero_state
    !> Remote processor to receive state from.
    integer, intent(in) :: remote_proc

#ifdef PMC_USE_MPI
    type(aero_state_t) :: aero_state_remote
    integer :: buffer_size, position, status(MPI_STATUS_SIZE), ierr
    character, allocatable :: buffer(:)

    call assert(206980035, pmc_mpi_rank() == 0)
    call assert(291452117, remote_proc /= 0)

    ! get buffer size
    call mpi_probe(remote_proc, TAG_OUTPUT_STATE_SINGLE, MPI_COMM_WORLD, &
         status, ierr)
    call pmc_mpi_check_ierr(ierr)
    call mpi_get_count(status, MPI_CHARACTER, buffer_size, ierr)

    ! get message
    allocate(buffer(buffer_size))
    call mpi_recv(buffer, buffer_size, MPI_CHARACTER, remote_proc, &
         TAG_OUTPUT_STATE_SINGLE, MPI_COMM_WORLD, status, ierr)
    call pmc_mpi_check_ierr(ierr)
    !>DEBUG
    !write(*,*) 'recv_output_state_single: received from proc ', status(MPI_SOURCE)
    !write(*,*) 'recv_output_state_single: buffer_size = ', buffer_size
    !write(*,*) 'proc 0, recv buffer(01) = ', ichar(buffer(01))
    !write(*,*) 'proc 0, recv buffer(02) = ', ichar(buffer(02))
    !write(*,*) 'proc 0, recv buffer(03) = ', ichar(buffer(03))
    !write(*,*) 'proc 0, recv buffer(04) = ', ichar(buffer(04))
    !write(*,*) 'proc 0, recv buffer(05) = ', ichar(buffer(05))
    !write(*,*) 'proc 0, recv buffer(06) = ', ichar(buffer(06))
    !write(*,*) 'proc 0, recv buffer(07) = ', ichar(buffer(07))
    !write(*,*) 'proc 0, recv buffer(08) = ', ichar(buffer(08))
    !write(*,*) 'proc 0, recv buffer(09) = ', ichar(buffer(09))
    !write(*,*) 'proc 0, recv buffer(10) = ', ichar(buffer(10))
    !write(*,*) 'proc 0, recv buffer(11) = ', ichar(buffer(11))
    !write(*,*) 'proc 0, recv buffer(12) = ', ichar(buffer(12))
    !<DEBUG

    ! unpack message
    position = 0
    call aero_state_allocate(aero_state_remote)
    call pmc_mpi_unpack_aero_state(buffer, position, aero_state_remote)
    call assert(466510584, position == buffer_size)
    deallocate(buffer)

    !>DEBUG
    !write(*,*) 'recv_output_state_single: received from proc ', status(MPI_SOURCE)
    !write(*,*) 'recv_output_state_single: n_part_remote ', aero_state_remote%n_part
    !write(*,*) 'recv_output_state_single: calling aero_state_add'
    !<DEBUG
    call aero_state_add(aero_state, aero_state_remote)
    
    call aero_state_deallocate(aero_state_remote)
#endif

  end subroutine recv_output_state_single

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read the current state.
  subroutine input_state(filename, bin_grid, aero_data, &
       aero_weight, aero_state, gas_data, gas_state, env_state, index, &
       time, del_t, i_loop)

    !> Prefix of state file.
    character(len=*), intent(in) :: filename
    !> Bin grid.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Aerosol data.
    type(aero_data_t), intent(out) :: aero_data
    !> Aerosol weight.
    type(aero_weight_t), intent(in) :: aero_weight
    !> Aerosol state.
    type(aero_state_t), intent(out) :: aero_state
    !> Gas data.
    type(gas_data_t), intent(out) :: gas_data
    !> Gas state.
    type(gas_state_t), intent(out) :: gas_state
    !> Environment state.
    type(env_state_t), intent(out) :: env_state
    !> Filename index.
    integer, intent(out) :: index
    !> Current time (s).
    real(kind=dp), intent(out) :: time
    !> Current timestep (s).
    real(kind=dp), intent(out) :: del_t
    !> Current loop number.
    integer, intent(out) :: i_loop
    
    integer :: ncid
    character(len=1000) :: unit

    ! only root node actually reads from the file
    if (pmc_mpi_rank() == 0) then
       call pmc_nc_open_read(filename, ncid)

       call pmc_nc_read_real(ncid, time, "time", unit)
       call pmc_nc_read_real(ncid, del_t, "timestep", unit)
       call pmc_nc_read_integer(ncid, i_loop, "loop", unit)
       call pmc_nc_read_integer(ncid, index, "timestep_index", unit)

       call env_state_input_netcdf(env_state, ncid)
       call gas_data_input_netcdf(gas_data, ncid)
       call gas_state_input_netcdf(gas_state, ncid, gas_data)
       call aero_data_input_netcdf(aero_data, ncid)
       call aero_state_input_netcdf(aero_state, ncid, bin_grid, &
            aero_data, aero_weight)

       call pmc_nc_close(ncid)
    end if
    
  end subroutine input_state
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Write the current sectional data.
  subroutine output_sectional(prefix, bin_grid, aero_data, &
       aero_binned, gas_data, gas_state, env_state, index, time, del_t)
    
    !> Prefix of filename to write
    character(len=*), intent(in) :: prefix
    !> Bin grid.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Binned aerosol data.
    type(aero_binned_t), intent(in) :: aero_binned
    !> Gas data.
    type(gas_data_t), intent(in) :: gas_data
    !> Gas state.
    type(gas_state_t), intent(in) :: gas_state
    !> Environment state.
    type(env_state_t), intent(in) :: env_state
    !> Filename index.
    integer, intent(in) :: index
    !> Current time (s).
    real(kind=dp), intent(in) :: time
    !> Current output time-step (s).
    real(kind=dp), intent(in) :: del_t

    integer :: ncid
    character(len=len(prefix)+100) :: filename

    write(filename, '(a,a,i8.8,a)') trim(prefix), &
         '_', index, '.nc'
    call pmc_nc_check_msg(nf90_create(filename, NF90_CLOBBER, ncid), &
         "opening " // trim(filename))

    ! write header attributes
    call pmc_nc_check(nf90_put_att(ncid, NF90_GLOBAL, "title", &
         "PartMC sectional output file"))
    call pmc_nc_check(nf90_enddef(ncid))

    call write_header_and_time(ncid, time, del_t, index)

    ! write data
    call env_state_output_netcdf(env_state, ncid)
    call gas_data_output_netcdf(gas_data, ncid)
    call gas_state_output_netcdf(gas_state, ncid, gas_data)
    call aero_data_output_netcdf(aero_data, ncid)
    call aero_binned_output_netcdf(aero_binned, ncid, bin_grid, &
         aero_data)

    call pmc_nc_check(nf90_close(ncid))

  end subroutine output_sectional
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
end module pmc_output
