! Copyright (C) 2005-2022 Nicole Riemer and Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_output module.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> \page output_format Output File Format
!!
!! PartMC output files are in the <a
!! href="http://www.unidata.ucar.edu/software/netcdf/">NetCDF Classic
!! Format</a> (also known as NetCDF-3 format). The dimensions and
!! variables in the files will depend on the type of run (particle,
!! analytical solution, etc), and options in the spec file (e.g. \c
!! record_removals and \c do_optical).
!!
!! The state of the simulation is periodically output during the run,
!! with frequency determined by the \c t_output input parameter. Each
!! output file has a filename of the form \c PREFIX_RRRR_SSSSSSSS.nc,
!! where \c PREFIX is given by the \c output_prefix input parameter,
!! \c RRRR is the four-digit repeat number (starting from 1), and \c
!! SSSSSSSS is the eight-digit output index (starting at 1 and
!! incremented each time the state is output). For exact and sectional
!! simulations all repeats would be identical so there is no support
!! for repeating and the filename is of the format \c
!! PREFIX_SSSSSSSS.nc.
!!
!! If run in parallel and \c output_type is \c central or \c dist,
!! then the output files have names like \c
!! PREFIX_RRRR_PPPP_SSSSSSSS.nc, where \c PPPP is a four-digit
!! process number (starting from 1) and the other variables are as
!! above. If \c output_type is \c single then the output file naming
!! scheme as the same as for serial runs.
!!
!! The data in each output file comes in several different groups, as
!! follows:
!!
!! \subpage output_format_general "General Information"
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
!! (only for particle-resolved simulations)
!!
!! \subpage output_format_aero_removed "Aerosol Particle Removal Information"
!! (only for particle-resolved simulations, if \c record_removals is \c yes)
!!
!! \subpage output_format_aero_weight_array "Aerosol Weighting Function"
!! (only for particle-resolved simulations)
!!
!! \subpage output_format_diam_bin_grid "Diameter Bin Grid Data"
!! (only for exact and sectional simulations)
!!
!! \subpage output_format_aero_binned "Aerosol Binned Sectional State"
!! (only for exact and sectional simulations)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> Write data in NetCDF format.
module pmc_output

  use pmc_bin_grid
  use pmc_aero_data
  use pmc_aero_state
  use pmc_aero_binned
  use pmc_netcdf
  use pmc_gas_state
  use pmc_env_state
  use pmc_util
  use pmc_gas_data
  use pmc_mpi
#ifdef PMC_USE_MPI
  use mpi
#endif

  !> PartMC verson number.
  character(len=100), parameter :: PARTMC_VERSION = "2.6.1"

  !> Type code for undefined or invalid output.
  integer, parameter :: OUTPUT_TYPE_INVALID = 0
  !> Type code for centralized output (one file per process, all written
  !> by process 0).
  integer, parameter :: OUTPUT_TYPE_CENTRAL = 1
  !> Type code for distributed output (one file per process, written by
  !> each process).
  integer, parameter :: OUTPUT_TYPE_DIST    = 2
  !> Type code for single output (one file for all processes, written by
  !> process 0).
  integer, parameter :: OUTPUT_TYPE_SINGLE  = 3

  !> Internal-use variable only.
  integer, parameter :: TAG_OUTPUT_STATE_CENTRAL = 4341
  !> Internal-use variable only.
  integer, parameter :: TAG_OUTPUT_STATE_SINGLE  = 4342

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Write the current state.
  subroutine output_state(prefix, output_type, aero_data, aero_state, &
       gas_data, gas_state, env_state, index, time, del_t, i_repeat, &
       record_removals, record_optical, uuid)

    !> Prefix of state file.
    character(len=*), intent(in) :: prefix
    !> Output type for parallel runs (see module constants).
    integer, intent(in) :: output_type
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
    real(kind=dp), intent(in) :: time
    !> Current timestep (s).
    real(kind=dp), intent(in) :: del_t
    !> Current repeat number.
    integer, intent(in) :: i_repeat
    !> Whether to output particle removal info.
    logical, intent(in) :: record_removals
    !> Whether to output aerosol optical properties.
    logical, intent(in) :: record_optical
    !> UUID of the simulation.
    character(len=PMC_UUID_LEN), intent(in) :: uuid

    integer :: rank, n_proc
#ifdef PMC_USE_MPI
    type(env_state_t) :: env_state_write
    type(gas_state_t) :: gas_state_write
    type(aero_state_t) :: aero_state_write
    integer :: ierr, status(MPI_STATUS_SIZE), i_proc, position
    character, allocatable :: buffer(:)
#endif

    rank = pmc_mpi_rank()
    n_proc = pmc_mpi_size()
    if (output_type == OUTPUT_TYPE_CENTRAL) then
       ! write per-process data to separate files, but do it by
       ! transferring data to process 0 and having it do the writes
       if (rank == 0) then
          call output_state_to_file(prefix, aero_data, aero_state, gas_data, &
               gas_state, env_state, index, time, del_t, i_repeat, &
               record_removals, record_optical, uuid, rank, n_proc)
#ifdef PMC_USE_MPI
          do i_proc = 1,(n_proc - 1)
             call recv_output_state_central(prefix, aero_data, gas_data, &
                  index, time, del_t, i_repeat, record_removals, &
                  record_optical, uuid, i_proc)
          end do
#endif
       else ! rank /= 0
          call send_output_state_central(aero_state, gas_state, env_state)
       end if
    elseif (output_type == OUTPUT_TYPE_DIST) then
       ! have each process write its own data directly
       call output_state_to_file(prefix, aero_data, aero_state, gas_data, &
            gas_state, env_state, index, time, del_t, i_repeat, &
            record_removals, record_optical, uuid, rank, n_proc)
    elseif (output_type == OUTPUT_TYPE_SINGLE) then
       if (n_proc == 1) then
          call output_state_to_file(prefix, aero_data, aero_state, gas_data, &
               gas_state, env_state, index, time, del_t, i_repeat, &
               record_removals, record_optical, uuid, rank, n_proc)
       else
#ifdef PMC_USE_MPI
          ! collect all data onto process 0 and then write it to a
          ! single file
          env_state_write = env_state
          gas_state_write = gas_state
          call env_state_reduce_avg(env_state_write)
          call gas_state_reduce_avg(gas_state_write)
          call aero_state_mpi_gather(aero_state, aero_state_write, aero_data)
          if (rank == 0) then
             call output_state_to_file(prefix, aero_data, aero_state_write, &
                  gas_data, gas_state_write, env_state_write, index, time, &
                  del_t, i_repeat, record_removals, record_optical, uuid, &
                  rank, 1)
          end if
#endif
       end if
    else
       call die_msg(626743323, "Unknown output_type: " &
            // trim(integer_to_string(output_type)))
    end if

  end subroutine output_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Make a filename from a given prefix and other information.
  subroutine make_filename(filename, prefix, suffix, index, i_repeat, &
       write_rank, write_n_proc)

    !> Filename to create.
    character(len=*), intent(out) :: filename
    !> Filename prefix.
    character(len=*), intent(in) :: prefix
    !> Filename suffix.
    character(len=*), intent(in) :: suffix
    !> Filename index.
    integer, intent(in), optional :: index
    !> Current repeat number.
    integer, intent(in), optional :: i_repeat
    !> Rank to write into file.
    integer, intent(in), optional :: write_rank
    !> Number of processes to write into file.
    integer, intent(in), optional :: write_n_proc

    integer :: ncid, use_rank, use_n_proc
    character(len=100) :: proc_string, index_string, repeat_string

    if (present(write_rank)) then
       use_rank = write_rank
    else
       use_rank = pmc_mpi_rank()
    end if
    if (present(write_n_proc)) then
       use_n_proc = write_n_proc
    else
       use_n_proc = pmc_mpi_size()
    end if

    repeat_string = ""
    proc_string = ""
    index_string = ""
    if (present(i_repeat)) write(repeat_string, '(a,i4.4)') '_', i_repeat
    if (use_n_proc > 1) write(proc_string, '(a,i4.4)') '_', (use_rank + 1)
    if (present(index)) write(index_string, '(a,i8.8)') '_', index
    write(filename, '(a,a,a,a,a)') trim(prefix), trim(repeat_string), &
         trim(proc_string), trim(index_string), trim(suffix)

  end subroutine make_filename

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Helper routine to write time variables. Do not call directly.
  subroutine write_time(ncid, time, del_t, index)

    !> NetCDF file ID, in data mode.
    integer, intent(in) :: ncid
    !> Current time (s).
    real(kind=dp), intent(in) :: time
    !> Current timestep (s).
    real(kind=dp), intent(in) :: del_t
    !> Filename index.
    integer, intent(in) :: index

    call pmc_nc_write_real(ncid, time, "time", unit="s", &
         description="time elapsed since simulation start")
    call pmc_nc_write_real(ncid, del_t, "timestep", unit="s", &
         description="current timestep size")
    call pmc_nc_write_integer(ncid, index, "timestep_index", &
         description="an integer that is 1 on the first timestep, " &
         // "2 on the second timestep, etc.")

  end subroutine write_time

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Write the current state for a single process. Do not call this
  !> subroutine directly, but rather call output_state().
  subroutine output_state_to_file(prefix, aero_data, aero_state, gas_data, &
       gas_state, env_state, index, time, del_t, i_repeat, record_removals, &
       record_optical, uuid, write_rank, write_n_proc)

    !> Prefix of state file.
    character(len=*), intent(in) :: prefix
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
    real(kind=dp), intent(in) :: time
    !> Current timestep (s).
    real(kind=dp), intent(in) :: del_t
    !> Current repeat number.
    integer, intent(in) :: i_repeat
    !> Whether to output particle removal info.
    logical, intent(in) :: record_removals
    !> Whether to output aerosol optical properties.
    logical, intent(in) :: record_optical
    !> UUID of the simulation.
    character(len=PMC_UUID_LEN), intent(in) :: uuid
    !> Rank to write into file.
    integer, intent(in), optional :: write_rank
    !> Number of processes to write into file.
    integer, intent(in), optional :: write_n_proc

    character(len=len(prefix)+100) :: filename
    integer :: ncid

    !> \page output_format_general Output File Format: General Information
    !!
    !! The general information global NetCDF attributes are:
    !!   - \b title: always set to the string "PartMC version V.V.V
    !!     output file" where V.V.V is the PartMC version that created
    !!     the file
    !!   - \b source: set to the string "PartMC version V.V.V"
    !!   - \b UUID: a string of the form F47AC10B-58CC-4372-A567-0E02B2C3D479
    !!     which is the same for all files generated by a single call of
    !!     PartMC.
    !!   - \b Conventions: set to the string "CF-1.4", indicating
    !!     compliance with the <a
    !!     href="http://cf-pcmdi.llnl.gov/documents/cf-conventions/1.4">CF
    !!     convention format</a>
    !!   - \b history: set to the string "YYYY-MM-DDThh:mm:ss[+-]ZZ:zz
    !!     created by PartMC version V.V.V" where the first term is
    !!     the file creation time in the <a
    !!     href="http://en.wikipedia.org/wiki/ISO_8601">ISO 8601
    !!     format</a>. For example, noon Pacific Standard Time (PST)
    !!     on February 1st, 2000 would be written
    !!     2000-02-01T12:00:00-08:00.  The date and time variables
    !!     are:
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
    !! The general information NetCDF variables are:
    !!   - \b time (unit s): time elapsed since the simulation start time,
    !!     as specified in the \ref output_format_env_state section
    !!   - \b timestep (unit s): the current timestep size
    !!   - \b repeat: the repeat number of this simulation (starting from 1)
    !!   - \b timestep_index: an integer that is 1 on the first timestep, 2
    !!     on the second timestep, etc.
    !!   - \b process (MPI only): the process number (starting from 1)
    !!     that output this data file
    !!   - \b total_processes (MPI only): the total number of processes
    !!     involved in writing data (may be less than the total number of
    !!     processes that computed the data)

    call make_filename(filename, prefix, ".nc", index, i_repeat, write_rank, &
         write_n_proc)
    call pmc_nc_open_write(filename, ncid)
    call pmc_nc_write_info(ncid, uuid, &
         "PartMC version " // trim(PARTMC_VERSION), write_rank, write_n_proc)
    call write_time(ncid, time, del_t, index)
    call pmc_nc_write_integer(ncid, i_repeat, "repeat", &
         description="repeat number of this simulation (starting from 1)")

    call env_state_output_netcdf(env_state, ncid)
    call gas_data_output_netcdf(gas_data, ncid)
    call gas_state_output_netcdf(gas_state, ncid, gas_data)
    call aero_data_output_netcdf(aero_data, ncid)
    call aero_state_output_netcdf(aero_state, ncid, aero_data, &
         record_removals, record_optical)

    call pmc_nc_check(nf90_close(ncid))

  end subroutine output_state_to_file

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Send the state for the "central" output method to the root process.
  subroutine send_output_state_central(aero_state, gas_state, env_state)

    !> Aerosol state.
    type(aero_state_t), intent(in) :: aero_state
    !> Gas state.
    type(gas_state_t), intent(in) :: gas_state
    !> Environment state.
    type(env_state_t), intent(in) :: env_state

#ifdef PMC_USE_MPI
    integer :: buffer_size, max_buffer_size, position, ierr
    character, allocatable :: buffer(:)

    call assert(645797304, pmc_mpi_rank() /= 0)

    max_buffer_size = 0
    max_buffer_size = max_buffer_size &
         + pmc_mpi_pack_size_env_state(env_state)
    max_buffer_size = max_buffer_size &
         + pmc_mpi_pack_size_gas_state(gas_state)
    max_buffer_size = max_buffer_size &
         + pmc_mpi_pack_size_aero_state(aero_state)
    allocate(buffer(max_buffer_size))
    position = 0
    call pmc_mpi_pack_env_state(buffer, position, env_state)
    call pmc_mpi_pack_gas_state(buffer, position, gas_state)
    call pmc_mpi_pack_aero_state(buffer, position, aero_state)
    call assert(839343839, position <= max_buffer_size)
    buffer_size = position
    call mpi_send(buffer, buffer_size, MPI_CHARACTER, 0, &
         TAG_OUTPUT_STATE_CENTRAL, MPI_COMM_WORLD, ierr)
    call pmc_mpi_check_ierr(ierr)
    deallocate(buffer)
#endif

  end subroutine send_output_state_central

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Receive the state for the "central" output method on the root
  !> process.
  subroutine recv_output_state_central(prefix, aero_data, gas_data, index, &
       time, del_t, i_repeat, record_removals, record_optical, uuid, &
       remote_proc)

    !> Prefix of state file.
    character(len=*), intent(in) :: prefix
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Gas data.
    type(gas_data_t), intent(in) :: gas_data
    !> Filename index.
    integer, intent(in) :: index
    !> Current time (s).
    real(kind=dp), intent(in) :: time
    !> Current timestep (s).
    real(kind=dp), intent(in) :: del_t
    !> Current repeat number.
    integer, intent(in) :: i_repeat
    !> Whether to output particle removal info.
    logical, intent(in) :: record_removals
    !> Whether to output aerosol_optical_properties.
    logical, intent(in) :: record_optical
    !> UUID of the simulation.
    character(len=PMC_UUID_LEN), intent(in) :: uuid
    !> Process number to receive from.
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
    call pmc_mpi_unpack_env_state(buffer, position, env_state)
    call pmc_mpi_unpack_gas_state(buffer, position, gas_state)
    call pmc_mpi_unpack_aero_state(buffer, position, aero_state)
    call assert(279581330, position == buffer_size)
    deallocate(buffer)

    call output_state_to_file(prefix, aero_data, aero_state, gas_data, &
         gas_state, env_state, index, time, del_t, i_repeat, &
         record_removals, record_optical, uuid, remote_proc, n_proc)
#endif

  end subroutine recv_output_state_central

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read the current state.
  subroutine input_state(filename, index, time, del_t, i_repeat, uuid, &
       aero_data, aero_state, gas_data, gas_state, env_state)

    !> Prefix of state file.
    character(len=*), intent(in) :: filename
    !> Filename index.
    integer, intent(out) :: index
    !> Current time (s).
    real(kind=dp), intent(out) :: time
    !> Current timestep (s).
    real(kind=dp), intent(out) :: del_t
    !> Current repeat number.
    integer, intent(out) :: i_repeat
    !> UUID of the simulation.
    character(len=PMC_UUID_LEN), intent(out) :: uuid
    !> Aerosol data.
    type(aero_data_t), optional, intent(inout) :: aero_data
    !> Aerosol state.
    type(aero_state_t), optional, intent(inout) :: aero_state
    !> Gas data.
    type(gas_data_t), optional, intent(inout) :: gas_data
    !> Gas state.
    type(gas_state_t), optional, intent(inout) :: gas_state
    !> Environment state.
    type(env_state_t), optional, intent(inout) :: env_state

    integer :: ncid

    call assert_msg(819739354, pmc_mpi_rank() == 0, &
         "can only call from process 0")

    call pmc_nc_open_read(filename, ncid)

    call pmc_nc_check(nf90_get_att(ncid, NF90_GLOBAL, "UUID", uuid))

    call pmc_nc_read_real(ncid, time, "time")
    call pmc_nc_read_real(ncid, del_t, "timestep")
    call pmc_nc_read_integer(ncid, i_repeat, "repeat")
    call pmc_nc_read_integer(ncid, index, "timestep_index")

    if (present(aero_data)) then
       call aero_data_input_netcdf(aero_data, ncid)
       if (present(aero_state)) then
          call aero_state_input_netcdf(aero_state, ncid, aero_data)
       end if
    else
       call assert_msg(289621231, present(aero_state) .eqv. .false., &
            "cannot input aero_state without aero_data")
    end if

    if (present(gas_data)) then
       call gas_data_input_netcdf(gas_data, ncid)
       if (present(gas_state)) then
          call gas_state_input_netcdf(gas_state, ncid, gas_data)
       end if
    else
       call assert_msg(874298496, present(gas_state) .eqv. .false., &
            "cannot input gas_state without gas_data")
    end if

    if (present(env_state)) then
       call env_state_input_netcdf(env_state, ncid)
    end if

    call pmc_nc_close(ncid)

  end subroutine input_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Find all NetCDF (.nc) filenames that match the given prefix.
  subroutine input_filename_list(prefix, filename_list)

    !> Filename prefix to search for.
    character(len=*), intent(in) :: prefix
    !> Filename list.
    character(len=*), intent(inout), allocatable :: filename_list(:)

    integer :: n_file, index, unit, ios
    character(len=len(prefix)+100) :: filename
    logical :: done

    call assert_msg(277193351, pmc_mpi_rank() == 0, &
         "can only call from process 0")

    index = 1
    done = .false.
    unit = get_unit()
    do while (.not. done)
       write(filename, '(a,a,i8.8,a)') trim(prefix), '_', index, '.nc'
       open(unit=unit, file=filename, status='old', iostat=ios)
       if (ios /= 0) then
          done = .true.
       else
          index = index + 1
          close(unit)
       end if
    end do
    call free_unit(unit)

    n_file = index - 1
    call ensure_string_array_size(filename_list, n_file)
    do index = 1,n_file
       write(filename, '(a,a,i8.8,a)') trim(prefix), '_', index, '.nc'
       filename_list(index) = filename
    end do

  end subroutine input_filename_list

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Find the number of repeats and indices for the given prefix.
  subroutine input_n_files(prefix, n_repeat, n_index)

    !> Filename prefix to search for.
    character(len=*), intent(in) :: prefix
    !> Number of repeats found.
    integer, intent(out) :: n_repeat
    !> Number of indices found.
    integer, intent(out) :: n_index

    integer :: repeat, index, unit, ios
    character(len=len(prefix)+100) :: filename
    logical :: done

    call assert_msg(711223711, pmc_mpi_rank() == 0, &
         "can only call from process 0")

    unit = get_unit()

    repeat = 1
    index = 1
    done = .false.
    do while (.not. done)
       call make_filename(filename, prefix, ".nc", index, repeat)
       open(unit=unit, file=filename, status='old', iostat=ios)
       if (ios /= 0) then
          done = .true.
       else
          repeat = repeat + 1
          close(unit)
       end if
    end do
    n_repeat = repeat - 1
    call assert_msg(252703940, n_repeat >= 1, &
         "no files found with prefix: " // trim(prefix))

    repeat = 1
    index = 1
    done = .false.
    do while (.not. done)
       call make_filename(filename, prefix, ".nc", index, repeat)
       open(unit=unit, file=filename, status='old', iostat=ios)
       if (ios /= 0) then
          done = .true.
       else
          index = index + 1
          close(unit)
       end if
    end do
    n_index = index - 1

  end subroutine input_n_files

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Write the current sectional data.
  subroutine output_sectional(prefix, bin_grid, aero_data, aero_binned, &
       gas_data, gas_state, env_state, index, time, del_t, uuid)

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
    !> UUID of the simulation.
    character(len=PMC_UUID_LEN), intent(in) :: uuid

    integer :: ncid
    character(len=len(prefix)+100) :: filename

    write(filename, '(a,a,i8.8,a)') trim(prefix), &
         '_', index, '.nc'
    call pmc_nc_open_write(filename, ncid)
    call pmc_nc_write_info(ncid, uuid, &
         "PartMC version " // trim(PARTMC_VERSION))
    call write_time(ncid, time, del_t, index)

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

  !> Input sectional data.
  subroutine input_sectional(filename, index, time, del_t, uuid, bin_grid, &
       aero_data, aero_binned, gas_data, gas_state, env_state)

    !> Filename to read.
    character(len=*), intent(in) :: filename
    !> Filename index.
    integer, intent(out) :: index
    !> Current time (s).
    real(kind=dp), intent(out) :: time
    !> Current output time-step (s).
    real(kind=dp), intent(out) :: del_t
    !> UUID of the simulation.
    character(len=PMC_UUID_LEN), intent(out) :: uuid
    !> Bin grid.
    type(bin_grid_t), optional, intent(inout) :: bin_grid
    !> Aerosol data.
    type(aero_data_t), optional, intent(inout) :: aero_data
    !> Binned aerosol data.
    type(aero_binned_t), optional, intent(inout) :: aero_binned
    !> Gas data.
    type(gas_data_t), optional, intent(inout) :: gas_data
    !> Gas state.
    type(gas_state_t), optional, intent(inout) :: gas_state
    !> Environment state.
    type(env_state_t), optional, intent(inout) :: env_state

    integer :: ncid

    call assert_msg(559676785, pmc_mpi_rank() == 0, &
         "can only call from process 0")

    call pmc_nc_open_read(filename, ncid)

    call pmc_nc_check(nf90_get_att(ncid, NF90_GLOBAL, "UUID", uuid))

    call pmc_nc_read_real(ncid, time, "time")
    call pmc_nc_read_real(ncid, del_t, "timestep")
    call pmc_nc_read_integer(ncid, index, "timestep_index")

    if (present(bin_grid)) then
       call bin_grid_input_netcdf(bin_grid, ncid, "aero_diam", scale=0.5d0)
    end if
    if (present(aero_data)) then
       call aero_data_input_netcdf(aero_data, ncid)
    end if
    if (present(aero_binned)) then
       call assert_msg(585353528, &
            present(bin_grid) .and. present(aero_data), &
            "cannot input aero_binned without bin_grid and aero_data")
       call aero_binned_input_netcdf(aero_binned, ncid, bin_grid, &
            aero_data)
    end if

    if (present(gas_data)) then
       call gas_data_input_netcdf(gas_data, ncid)
       if (present(gas_state)) then
          call gas_state_input_netcdf(gas_state, ncid, gas_data)
       end if
    else
       call assert_msg(214545112, present(gas_state) .eqv. .false., &
            "cannot input gas_state without gas_data")
    end if

    if (present(env_state)) then
       call env_state_input_netcdf(env_state, ncid)
    end if

    call pmc_nc_close(ncid)

  end subroutine input_sectional

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read the specification for an output type from a spec file and
  !> generate it.
  subroutine spec_file_read_output_type(file, output_type)

    !> Spec file.
    type(spec_file_t), intent(inout) :: file
    !> Kernel type.
    integer, intent(out) :: output_type

    character(len=SPEC_LINE_MAX_VAR_LEN) :: output_type_name

    !> \page input_format_output Input File Format: Output Type
    !!
    !! The output type is specified by the parameter:
    !!   - \b output_type (string): type of disk output --- must be
    !!     one of: \c central to write one file per process, but all
    !!     written by process 0; \c dist for every process to
    !!     write its own state file; or \c single to transfer all data
    !!     to process 0 and write a single unified output file
    !!
    !! See also:
    !!   - \ref spec_file_format --- the input file text format

    call spec_file_read_string(file, 'output_type', output_type_name)
    if (trim(output_type_name) == 'central') then
       output_type = OUTPUT_TYPE_CENTRAL
    elseif (trim(output_type_name) == 'dist') then
       output_type = OUTPUT_TYPE_DIST
    elseif (trim(output_type_name) == 'single') then
       output_type = OUTPUT_TYPE_SINGLE
    else
       call spec_file_die_msg(392313600, file, &
            "Unknown output type: " // trim(output_type_name))
    end if

  end subroutine spec_file_read_output_type

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_output
