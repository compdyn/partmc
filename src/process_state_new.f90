! Copyright (C) 2005-2007 Nicole Riemer and Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.
!
! Process the saved state files to obtain summary data.

program process_state_new

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
  use pmc_process_spec
  use pmc_util

  character(len=300) :: filename
  type(process_spec_t), pointer :: process_spec_list(:)
  integer :: i

  call pmc_mpi_init()

  if (iargc() < 2) then
     call print_usage()
     call exit(1)
  end if

  call getarg(1, filename)
  call inout_read_process_spec_list(filename, process_spec_list)

  do i = 2,iargc()
     call getarg(i, filename)
     call process_state_file(filename, process_spec_list)
  end do

  do i = 1,size(process_spec_list)
     call process_spec_free(process_spec_list(i))
  end do
  deallocate(process_spec_list)

  call pmc_mpi_finalize()

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine print_usage()

    write(*,*) 'Usage: process_state <process.spec> <state_filenames...>'
    
  end subroutine print_usage

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine process_state_file(filename, process_spec_list)

    character(len=*), intent(in) :: filename ! state filename
    type(process_spec_t), intent(in) :: process_spec_list(:) ! process specs

    character(len=len(filename)) :: basename
    type(bin_grid_t) :: bin_grid
    type(aero_data_t) :: aero_data
    type(aero_state_t) :: aero_state
    type(gas_data_t) :: gas_data
    type(gas_state_t) :: gas_state
    type(env_t) :: env
    real*8 :: time
    integer :: index, i

    call get_basename(filename, basename)

    call inout_read_state(filename, bin_grid, aero_data, aero_state, &
         gas_data, gas_state, env, time, index)
    write(*,'(a,e20.10)') 'time (s) = ', time

    call process_state_spec_list(basename, process_spec_list, &
         bin_grid, aero_data, aero_state, gas_data, gas_state, &
         env, time, index)

  end subroutine process_state_file

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine process_state_spec_list(basename, process_spec_list, &
       bin_grid, aero_data, aero_state, gas_data, gas_state, env, &
       time, index)

    character(len=*), intent(in) :: basename  ! basename of the output
    type(process_spec_t), intent(in) :: process_spec_list(:) ! process specs
    type(bin_grid_t), intent(in) :: bin_grid ! bin grid
    type(aero_data_t), intent(in) :: aero_data ! aerosol data
    type(aero_state_t), intent(in) :: aero_state ! aerosol state
    type(gas_data_t), intent(in) :: gas_data ! gas data
    type(gas_state_t), intent(in) :: gas_state ! gas state
    type(env_t), intent(in) :: env      ! environment state
    real*8, intent(in) :: time          ! current time (s)
    integer, intent(in) :: index        ! current index

    integer :: i

    do i = 1,size(process_spec_list)
       if (process_spec_list(i)%type == "env") then
          call process_env(basename, process_spec_list(i)%suffix, &
               time, index, env)
       elseif (process_spec_list(i)%type == "gas") then
          call process_gas(basename, process_spec_list(i)%suffix, &
               time, index, gas_data, gas_state)
       elseif (process_spec_list(i)%type == "aero") then
          call process_aero(basename, process_spec_list(i)%suffix, &
               time, index, bin_grid, aero_data, aero_state)
       elseif ((process_spec_list(i)%type == "kappa") &
          .or. (process_spec_list(i)%type == "comp") &
          .or. (process_spec_list(i)%type == "n_orig_part") &
          .or. (process_spec_list(i)%type == "optic_absorb") &
          .or. (process_spec_list(i)%type == "optic_scatter") &
          .or. (process_spec_list(i)%type == "optic_extinct")) then
          call process_hist_new(basename, time, index, bin_grid, &
               env, aero_data, aero_state, process_spec_list(i))
       else
          call die(450985234)
       end if
    end do

  end subroutine process_state_spec_list

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine process_state_open_output(file, basename, suffix)

    ! Open the file basename + suffix + ".dat" for writing.

    type(inout_file_t), intent(out) :: file ! file to open
    character(len=*), intent(in) :: basename ! basename of the file
    character(len=*), intent(in) :: suffix ! suffix of the file
    
    character(len=(len(basename)+len(suffix)+10)) :: filename

    filename = ""
    filename = basename
    filename((len_trim(filename)+1):) = "_"
    filename((len_trim(filename)+1):) = suffix
    filename((len_trim(filename)+1):) = ".dat"
    call inout_open_write(filename, file)
    
  end subroutine process_state_open_output

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine process_env(basename, suffix, time, index, env)

    character(len=*), intent(in) :: basename ! basename of the file
    character(len=*), intent(in) :: suffix ! suffix of the file
    real*8, intent(in) :: time          ! current time (s)
    integer, intent(in) :: index        ! current index
    type(env_t), intent(in) :: env      ! environment state

    type(inout_file_t) :: file

    call process_state_open_output(file, basename, suffix)
    call inout_write_string(file, 'name', 'env')
    call inout_write_real(file, 'time', time)
    call inout_write_integer(file, 'index', index)
    call inout_write_integer(file, 'n_dim', 1)
    call inout_write_integer(file, 'unit_dim', 1)

    call inout_write_integer(file, 'dim', 1)
    call inout_write_string(file, 'name', 'env_quantity')
    call inout_write_string(file, 'grid_type', 'center')
    call inout_write_string(file, 'data_type', 'string')
    call inout_write_string(file, 'unit', 'none')
    call inout_write_integer(file, 'length', 4)

    call inout_write_indexed_string(file, 'grid', 1, 'temp')
    call inout_write_indexed_real(file, 'grid_width', 1, 1d0)
    call inout_write_indexed_string(file, 'unit', 1, 'K')

    call inout_write_indexed_string(file, 'grid', 2, 'rel_humid')
    call inout_write_indexed_real(file, 'grid_width', 2, 1d0)
    call inout_write_indexed_string(file, 'unit', 2, '1')

    call inout_write_indexed_string(file, 'grid', 3, 'pressure')
    call inout_write_indexed_real(file, 'grid_width', 3, 1d0)
    call inout_write_indexed_string(file, 'unit', 3, 'Pa')

    call inout_write_indexed_string(file, 'grid', 4, 'height')
    call inout_write_indexed_real(file, 'grid_width', 4, 1d0)
    call inout_write_indexed_string(file, 'unit', 4, 'm')

    call inout_write_comment(file, 'data values follow, row major order')
    call inout_write_unnamed_real(file, env%temp)
    call inout_write_unnamed_real(file, env%rel_humid)
    call inout_write_unnamed_real(file, env%pressure)
    call inout_write_unnamed_real(file, env%height)

  end subroutine process_env

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine process_gas(basename, suffix, time, index, &
       gas_data, gas_state)

    character(len=*), intent(in) :: basename ! basename of the file
    character(len=*), intent(in) :: suffix ! suffix of the file
    real*8, intent(in) :: time          ! current time (s)
    integer, intent(in) :: index        ! current index
    type(gas_data_t), intent(in) :: gas_data ! gas data
    type(gas_state_t), intent(in) :: gas_state ! gas state

    type(inout_file_t) :: file
    integer :: i_spec

    call process_state_open_output(file, basename, suffix)
    call inout_write_string(file, 'name', 'gas')
    call inout_write_real(file, 'time', time)
    call inout_write_integer(file, 'index', index)
    call inout_write_integer(file, 'n_dim', 1)
    call inout_write_integer(file, 'unit_dim', 1)

    call inout_write_integer(file, 'dim', 1)
    call inout_write_string(file, 'name', 'species')
    call inout_write_string(file, 'grid_type', 'center')
    call inout_write_string(file, 'data_type', 'string')
    call inout_write_string(file, 'unit', 'none')
    call inout_write_integer(file, 'length', gas_data%n_spec)

    do i_spec = 1,gas_data%n_spec
       call inout_write_indexed_string(file, 'grid', i_spec, &
            gas_data%name(i_spec))
       call inout_write_indexed_real(file, 'grid_width', 1, 1d0)
       call inout_write_indexed_string(file, 'unit', 1, 'ppb')
    end do

    call inout_write_comment(file, 'data values follow, row major order')
    do i_spec = 1,gas_data%n_spec
       call inout_write_unnamed_real(file, gas_state%conc(i_spec))
    end do

  end subroutine process_gas

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine process_aero(basename, suffix, time, index, bin_grid, &
       aero_data, aero_state)

    use pmc_util
    use pmc_aero_binned

    character(len=*), intent(in) :: basename ! basename of the output filename
    character(len=*), intent(in) :: suffix ! suffix for the output filename
    real*8, intent(in) :: time          ! current time (s)
    integer, intent(in) :: index        ! current index
    type(bin_grid_t), intent(in) :: bin_grid ! bin_grid structure
    type(aero_data_t), intent(in) :: aero_data ! aero_data structure
    type(aero_state_t), intent(in) :: aero_state ! aero_state structure

    type(aero_binned_t) :: aero_binned
    type(inout_file_t) :: file
    integer :: i_bin, i_spec

    call aero_binned_alloc(aero_binned, bin_grid%n_bin, aero_data%n_spec)
    call aero_state_to_binned(bin_grid, aero_data, aero_state, aero_binned)

    call process_state_open_output(file, basename, suffix)
    call inout_write_string(file, 'name', 'aero')
    call inout_write_real(file, 'time', time)
    call inout_write_integer(file, 'index', index)
    call inout_write_integer(file, 'n_dim', 3)
    call inout_write_integer(file, 'unit_dim', 3)

    ! dim 1: bin
    call inout_write_integer(file, 'dim', 1)
    call inout_write_string(file, 'name', 'bin')
    call inout_write_string(file, 'grid_type', 'edge')
    call inout_write_string(file, 'data_type', 'real')
    call inout_write_string(file, 'unit', 'm')
    call inout_write_integer(file, 'length', bin_grid%n_bin)
    do i_bin = 1,(bin_grid%n_bin + 1)
       call inout_write_indexed_real(file, 'grid', i_bin, &
            bin_edge(bin_grid, i_bin))
       call inout_write_indexed_real(file, 'grid_width', i_bin, bin_grid%dlnr)
    end do

    ! dim 2: species
    call inout_write_integer(file, 'dim', 2)
    call inout_write_string(file, 'name', 'species')
    call inout_write_string(file, 'grid_type', 'center')
    call inout_write_string(file, 'data_type', 'string')
    call inout_write_string(file, 'unit', 'none')
    call inout_write_integer(file, 'length', aero_data%n_spec)
    do i_spec = 1,aero_data%n_spec
       call inout_write_indexed_string(file, 'grid', i_spec, &
            aero_data%name(i_spec))
       call inout_write_indexed_real(file, 'grid_width', i_spec, 1d0)
    end do

    ! dim 3: unit
    call inout_write_integer(file, 'dim', 3)
    call inout_write_string(file, 'name', 'unit')
    call inout_write_string(file, 'grid_type', 'center')
    call inout_write_string(file, 'data_type', 'string')
    call inout_write_string(file, 'unit', 'none')
    call inout_write_integer(file, 'length', 4)
    call inout_write_indexed_string(file, 'grid', 1, 'num_den')
    call inout_write_indexed_real(file, 'grid_width', 1, 1d0)
    call inout_write_indexed_string(file, 'unit', 1, '#/m^3')
    call inout_write_indexed_string(file, 'grid', 2, 'vol_den')
    call inout_write_indexed_real(file, 'grid_width', 2, 1d0)
    call inout_write_indexed_string(file, 'unit', 2, 'm^3/m^3')
    call inout_write_indexed_string(file, 'grid', 3, 'mass_den')
    call inout_write_indexed_real(file, 'grid_width', 3, 1d0)
    call inout_write_indexed_string(file, 'unit', 3, 'kg/m^3')
    call inout_write_indexed_string(file, 'grid', 4, 'mole_den')
    call inout_write_indexed_real(file, 'grid_width', 4, 1d0)
    call inout_write_indexed_string(file, 'unit', 4, 'moles/m^3')

    call inout_write_comment(file, 'data values follow, row major order')
    do i_bin = 1,bin_grid%n_bin
       do i_spec = 1,aero_data%n_spec
          call inout_write_unnamed_real(file, &
               aero_binned%num_den(i_bin) / dble(aero_data%n_spec))
          call inout_write_unnamed_real(file, &
               aero_binned%vol_den(i_bin, i_spec))
          call inout_write_unnamed_real(file, &
               aero_binned%vol_den(i_bin, i_spec) &
               * aero_data%density(i_spec))
          call inout_write_unnamed_real(file, &
               aero_binned%vol_den(i_bin, i_spec) &
               * aero_data%density(i_spec) &
               / aero_data%molec_weight(i_spec))
       end do
    end do

    call aero_binned_free(aero_binned)

  end subroutine process_aero

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine process_hist_new(basename, time, index, bin_grid, &
       env, aero_data, aero_state, process_spec)

    character(len=*), intent(in) :: basename  ! basename of the output
    real*8, intent(in) :: time          ! current time (s)
    integer, intent(in) :: index        ! current index
    type(bin_grid_t), intent(in) :: bin_grid ! bin grid
    type(env_t), intent(in) :: env      ! environment state
    type(aero_data_t), intent(in) :: aero_data ! aerosol data
    type(aero_state_t), intent(in) :: aero_state ! aerosol state
    type(process_spec_t), intent(in) :: process_spec ! process spec

  end subroutine process_hist_new

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
end program process_state_new
