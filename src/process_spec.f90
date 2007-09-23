! Copyright (C) 2007 Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.
!
! Datatype for process_state.

module pmc_process_spec

  use pmc_inout
  use pmc_bin_grid
  use pmc_aero_data
  use pmc_aero_state
  use pmc_gas_data
  use pmc_gas_state
  use pmc_env
  use pmc_util

  integer, parameter :: PROCESS_SPEC_TYPE_LEN = 40
  integer, parameter :: PROCESS_SPEC_FILE_LEN = 200

  type process_spec_t
     character(len=PROCESS_SPEC_TYPE_LEN) :: type ! processing type
     character(len=PROCESS_SPEC_FILE_LEN) :: suffix ! output file suffix
     integer :: n_step                  ! number of steps for histograms
     real*8 :: min_real                 ! minimum histogram value (real)
     real*8 :: max_real                 ! maximum histogram value (real)
     integer :: min_int                 ! minimum histogram value (int)
     integer :: max_int                 ! maximum histogram value (int)
     logical :: log_scale               ! use a log-scale?
     character(len=AERO_NAME_LEN), pointer :: a_species(:) ! comp A species
     character(len=AERO_NAME_LEN), pointer :: b_species(:) ! comp B species
  end type process_spec_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine process_spec_alloc(process_spec)

    ! Allocate and initialize storage.

    type(process_spec_t), intent(out) :: process_spec ! structure to alloc

    process_spec%type = "none"
    process_spec%suffix = "none"
    process_spec%n_step = 0
    process_spec%min_real = 0d0
    process_spec%max_real = 0d0
    process_spec%min_int = 0
    process_spec%max_int = 0
    process_spec%log_scale = .false.
    allocate(process_spec%a_species(0))
    allocate(process_spec%b_species(0))

  end subroutine process_spec_alloc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine process_spec_free(process_spec)

    ! Free storage.

    type(process_spec_t), intent(inout) :: process_spec ! structure to free

    deallocate(process_spec%a_species)
    deallocate(process_spec%b_species)

  end subroutine process_spec_free

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine process_spec_list_free(process_spec_list)

    ! Free storage.

    type(process_spec_t), pointer :: process_spec_list(:) ! structure to free

    integer :: i

    do i = 1,size(process_spec_list)
       call process_spec_free(process_spec_list(i))
    end do
    deallocate(process_spec_list)

  end subroutine process_spec_list_free

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine process_spec_copy(process_spec_from, process_spec_to)

    ! Copy all data.

    type(process_spec_t), intent(in) :: process_spec_from ! source of data
    type(process_spec_t), intent(inout) :: process_spec_to ! dest of data

    process_spec_to%type = process_spec_from%type
    process_spec_to%suffix = process_spec_from%suffix
    process_spec_to%n_step = process_spec_from%n_step
    process_spec_to%min_real = process_spec_from%min_real
    process_spec_to%max_real = process_spec_from%max_real
    process_spec_to%min_int = process_spec_from%min_int
    process_spec_to%max_int = process_spec_from%max_int
    process_spec_to%log_scale = process_spec_from%log_scale
    allocate(process_spec_to%a_species(size(process_spec_from%a_species)))
    allocate(process_spec_to%b_species(size(process_spec_from%b_species)))
    process_spec_to%a_species = process_spec_from%a_species
    process_spec_to%b_species = process_spec_from%b_species

  end subroutine process_spec_copy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine spec_read_process_spec_list_filename(file, name, &
       process_spec_list)

    ! Read process_spec_list from filename on line in file.

    type(inout_file_t), intent(inout) :: file ! inout file
    character(len=*), intent(in) :: name ! name of data line for filename
    type(process_spec_t), pointer :: process_spec_list(:) ! data to read

    character(len=MAX_CHAR_LEN) :: read_name
    type(inout_file_t) :: read_file

    call inout_read_string(file, name, read_name)
    call inout_open_read(read_name, read_file)
    call inout_read_process_spec_list(read_file, process_spec_list)
    call inout_close(read_file)

  end subroutine spec_read_process_spec_list_filename
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine inout_read_process_spec_list(file, process_spec_list)

    ! Read a list of process_specs from the given file.

    type(inout_file_t), intent(inout) :: file ! file to read
    type(process_spec_t), pointer :: process_spec_list(:) ! data to read

    type(process_spec_t), pointer :: new_process_spec_list(:)
    type(process_spec_t) :: process_spec
    integer :: i, n_process_spec
    logical :: eof
    
    n_process_spec = 0
    allocate(process_spec_list(n_process_spec))
    call inout_read_process_spec(file, process_spec, eof)
    do while (.not. eof)
       n_process_spec = n_process_spec + 1
       allocate(new_process_spec_list(n_process_spec))
       call process_spec_copy(process_spec, &
            new_process_spec_list(n_process_spec))
       call process_spec_free(process_spec)
       do i = 1,(n_process_spec - 1)
          call process_spec_copy(process_spec_list(i), new_process_spec_list(i))
          call process_spec_free(process_spec_list(i))
       end do
       deallocate(process_spec_list)
       process_spec_list => new_process_spec_list
       nullify(new_process_spec_list)
       call inout_read_process_spec(file, process_spec, eof)
    end do
    
  end subroutine inout_read_process_spec_list

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine inout_read_process_spec(file, process_spec, eof)

    ! Read from an inout_file.

    type(inout_file_t), intent(inout) :: file ! inout file
    type(process_spec_t), intent(out) :: process_spec ! data to read
    logical :: eof                      ! if eof instead of reading data

    type(inout_line_t) :: line

    call inout_read_line(file, line, eof)
    if (.not. eof) then
       call inout_check_line_name(file, line, "process")
       call inout_check_line_length(file, line, 1)
       call process_spec_alloc(process_spec)
       process_spec%type = line%data(1)(1:PROCESS_SPEC_TYPE_LEN)
       call inout_read_string(file, "suffix", process_spec%suffix)
       if (process_spec%type == "kappa") then
          call inout_read_process_spec_kappa(file, process_spec)
       elseif (process_spec%type == "comp") then
          call inout_read_process_spec_comp(file, process_spec)
       elseif (process_spec%type == "n_orig_part") then
          call inout_read_process_spec_n_orig_part(file, process_spec)
       elseif (process_spec%type == "optic_absorb") then
          call inout_read_process_spec_optic(file, process_spec)
       elseif (process_spec%type == "optic_scatter") then
          call inout_read_process_spec_optic(file, process_spec)
       elseif (process_spec%type == "optic_extinct") then
          call inout_read_process_spec_optic(file, process_spec)
       elseif ((process_spec%type /= "env") &
            .and. (process_spec%type /= "gas") &
            .and. (process_spec%type /= "aero")) then
          write(0,*) 'ERROR: unknown process type on line ', &
               file%line_num, ' of file ', trim(file%name)
          call exit(1)
       end if
       call inout_line_free(line)
    end if

  end subroutine inout_read_process_spec

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine inout_read_process_spec_kappa(file, process_spec)

    ! Read kappa spec from an inout_file.

    type(inout_file_t), intent(inout) :: file ! inout file
    type(process_spec_t), intent(out) :: process_spec ! data to read

    call inout_read_integer(file, "n_step", process_spec%n_step)
    call inout_read_real(file, "min", process_spec%min_real)
    call inout_read_real(file, "max", process_spec%max_real)

  end subroutine inout_read_process_spec_kappa

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine inout_read_process_spec_comp(file, process_spec)

    ! Read comp spec from an inout_file.

    type(inout_file_t), intent(inout) :: file ! inout file
    type(process_spec_t), intent(out) :: process_spec ! data to read

    type(inout_line_t) :: line
    integer :: i

    call inout_read_integer(file, "n_step", process_spec%n_step)
    call inout_read_real(file, "min", process_spec%min_real)
    call inout_read_real(file, "max", process_spec%max_real)

    call inout_read_line_no_eof(file, line)
    call inout_check_line_name(file, line, "a_species")
    allocate(process_spec%a_species(size(line%data)))
    do i = 1,size(line%data)
       process_spec%a_species(i) = line%data(i)(1:AERO_NAME_LEN)
    end do
    call inout_line_free(line)

    call inout_read_line_no_eof(file, line)
    call inout_check_line_name(file, line, "b_species")
    allocate(process_spec%b_species(size(line%data)))
    do i = 1,size(line%data)
       process_spec%b_species(i) = line%data(i)(1:AERO_NAME_LEN)
    end do
    call inout_line_free(line)

  end subroutine inout_read_process_spec_comp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine inout_read_process_spec_n_orig_part(file, process_spec)

    ! Read n_orig_part spec from an inout_file.

    type(inout_file_t), intent(inout) :: file ! inout file
    type(process_spec_t), intent(out) :: process_spec ! data to read

    call inout_read_integer(file, "min", process_spec%min_int)
    call inout_read_integer(file, "max", process_spec%max_int)

  end subroutine inout_read_process_spec_n_orig_part

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine inout_read_process_spec_optic(file, process_spec)

    ! Read optic spec from an inout_file.

    type(inout_file_t), intent(inout) :: file ! inout file
    type(process_spec_t), intent(out) :: process_spec ! data to read

    call inout_read_integer(file, "n_step", process_spec%n_step)
    call inout_read_real(file, "min", process_spec%min_real)
    call inout_read_real(file, "max", process_spec%max_real)

  end subroutine inout_read_process_spec_optic

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
    call inout_write_real(file, 'time', time)
    call inout_write_integer(file, 'index', index)
    call inout_write_string(file, 'name', 'env')
    call inout_write_integer(file, 'n_dim', 1)

    call inout_write_integer(file, 'dim', 1)
    call inout_write_string(file, 'name', 'env_quantity')
    call inout_write_string(file, 'grid_type', 'center')
    call inout_write_string(file, 'data_type', 'string')
    call inout_write_string(file, 'unit', 'none')
    call inout_write_string(file, 'have_grid_units', 'yes')
    call inout_write_integer(file, 'length', 4)

    call inout_write_indexed_string(file, 'grid_center', 1, 'temp')
    call inout_write_indexed_string(file, 'grid_center', 2, 'rel_humid')
    call inout_write_indexed_string(file, 'grid_center', 3, 'pressure')
    call inout_write_indexed_string(file, 'grid_center', 4, 'height')
    call inout_write_indexed_real(file, 'grid_width', 1, 1d0)
    call inout_write_indexed_real(file, 'grid_width', 2, 1d0)
    call inout_write_indexed_real(file, 'grid_width', 3, 1d0)
    call inout_write_indexed_real(file, 'grid_width', 4, 1d0)
    call inout_write_indexed_string(file, 'grid_unit', 1, 'K')
    call inout_write_indexed_string(file, 'grid_unit', 2, '1')
    call inout_write_indexed_string(file, 'grid_unit', 3, 'Pa')
    call inout_write_indexed_string(file, 'grid_unit', 4, 'm')

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
    call inout_write_real(file, 'time', time)
    call inout_write_integer(file, 'index', index)
    call inout_write_string(file, 'name', 'gas')
    call inout_write_integer(file, 'n_dim', 1)

    call inout_write_integer(file, 'dim', 1)
    call inout_write_string(file, 'name', 'species')
    call inout_write_string(file, 'grid_type', 'center')
    call inout_write_string(file, 'data_type', 'string')
    call inout_write_string(file, 'unit', 'none')
    call inout_write_string(file, 'have_grid_units', 'yes')
    call inout_write_integer(file, 'length', gas_data%n_spec)

    do i_spec = 1,gas_data%n_spec
       call inout_write_indexed_string(file, 'grid_center', i_spec, &
            gas_data%name(i_spec))
    end do
    do i_spec = 1,gas_data%n_spec
       call inout_write_indexed_real(file, 'grid_width', i_spec, 1d0)
    end do
    do i_spec = 1,gas_data%n_spec
       call inout_write_indexed_string(file, 'grid_unit', i_spec, 'ppb')
    end do

    call inout_write_comment(file, 'data values follow, row major order')
    do i_spec = 1,gas_data%n_spec
       call inout_write_unnamed_real(file, gas_state%conc(i_spec))
    end do

  end subroutine process_gas

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine output_aero(basename, suffix, time, index, bin_grid, &
       aero_data, aero_binned)

    use pmc_util
    use pmc_aero_binned

    character(len=*), intent(in) :: basename ! basename of the output filename
    character(len=*), intent(in) :: suffix ! suffix for the output filename
    real*8, intent(in) :: time          ! current time (s)
    integer, intent(in) :: index        ! current index
    type(bin_grid_t), intent(in) :: bin_grid ! bin_grid structure
    type(aero_data_t), intent(in) :: aero_data ! aero_data structure
    type(aero_binned_t), intent(in) :: aero_binned ! aero_binned structure

    type(inout_file_t) :: file
    integer :: i_bin, i_spec

    call process_state_open_output(file, basename, suffix)
    call inout_write_real(file, 'time', time)
    call inout_write_integer(file, 'index', index)
    call inout_write_string(file, 'name', 'aero')
    call inout_write_integer(file, 'n_dim', 3)

    ! dim 1: bin
    call inout_write_integer(file, 'dim', 1)
    call inout_write_string(file, 'name', 'radius')
    call inout_write_string(file, 'grid_type', 'center_edge')
    call inout_write_string(file, 'data_type', 'real')
    call inout_write_string(file, 'unit', 'm')
    call inout_write_string(file, 'have_grid_units', 'no')
    call inout_write_integer(file, 'length', bin_grid%n_bin)
    do i_bin = 1,bin_grid%n_bin
       call inout_write_indexed_real(file, 'grid_center', i_bin, &
            vol2rad(bin_grid%v(i_bin)))
    end do
    do i_bin = 1,(bin_grid%n_bin + 1)
       call inout_write_indexed_real(file, 'grid_edge', i_bin, &
            vol2rad(bin_edge(bin_grid, i_bin)))
    end do
    do i_bin = 1,bin_grid%n_bin
       call inout_write_indexed_real(file, 'grid_width', i_bin, bin_grid%dlnr)
    end do

    ! dim 2: species
    call inout_write_integer(file, 'dim', 2)
    call inout_write_string(file, 'name', 'species')
    call inout_write_string(file, 'grid_type', 'center')
    call inout_write_string(file, 'data_type', 'string')
    call inout_write_string(file, 'unit', 'none')
    call inout_write_string(file, 'have_grid_units', 'no')
    call inout_write_integer(file, 'length', aero_data%n_spec)
    do i_spec = 1,aero_data%n_spec
       call inout_write_indexed_string(file, 'grid_center', i_spec, &
            aero_data%name(i_spec))
    end do
    do i_spec = 1,aero_data%n_spec
       call inout_write_indexed_real(file, 'grid_width', i_spec, 1d0)
    end do

    ! dim 3: unit
    call inout_write_integer(file, 'dim', 3)
    call inout_write_string(file, 'name', 'unit')
    call inout_write_string(file, 'grid_type', 'center')
    call inout_write_string(file, 'data_type', 'string')
    call inout_write_string(file, 'unit', 'none')
    call inout_write_string(file, 'have_grid_units', 'yes')
    call inout_write_integer(file, 'length', 4)
    call inout_write_indexed_string(file, 'grid_center', 1, 'num_den')
    call inout_write_indexed_string(file, 'grid_center', 2, 'vol_den')
    call inout_write_indexed_string(file, 'grid_center', 3, 'mass_den')
    call inout_write_indexed_string(file, 'grid_center', 4, 'mole_den')
    call inout_write_indexed_real(file, 'grid_width', 1, 1d0)
    call inout_write_indexed_real(file, 'grid_width', 2, 1d0)
    call inout_write_indexed_real(file, 'grid_width', 3, 1d0)
    call inout_write_indexed_real(file, 'grid_width', 4, 1d0)
    call inout_write_indexed_string(file, 'grid_unit', 1, '#/m^3')
    call inout_write_indexed_string(file, 'grid_unit', 2, 'm^3/m^3')
    call inout_write_indexed_string(file, 'grid_unit', 3, 'kg/m^3')
    call inout_write_indexed_string(file, 'grid_unit', 4, 'moles/m^3')

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

  end subroutine output_aero

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

    call aero_binned_alloc(aero_binned, bin_grid%n_bin, aero_data%n_spec)
    call aero_state_to_binned(bin_grid, aero_data, aero_state, aero_binned)

    call output_aero(basename, suffix, time, index, bin_grid, &
         aero_data, aero_binned)

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
  
end module pmc_process_spec
