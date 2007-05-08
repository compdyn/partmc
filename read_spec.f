! -*- mode: f90; -*-
! Copyright (C) 2007 Nicole Riemer and Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.
!
! Routines to read data out of .spec files.
!
! .spec files have a very simple format, with one variable (scalar or
! 1D array) per line. Each line must begin with the variable
! name. Variables must be specified in a certain order, and they may
! not be skipped. Blank lines are ok. Comments are everything after a
! # character and are ignored.
!
! FIXME: lots of memory leaks in this code

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module mod_read_spec

  logical, parameter :: DEBUG_OUTPUT = .false. ! .true. for verbose output
  integer, parameter :: MAX_CHAR_LEN = 300 ! max size of a line or variable
  integer, parameter :: MAX_LINES = 500 ! max lines in an array

  type spec_file
     character(len=MAX_CHAR_LEN) :: name ! filename
     integer :: unit                    ! attached unit
     integer :: line_num                ! current line number
  end type spec_file

  type spec_line
     character(len=MAX_CHAR_LEN) :: name ! variable name
     character(len=MAX_CHAR_LEN), pointer :: data(:) ! array of data as strings
  end type spec_line

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine open_spec(spec, filename)

    ! Open a .spec file.

    use mod_util

    type(spec_file), intent(out) :: spec ! spec file
    character(len=*), intent(in) :: filename ! name of file to open

    integer :: ios, unit

    spec%name = trim(filename)
    spec%unit = get_unit()
    open(unit=spec%unit, status='old', file=spec%name, iostat=ios)
    if (ios /= 0) then
       write(0,*) 'ERROR: unable to open file ', spec%name, ': ', ios
       call exit(1)
    end if
    spec%line_num = 0

  end subroutine open_spec

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine close_spec(spec)

    ! Close a .spec file.

    use mod_util

    type(spec_file), intent(in) :: spec ! spec file

    close(spec%unit)
    call free_unit(spec%unit)

  end subroutine close_spec

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine copy_spec_line(from_line, to_line)

    ! Copies a spec_line.

    type(spec_line), intent(in) :: from_line ! original line
    type(spec_line), intent(out) :: to_line ! new line

    to_line%name = from_line%name
    allocate(to_line%data(size(from_line%data)))
    to_line%data = from_line%data

  end subroutine copy_spec_line

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine read_line(spec, line, eof)

    ! Read a single line from a .spec file, signaling if we have hit EOF.

    type(spec_file), intent(inout) :: spec ! spec file
    character(len=*), intent(out) :: line ! complete line read
    logical, intent(out) :: eof           ! true if at EOF

    integer ios

    read(unit=spec%unit, fmt='(a)', end=100, iostat=ios) line
    spec%line_num = spec%line_num + 1
    if (ios /= 0) then
       write(0,*) 'ERROR: reading from ', trim(spec%name), &
            ' at line ', spec%line_num, ': IOSTAT = ', ios
       call exit(1)
    end if
    eof = .false.
    return

100 line = ""
    eof = .true.
    
  end subroutine read_line

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine strip_comment(line)

    ! Strip the comments from a line. Comments are everything after
    ! the first # character.
    
    character(len=*), intent(inout) :: line ! complete input line
    
    integer hash_index
    
    hash_index = index(line, '#')
    if (hash_index > 0) then
       line = line(1:(hash_index - 1))
    end if
    
  end subroutine strip_comment
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine tabs_to_spaces(line)

    ! Expand all tabs in a line into single spaces (one tab makes one
    ! space).

    character(len=*), intent(inout) :: line ! complete input line

    integer i

    do i = 1,len(line)
       if (ichar(line(i:i)) == 9) then
          line(i:i) = ' '
       end if
    end do

  end subroutine tabs_to_spaces

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine strip_leading_spaces(line)

    ! Strip leading spaces from a string.

    character(len=*), intent(inout) :: line ! complete input line

    integer i

    if (len_trim(line) > 0) then
       i = verify(line, ' ') ! first character not a space
       line = line(i:)
    end if

  end subroutine strip_leading_spaces

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine read_next_data_line(spec, line, eof)

    ! Read the next line from the .spec file that contains useful data
    ! (stripping comments and blank lines).

    type(spec_file), intent(inout) :: spec ! spec file
    character(len=*), intent(out) :: line ! complete line read
    logical, intent(out) :: eof           ! true if EOF encountered

    logical :: done

    done = .false.
    do while (.not. done)
       call read_line(spec, line, eof)
       if (eof) then
          done = .true.
       else
          call strip_comment(line)
          call tabs_to_spaces(line)
          call strip_leading_spaces(line)
          if (len_trim(line) > 0) then
             done = .true.
          end if
       end if
    end do

  end subroutine read_next_data_line

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine read_spec_line(spec, line, eof)

    ! Read a spec_line from the spec_file.

    type(spec_file), intent(inout) :: spec ! spec file
    type(spec_line), intent(out) :: line  ! spec line
    logical, intent(out) :: eof           ! true if EOF encountered

    character(len=MAX_CHAR_LEN) :: line_string, rest
    integer i, n_data
    logical done

    call read_next_data_line(spec, line_string, eof)
    if (eof) return

    ! strip off the name
    i = index(line_string, ' ') ! first space
    if (i == 0) then
       write(0,'(a,i3,a,a,a)') 'ERROR: line ', spec%line_num, &
            ' of input file ', trim(spec%name), &
            ' is malformed'
       call exit(1)
    end if
    line%name = line_string(1:(i-1))
    line_string = line_string(i:)
    call strip_leading_spaces(line_string)

    ! figure out how many data items we have (consecutive non-spaces)
    n_data = 0
    rest = line_string
    done = .false.
    do while (.not. done)
       i = verify(rest, ' ') ! first non-space
       if (i == 0) then ! only spaces left
          done = .true.
       else
          ! strip the data element
          n_data = n_data + 1
          rest = rest(i:)
          call strip_leading_spaces(rest)
       end if
    end do

    ! allocate the data and read out the data items
    allocate(line%data(n_data))
    n_data = 0
    rest = line_string
    done = .false.
    do while (.not. done)
       i = verify(rest, ' ') ! first non-space
       if (i == 0) then ! only spaces left
          done = .true.
       else
          ! strip the data element
          n_data = n_data + 1
          line%data(n_data) = rest(1:(i-1))
          rest = rest(i:)
          call strip_leading_spaces(rest)
       end if
    end do
    
  end subroutine read_spec_line

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine read_spec_line_no_eof(spec, line)

    ! Read a spec_line from the spec_file. This will always succeed or
    ! error out, so should only be called if we know there should be a
    ! valid line coming.

    type(spec_file), intent(inout) :: spec ! spec file
    type(spec_line), intent(out) :: line  ! spec line

    logical :: eof

    call read_spec_line(spec, line, eof)
    if (eof) then
       write(0,*) 'ERROR: end of file ', trim(spec%name), &
            ' unexpected at line ', spec%line_num
       call exit(1)
    end if

  end subroutine read_spec_line_no_eof

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine read_spec_line_list(spec, max_lines, line_list)

    ! Read a list of spec_lines from a file, stopping at max_lines
    ! or EOF.

    type(spec_file), intent(inout) :: spec ! spec file
    integer, intent(in) :: max_lines      ! max lines to read (0 = no max)
    type(spec_line), pointer :: line_list(:) ! list of spec_lines

    logical :: eof
    integer :: i, num_lines
    type(spec_line) :: temp_line_list(MAX_LINES)

    ! read file, working out how many lines we have
    num_lines = 0
    eof = .false.
    call read_spec_line(spec, temp_line_list(num_lines + 1), eof)
    do while (.not. eof)
       num_lines = num_lines + 1
       if (max_lines > 0) then
          if (num_lines >= max_lines) then
             eof = .true.
          end if
       end if
       if (.not. eof) then
          call read_spec_line(spec, temp_line_list(num_lines + 1), eof)
       end if
    end do

    ! allocate actual list
    allocate(line_list(num_lines))

    ! copy data to actual list
    do i = 1,num_lines
       call copy_spec_line(temp_line_list(i), line_list(i))
    end do

  end subroutine read_spec_line_list

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine read_spec_line_array(spec, max_lines, line_array)

    ! Read an array of spec_lines from a file, stopping at max_lines
    ! or EOF. All lines must have the same length.

    type(spec_file), intent(inout) :: spec ! spec file
    integer, intent(in) :: max_lines      ! max lines to read (0 = no max)
    type(spec_line), pointer :: line_array(:) ! array of spec_lines

    integer i, line_length

    call read_spec_line_list(spec, max_lines, line_array)
    if (size(line_array) > 0) then
       line_length = size(line_array(i)%data)
       do i = 2,size(line_array)
          if (size(line_array(i)%data) /= line_length) then
             write(0,'(a,a,i3,a,a,a)') 'ERROR: tried to read ', &
                  'array before line ', spec%line_num, &
                  ' of input file ', trim(spec%name), &
                  ' but lines contain varying numbers of elements'
             call exit(1)
          end if
       end do
    end if
    
  end subroutine read_spec_line_array

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine check_line_name(spec, line, name)

    ! Check that the name of the line data is as given.

    type(spec_file), intent(in) :: spec ! spec file
    type(spec_line), intent(in) :: line ! spec line
    character(len=*), intent(in) :: name ! expected line name

    if (line%name /= name) then
       write(0,'(a,i3,a,a,a,a)') 'ERROR: line ', spec%line_num, &
            ' of input file ', trim(spec%name), &
            ' must begin with: ', trim(name)
       call exit(1)
    end if

  end subroutine check_line_name
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine check_name(spec, name, read_name)

    ! Checks that the read_name is the same as name.
    
    type(spec_file), intent(inout) :: spec ! spec file
    character(len=*), intent(in) :: name   ! name that we should have
    character(len=*), intent(in) :: read_name ! name that we do have
    
    integer name_len, read_name_len

    if (name /= read_name) then
       write(0,'(a,i3,a,a,a,a)') 'ERROR: line ', spec%line_num, &
            ' of input file ', trim(spec%name), &
            ' must begin with: ', trim(name)
       call exit(1)
    end if
    
  end subroutine check_name
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine check_line_length(spec, line, length)

    ! Check that the length of the line data is as given.

    type(spec_file), intent(in) :: spec ! spec file
    type(spec_line), intent(in) :: line ! spec line
    integer, intent(in) :: length       ! expected data length

    if (size(line%data) /= length) then
       write(0,*) 'ERROR: expected ', length, ' data items on line ', &
            spec%line_num, ' of file ', trim(spec%name)
       call exit(1)
    end if

  end subroutine check_line_length
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine check_read_iostat(spec, ios, type)

    ! Check the IOSTAT and error if it is bad.

    type(spec_file), intent(in) :: spec ! spec file
    integer, intent(in) :: ios          ! iostat result
    character(len=*), intent(in) :: type ! type being read during error

    if (ios /= 0) then
       write(0,'(a,a,a,a,i3,a,i4)') 'ERROR: reading ', trim(type), &
            ' from file ', trim(spec%name), ' at line ', spec%line_num, &
            ': IOSTAT = ', ios
       call exit(1)
    end if

  end subroutine check_read_iostat

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function string_to_integer(spec, string)

    ! Convert a string to an integer.

    type(spec_file), intent(in) :: spec ! spec file
    character(len=*), intent(in) :: string ! string to convert
    
    integer :: val
    integer :: ios

    read(string, '(i20)', iostat=ios) val
    call check_read_iostat(spec, ios, "integer")
    string_to_integer = val

  end function string_to_integer

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real*8 function string_to_real(spec, string)

    ! Convert a string to an real.

    type(spec_file), intent(in) :: spec ! spec file
    character(len=*), intent(in) :: string ! string to convert
    
    real*8 :: val
    integer :: ios

    read(string, '(f20.0)', iostat=ios) val
    call check_read_iostat(spec, ios, "real")
    string_to_real = val

  end function string_to_real

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  logical function string_to_logical(spec, string)

    ! Convert a string to an logical.

    type(spec_file), intent(in) :: spec ! spec file
    character(len=*), intent(in) :: string ! string to convert
    
    logical :: val
    integer :: ios

    if ((trim(string) == 'yes') &
         .or. (trim(string) == 'y') &
         .or. (trim(string) == 'true') &
         .or. (trim(string) == 't') &
         .or. (trim(string) == '1')) then
       val = .true.
    elseif ((trim(string) == 'no') &
         .or. (trim(string) == 'n') &
         .or. (trim(string) == 'false') &
         .or. (trim(string) == 'f') &
         .or. (trim(string) == '0')) then
       val = .false.
    else
       call check_read_iostat(spec, 1, "logical")
    end if
    string_to_logical = val

  end function string_to_logical

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine read_integer(spec, name, var)

    ! Read an integer from a .spec file that must have the given name.

    type(spec_file), intent(inout) :: spec ! spec file
    character(len=*), intent(in) :: name  ! name
    integer, intent(out) :: var           ! variable to store data

    type(spec_line) :: line

    call read_spec_line_no_eof(spec, line)
    call check_line_name(spec, line, name)
    call check_line_length(spec, line, 1)
    var = string_to_integer(spec, line%data(1))

  end subroutine read_integer
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine read_real(spec, name, var)

    ! Read a real number from a .spec file that must have the given
    ! name.

    type(spec_file), intent(inout) :: spec ! spec file
    character(len=*), intent(in) :: name  ! name
    real*8, intent(out) :: var            ! variable to store data

    type(spec_line) :: line

    call read_spec_line_no_eof(spec, line)
    call check_line_name(spec, line, name)
    call check_line_length(spec, line, 1)
    var = string_to_real(spec, line%data(1))

  end subroutine read_real

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine read_logical(spec, name, var)

    ! Read a logical from a .spec file that must have a given name.

    type(spec_file), intent(inout) :: spec ! spec file
    character(len=*), intent(in) :: name  ! name
    logical, intent(out) :: var           ! variable to store data

    type(spec_line) :: line

    call read_spec_line_no_eof(spec, line)
    call check_line_name(spec, line, name)
    call check_line_length(spec, line, 1)
    var = string_to_logical(spec, line%data(1))

  end subroutine read_logical

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine read_string(spec, name, var)

    ! Read a string from a .spec file that must have a given name.

    type(spec_file), intent(inout) :: spec ! spec file
    character(len=*), intent(in) :: name  ! name
    character(len=*), intent(out) :: var  ! variable to store data

    type(spec_line) :: line

    call read_spec_line_no_eof(spec, line)
    call check_line_name(spec, line, name)
    call check_line_length(spec, line, 1)
    var = line%data(1)

  end subroutine read_string

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine read_real_array(spec, max_lines, names, vals)

    ! Read an array of lines with real data. All lines must have the
    ! same number of data elements.

    type(spec_file), intent(inout) :: spec ! spec file
    integer, intent(in) :: max_lines      ! max lines to read (0 = no max)
    character(len=MAX_CHAR_LEN), pointer :: names(:) ! names of lines
    real*8, pointer :: vals(:,:)          ! data values

    type(spec_line), pointer :: line_array(:)
    integer :: num_lines, line_length, i, j

    call read_spec_line_array(spec, max_lines, line_array)
    num_lines = size(line_array)
    if (num_lines > 0) then
       line_length = size(line_array(1)%data)
       allocate(names(num_lines))
       allocate(vals(num_lines, line_length))
       do i = 1,num_lines
          names(i) = line_array(i)%name
          do j = 1,line_length
             vals(i,j) = string_to_real(spec, line_array(i)%data(j))
          end do
       end do
    else
       allocate(names(0))
       allocate(vals(0,0))
    end if

  end subroutine read_real_array

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine read_init_dist(spec, dist_type, dist_args)

    ! Read initial distribution specification from a .spec file.

    type(spec_file), intent(inout) :: spec ! spec file
    character(len=*), intent(out) :: dist_type ! type of init distribution
    real*8, intent(out) :: dist_args(:) ! distribution parameters

    call read_string(spec, 'dist_type', dist_type)
    if (trim(dist_type) == 'log_normal') then
       call read_real(spec, 'dist_mean_diam', dist_args(1))
       call read_real(spec, 'dist_std_dev', dist_args(2))
    elseif (trim(dist_type) == 'exp') then
       call read_real(spec, 'dist_mean_vol', dist_args(1))
    elseif (trim(dist_type) == 'bidisperse') then
       call read_real(spec, 'dist_small_vol', dist_args(1))
       call read_real(spec, 'dist_big_vol', dist_args(2))
       call read_real(spec, 'dist_big_num', dist_args(3))
    else
       write(0,'(a,a,a,a,a,i3)') 'ERROR: Unknown distribution type ', &
            trim(dist_type), ' in file ', trim(spec%name), &
            ' at line ', spec%line_num
       call exit(1)
    end if

  end subroutine read_init_dist

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine read_gas(spec, gas)

    ! Read gas concentrations from a .spec file.

    use mod_gas

    type(spec_file), intent(inout) :: spec ! spec file
    type(gas_chem), intent(out) :: gas     ! gas data

    integer :: n_species, species, i
    character(len=MAX_CHAR_LEN) :: gas_name
    type(spec_file) :: gas_spec
    character(len=MAX_CHAR_LEN), pointer :: species_name(:)
    real*8, pointer :: species_conc(:,:)
    integer :: species_conc_shape(2)

    ! read the gas data from the specified file
    call read_string(spec, 'gas_init_conc', gas_name)
    call open_spec(gas_spec, gas_name)
    call read_real_array(gas_spec, 0, species_name, species_conc)
    call close_spec(gas_spec)

    ! check the data size
    species_conc_shape = shape(species_conc)
    if (species_conc_shape(2) /= 1) then
       write(0,*) 'ERROR: each line in ', trim(gas_name), &
            ' should only contain one value'
       call exit(1)
    end if

    ! allocate and copy over the data
    n_species = species_conc_shape(1)
    call allocate_gas(gas, n_species)
    do i = 1,n_species
       gas%name(i) = species_name(i)
       gas%conc(i) = species_conc(i,1)
    end do
    call set_gas_mosaic_map(gas)

  end subroutine read_gas

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine read_material(spec, mat)

    ! Read material specification from a .spec file.

    use mod_material

    type(spec_file), intent(inout) :: spec ! spec file
    type(material), intent(out) :: mat  ! material data

    integer :: n_species, species, i
    character(len=MAX_CHAR_LEN) :: aero_name
    type(spec_file) :: aero_spec
    character(len=MAX_CHAR_LEN), pointer :: species_name(:)
    real*8, pointer :: species_data(:,:)
    integer :: species_data_shape(2)

    ! read the aerosol data from the specified file
    call read_string(spec, 'aerosol_data', aero_name)
    call open_spec(aero_spec, aero_name)
    call read_real_array(aero_spec, 0, species_name, species_data)
    call close_spec(aero_spec)

    ! check the data size
    species_data_shape = shape(species_data)
    if (species_data_shape(2) /= 4) then
       write(0,*) 'ERROR: each line in ', trim(aero_name), &
            ' should contain exactly 4 values'
       call exit(1)
    end if

    ! allocate and copy over the data
    n_species = species_data_shape(1)
    call allocate_material(mat, n_species)
    do i = 1,n_species
       mat%name(i) = species_name(i)
       if (species_name(i) == "H2O") then
          mat%i_water = i
       end if
       mat%rho(i) = species_data(i,1)
       mat%nu(i) = int(species_data(i,2))
       mat%eps(i) = species_data(i,3)
       mat%M_w(i) = species_data(i,4)
    end do
    call set_material_mosaic_map(mat)

  end subroutine read_material

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine read_environ(spec, env)

    ! Read environment specification from a .spec file.

    use mod_environ

    type(spec_file), intent(inout) :: spec ! spec file
    type(environ), intent(out) :: env   ! environment data

    integer :: n_temps
    character(len=MAX_CHAR_LEN) :: temp_name
    type(spec_file) :: temp_spec
    character(len=MAX_CHAR_LEN), pointer :: times_name(:), temps_name(:)
    real*8, pointer :: times_data(:,:), temps_data(:,:)
    integer :: times_data_shape(2), temps_data_shape(2)

    ! read the tempurature data from the specified file
    call read_string(spec, 'temp_profile', temp_name)
    call open_spec(temp_spec, temp_name)
    call read_real_array(temp_spec, 1, times_name, times_data)
    call check_name(temp_spec, "time", times_name(1))
    ! FIXME: add a min_lines arg to read_real_array to ensure that
    ! really got one line here
    call read_real_array(temp_spec, 1, temps_name, temps_data)
    call check_name(temp_spec, "temp", temps_name(1))
    call close_spec(temp_spec)

    ! check the data size
    times_data_shape = shape(times_data)
    temps_data_shape = shape(temps_data)
    if (times_data_shape(2) /= temps_data_shape(2)) then
       write(0,*) 'ERROR: file ', trim(temp_name), &
            ' should contain exactly two lines with equal numbers of values'
       call exit(1)
    end if
    n_temps = temps_data_shape(2)

    call allocate_environ_temps(env, n_temps)
    env%temp_times = times_data(1,:)
    env%temps = temps_data(1,:)
    call read_real(spec, 'RH', env%RH)
    call read_real(spec, 'pressure', env%p)
    call read_real(spec, 'rho_a', env%rho_a)
    call read_real(spec, 'latitude', env%latitude)
    call read_real(spec, 'longitude', env%longitude)
    call read_real(spec, 'altitude', env%altitude)
    call read_real(spec, 'start_time', env%start_time)
    call read_integer(spec, 'start_day', env%start_day)

  end subroutine read_environ

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine read_vol_frac(spec, mat, vol_frac)

    ! Read volume fractions from a data file.

    use mod_material

    type(spec_file), intent(inout) :: spec ! spec file
    type(material), intent(in) :: mat   ! material data
    real*8, intent(out) :: vol_frac(:)  ! aerosol species volume fractions

    integer :: n_species, species, i
    character(len=MAX_CHAR_LEN) :: aero_name
    type(spec_file) :: aero_spec
    character(len=MAX_CHAR_LEN), pointer :: species_name(:)
    real*8, pointer :: species_data(:,:)
    integer :: species_data_shape(2)

    ! read the aerosol data from the specified file
    call read_string(spec, 'aerosol_data', aero_name)
    call open_spec(aero_spec, aero_name)
    call read_real_array(aero_spec, 0, species_name, species_data)
    call close_spec(aero_spec)

    ! check the data size
    species_data_shape = shape(species_data)
    if (species_data_shape(2) /= 1) then
       write(0,*) 'ERROR: each line in ', trim(aero_name), &
            ' should contain exactly one data value'
       call exit(1)
    end if

    ! copy over the data
    vol_frac = 0d0
    do i = 1,n_species
       species = material_spec_by_name(mat, species_name(i))
       if (species == 0) then
          write(0,*) 'ERROR: unknown species ', trim(species_name(i)), &
               ' in file ', trim(aero_name)
          call exit(1)
       end if
       vol_frac(species) = species_data(i,1)
    end do

  end subroutine read_vol_frac

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module mod_read_spec
