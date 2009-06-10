! Copyright (C) 2007, 2008, 2009 Nicole Riemer and Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_spec_read module.

!> Reading formatted text input.
module pmc_spec_read

  use pmc_util
  
  !> Maximum size of a single line.
  integer, parameter :: MAX_LINE_LEN = 10000
  !> Maximum size of a variable.
  integer, parameter :: MAX_VAR_LEN = 300
  !> Maximum number of lines in an array.
  integer, parameter :: MAX_LIST_LINES = 1000

  !> An input file with extra data for printing messages.
  !!
  !! An spec_file_t is just a simple wrapper around a Fortran unit
  !! number together with the filename and current line number. The
  !! line number is updated manually by the various \c spec_read_*()
  !! subroutine. To maintain its validity all file accesses must be
  !! done via the \c spec_read_*() subroutines, and no data should be
  !! accessed directly via \c spec_file%%unit.
  type spec_file_t
     !> Filename.
     character(len=MAX_VAR_LEN) :: name
     !> Attached unit.
     integer :: unit
     !> Current line number.
     integer :: line_num
  end type spec_file_t

  !> A single line of input data, split at whitespace.
  !!
  !! Input lines are assumed to be in the format
  !! <pre>
  !! &lt;name&gt; &lt;whitespace&gt; &lt;data1&gt; &lt;whitespace&gt; &lt;data2&gt; ... # optional comment
  !! </pre>
  !! An spec_line_t structure stores the name and data split at
  !! whitespace.
  type spec_line_t
     !> Variable name.
     character(len=MAX_VAR_LEN) :: name
     !> Array of data as strings.
     character(len=MAX_VAR_LEN), pointer :: data(:)
  end type spec_line_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Exit with an error message containing filename and line number.
  subroutine spec_read_die_msg(code, file, msg)

    !> Failure status code.
    integer, intent(in) :: code
    !> Spec file.
    type(spec_file_t), intent(in) :: file
    !> Error message.
    character(len=*), intent(in) :: msg
    
    character(len=MAX_LINE_LEN) :: error_msg

    write(error_msg, *) "file ", trim(file%name), &
         "line ", file%line_num, ": ", msg
    call die_msg(code, error_msg)

  end subroutine spec_read_die_msg

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Open a spec file for reading.
  subroutine spec_read_open(filename, file)

    !> Name of file to open.
    character(len=*), intent(in) :: filename
    !> Spec file.
    type(spec_file_t), intent(out) :: file

    integer :: ios, unit
    character(len=MAX_LINE_LEN) :: error_msg

    file%name = trim(filename)
    file%unit = get_unit()
    open(unit=file%unit, status='old', file=file%name, iostat=ios)
    if (ios /= 0) then
       write(error_msg, *) 'unable to open file ', &
            trim(file%name), ' for reading: IOSTAT = ', ios
       call die_msg(173932734, error_msg)
    end if
    file%line_num = 0

  end subroutine spec_read_open

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Close a spec file.
  subroutine spec_read_close(file)

    !> Spec file.
    type(spec_file_t), intent(in) :: file

    close(file%unit)
    call free_unit(file%unit)

  end subroutine spec_read_close

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Allocates memory for a spec_line.
  subroutine spec_line_alloc(n_data, spec_line)

    !> Number of data items.
    integer, intent(in) :: n_data
    !> Struct to alloc.
    type(spec_line_t), intent(inout) :: spec_line

    allocate(spec_line%data(n_data))

  end subroutine spec_line_alloc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Frees all storage.
  subroutine spec_line_free(spec_line)

    !> Struct to free.
    type(spec_line_t), intent(inout) :: spec_line

    deallocate(spec_line%data)

  end subroutine spec_line_free

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Copies a spec_line.
  subroutine spec_line_copy(from_line, to_line)

    !> Original line.
    type(spec_line_t), intent(in) :: from_line
    !> Destination, already alloced.
    type(spec_line_t), intent(out) :: to_line

    to_line%name = from_line%name
    if (size(to_line%data) /= size(from_line%data)) then
       deallocate(to_line%data)
    end if
    allocate(to_line%data(size(from_line%data)))
    to_line%data = from_line%data

  end subroutine spec_line_copy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read a single line from a spec file, signaling if we have hit EOF.
  subroutine spec_read_line_raw(file, line, eof)

    !> Spec file.
    type(spec_file_t), intent(inout) :: file
    !> Complete line read.
    character(len=*), intent(out) :: line
    !> True if at EOF.
    logical, intent(out) :: eof

    integer :: ios, n_read
    character(len=MAX_LINE_LEN) :: error_msg

    file%line_num = file%line_num + 1
    eof = .false.
    line = "" ! needed for pgf95 for reading blank lines
    read(unit=file%unit, fmt='(a)', advance='no', end=100, eor=110, &
         iostat=ios) line
    if (ios /= 0) then
       write(error_msg, *) 'error reading: IOSTAT = ', ios
       call spec_read_die_msg(869855853, file, error_msg)
    end if
    ! only reach here if we didn't hit end-of-record (end-of-line) in
    ! the above read, meaning the line was too long
    write(error_msg, *) 'line exceeds length: ', len(line)
    call spec_read_die_msg(468785871, file, error_msg)

100 line = "" ! goto here if end-of-file was encountered immediately
    eof = .true.

110 return ! goto here if end-of-record, meaning everything is ok
    
  end subroutine spec_read_line_raw

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Strip the comments from a line. Comments are everything after
  !> the first # character.
  subroutine spec_read_strip_comment(line)
    
    !> Complete input line.
    character(len=*), intent(inout) :: line
    
    integer hash_index
    
    hash_index = index(line, '#')
    if (hash_index > 0) then
       line = line(1:(hash_index - 1))
    end if
    
  end subroutine spec_read_strip_comment
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Expand all tabs in a line into single spaces (one tab makes one
  !> space).
  subroutine spec_read_tabs_to_spaces(line)

    !> Complete input line.
    character(len=*), intent(inout) :: line

    integer i

    do i = 1,len(line)
       if (ichar(line(i:i)) == 9) then
          line(i:i) = ' '
       end if
    end do

  end subroutine spec_read_tabs_to_spaces

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Strip leading spaces from a string.
  subroutine spec_read_strip_leading_spaces(line)

    !> Complete input line.
    character(len=*), intent(inout) :: line

    integer i

    if (len_trim(line) > 0) then
       i = verify(line, ' ') ! first character not a space
       line = line(i:)
    end if

  end subroutine spec_read_strip_leading_spaces

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read the next line from the spec file that contains useful data
  !> (stripping comments and blank lines).
  subroutine spec_read_next_data_line(file, line, eof)

    !> Spec file.
    type(spec_file_t), intent(inout) :: file
    !> Complete line read.
    character(len=*), intent(out) :: line
    !> True if EOF encountered.
    logical, intent(out) :: eof

    logical :: done

    done = .false.
    do while (.not. done)
       call spec_read_line_raw(file, line, eof)
       if (eof) then
          done = .true.
       else
          call spec_read_strip_comment(line)
          call spec_read_tabs_to_spaces(line)
          call spec_read_strip_leading_spaces(line)
          if (len_trim(line) > 0) then
             done = .true.
          end if
       end if
    end do

  end subroutine spec_read_next_data_line

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read a spec_line from the spec_file.
  subroutine spec_read_line(file, line, eof)

    !> Spec file.
    type(spec_file_t), intent(inout) :: file
    !> Spec line.
    type(spec_line_t), intent(out) :: line
    !> True if EOF encountered.
    logical, intent(out) :: eof

    character(len=MAX_LINE_LEN) :: line_string, rest
    integer i, n_data
    logical done
    character(len=MAX_LINE_LEN) :: error_msg

    call spec_read_next_data_line(file, line_string, eof)
    if (eof) return

    ! strip off the name
    i = index(line_string, ' ') ! first space
    if (i == 0) then
       call spec_read_die_msg(117442928, file, 'line contains no whitespace')
    end if
    if (i == 1) then
       call spec_read_die_msg(650916702, file, 'line starts with whitespace')
    end if
    if (i >= MAX_VAR_LEN) then
       write(error_msg, *) 'line name longer than: ', MAX_VAR_LEN
       call spec_read_die_msg(170403881, file, error_msg)
    end if
    line%name = line_string(1:(i-1))
    line_string = line_string(i:)
    call spec_read_strip_leading_spaces(line_string)

    ! figure out how many data items we have (consecutive non-spaces)
    n_data = 0
    rest = line_string
    done = .false.
    do while (.not. done)
       if (len_trim(rest) == 0) then ! only spaces left
          done = .true.
       else
          ! strip the data element
          n_data = n_data + 1
          i = index(rest, ' ') ! first space
          rest = rest(i:)
          call spec_read_strip_leading_spaces(rest)
       end if
    end do

    ! allocate the data and read out the data items
    allocate(line%data(n_data))
    n_data = 0
    rest = line_string
    done = .false.
    do while (.not. done)
       if (len_trim(rest) == 0) then ! only spaces left
          done = .true.
       else
          ! strip the data element
          n_data = n_data + 1
          i = index(rest, ' ') ! first space
          if (i <= 1) then
             call spec_read_die_msg(332939443, file, 'internal processing error')
          end if
          if (i >= MAX_VAR_LEN) then
             write(error_msg, *) 'data element ', n_data, ' longer than: ', &
                  MAX_VAR_LEN
             call spec_read_die_msg(145508629, file, error_msg)
          end if
          line%data(n_data) = rest(1:(i-1))
          rest = rest(i:)
          call spec_read_strip_leading_spaces(rest)
       end if
    end do
    
  end subroutine spec_read_line

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read a spec_line from the spec_file. This will always succeed or
  !> error out, so should only be called if we know there should be a
  !> valid line coming.
  subroutine spec_read_line_no_eof(file, line)

    !> Spec file.
    type(spec_file_t), intent(inout) :: file
    !> Spec line.
    type(spec_line_t), intent(out) :: line

    logical :: eof

    call spec_read_line(file, line, eof)
    if (eof) then
       call spec_read_die_msg(358475502, file, 'unexpected end of file')
    end if

  end subroutine spec_read_line_no_eof

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read a list of spec_lines from a file, stopping at max_lines
  !> or EOF, whichever comes first.
  subroutine spec_read_line_list(file, max_lines, line_list)

    !> Spec file.
    type(spec_file_t), intent(inout) :: file
    !> Max lines to read (0 = no max).
    integer, intent(in) :: max_lines
    !> List of spec_lines,.
    type(spec_line_t), pointer :: line_list(:)
                                                ! will be allocated

    logical :: eof
    integer :: i, num_lines
    type(spec_line_t) :: temp_line_list(MAX_LIST_LINES)

    ! read file, working out how many lines we have
    num_lines = 0
    eof = .false.
    call spec_read_line(file, temp_line_list(num_lines + 1), eof)
    do while (.not. eof)
       num_lines = num_lines + 1
       if (num_lines > MAX_LIST_LINES) then
          call spec_read_die_msg(450564159, file, 'maximum number of lines exceeded')
       end if
       if (max_lines > 0) then
          if (num_lines >= max_lines) then
             eof = .true.
          end if
       end if
       if (.not. eof) then
          call spec_read_line(file, temp_line_list(num_lines + 1), eof)
       end if
    end do

    ! copy data to actual list
    allocate(line_list(num_lines))
    do i = 1,num_lines
       call spec_line_alloc(0, line_list(i))
       call spec_line_copy(temp_line_list(i), line_list(i))
       call spec_line_free(temp_line_list(i))
    end do

  end subroutine spec_read_line_list

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read an array of spec_lines from a file, stopping at max_lines
  !> or EOF. All lines must have the same number of elements.
  subroutine spec_read_line_array(file, max_lines, line_array)

    !> Spec file.
    type(spec_file_t), intent(inout) :: file
    !> Max lines to read (0 = no max).
    integer, intent(in) :: max_lines
    !> Array of spec_lines,.
    type(spec_line_t), pointer :: line_array(:)
                                                 ! will be allocated

    integer :: i, line_length

    call spec_read_line_list(file, max_lines, line_array)
    if (size(line_array) > 0) then
       line_length = size(line_array(1)%data)
       do i = 2,size(line_array)
          if (size(line_array(i)%data) /= line_length) then
             call spec_read_die_msg(298076484, file, &
                  'lines have unequal numbers of entries for array')
          end if
       end do
    end if
    
  end subroutine spec_read_line_array

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Check that the name of the line data is as given.
  subroutine spec_read_check_line_name(file, line, name)

    !> Spec file.
    type(spec_file_t), intent(in) :: file
    !> Spec line.
    type(spec_line_t), intent(in) :: line
    !> Expected line name.
    character(len=*), intent(in) :: name

    if (line%name /= name) then
       call spec_read_die_msg(462932478, file, &
            'line must begin with: ' // trim(name) &
            // ' not: ' // trim(line%name))
    end if

  end subroutine spec_read_check_line_name
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Checks that the read_name is the same as name.
  subroutine spec_read_check_name(file, name, read_name)
    
    !> Spec file.
    type(spec_file_t), intent(inout) :: file
    !> Name that we should have.
    character(len=*), intent(in) :: name
    !> Name that we do have.
    character(len=*), intent(in) :: read_name
    
    integer name_len, read_name_len

    if (name /= read_name) then
       call spec_read_die_msg(683719069, file, &
            'line must begin with: ' // trim(name) &
            // ' not: ' // trim(read_name))
    end if
    
  end subroutine spec_read_check_name
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Check that the length of the line data is as given.
  subroutine spec_read_check_line_length(file, line, length)

    !> Spec file.
    type(spec_file_t), intent(in) :: file
    !> Spec line.
    type(spec_line_t), intent(in) :: line
    !> Expected data length.
    integer, intent(in) :: length

    character(len=MAX_LINE_LEN) :: error_msg

    if (size(line%data) /= length) then
       write(error_msg, *) 'expected ', length, ' data items on line'
       call spec_read_die_msg(189339129, file, error_msg)
    end if

  end subroutine spec_read_check_line_length
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Check the IOSTAT and error if it is bad.
  subroutine spec_read_check_read_iostat(file, ios, type)

    !> Spec file.
    type(spec_file_t), intent(in) :: file
    !> Iostat result.
    integer, intent(in) :: ios
    !> Type being read during error.
    character(len=*), intent(in) :: type

    character(len=MAX_LINE_LEN) :: error_msg

    if (ios /= 0) then
       write(error_msg, *) 'error reading: IOSTAT = ', ios
       call spec_read_die_msg(704342497, file, error_msg)
    end if

  end subroutine spec_read_check_read_iostat

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Convert a string to an integer.
  integer function spec_read_string_to_integer(file, string)

    !> Spec file.
    type(spec_file_t), intent(in) :: file
    !> String to convert.
    character(len=*), intent(in) :: string
    
    integer :: val
    integer :: ios

    read(string, '(i20)', iostat=ios) val
    call spec_read_check_read_iostat(file, ios, "integer")
    spec_read_string_to_integer = val

  end function spec_read_string_to_integer

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Convert a string to an real.
  real*8 function spec_read_string_to_real(file, string)

    !> Spec file.
    type(spec_file_t), intent(in) :: file
    !> String to convert.
    character(len=*), intent(in) :: string
    
    real*8 :: val
    integer :: ios

    read(string, '(f30.0)', iostat=ios) val
    call spec_read_check_read_iostat(file, ios, "real")
    spec_read_string_to_real = val

  end function spec_read_string_to_real

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Convert a string to an logical.
  logical function spec_read_string_to_logical(file, string)

    !> Spec file.
    type(spec_file_t), intent(in) :: file
    !> String to convert.
    character(len=*), intent(in) :: string
    
    logical :: val
    integer :: ios

    val = .false.
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
       call spec_read_check_read_iostat(file, 1, "logical")
    end if
    spec_read_string_to_logical = val

  end function spec_read_string_to_logical

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read an integer from a spec file that must have the given name.
  subroutine spec_read_integer(file, name, var)

    !> Spec file.
    type(spec_file_t), intent(inout) :: file
    !> Name.
    character(len=*), intent(in) :: name
    !> Variable to store data.
    integer, intent(out) :: var

    type(spec_line_t) :: line

    call spec_read_line_no_eof(file, line)
    call spec_read_check_line_name(file, line, name)
    call spec_read_check_line_length(file, line, 1)
    var = spec_read_string_to_integer(file, line%data(1))
    call spec_line_free(line)

  end subroutine spec_read_integer
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read a real number from a spec file that must have the given
  !> name.
  subroutine spec_read_real(file, name, var)

    !> Spec file.
    type(spec_file_t), intent(inout) :: file
    !> Name.
    character(len=*), intent(in) :: name
    !> Variable to store data.
    real*8, intent(out) :: var

    type(spec_line_t) :: line

    call spec_read_line_no_eof(file, line)
    call spec_read_check_line_name(file, line, name)
    call spec_read_check_line_length(file, line, 1)
    var = spec_read_string_to_real(file, line%data(1))
    call spec_line_free(line)

  end subroutine spec_read_real

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read a logical from a spec file that must have a given name.
  subroutine spec_read_logical(file, name, var)

    !> Spec file.
    type(spec_file_t), intent(inout) :: file
    !> Name.
    character(len=*), intent(in) :: name
    !> Variable to store data.
    logical, intent(out) :: var

    type(spec_line_t) :: line

    call spec_read_line_no_eof(file, line)
    call spec_read_check_line_name(file, line, name)
    call spec_read_check_line_length(file, line, 1)
    var = spec_read_string_to_logical(file, line%data(1))
    call spec_line_free(line)

  end subroutine spec_read_logical

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read a string from a spec file that must have a given name.
  subroutine spec_read_string(file, name, var)

    !> Spec file.
    type(spec_file_t), intent(inout) :: file
    !> Name.
    character(len=*), intent(in) :: name
    !> Variable to store data.
    character(len=*), intent(out) :: var

    type(spec_line_t) :: line

    call spec_read_line_no_eof(file, line)
    call spec_read_check_line_name(file, line, name)
    call spec_read_check_line_length(file, line, 1)
    var = line%data(1)
    call spec_line_free(line)

  end subroutine spec_read_string

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read a complex number from a spec file that must have the given
  !> name.
  subroutine spec_read_complex(file, name, var)

    !> Spec file.
    type(spec_file_t), intent(inout) :: file
    !> Name.
    character(len=*), intent(in) :: name
    !> Variable to store data.
    complex*16, intent(out) :: var

    type(spec_line_t) :: line

    call spec_read_line_no_eof(file, line)
    call spec_read_check_line_name(file, line, name)
    call spec_read_check_line_length(file, line, 2)
    var = cmplx(spec_read_string_to_real(file, line%data(1)), &
         spec_read_string_to_real(file, line%data(2)), 8)
    call spec_line_free(line)

  end subroutine spec_read_complex

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read an array of named lines with real data. All lines must have
  !> the same number of data elements.
  subroutine spec_read_real_named_array(file, max_lines, names, vals)

    !> Spec file.
    type(spec_file_t), intent(inout) :: file
    !> Max lines to read (0 = no max).
    integer, intent(in) :: max_lines
    !> Names of lines.
    character(len=MAX_VAR_LEN), pointer :: names(:)
    !> Data values.
    real*8, pointer :: vals(:,:)

    type(spec_line_t), pointer :: line_array(:)
    integer :: num_lines, line_length, i, j

    call spec_read_line_array(file, max_lines, line_array)
    num_lines = size(line_array)
    if (num_lines > 0) then
       line_length = size(line_array(1)%data)
       allocate(names(num_lines))
       allocate(vals(num_lines, line_length))
       do i = 1,num_lines
          names(i) = line_array(i)%name
          do j = 1,line_length
             vals(i,j) = spec_read_string_to_real(file, line_array(i)%data(j))
          end do
       end do
    else
       allocate(names(0))
       allocate(vals(0,0))
    end if
    do i = 1,num_lines
       call spec_line_free(line_array(i))
    end do
    deallocate(line_array)

  end subroutine spec_read_real_named_array

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read an a time-indexed array of real data.
  subroutine spec_read_timed_real_array(file, line_name, name, times, vals)

    !> Spec file.
    type(spec_file_t), intent(inout) :: file
    !> Name of line for filename.
    character(len=*), intent(in) :: line_name
    !> Variable name.
    character(len=*), intent(in) :: name
    !> Names of lines.
    real*8, pointer :: times(:)
    !> Data values.
    real*8, pointer :: vals(:)
    
    integer :: n_lines, n_times
    character(len=MAX_VAR_LEN) :: read_name
    type(spec_file_t) :: read_file
    character(len=MAX_VAR_LEN), pointer :: read_names(:)
    real*8, pointer :: read_data(:,:)
    character(len=MAX_LINE_LEN) :: error_msg

    call spec_read_string(file, line_name, read_name)
    call spec_read_open(read_name, read_file)
    call spec_read_real_named_array(read_file, 0, read_names, read_data)
    call spec_read_close(read_file)

    n_lines = size(read_names)
    if (n_lines /= 2) then
       call die_msg(694159200, 'must have exactly two data lines in file ' &
            // trim(read_name))
    end if
    n_times = size(read_data,2)
    if (n_times < 1) then
       call die_msg(925956383, 'must have at least one data poin in file ' &
            // trim(read_name))
    end if
    if (trim(read_names(1)) /= "time") then
       call die_msg(692842968, 'first data line in ' // trim(read_name) &
            // ' must start with: time not: ' // trim(read_names(1)))
    end if
    if (trim(read_names(2)) /= name) then
       call die_msg(692842968, 'second data line in ' // trim(read_name) &
            // ' must start with: ' // trim(name) &
            // ' not: ' // trim(read_names(2)))
    end if
    
    allocate(times(n_times))
    allocate(vals(n_times))
    times = read_data(1,:)
    vals = read_data(2,:)
    deallocate(read_names)
    deallocate(read_data)

  end subroutine spec_read_timed_real_array

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_spec_read
