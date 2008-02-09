! Copyright (C) 2007, 2008 Nicole Riemer and Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_inout module.

!> Reading formatted text input.
module pmc_inout

  use pmc_util
  
  !> Maximum size of a single line.
  integer, parameter :: MAX_LINE_LEN = 10000
  !> Maximum size of a variable.
  integer, parameter :: MAX_VAR_LEN = 300
  !> Maximum number of lines in an array.
  integer, parameter :: MAX_LIST_LINES = 1000

  !> An input file with extra data for printing messages.
  !!
  !! An inout_file_t is just a simple wrapper around a Fortran unit
  !! number together with the filename and current line number. The
  !! line number is updated manually by the various \c inout_*()
  !! subroutine. To maintain its validity all file accesses must be
  !! done via the \c inout_*() subroutines, and no data should be
  !! accessed directly via \c inout_file%%unit.
  type inout_file_t
     !> Filename.
     character(len=MAX_VAR_LEN) :: name
     !> Attached unit.
     integer :: unit
     !> Current line number.
     integer :: line_num
  end type inout_file_t

  !> A single line of input data, split at whitespace.
  !!
  !! Input lines are assumed to be in the format
  !! <pre>
  !! <name> <whitespace> <data1> <whitespace> <data2> ... # optional comment
  !! </pre>
  !! An inout_line_t structure stores the name and data split at
  !! whitespace.
  type inout_line_t
     !> Variable name.
     character(len=MAX_VAR_LEN) :: name
     !> Array of data as strings.
     character(len=MAX_VAR_LEN), pointer :: data(:)
  end type inout_line_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Open an inout file for reading.
  subroutine inout_open_read(filename, file)

    !> Name of file to open.
    character(len=*), intent(in) :: filename
    !> Inout file.
    type(inout_file_t), intent(out) :: file

    integer :: ios, unit

    file%name = trim(filename)
    file%unit = get_unit()
    open(unit=file%unit, status='old', file=file%name, iostat=ios)
    if (ios /= 0) then
       write(0,'(a,a,a,i4)') 'ERROR: unable to open file ', &
            trim(file%name), ' for reading: ', ios
       call exit(1)
    end if
    file%line_num = 0

  end subroutine inout_open_read

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Open an inout file for writing.
  subroutine inout_open_write(filename, file)

    !> Name of file to open.
    character(len=*), intent(in) :: filename
    !> Inout file.
    type(inout_file_t), intent(out) :: file

    integer :: ios, unit

    file%name = trim(filename)
    file%unit = get_unit()
    open(unit=file%unit, file=file%name, iostat=ios)
    if (ios /= 0) then
       write(0,'(a,a,a,i4)') 'ERROR: unable to open file ', &
            trim(file%name), ' for writing: ', ios
       call exit(1)
    end if
    file%line_num = 0

  end subroutine inout_open_write

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Close an inout file.
  subroutine inout_close(file)

    !> Inout file.
    type(inout_file_t), intent(in) :: file

    close(file%unit)
    call free_unit(file%unit)

  end subroutine inout_close

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Allocates memory for an inout_line.
  subroutine inout_line_alloc(n_data, inout_line)

    !> Number of data items.
    integer, intent(in) :: n_data
    !> Struct to alloc.
    type(inout_line_t), intent(inout) :: inout_line

    allocate(inout_line%data(n_data))

  end subroutine inout_line_alloc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Frees all storage.
  subroutine inout_line_free(inout_line)

    !> Struct to free.
    type(inout_line_t), intent(inout) :: inout_line

    deallocate(inout_line%data)

  end subroutine inout_line_free

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Copies a inout_line.
  subroutine inout_line_copy(from_line, to_line)

    !> Original line.
    type(inout_line_t), intent(in) :: from_line
    !> Destination, already alloced.
    type(inout_line_t), intent(out) :: to_line

    to_line%name = from_line%name
    if (size(to_line%data) /= size(from_line%data)) then
       deallocate(to_line%data)
    end if
    allocate(to_line%data(size(from_line%data)))
    to_line%data = from_line%data

  end subroutine inout_line_copy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read a single line from a inout file, signaling if we have hit EOF.
  subroutine inout_read_line_raw(file, line, eof)

    !> Inout file.
    type(inout_file_t), intent(inout) :: file
    !> Complete line read.
    character(len=*), intent(out) :: line
    !> True if at EOF.
    logical, intent(out) :: eof

    integer :: ios, n_read

    file%line_num = file%line_num + 1
    eof = .false.
    read(unit=file%unit, fmt='(a)', advance='no', end=100, eor=110, &
         iostat=ios) line
    if (ios /= 0) then
       write(0,*) 'ERROR: reading from ', trim(file%name), &
            ' at line ', file%line_num, ': IOSTAT = ', ios
       call exit(1)
    end if
    ! only reach here if we didn't hit end-of-record (end-of-line) in
    ! the above read, meaning the line was too long
    write(0,*) 'ERROR: reading from ', trim(file%name), &
         ' at line ', file%line_num, ': line exceeds ', len(line), &
         ' characters'
    call exit(1)

100 line = "" ! goto here if end-of-file was encountered immediately
    eof = .true.

110 return ! goto here if end-of-record, meaning everything is ok
    
  end subroutine inout_read_line_raw

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Strip the comments from a line. Comments are everything after
  !> the first # character.
  subroutine inout_strip_comment(line)
    
    !> Complete input line.
    character(len=*), intent(inout) :: line
    
    integer hash_index
    
    hash_index = index(line, '#')
    if (hash_index > 0) then
       line = line(1:(hash_index - 1))
    end if
    
  end subroutine inout_strip_comment
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Expand all tabs in a line into single spaces (one tab makes one
  !> space).
  subroutine inout_tabs_to_spaces(line)

    !> Complete input line.
    character(len=*), intent(inout) :: line

    integer i

    do i = 1,len(line)
       if (ichar(line(i:i)) == 9) then
          line(i:i) = ' '
       end if
    end do

  end subroutine inout_tabs_to_spaces

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Strip leading spaces from a string.
  subroutine inout_strip_leading_spaces(line)

    !> Complete input line.
    character(len=*), intent(inout) :: line

    integer i

    if (len_trim(line) > 0) then
       i = verify(line, ' ') ! first character not a space
       line = line(i:)
    end if

  end subroutine inout_strip_leading_spaces

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read the next line from the inout file that contains useful data
  !> (stripping comments and blank lines).
  subroutine inout_read_next_data_line(file, line, eof)

    !> Inout file.
    type(inout_file_t), intent(inout) :: file
    !> Complete line read.
    character(len=*), intent(out) :: line
    !> True if EOF encountered.
    logical, intent(out) :: eof

    logical :: done

    done = .false.
    do while (.not. done)
       call inout_read_line_raw(file, line, eof)
       if (eof) then
          done = .true.
       else
          call inout_strip_comment(line)
          call inout_tabs_to_spaces(line)
          call inout_strip_leading_spaces(line)
          if (len_trim(line) > 0) then
             done = .true.
          end if
       end if
    end do

  end subroutine inout_read_next_data_line

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read a inout_line from the inout_file.
  subroutine inout_read_line(file, line, eof)

    !> Inout file.
    type(inout_file_t), intent(inout) :: file
    !> Inout line.
    type(inout_line_t), intent(out) :: line
    !> True if EOF encountered.
    logical, intent(out) :: eof

    character(len=MAX_LINE_LEN) :: line_string, rest
    integer i, n_data
    logical done

    call inout_read_next_data_line(file, line_string, eof)
    if (eof) return

    ! strip off the name
    i = index(line_string, ' ') ! first space
    if (i == 0) then
       write(0,'(a,i3,a,a,a)') 'ERROR: line ', file%line_num, &
            ' of input file ', trim(file%name), &
            ' contains no whitespace'
       call exit(1)
    end if
    if (i == 1) then
       write(0,'(a,i3,a,a,a)') 'ERROR: line ', file%line_num, &
            ' of input file ', trim(file%name), &
            ' starts with whitespace'
       call exit(1)
    end if
    if (i >= MAX_VAR_LEN) then
       write(0,'(a,i3,a,a,a,i6,a)') 'ERROR: line ', file%line_num, &
            ' of input file ', trim(file%name), &
            ' has a name longer than ', MAX_VAR_LEN, ' characters'
    end if
    line%name = line_string(1:(i-1))
    line_string = line_string(i:)
    call inout_strip_leading_spaces(line_string)

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
          call inout_strip_leading_spaces(rest)
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
             write(0,'(a,i3,a,a,a)') 'ERROR: line ', file%line_num, &
                  ' of input file ', trim(file%name), &
                  ' internal processing error'
             call exit(1)
          end if
          if (i >= MAX_VAR_LEN) then
             write(0,'(a,i3,a,a,a,i6,a,i6,a)') 'ERROR: line ', file%line_num, &
                  ' of input file ', trim(file%name), &
                  ' has data element ', n_data, &
                  ' longer than ', MAX_VAR_LEN, ' characters'
          end if
          line%data(n_data) = rest(1:(i-1))
          rest = rest(i:)
          call inout_strip_leading_spaces(rest)
       end if
    end do
    
  end subroutine inout_read_line

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read a inout_line from the inout_file. This will always succeed or
  !> error out, so should only be called if we know there should be a
  !> valid line coming.
  subroutine inout_read_line_no_eof(file, line)

    !> Inout file.
    type(inout_file_t), intent(inout) :: file
    !> Inout line.
    type(inout_line_t), intent(out) :: line

    logical :: eof

    call inout_read_line(file, line, eof)
    if (eof) then
       write(0,*) 'ERROR: end of file ', trim(file%name), &
            ' unexpected at line ', file%line_num
       call exit(1)
    end if

  end subroutine inout_read_line_no_eof

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read a list of inout_lines from a file, stopping at max_lines
  !> or EOF, whichever comes first.
  subroutine inout_read_line_list(file, max_lines, line_list)

    !> Inout file.
    type(inout_file_t), intent(inout) :: file
    !> Max lines to read (0 = no max).
    integer, intent(in) :: max_lines
    !> List of inout_lines,.
    type(inout_line_t), pointer :: line_list(:)
                                                ! will be allocated

    logical :: eof
    integer :: i, num_lines
    type(inout_line_t) :: temp_line_list(MAX_LIST_LINES)

    ! read file, working out how many lines we have
    num_lines = 0
    eof = .false.
    call inout_read_line(file, temp_line_list(num_lines + 1), eof)
    do while (.not. eof)
       num_lines = num_lines + 1
       if (num_lines > MAX_LIST_LINES) then
          write(0,*) 'ERROR: maximum number of lines exceeded in file ', &
               trim(file%name), ' at line ', file%line_num
          call exit(1)
       end if
       if (max_lines > 0) then
          if (num_lines >= max_lines) then
             eof = .true.
          end if
       end if
       if (.not. eof) then
          call inout_read_line(file, temp_line_list(num_lines + 1), eof)
       end if
    end do

    ! copy data to actual list
    allocate(line_list(num_lines))
    do i = 1,num_lines
       call inout_line_alloc(0, line_list(i))
       call inout_line_copy(temp_line_list(i), line_list(i))
       call inout_line_free(temp_line_list(i))
    end do

  end subroutine inout_read_line_list

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read an array of inout_lines from a file, stopping at max_lines
  !> or EOF. All lines must have the same number of elements.
  subroutine inout_read_line_array(file, max_lines, line_array)

    !> Inout file.
    type(inout_file_t), intent(inout) :: file
    !> Max lines to read (0 = no max).
    integer, intent(in) :: max_lines
    !> Array of inout_lines,.
    type(inout_line_t), pointer :: line_array(:)
                                                 ! will be allocated

    integer :: i, line_length

    call inout_read_line_list(file, max_lines, line_array)
    if (size(line_array) > 0) then
       line_length = size(line_array(1)%data)
       do i = 2,size(line_array)
          if (size(line_array(i)%data) /= line_length) then
             write(0,'(a,a,i3,a,a,a,i3,a,a)') 'ERROR: tried to read ', &
                  'array before line ', file%line_num, &
                  ' of input file ', trim(file%name), &
                  ' but line ', i,' contains a different number of', &
                  ' elements to line 1'
             call exit(1)
          end if
       end do
    end if
    
  end subroutine inout_read_line_array

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Check that the name of the line data is as given.
  subroutine inout_check_line_name(file, line, name)

    !> Inout file.
    type(inout_file_t), intent(in) :: file
    !> Inout line.
    type(inout_line_t), intent(in) :: line
    !> Expected line name.
    character(len=*), intent(in) :: name

    if (line%name /= name) then
       write(0,'(a,i3,a,a,a,a,a,a)') 'ERROR: line ', file%line_num, &
            ' of input file ', trim(file%name), &
            ' must begin with: ', trim(name), ', not: ', trim(line%name)
       call exit(1)
    end if

  end subroutine inout_check_line_name
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Checks that the read_name is the same as name.
  subroutine inout_check_name(file, name, read_name)
    
    !> Inout file.
    type(inout_file_t), intent(inout) :: file
    !> Name that we should have.
    character(len=*), intent(in) :: name
    !> Name that we do have.
    character(len=*), intent(in) :: read_name
    
    integer name_len, read_name_len

    if (name /= read_name) then
       write(0,'(a,i3,a,a,a,a,a,a)') 'ERROR: line ', file%line_num, &
            ' of input file ', trim(file%name), &
            ' must begin with: ', trim(name), ', not: ', trim(read_name)
       call exit(1)
    end if
    
  end subroutine inout_check_name
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Check that the length of the line data is as given.
  subroutine inout_check_line_length(file, line, length)

    !> Inout file.
    type(inout_file_t), intent(in) :: file
    !> Inout line.
    type(inout_line_t), intent(in) :: line
    !> Expected data length.
    integer, intent(in) :: length

    if (size(line%data) /= length) then
       write(0,*) 'ERROR: expected ', length, ' data items on line ', &
            file%line_num, ' of file ', trim(file%name)
       call exit(1)
    end if

  end subroutine inout_check_line_length
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Check the IOSTAT and error if it is bad.
  subroutine inout_check_read_iostat(file, ios, type)

    !> Inout file.
    type(inout_file_t), intent(in) :: file
    !> Iostat result.
    integer, intent(in) :: ios
    !> Type being read during error.
    character(len=*), intent(in) :: type

    if (ios /= 0) then
       write(0,'(a,a,a,a,a,i8,a,i4)') 'ERROR: reading ', trim(type), &
            ' from file ', trim(file%name), ' at line ', file%line_num, &
            ': IOSTAT = ', ios
       call exit(1)
    end if

  end subroutine inout_check_read_iostat

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Convert a string to an integer.
  integer function inout_string_to_integer(file, string)

    !> Inout file.
    type(inout_file_t), intent(in) :: file
    !> String to convert.
    character(len=*), intent(in) :: string
    
    integer :: val
    integer :: ios

    read(string, '(i20)', iostat=ios) val
    call inout_check_read_iostat(file, ios, "integer")
    inout_string_to_integer = val

  end function inout_string_to_integer

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Convert a string to an real.
  real*8 function inout_string_to_real(file, string)

    !> Inout file.
    type(inout_file_t), intent(in) :: file
    !> String to convert.
    character(len=*), intent(in) :: string
    
    real*8 :: val
    integer :: ios

    read(string, '(f30.0)', iostat=ios) val
    call inout_check_read_iostat(file, ios, "real")
    inout_string_to_real = val

  end function inout_string_to_real

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Convert a string to an logical.
  logical function inout_string_to_logical(file, string)

    !> Inout file.
    type(inout_file_t), intent(in) :: file
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
       call inout_check_read_iostat(file, 1, "logical")
    end if
    inout_string_to_logical = val

  end function inout_string_to_logical

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read an integer from a inout file that must have the given name.
  subroutine inout_read_integer(file, name, var)

    !> Inout file.
    type(inout_file_t), intent(inout) :: file
    !> Name.
    character(len=*), intent(in) :: name
    !> Variable to store data.
    integer, intent(out) :: var

    type(inout_line_t) :: line

    call inout_read_line_no_eof(file, line)
    call inout_check_line_name(file, line, name)
    call inout_check_line_length(file, line, 1)
    var = inout_string_to_integer(file, line%data(1))
    call inout_line_free(line)

  end subroutine inout_read_integer
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read a real number from a inout file that must have the given
  !> name.
  subroutine inout_read_real(file, name, var)

    !> Inout file.
    type(inout_file_t), intent(inout) :: file
    !> Name.
    character(len=*), intent(in) :: name
    !> Variable to store data.
    real*8, intent(out) :: var

    type(inout_line_t) :: line

    call inout_read_line_no_eof(file, line)
    call inout_check_line_name(file, line, name)
    call inout_check_line_length(file, line, 1)
    var = inout_string_to_real(file, line%data(1))
    call inout_line_free(line)

  end subroutine inout_read_real

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read a logical from a inout file that must have a given name.
  subroutine inout_read_logical(file, name, var)

    !> Inout file.
    type(inout_file_t), intent(inout) :: file
    !> Name.
    character(len=*), intent(in) :: name
    !> Variable to store data.
    logical, intent(out) :: var

    type(inout_line_t) :: line

    call inout_read_line_no_eof(file, line)
    call inout_check_line_name(file, line, name)
    call inout_check_line_length(file, line, 1)
    var = inout_string_to_logical(file, line%data(1))
    call inout_line_free(line)

  end subroutine inout_read_logical

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read a string from a inout file that must have a given name.
  subroutine inout_read_string(file, name, var)

    !> Inout file.
    type(inout_file_t), intent(inout) :: file
    !> Name.
    character(len=*), intent(in) :: name
    !> Variable to store data.
    character(len=*), intent(out) :: var

    type(inout_line_t) :: line

    call inout_read_line_no_eof(file, line)
    call inout_check_line_name(file, line, name)
    call inout_check_line_length(file, line, 1)
    var = line%data(1)
    call inout_line_free(line)

  end subroutine inout_read_string

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read a complex number from a inout file that must have the given
  !> name.
  subroutine inout_read_complex(file, name, var)

    !> Inout file.
    type(inout_file_t), intent(inout) :: file
    !> Name.
    character(len=*), intent(in) :: name
    !> Variable to store data.
    complex*16, intent(out) :: var

    type(inout_line_t) :: line

    call inout_read_line_no_eof(file, line)
    call inout_check_line_name(file, line, name)
    call inout_check_line_length(file, line, 2)
    var = cmplx(inout_string_to_real(file, line%data(1)), &
         inout_string_to_real(file, line%data(2)), 8)
    call inout_line_free(line)

  end subroutine inout_read_complex

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Check that indices are equal as expected.
  subroutine inout_check_index(file, index, check_index)

    !> Inout file read from.
    type(inout_file_t), intent(in) :: file
    !> Index that was read.
    integer, intent(in) :: index
    !> Expected index.
    integer, intent(in) :: check_index

    if (index /= check_index) then
       write(0,*) 'ERROR: reading from ', trim(file%name), &
            ':', file%line_num, ' expected index ', index, &
            ' but got ', check_index
       call exit(1)
    end if

  end subroutine inout_check_index

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read an integer from a inout file that must have the given index.
  subroutine inout_read_indexed_integer(file, index, var)

    !> Inout file.
    type(inout_file_t), intent(inout) :: file
    !> Expected index.
    integer, intent(in) :: index
    !> Variable to store data.
    integer, intent(out) :: var

    type(inout_line_t) :: line
    integer :: check_index

    call inout_read_line_no_eof(file, line)
    check_index = inout_string_to_integer(file, line%name)
    call inout_check_index(file, index, check_index)
    call inout_check_line_length(file, line, 1)
    var = inout_string_to_integer(file, line%data(1))
    call inout_line_free(line)

  end subroutine inout_read_indexed_integer
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read a real from a inout file that must have the given index.
  subroutine inout_read_indexed_real(file, index, var)

    !> Inout file.
    type(inout_file_t), intent(inout) :: file
    !> Expected index.
    integer, intent(in) :: index
    !> Variable to store data.
    real*8, intent(out) :: var

    type(inout_line_t) :: line
    integer :: check_index

    call inout_read_line_no_eof(file, line)
    check_index = inout_string_to_integer(file, line%name)
    call inout_check_index(file, index, check_index)
    call inout_check_line_length(file, line, 1)
    var = inout_string_to_real(file, line%data(1))
    call inout_line_free(line)

  end subroutine inout_read_indexed_real

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read a logical from a inout file that must have a given index.
  subroutine inout_read_indexed_logical(file, index, var)

    !> Inout file.
    type(inout_file_t), intent(inout) :: file
    !> Expected index.
    integer, intent(in) :: index
    !> Variable to store data.
    logical, intent(out) :: var

    type(inout_line_t) :: line
    integer :: check_index

    call inout_read_line_no_eof(file, line)
    check_index = inout_string_to_integer(file, line%name)
    call inout_check_index(file, index, check_index)
    call inout_check_line_length(file, line, 1)
    var = inout_string_to_logical(file, line%data(1))
    call inout_line_free(line)

  end subroutine inout_read_indexed_logical

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read a string from a inout file that must have a given index.
  subroutine inout_read_indexed_string(file, index, var)

    !> Inout file.
    type(inout_file_t), intent(inout) :: file
    !> Expected index.
    integer, intent(in) :: index
    !> Variable to store data.
    character(len=*), intent(out) :: var

    type(inout_line_t) :: line
    integer :: check_index

    call inout_read_line_no_eof(file, line)
    check_index = inout_string_to_integer(file, line%name)
    call inout_check_index(file, index, check_index)
    call inout_check_line_length(file, line, 1)
    var = line%data(1)
    call inout_line_free(line)

  end subroutine inout_read_indexed_string

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read a real from a inout file that must have the given indices.
  subroutine inout_read_indexed_2d_real(file, index1, index2, var)

    !> Inout file.
    type(inout_file_t), intent(inout) :: file
    !> Expected first index.
    integer, intent(in) :: index1
    !> Expected second index.
    integer, intent(in) :: index2
    !> Variable to store data.
    real*8, intent(out) :: var

    type(inout_line_t) :: line
    integer :: check_index1, check_index2

    call inout_read_line_no_eof(file, line)
    call inout_check_line_length(file, line, 2)
    check_index1 = inout_string_to_integer(file, line%name)
    check_index2 = inout_string_to_integer(file, line%data(1))
    call inout_check_index(file, index1, check_index1)
    call inout_check_index(file, index2, check_index2)
    var = inout_string_to_real(file, line%data(2))
    call inout_line_free(line)

  end subroutine inout_read_indexed_2d_real

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read an integer array.
  subroutine inout_read_integer_array(file, name, vals)

    !> Inout file.
    type(inout_file_t), intent(inout) :: file
    !> Name.
    character(len=*), intent(in) :: name
    !> Data values.
    integer, pointer :: vals(:)

    integer :: length, i

    call inout_read_integer(file, name, length)
    allocate(vals(length))
    do i = 1,length
       call inout_read_indexed_integer(file, i, vals(i))
    end do

  end subroutine inout_read_integer_array

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read a real array.
  subroutine inout_read_real_array(file, name, vals)

    !> Inout file.
    type(inout_file_t), intent(inout) :: file
    !> Name.
    character(len=*), intent(in) :: name
    !> Data values.
    real*8, pointer :: vals(:)

    integer :: length, i

    call inout_read_integer(file, name, length)
    allocate(vals(length))
    do i = 1,length
       call inout_read_indexed_real(file, i, vals(i))
    end do

  end subroutine inout_read_real_array

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read a logical array.
  subroutine inout_read_logical_array(file, name, vals)

    !> Inout file.
    type(inout_file_t), intent(inout) :: file
    !> Name.
    character(len=*), intent(in) :: name
    !> Data values.
    logical, pointer :: vals(:)

    integer :: length, i

    call inout_read_integer(file, name, length)
    allocate(vals(length))
    do i = 1,length
       call inout_read_indexed_logical(file, i, vals(i))
    end do

  end subroutine inout_read_logical_array

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read a string array.
  subroutine inout_read_string_array(file, name, vals)

    !> Inout file.
    type(inout_file_t), intent(inout) :: file
    !> Name.
    character(len=*), intent(in) :: name
    !> Data values.
    character(len=*), pointer :: vals(:)

    integer :: length, i

    call inout_read_integer(file, name, length)
    allocate(vals(length))
    do i = 1,length
       call inout_read_indexed_string(file, i, vals(i))
    end do

  end subroutine inout_read_string_array

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read a real 2d array.
  subroutine inout_read_real_array_2d(file, name, vals)

    !> Inout file.
    type(inout_file_t), intent(inout) :: file
    !> Name.
    character(len=*), intent(in) :: name
    !> Data values.
    real*8, pointer :: vals(:,:)

    integer :: size1, size2, i, j
    type(inout_line_t) :: line

    call inout_read_line_no_eof(file, line)
    call inout_check_line_name(file, line, name)
    call inout_check_line_length(file, line, 2)
    size1 = inout_string_to_integer(file, line%data(1))
    size2 = inout_string_to_integer(file, line%data(2))
    call inout_line_free(line)
    allocate(vals(size1,size2))
    do i = 1,size1
       do j = 1,size2
          call inout_read_indexed_2d_real(file, i, j, vals(i,j))
       end do
    end do

  end subroutine inout_read_real_array_2d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read an array of named lines with real data. All lines must have
  !> the same number of data elements.
  subroutine inout_read_real_named_array(file, max_lines, names, vals)

    !> Inout file.
    type(inout_file_t), intent(inout) :: file
    !> Max lines to read (0 = no max).
    integer, intent(in) :: max_lines
    !> Names of lines.
    character(len=MAX_VAR_LEN), pointer :: names(:)
    !> Data values.
    real*8, pointer :: vals(:,:)

    type(inout_line_t), pointer :: line_array(:)
    integer :: num_lines, line_length, i, j

    call inout_read_line_array(file, max_lines, line_array)
    num_lines = size(line_array)
    if (num_lines > 0) then
       line_length = size(line_array(1)%data)
       allocate(names(num_lines))
       allocate(vals(num_lines, line_length))
       do i = 1,num_lines
          names(i) = line_array(i)%name
          do j = 1,line_length
             vals(i,j) = inout_string_to_real(file, line_array(i)%data(j))
          end do
       end do
    else
       allocate(names(0))
       allocate(vals(0,0))
    end if
    do i = 1,num_lines
       call inout_line_free(line_array(i))
    end do
    deallocate(line_array)

  end subroutine inout_read_real_named_array

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read an a time-indexed array of real data.
  subroutine inout_read_timed_real_array(file, line_name, name, times, vals)

    !> Inout file.
    type(inout_file_t), intent(inout) :: file
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
    type(inout_file_t) :: read_file
    character(len=MAX_VAR_LEN), pointer :: read_names(:)
    real*8, pointer :: read_data(:,:)

    call inout_read_string(file, line_name, read_name)
    call inout_open_read(read_name, read_file)
    call inout_read_real_named_array(read_file, 0, read_names, read_data)
    call inout_close(read_file)

    n_lines = size(read_names)
    if (n_lines /= 2) then
       write(0,*) 'ERROR: must have exact two data lines in file ', &
            trim(read_name)
       call exit(1)
    end if
    n_times = size(read_data,2)
    if (n_times < 1) then
       write(0,*) 'ERROR: must have at least one data point in file ', &
            trim(read_name)
       call exit(1)
    end if
    if (trim(read_names(1)) /= "time") then
       write(0,*) 'ERROR: first data line in ', trim(read_name), &
            ' must start with: time'
       call exit(1)
    end if
    if (trim(read_names(2)) /= name) then
       write(0,*) 'ERROR: second data line in ', trim(read_name), &
            ' must start with: ', trim(name), ', not: ', trim(read_names(2))
       call exit(1)
    end if
    
    allocate(times(n_times))
    allocate(vals(n_times))
    times = read_data(1,:)
    vals = read_data(2,:)
    deallocate(read_names)
    deallocate(read_data)

  end subroutine inout_read_timed_real_array

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Write a comment string to an inout file.
  subroutine inout_write_comment(file, comment)

    !> Inout file.
    type(inout_file_t), intent(inout) :: file
    !> Comment string.
    character(len=*), intent(in) :: comment

    write(file%unit, '(a,a)') '# ', trim(comment)
    file%line_num = file%line_num + 1

  end subroutine inout_write_comment

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Write a comment string to an inout file.
  subroutine inout_check_comment(file, comment)

    !> Inout file.
    type(inout_file_t), intent(inout) :: file
    !> Comment string.
    character(len=*), intent(in) :: comment

    character(len=MAX_LINE_LEN) :: line
    logical :: eof

    call inout_read_line_raw(file, line, eof)
    if (eof) then
       write(0,*) "ERROR: EOF encountered at line ", file%line_num, &
            " of file ", trim(file%name), &
            " but expected comment: ", trim(comment)
    end if
    if ((line(1:2) /= "# ") .or. (line(3:) /= comment)) then
       write(0,*) "ERROR: at line ", file%line_num, &
            " of file ", trim(file%name), &
            " expected comment: ", trim(comment)
    end if

  end subroutine inout_check_comment

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Write an integer to an inout file.
  subroutine inout_write_integer(file, name, var)

    !> Inout file.
    type(inout_file_t), intent(inout) :: file
    !> Name.
    character(len=*), intent(in) :: name
    !> Variable to write.
    integer, intent(in) :: var

    write(file%unit, '(a20,i20)') trim(name), var
    file%line_num = file%line_num + 1

  end subroutine inout_write_integer
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Write a real to an inout file.
  subroutine inout_write_real(file, name, var)

    !> Inout file.
    type(inout_file_t), intent(inout) :: file
    !> Name.
    character(len=*), intent(in) :: name
    !> Variable to write.
    real*8, intent(in) :: var

    write(file%unit, '(a20,e30.20)') trim(name), var
    file%line_num = file%line_num + 1

  end subroutine inout_write_real
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Write a logical to an inout file.
  subroutine inout_write_logical(file, name, var)

    !> Inout file.
    type(inout_file_t), intent(inout) :: file
    !> Name.
    character(len=*), intent(in) :: name
    !> Variable to write.
    logical, intent(in) :: var

    write(file%unit, '(a20,l20)') trim(name), var
    file%line_num = file%line_num + 1

  end subroutine inout_write_logical
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Write a string to an inout file.
  subroutine inout_write_string(file, name, var)

    !> Inout file.
    type(inout_file_t), intent(inout) :: file
    !> Name.
    character(len=*), intent(in) :: name
    !> Variable to write.
    character(len=*), intent(in) :: var

    write(file%unit, '(a20,a,a)') trim(name), ' ', trim(var)
    file%line_num = file%line_num + 1

  end subroutine inout_write_string
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Write a complex to an inout file.
  subroutine inout_write_complex(file, name, var)

    !> Inout file.
    type(inout_file_t), intent(inout) :: file
    !> Name.
    character(len=*), intent(in) :: name
    !> Variable to write.
    complex*16, intent(in) :: var

    write(file%unit, '(a20,e30.20,e30.20)') trim(name), var
    file%line_num = file%line_num + 1

  end subroutine inout_write_complex
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Write an integer to an inout file.
  subroutine inout_write_indexed_integer(file, name, index, var)

    !> Inout file.
    type(inout_file_t), intent(inout) :: file
    !> Name.
    character(len=*), intent(in) :: name
    !> Index to write.
    integer, intent(in) :: index
    !> Variable to write.
    integer, intent(in) :: var

    write(file%unit, '(a20,i20,i20)') trim(name), index, var
    file%line_num = file%line_num + 1

  end subroutine inout_write_indexed_integer
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Write a real to an inout file.
  subroutine inout_write_indexed_real(file, name, index, var)

    !> Inout file.
    type(inout_file_t), intent(inout) :: file
    !> Name.
    character(len=*), intent(in) :: name
    !> Index to write.
    integer, intent(in) :: index
    !> Variable to write.
    real*8, intent(in) :: var

    write(file%unit, '(a20,i20,e30.20)') trim(name), index, var
    file%line_num = file%line_num + 1

  end subroutine inout_write_indexed_real
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Write a logical to an inout file.
  subroutine inout_write_indexed_logical(file, name, index, var)

    !> Inout file.
    type(inout_file_t), intent(inout) :: file
    !> Name.
    character(len=*), intent(in) :: name
    !> Index to write.
    integer, intent(in) :: index
    !> Variable to write.
    logical, intent(in) :: var

    write(file%unit, '(a20,i20,l20)') trim(name), index, var
    file%line_num = file%line_num + 1

  end subroutine inout_write_indexed_logical
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Write a string to an inout file.
  subroutine inout_write_indexed_string(file, name, index, var)

    !> Inout file.
    type(inout_file_t), intent(inout) :: file
    !> Name.
    character(len=*), intent(in) :: name
    !> Index to write.
    integer, intent(in) :: index
    !> Variable to write.
    character(len=*), intent(in) :: var

    write(file%unit, '(a20,i20,a20)') trim(name), index, var
    file%line_num = file%line_num + 1

  end subroutine inout_write_indexed_string
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Write a complex to an inout file.
  subroutine inout_write_indexed_complex(file, name, index, var)

    !> Inout file.
    type(inout_file_t), intent(inout) :: file
    !> Name.
    character(len=*), intent(in) :: name
    !> Index to write.
    integer, intent(in) :: index
    !> Variable to write.
    complex*16, intent(in) :: var

    write(file%unit, '(a20,i20,e30.20,e30.20)') trim(name), index, var
    file%line_num = file%line_num + 1

  end subroutine inout_write_indexed_complex
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Write a real to an inout file.
  subroutine inout_write_unnamed_real(file, var)

    !> Inout file.
    type(inout_file_t), intent(inout) :: file
    !> Variable to write.
    real*8, intent(in) :: var

    write(file%unit, '(e30.20)') var
    file%line_num = file%line_num + 1

  end subroutine inout_write_unnamed_real
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Write an integer array.
  subroutine inout_write_integer_array(file, name, vals)

    !> Inout file.
    type(inout_file_t), intent(inout) :: file
    !> Variable name.
    character(len=*), intent(in) :: name
    !> Values to write.
    integer, intent(in) :: vals(:)

    integer :: length, i

    length = size(vals)
    write(file%unit, '(a20,i20)'), trim(name), length
    do i = 1,length
       write(file%unit, '(i20,i20)'), i, vals(i)
    end do

  end subroutine inout_write_integer_array

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Write a real array.
  subroutine inout_write_real_array(file, name, vals)

    !> Inout file.
    type(inout_file_t), intent(inout) :: file
    !> Variable name.
    character(len=*), intent(in) :: name
    !> Values to write.
    real*8, intent(in) :: vals(:)

    integer :: length, i

    length = size(vals)
    write(file%unit, '(a20,i20)'), trim(name), length
    do i = 1,length
       write(file%unit, '(i20,e30.20)'), i, vals(i)
    end do

  end subroutine inout_write_real_array

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Write a logical array.
  subroutine inout_write_logical_array(file, name, vals)

    !> Inout file.
    type(inout_file_t), intent(inout) :: file
    !> Variable name.
    character(len=*), intent(in) :: name
    !> Values to write.
    logical, intent(in) :: vals(:)

    integer :: length, i

    length = size(vals)
    write(file%unit, '(a20,i20)'), trim(name), length
    do i = 1,length
       write(file%unit, '(i20,l20)'), i, vals(i)
    end do

  end subroutine inout_write_logical_array

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Write a character array.
  subroutine inout_write_string_array(file, name, vals)

    !> Inout file.
    type(inout_file_t), intent(inout) :: file
    !> Variable name.
    character(len=*), intent(in) :: name
    !> Values to write.
    character(len=*), intent(in) :: vals(:)

    integer :: length, i

    length = size(vals)
    write(file%unit, '(a20,i20)'), trim(name), length
    do i = 1,length
       write(file%unit, '(i20,a,a)'), i, ' ', trim(vals(i))
    end do

  end subroutine inout_write_string_array

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Write a real 2d array.
  subroutine inout_write_real_array_2d(file, name, vals)

    !> Inout file.
    type(inout_file_t), intent(inout) :: file
    !> Variable name.
    character(len=*), intent(in) :: name
    !> Values to write.
    real*8, intent(in) :: vals(:,:)

    integer :: size1, size2, i, j

    size1 = size(vals, 1)
    size2 = size(vals, 2)
    write(file%unit, '(a20,i20,i20)'), trim(name), size1, size2
    do i = 1,size1
       do j = 1,size2
          write(file%unit, '(i20,i20,e30.20)'), i, j, vals(i,j)
       end do
    end do

  end subroutine inout_write_real_array_2d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_inout
