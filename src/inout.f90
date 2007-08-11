! Copyright (C) 2007 Nicole Riemer and Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.
!
! Generic input and output routines.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module pmc_inout

  integer, parameter :: MAX_CHAR_LEN = 300 ! max size of line or variable
  integer, parameter :: MAX_LIST_LINES = 500 ! max lines in an array

  type inout_file_t
     character(len=MAX_CHAR_LEN) :: name ! filename
     integer :: unit                    ! attached unit
     integer :: line_num                ! current line number
  end type inout_file_t

  type inout_line_t
     character(len=MAX_CHAR_LEN) :: name ! variable name
     character(len=MAX_CHAR_LEN), pointer :: data(:) ! array of data as strings
  end type inout_line_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine inout_open_read(filename, file)

    ! Open an inout file for reading.

    use pmc_util

    character(len=*), intent(in) :: filename ! name of file to open
    type(inout_file_t), intent(out) :: file  ! inout file

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

  subroutine inout_open_write(filename, file)

    ! Open an inout file for writing.

    use pmc_util

    character(len=*), intent(in) :: filename ! name of file to open
    type(inout_file_t), intent(out) :: file  ! inout file

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

  subroutine inout_close(file)

    ! Close an inout file.

    use pmc_util

    type(inout_file_t), intent(in) :: file ! inout file

    close(file%unit)
    call free_unit(file%unit)

  end subroutine inout_close

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine inout_line_alloc(n_data, inout_line)

    ! Allocates memory for an inout_line.

    integer, intent(in) :: n_data       ! number of data items
    type(inout_line_t), intent(inout) :: inout_line ! struct to alloc

    allocate(inout_line%data(n_data))

  end subroutine inout_line_alloc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine inout_line_free(inout_line)

    ! Frees all storage.

    type(inout_line_t), intent(inout) :: inout_line ! struct to free

    deallocate(inout_line%data)

  end subroutine inout_line_free

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine inout_line_copy(from_line, to_line)

    ! Copies a inout_line.

    type(inout_line_t), intent(in) :: from_line ! original line
    type(inout_line_t), intent(out) :: to_line ! destination, already alloced

    to_line%name = from_line%name
    if (size(to_line%data) /= size(from_line%data)) then
       deallocate(to_line%data)
    end if
    allocate(to_line%data(size(from_line%data)))
    to_line%data = from_line%data

  end subroutine inout_line_copy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine inout_read_line_raw(file, line, eof)

    ! Read a single line from a inout file, signaling if we have hit EOF.

    type(inout_file_t), intent(inout) :: file ! inout file
    character(len=*), intent(out) :: line ! complete line read
    logical, intent(out) :: eof           ! true if at EOF

    integer ios

    read(unit=file%unit, fmt='(a)', end=100, iostat=ios) line
    file%line_num = file%line_num + 1
    if (ios /= 0) then
       write(0,*) 'ERROR: reading from ', trim(file%name), &
            ' at line ', file%line_num, ': IOSTAT = ', ios
       call exit(1)
    end if
    eof = .false.
    return

100 line = ""
    eof = .true.
    
  end subroutine inout_read_line_raw

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine inout_strip_comment(line)

    ! Strip the comments from a line. Comments are everything after
    ! the first # character.
    
    character(len=*), intent(inout) :: line ! complete input line
    
    integer hash_index
    
    hash_index = index(line, '#')
    if (hash_index > 0) then
       line = line(1:(hash_index - 1))
    end if
    
  end subroutine inout_strip_comment
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine inout_tabs_to_spaces(line)

    ! Expand all tabs in a line into single spaces (one tab makes one
    ! space).

    character(len=*), intent(inout) :: line ! complete input line

    integer i

    do i = 1,len(line)
       if (ichar(line(i:i)) == 9) then
          line(i:i) = ' '
       end if
    end do

  end subroutine inout_tabs_to_spaces

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine inout_strip_leading_spaces(line)

    ! Strip leading spaces from a string.

    character(len=*), intent(inout) :: line ! complete input line

    integer i

    if (len_trim(line) > 0) then
       i = verify(line, ' ') ! first character not a space
       line = line(i:)
    end if

  end subroutine inout_strip_leading_spaces

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine inout_read_next_data_line(file, line, eof)

    ! Read the next line from the inout file that contains useful data
    ! (stripping comments and blank lines).

    type(inout_file_t), intent(inout) :: file ! inout file
    character(len=*), intent(out) :: line ! complete line read
    logical, intent(out) :: eof           ! true if EOF encountered

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

  subroutine inout_read_line(file, line, eof)

    ! Read a inout_line from the inout_file.

    type(inout_file_t), intent(inout) :: file ! inout file
    type(inout_line_t), intent(out) :: line  ! inout line
    logical, intent(out) :: eof           ! true if EOF encountered

    character(len=MAX_CHAR_LEN) :: line_string, rest
    integer i, n_data
    logical done

    call inout_read_next_data_line(file, line_string, eof)
    if (eof) return

    ! strip off the name
    i = index(line_string, ' ') ! first space
    if (i == 0) then
       write(0,'(a,i3,a,a,a)') 'ERROR: line ', file%line_num, &
            ' of input file ', trim(file%name), &
            ' is malformed'
       call exit(1)
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
          line%data(n_data) = rest(1:(i-1))
          rest = rest(i:)
          call inout_strip_leading_spaces(rest)
       end if
    end do
    
  end subroutine inout_read_line

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine inout_read_line_no_eof(file, line)

    ! Read a inout_line from the inout_file. This will always succeed or
    ! error out, so should only be called if we know there should be a
    ! valid line coming.

    type(inout_file_t), intent(inout) :: file ! inout file
    type(inout_line_t), intent(out) :: line  ! inout line

    logical :: eof

    call inout_read_line(file, line, eof)
    if (eof) then
       write(0,*) 'ERROR: end of file ', trim(file%name), &
            ' unexpected at line ', file%line_num
       call exit(1)
    end if

  end subroutine inout_read_line_no_eof

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine inout_read_line_list(file, max_lines, line_list)

    ! Read a list of inout_lines from a file, stopping at max_lines
    ! or EOF, whichever comes first.

    type(inout_file_t), intent(inout) :: file ! inout file
    integer, intent(in) :: max_lines      ! max lines to read (0 = no max)
    type(inout_line_t), pointer :: line_list(:) ! list of inout_lines,
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

  subroutine inout_read_line_array(file, max_lines, line_array)

    ! Read an array of inout_lines from a file, stopping at max_lines
    ! or EOF. All lines must have the same number of elements.

    type(inout_file_t), intent(inout) :: file ! inout file
    integer, intent(in) :: max_lines      ! max lines to read (0 = no max)
    type(inout_line_t), pointer :: line_array(:) ! array of inout_lines,
                                                 ! will be allocated

    integer i, line_length

    call inout_read_line_list(file, max_lines, line_array)
    if (size(line_array) > 0) then
       line_length = size(line_array(1)%data)
       do i = 2,size(line_array)
          if (size(line_array(i)%data) /= line_length) then
             write(0,'(a,a,i3,a,a,a)') 'ERROR: tried to read ', &
                  'array before line ', file%line_num, &
                  ' of input file ', trim(file%name), &
                  ' but lines contain varying numbers of elements'
             call exit(1)
          end if
       end do
    end if
    
  end subroutine inout_read_line_array

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine inout_check_line_name(file, line, name)

    ! Check that the name of the line data is as given.

    type(inout_file_t), intent(in) :: file ! inout file
    type(inout_line_t), intent(in) :: line ! inout line
    character(len=*), intent(in) :: name ! expected line name

    if (line%name /= name) then
       write(0,'(a,i3,a,a,a,a,a,a)') 'ERROR: line ', file%line_num, &
            ' of input file ', trim(file%name), &
            ' must begin with: ', trim(name), ', not: ', trim(line%name)
       call exit(1)
    end if

  end subroutine inout_check_line_name
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine inout_check_name(file, name, read_name)

    ! Checks that the read_name is the same as name.
    
    type(inout_file_t), intent(inout) :: file ! inout file
    character(len=*), intent(in) :: name   ! name that we should have
    character(len=*), intent(in) :: read_name ! name that we do have
    
    integer name_len, read_name_len

    if (name /= read_name) then
       write(0,'(a,i3,a,a,a,a,a,a)') 'ERROR: line ', file%line_num, &
            ' of input file ', trim(file%name), &
            ' must begin with: ', trim(name), ', not: ', trim(read_name)
       call exit(1)
    end if
    
  end subroutine inout_check_name
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine inout_check_line_length(file, line, length)

    ! Check that the length of the line data is as given.

    type(inout_file_t), intent(in) :: file ! inout file
    type(inout_line_t), intent(in) :: line ! inout line
    integer, intent(in) :: length       ! expected data length

    if (size(line%data) /= length) then
       write(0,*) 'ERROR: expected ', length, ' data items on line ', &
            file%line_num, ' of file ', trim(file%name)
       call exit(1)
    end if

  end subroutine inout_check_line_length
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine inout_check_read_iostat(file, ios, type)

    ! Check the IOSTAT and error if it is bad.

    type(inout_file_t), intent(in) :: file ! inout file
    integer, intent(in) :: ios          ! iostat result
    character(len=*), intent(in) :: type ! type being read during error

    if (ios /= 0) then
       write(0,'(a,a,a,a,a,i8,a,i4)') 'ERROR: reading ', trim(type), &
            ' from file ', trim(file%name), ' at line ', file%line_num, &
            ': IOSTAT = ', ios
       call exit(1)
    end if

  end subroutine inout_check_read_iostat

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function inout_string_to_integer(file, string)

    ! Convert a string to an integer.

    type(inout_file_t), intent(in) :: file ! inout file
    character(len=*), intent(in) :: string ! string to convert
    
    integer :: val
    integer :: ios

    read(string, '(i20)', iostat=ios) val
    call inout_check_read_iostat(file, ios, "integer")
    inout_string_to_integer = val

  end function inout_string_to_integer

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real*8 function inout_string_to_real(file, string)

    ! Convert a string to an real.

    type(inout_file_t), intent(in) :: file ! inout file
    character(len=*), intent(in) :: string ! string to convert
    
    real*8 :: val
    integer :: ios

    read(string, '(f30.0)', iostat=ios) val
    call inout_check_read_iostat(file, ios, "real")
    inout_string_to_real = val

  end function inout_string_to_real

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  logical function inout_string_to_logical(file, string)

    ! Convert a string to an logical.

    type(inout_file_t), intent(in) :: file ! inout file
    character(len=*), intent(in) :: string ! string to convert
    
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

  subroutine inout_read_integer(file, name, var)

    ! Read an integer from a inout file that must have the given name.

    type(inout_file_t), intent(inout) :: file ! inout file
    character(len=*), intent(in) :: name  ! name
    integer, intent(out) :: var           ! variable to store data

    type(inout_line_t) :: line

    call inout_read_line_no_eof(file, line)
    call inout_check_line_name(file, line, name)
    call inout_check_line_length(file, line, 1)
    var = inout_string_to_integer(file, line%data(1))
    call inout_line_free(line)

  end subroutine inout_read_integer
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine inout_read_real(file, name, var)

    ! Read a real number from a inout file that must have the given
    ! name.

    type(inout_file_t), intent(inout) :: file ! inout file
    character(len=*), intent(in) :: name  ! name
    real*8, intent(out) :: var            ! variable to store data

    type(inout_line_t) :: line

    call inout_read_line_no_eof(file, line)
    call inout_check_line_name(file, line, name)
    call inout_check_line_length(file, line, 1)
    var = inout_string_to_real(file, line%data(1))
    call inout_line_free(line)

  end subroutine inout_read_real

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine inout_read_logical(file, name, var)

    ! Read a logical from a inout file that must have a given name.

    type(inout_file_t), intent(inout) :: file ! inout file
    character(len=*), intent(in) :: name  ! name
    logical, intent(out) :: var           ! variable to store data

    type(inout_line_t) :: line

    call inout_read_line_no_eof(file, line)
    call inout_check_line_name(file, line, name)
    call inout_check_line_length(file, line, 1)
    var = inout_string_to_logical(file, line%data(1))
    call inout_line_free(line)

  end subroutine inout_read_logical

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine inout_read_string(file, name, var)

    ! Read a string from a inout file that must have a given name.

    type(inout_file_t), intent(inout) :: file ! inout file
    character(len=*), intent(in) :: name  ! name
    character(len=*), intent(out) :: var  ! variable to store data

    type(inout_line_t) :: line

    call inout_read_line_no_eof(file, line)
    call inout_check_line_name(file, line, name)
    call inout_check_line_length(file, line, 1)
    var = line%data(1)
    call inout_line_free(line)

  end subroutine inout_read_string

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine inout_check_index(file, index, check_index)

    ! Check that indices are equal as expected.

    type(inout_file_t), intent(in) :: file ! inout file read from
    integer, intent(in) :: index        ! index that was read
    integer, intent(in) :: check_index  ! expected index

    if (index /= check_index) then
       write(0,*) 'ERROR: reading from ', trim(file%name), &
            ':', file%line_num, ' expected index ', index, &
            ' but got ', check_index
       call exit(1)
    end if

  end subroutine inout_check_index

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine inout_read_indexed_integer(file, index, var)

    ! Read an integer from a inout file that must have the given index.

    type(inout_file_t), intent(inout) :: file ! inout file
    integer, intent(in) :: index          ! expected index
    integer, intent(out) :: var           ! variable to store data

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

  subroutine inout_read_indexed_real(file, index, var)

    ! Read a real from a inout file that must have the given index.

    type(inout_file_t), intent(inout) :: file ! inout file
    integer, intent(in) :: index          ! expected index
    real*8, intent(out) :: var            ! variable to store data

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

  subroutine inout_read_indexed_logical(file, index, var)

    ! Read a logical from a inout file that must have a given index.

    type(inout_file_t), intent(inout) :: file ! inout file
    integer, intent(in) :: index          ! expected index
    logical, intent(out) :: var           ! variable to store data

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

  subroutine inout_read_indexed_string(file, index, var)

    ! Read a string from a inout file that must have a given index.

    type(inout_file_t), intent(inout) :: file ! inout file
    integer, intent(in) :: index          ! expected index
    character(len=*), intent(out) :: var  ! variable to store data

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

  subroutine inout_read_indexed_2d_real(file, index1, index2, var)

    ! Read a real from a inout file that must have the given indices.

    type(inout_file_t), intent(inout) :: file ! inout file
    integer, intent(in) :: index1         ! expected first index
    integer, intent(in) :: index2         ! expected second index
    real*8, intent(out) :: var            ! variable to store data

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

  subroutine inout_read_integer_array(file, name, vals)

    ! Read an integer array.

    type(inout_file_t), intent(inout) :: file ! inout file
    character(len=*), intent(in) :: name  ! name
    integer, pointer :: vals(:)            ! data values

    integer :: length, i

    call inout_read_integer(file, name, length)
    allocate(vals(length))
    do i = 1,length
       call inout_read_indexed_integer(file, i, vals(i))
    end do

  end subroutine inout_read_integer_array

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine inout_read_real_array(file, name, vals)

    ! Read a real array.

    type(inout_file_t), intent(inout) :: file ! inout file
    character(len=*), intent(in) :: name  ! name
    real*8, pointer :: vals(:)            ! data values

    integer :: length, i

    call inout_read_integer(file, name, length)
    allocate(vals(length))
    do i = 1,length
       call inout_read_indexed_real(file, i, vals(i))
    end do

  end subroutine inout_read_real_array

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine inout_read_logical_array(file, name, vals)

    ! Read a logical array.

    type(inout_file_t), intent(inout) :: file ! inout file
    character(len=*), intent(in) :: name  ! name
    logical, pointer :: vals(:)            ! data values

    integer :: length, i

    call inout_read_integer(file, name, length)
    allocate(vals(length))
    do i = 1,length
       call inout_read_indexed_logical(file, i, vals(i))
    end do

  end subroutine inout_read_logical_array

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine inout_read_string_array(file, name, vals)

    ! Read a string array.

    type(inout_file_t), intent(inout) :: file ! inout file
    character(len=*), intent(in) :: name  ! name
    character(len=*), pointer :: vals(:) ! data values

    integer :: length, i

    call inout_read_integer(file, name, length)
    allocate(vals(length))
    do i = 1,length
       call inout_read_indexed_string(file, i, vals(i))
    end do

  end subroutine inout_read_string_array

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine inout_read_real_array_2d(file, name, vals)

    ! Read a real 2d array.

    type(inout_file_t), intent(inout) :: file ! inout file
    character(len=*), intent(in) :: name  ! name
    real*8, pointer :: vals(:,:)          ! data values

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

  subroutine inout_read_real_named_array(file, max_lines, names, vals)

    ! Read an array of named lines with real data. All lines must have
    ! the same number of data elements.

    type(inout_file_t), intent(inout) :: file ! spec file
    integer, intent(in) :: max_lines      ! max lines to read (0 = no max)
    character(len=MAX_CHAR_LEN), pointer :: names(:) ! names of lines
    real*8, pointer :: vals(:,:)          ! data values

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

  subroutine inout_read_timed_real_array(file, line_name, name, times, vals)

    ! Read an a time-indexed array of real data.

    type(inout_file_t), intent(inout) :: file ! spec file
    character(len=*), intent(in) :: line_name ! name of line for filename
    character(len=*), intent(in) :: name  ! variable name
    real*8, pointer :: times(:)         ! names of lines
    real*8, pointer :: vals(:)          ! data values
    
    integer :: n_lines, n_times
    character(len=MAX_CHAR_LEN) :: read_name
    type(inout_file_t) :: read_file
    character(len=MAX_CHAR_LEN), pointer :: read_names(:)
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

  subroutine inout_write_integer(file, name, var)

    ! Write an integer to an inout file.

    type(inout_file_t), intent(inout) :: file ! inout file
    character(len=*), intent(in) :: name ! name
    integer, intent(in) :: var          ! variable to write

    write(file%unit, '(a20,i20)') trim(name), var
    file%line_num = file%line_num + 1

  end subroutine inout_write_integer
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine inout_write_real(file, name, var)

    ! Write a real to an inout file.

    type(inout_file_t), intent(inout) :: file ! inout file
    character(len=*), intent(in) :: name ! name
    real*8, intent(in) :: var          ! variable to write

    write(file%unit, '(a20,e30.20)') trim(name), var
    file%line_num = file%line_num + 1

  end subroutine inout_write_real
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine inout_write_logical(file, name, var)

    ! Write a logical to an inout file.

    type(inout_file_t), intent(inout) :: file ! inout file
    character(len=*), intent(in) :: name ! name
    logical, intent(in) :: var          ! variable to write

    write(file%unit, '(a20,l20)') trim(name), var
    file%line_num = file%line_num + 1

  end subroutine inout_write_logical
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine inout_write_string(file, name, var)

    ! Write a string to an inout file.

    type(inout_file_t), intent(inout) :: file ! inout file
    character(len=*), intent(in) :: name ! name
    character(len=*), intent(in) :: var ! variable to write

    write(file%unit, '(a20,a20)') trim(name), var
    file%line_num = file%line_num + 1

  end subroutine inout_write_string
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine inout_write_integer_array(file, name, vals)

    ! Write an integer array.

    type(inout_file_t), intent(inout) :: file ! inout file
    character(len=*), intent(in) :: name ! variable name
    integer, intent(in) :: vals(:) ! values to write

    integer :: length, i

    length = size(vals)
    write(file%unit, '(a20,i20)'), trim(name), length
    do i = 1,length
       write(file%unit, '(i20,i20)'), i, vals(i)
    end do

  end subroutine inout_write_integer_array

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine inout_write_real_array(file, name, vals)

    ! Write a real array.

    type(inout_file_t), intent(inout) :: file ! inout file
    character(len=*), intent(in) :: name ! variable name
    real*8, intent(in) :: vals(:) ! values to write

    integer :: length, i

    length = size(vals)
    write(file%unit, '(a20,i20)'), trim(name), length
    do i = 1,length
       write(file%unit, '(i20,e30.20)'), i, vals(i)
    end do

  end subroutine inout_write_real_array

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine inout_write_logical_array(file, name, vals)

    ! Write a logical array.

    type(inout_file_t), intent(inout) :: file ! inout file
    character(len=*), intent(in) :: name ! variable name
    logical, intent(in) :: vals(:) ! values to write

    integer :: length, i

    length = size(vals)
    write(file%unit, '(a20,i20)'), trim(name), length
    do i = 1,length
       write(file%unit, '(i20,l20)'), i, vals(i)
    end do

  end subroutine inout_write_logical_array

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine inout_write_string_array(file, name, vals)

    ! Write a character array.

    type(inout_file_t), intent(inout) :: file ! inout file
    character(len=*), intent(in) :: name ! variable name
    character(len=*), intent(in) :: vals(:) ! values to write

    integer :: length, i

    length = size(vals)
    write(file%unit, '(a20,i20)'), trim(name), length
    do i = 1,length
       write(file%unit, '(i20,a20)'), i, vals(i)
    end do

  end subroutine inout_write_string_array

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine inout_write_real_array_2d(file, name, vals)

    ! Write a real 2d array.

    type(inout_file_t), intent(inout) :: file ! inout file
    character(len=*), intent(in) :: name ! variable name
    real*8, intent(in) :: vals(:,:) ! values to write

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
