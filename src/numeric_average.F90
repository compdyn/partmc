! Copyright (C) 2009, 2011 Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The numeric_average program.

!> Compute the mean of a sequence of files containing numerical arrays,
!> all of the same size.
program numeric_average

  integer, parameter :: dp = kind(0.d0)
  integer, parameter :: MAX_INPUT_FILES = 10000
  integer, parameter :: out_unit = 40
  integer, parameter :: in_unit_start = 41

  character(len=1000) :: filename
  integer :: ios
  character(len=1000) :: word1, word2
  logical :: eol1, eol2, eof1, eof2
  real(kind=dp) :: total
  integer :: row, col, i_file, n_file
  integer :: in_units(MAX_INPUT_FILES)

  ! process commandline arguments
  if (command_argument_count() < 2) then
     write(6,*) 'Usage: numeric_average <out_filename>' &
          // ' <in_filename_1> ... <in_filename_N>'
     stop 2
  endif
  n_file = command_argument_count() - 1
  if (n_file > MAX_INPUT_FILES) then
     write(0,*) 'ERROR: Too many input files'
     stop 1
  end if
  call get_command_argument(1, filename)
  write(*,*) "averaging output: ", trim(filename)
  open(unit=out_unit, file=filename, iostat=ios)
  if (ios /= 0) then
     write(0,'(a,a,a,i4)') 'ERROR: unable to open file ', &
          trim(filename), ' for writing: ', ios
     stop 1
  end if
  do i_file = 1,n_file
     call get_command_argument(i_file + 1, filename)
     in_units(i_file) = in_unit_start + i_file - 1
     write(*,*) "averaging input: ", trim(filename)
     open(unit=in_units(i_file), status='old', file=filename, iostat=ios)
     if (ios /= 0) then
        write(0,'(a,a,a,i4)') 'ERROR: unable to open file ', &
             trim(filename), ' for reading: ', ios
        stop 2
     end if
  end do

  ! read data and compute average
  eof1 = .false.
  row = 1
  col = 1
  do while (.not. eof1)
     total = 0d0
     call read_word_raw(in_units(1), word1, eol1, eof1)
     if (len(word1) > 0) then
        total = string_to_real(word1)
     end if
     do i_file = 2,n_file
        call read_word_raw(in_units(i_file), word2, eol2, eof2)
        if (((len(word1) > 0) .and. (len(word2) == 0)) &
             .or. ((len(word1) > 0) .and. (len(word2) == 0)) &
             .or. (eol1 .and. (.not. eol2)) &
             .or. ((.not. eol1) .and. eol2) &
             .or. (eof1 .and. (.not. eof2)) &
             .or. ((.not. eof1) .and. eof2)) then
           write(*,'(a,i8,i8,i8)') 'different shape at row/col/file:', &
                row, col, i_file
           stop 1
        end if
        if (len(word1) > 0) then
           total = total + string_to_real(word2)
        end if
     end do
     if (len(word1) > 0) then
        if (eol1) then
           row = row + 1
           col = 1
        else
           col = col + 1
        end if
        if (.not. eof1) then
           write(out_unit,'(e30.15e3)', advance='no') &
                (total / real(n_file, kind=dp))
           if (eol1) write(out_unit, '(a)') ''
        end if
     end if
  end do

  close(out_unit)
  do i_file = 1,n_file
     close(in_units(i_file))
  end do

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Convert a string to a real.
  real(kind=dp) function string_to_real(string)

    !> String to convert.
    character(len=*), intent(in) :: string

    real(kind=dp) :: val
    integer :: ios

    read(string, '(e40.0)', iostat=ios) val
    if (ios /= 0) then
       write(0,'(a,a,a,i3)') 'Error converting ', trim(string), &
            ' to real: IOSTAT = ', ios
       stop 2
    end if
    string_to_real = val

  end function string_to_real

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read a single character from a file, signaling if we have hit EOL
  !> or EOF. If EOL or EOF are true then the character value should be
  !> ignored. A file containing a single line with a single character
  !> on it will first return the character with EOL and EOF both
  !> false, then will return with EOL true but EOF false, and finally
  !> will return with EOL false and EOF true.
  subroutine read_char_raw(unit, char, eol, eof)

    !> Unit number to read from.
    integer, intent(in) :: unit
    !> Character read.
    character, intent(out) :: char
    !> True if at EOL (end of line).
    logical, intent(out) :: eol
    !> True if at EOF (end of file).
    logical, intent(out) :: eof

    integer :: ios, n_read
    character(len=1) :: read_char

    eol = .false.
    eof = .false.
    char = " " ! shut up uninitialized variable warnings
    read_char = "" ! needed for pgf95 for reading blank lines
    read(unit=unit, fmt='(a)', advance='no', end=100, eor=110, &
         iostat=ios) read_char
    if (ios /= 0) then
       write(0,*) 'ERROR: reading file: IOSTAT = ', ios
       stop 2
    end if
    ! only reach here if we didn't hit end-of-record (end-of-line) in
    ! the above read
    char = read_char
    goto 120

100 eof = .true. ! goto here if end-of-file was encountered immediately
    goto 120

110 eol = .true. ! goto here if end-of-record, meaning end-of-line

120 return

  end subroutine read_char_raw

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read a white-space delimited word from a file, signaling if we
  !> have EOL or EOF. If EOL or EOF are true then the word will still
  !> be meaningful data. If there was no data to be read then
  !> len(word) will be 0.
  subroutine read_word_raw(unit, word, eol, eof)

    !> Unit number to read from.
    integer, intent(in) :: unit
    !> Word read.
    character(len=*), intent(out) :: word
    !> True if at EOL (end of line).
    logical, intent(out) :: eol
    !> True if at EOF (end of file).
    logical, intent(out) :: eof

    integer :: i
    character :: char

    word = ""

    ! skip over spaces
    call read_char_raw(unit, char, eol, eof)
    do while (((ichar(char) == 9) .or. (ichar(char) == 32)) &
         .and. (.not. eol) .and. (.not. eof))
       call read_char_raw(unit, char, eol, eof)
    end do
    if (eol .or. eof) return

    ! char is now the first word character
    i = 1
    word(i:i) = char
    call read_char_raw(unit, char, eol, eof)
    do while ((ichar(char) /= 9) .and. (ichar(char) /= 32) &
         .and. (.not. eol) .and. (.not. eof))
       i = i + 1
       word(i:i) = char
       call read_char_raw(unit, char, eol, eof)
    end do

  end subroutine read_word_raw

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program numeric_average
