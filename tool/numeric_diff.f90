! Copyright (C) 2009 Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.
!
! Compare two files containing numerical arrays and check whether they
! are the same as each other, to within the specified tolerance.
!
! If the arrays in the two files are of different sizes then they are
! automatically different. Otherwise the are the same if
!       | A_1 - A_2 |_2 < abs_tol
! and
!       | A_1 - A_2 |_2 / (| A_1 |_2 + | A_2 |_2) < rel_tol
! and are otherwise different. Setting abs_tol or rel_tol to zero
! skips the corresponding test.
!
! If the files are the same then nothing is printed on stdout,
! otherwise "different" is printed, followed by the absolute and
! relative differences, as above, or a message describing the
! difference. The files will be reported as different if they have a
! different pattern of end-of-lines and end-of-files, or if they have
! whitespace in different places (amount of whitespace is irrelevant).
!
! The exit status of the program is:
!   0 if the files are the same
!   1 if the files are different
!   2 if an error occurred

program numeric_diff

  integer, parameter :: unit1 = 40
  integer, parameter :: unit2 = 41

  character(len=1000) :: filename1, filename2, tmp
  real*8 :: abs_tol, rel_tol
  integer :: ios

  character(len=1000) :: word1, word2
  logical :: eol1, eol2, eof1, eof2
  real*8 :: value1, value2, norm1, norm2, abs_error, rel_error
  integer :: row, col
  integer :: min_row, max_row, min_col, max_col

  ! process commandline arguments
  if ((iargc() < 2) .or. (iargc() > 8)) then
     write(6,*) 'Usage: numeric_diff <filename1> <filename2> [abs_tol]' &
          // ' [rel_tol] [min_row] [max_row] [min_col] [max_col]'
     write(6,*) 'Setting tolerances or min/max values to 0 disables' &
          // ' that check.'
     write(6,*) 'If both tolerances are 0 then just print the differences.'
     write(6,*) 'All parameters default to 0 if not specified.'
     call exit(2)
  endif
  call getarg(1, filename1)
  call getarg(2, filename2)
  abs_tol = 0d0
  if (iargc() >= 3) then
     call getarg(3, tmp)
     abs_tol = string_to_real(tmp)
  end if
  rel_tol = 0d0
  if (iargc() >= 3) then
     call getarg(4, tmp)
     rel_tol = string_to_real(tmp)
  end if
  min_row = 0
  if (iargc() >= 5) then
     call getarg(5, tmp)
     min_row = string_to_integer(tmp)
  end if
  max_row = 0
  if (iargc() >= 6) then
     call getarg(6, tmp)
     max_row = string_to_integer(tmp)
  end if
  min_col = 0
  if (iargc() >= 7) then
     call getarg(7, tmp)
     min_col = string_to_integer(tmp)
  end if
  max_col = 0
  if (iargc() >= 8) then
     call getarg(8, tmp)
     max_col = string_to_integer(tmp)
  end if

  ! open files
  open(unit=unit1, status='old', file=filename1, iostat=ios)
  if (ios /= 0) then
     write(0,'(a,a,a,i4)') 'ERROR: unable to open file ', &
          trim(filename1), ' for reading: ', ios
     call exit(2)
  end if

  open(unit=unit2, status='old', file=filename2, iostat=ios)
  if (ios /= 0) then
     write(0,'(a,a,a,i4)') 'ERROR: unable to open file ', &
          trim(filename2), ' for reading: ', ios
     call exit(2)
  end if

  ! read data and compute norms
  eof1 = .false.
  row = 1
  col = 1
  norm1 = 0d0
  norm2 = 0d0
  abs_error = 0d0
  do while (.not. eof1)
     call read_word_raw(unit1, word1, eol1, eof1)
     call read_word_raw(unit2, word2, eol2, eof2)
     if (((len(word1) > 0) .and. (len(word2) == 0)) &
          .or. ((len(word1) > 0) .and. (len(word2) == 0)) &
          .or. (eol1 .and. (.not. eol2)) &
          .or. ((.not. eol1) .and. eol2) &
          .or. (eof1 .and. (.not. eof2)) &
          .or. ((.not. eof1) .and. eof2)) then
        write(*,'(a,i8,i8)') 'different shape at', row, col
        call exit(1)
     end if
     if (len(word1) > 0) then
        value1 = string_to_real(word1)
        value2 = string_to_real(word2)
        if (((min_row == 0) .or. (row >= min_row)) &
             .and. ((max_row == 0) .or. (row <= max_row)) &
             .and. ((min_col == 0) .or. (col >= min_col)) &
             .and. ((max_col == 0) .or. (col <= max_col))) then
             norm1 = norm1 + value1**2
             norm2 = norm2 + value2**2
             abs_error = abs_error + (value1 - value2)**2
          end if
        if (eol1) then
           row = row + 1
           col = 1
        else
           col = col + 1
        end if
     end if
  end do
  norm1 = sqrt(norm1)
  norm2 = sqrt(norm2)
  abs_error = sqrt(abs_error)
  rel_error = 2d0 * abs_error / (norm1 + norm2)
  
  ! check equivalence
  if ((abs_tol == 0d0) .and. (rel_tol == 0d0)) then
     write(*,'(e12.3,e12.3)') abs_error, rel_error
     call exit(0)
  end if
  if (((abs_tol == 0d0) .or. (abs_error < abs_tol)) &
       .and. ((rel_tol == 0d0) .or. (rel_error < rel_tol))) then
     write(*,*) 'files match within the given tolerances'
     call exit(0)
  end if
  write(*,'(a,e12.3,e12.3)') 'files are different', abs_error, rel_error
  call exit(1)

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Convert a string to a real.
  real*8 function string_to_real(string)

    !> String to convert.
    character(len=*), intent(in) :: string
    
    real*8 :: val
    integer :: ios

    read(string, '(e40.0)', iostat=ios) val
    if (ios /= 0) then
       write(0,'(a,a,a,i3)') 'Error converting ', trim(string), &
            ' to real: IOSTAT = ', ios
       call exit(2)
    end if
    string_to_real = val

  end function string_to_real

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Convert a string to an integer.
  integer function string_to_integer(string)

    !> String to convert.
    character(len=*), intent(in) :: string
    
    integer :: val
    integer :: ios

    read(string, '(i20)', iostat=ios) val
    if (ios /= 0) then
       write(0,'(a,a,a,i3)') 'Error converting ', trim(string), &
            ' to integer: IOSTAT = ', ios
       call exit(1)
    end if
    string_to_integer = val

  end function string_to_integer

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
    read_char = "" ! needed for pgf95 for reading blank lines
    read(unit=unit, fmt='(a)', advance='no', end=100, eor=110, &
         iostat=ios) read_char
    if (ios /= 0) then
       write(0,*) 'ERROR: reading file: IOSTAT = ', ios
       call exit(2)
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

end program numeric_diff
