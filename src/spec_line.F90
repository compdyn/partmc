! Copyright (C) 2007-2010 Nicole Riemer and Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_spec_line module.

!> A single line of formatted test for input.
module pmc_spec_line

  use pmc_util

  !> Maximum size of a single line.
  integer, parameter :: SPEC_LINE_MAX_LEN = 10000
  !> Maximum size of a variable.
  integer, parameter :: SPEC_LINE_MAX_VAR_LEN = 300

  !> A single line of input data, split at whitespace.
  !!
  !! Input lines are assumed to be in the format
  !! <pre>
  !! # stand-alone comment
  !! &lt;name&gt; &lt;whitespace&gt; &lt;data1&gt; &lt;whitespace&gt; &lt;data2&gt; ... # optional comment
  !! </pre>
  !! An spec_line_t structure stores the name and data split at
  !! whitespace.
  type spec_line_t
     !> Variable name.
     character(len=SPEC_LINE_MAX_VAR_LEN) :: name
     !> Array of data as strings.
     character(len=SPEC_LINE_MAX_VAR_LEN), pointer :: data(:)
  end type spec_line_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Allocates memory for a spec_line.
  subroutine spec_line_allocate(spec_line)

    !> Struct to alloc.
    type(spec_line_t), intent(out) :: spec_line

    allocate(spec_line%data(0))

  end subroutine spec_line_allocate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Allocates memory for a spec_line of the given size.
  subroutine spec_line_allocate_size(spec_line, n_data)

    !> Struct to alloc.
    type(spec_line_t), intent(out) :: spec_line
    !> Number of data items.
    integer, intent(in) :: n_data

    allocate(spec_line%data(n_data))

  end subroutine spec_line_allocate_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Frees all storage.
  subroutine spec_line_deallocate(spec_line)

    !> Struct to free.
    type(spec_line_t), intent(inout) :: spec_line

    deallocate(spec_line%data)

  end subroutine spec_line_deallocate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Copies a spec_line.
  subroutine spec_line_copy(from_spec_line, to_spec_line)

    !> Original spec_line.
    type(spec_line_t), intent(in) :: from_spec_line
    !> Destination, already alloced.
    type(spec_line_t), intent(inout) :: to_spec_line

    if (size(to_spec_line%data) /= size(from_spec_line%data)) then
       call spec_line_deallocate(to_spec_line)
       call spec_line_allocate_size(to_spec_line, &
            size(from_spec_line%data))
    end if
    to_spec_line%name = from_spec_line%name
    to_spec_line%data = from_spec_line%data

  end subroutine spec_line_copy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Strip the comments from a string. Comments are everything after
  !> the first # character.
  subroutine spec_line_strip_comment(string)

    !> Complete input string.
    character(len=*), intent(inout) :: string

    integer :: hash_index

    hash_index = index(string, '#')
    if (hash_index > 0) then
       string = string(1:(hash_index - 1))
    end if

  end subroutine spec_line_strip_comment

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Expand all tabs in a string into single spaces (one tab makes one
  !> space).
  subroutine spec_line_tabs_to_spaces(string)

    !> Complete input string.
    character(len=*), intent(inout) :: string

    integer :: i

    do i = 1,len(string)
       if (ichar(string(i:i)) == 9) then
          string(i:i) = ' '
       end if
    end do

  end subroutine spec_line_tabs_to_spaces

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Strip leading spaces from a string.
  subroutine spec_line_strip_leading_spaces(string)

    !> Complete input string.
    character(len=*), intent(inout) :: string

    integer :: i

    if (len_trim(string) > 0) then
       i = verify(string, ' ') ! first character not a space
       string = string(i:)
    end if

  end subroutine spec_line_strip_leading_spaces

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_spec_line
