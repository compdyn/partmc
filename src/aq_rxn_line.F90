! Copyright (C) 2007-2015 Nicole Riemer and Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_aq_rxn_line module.

!> A single line of formatted text for input.
module pmc_aq_rxn_line

  use pmc_util

  implicit none

  !> Maximum size of a single line.
  integer, parameter :: AQ_RXN_LINE_MAX_LEN = 10000
  !> Maximum size of a variable.
  integer, parameter :: AQ_RXN_LINE_MAX_VAR_LEN = 300

  !> A single line of input data, split at whitespace into a variable number of terms
  !!
  !! An aq_rxn_line_t structure stores the terms from a single line
  type aq_rxn_line_t
     !> Array of data as strings.
     character(len=AQ_RXN_LINE_MAX_VAR_LEN), pointer :: data(:)
  end type aq_rxn_line_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Allocates memory for a aq_rxn_line.
  subroutine aq_rxn_line_allocate(aq_rxn_line)

    !> Struct to alloc.
    type(aq_rxn_line_t), intent(out) :: aq_rxn_line

    allocate(aq_rxn_line%data(0))

  end subroutine aq_rxn_line_allocate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Allocates memory for a aq_rxn_line of the given size.
  subroutine aq_rxn_line_allocate_size(aq_rxn_line, n_data)

    !> Struct to alloc.
    type(aq_rxn_line_t), intent(out) :: aq_rxn_line
    !> Number of data items.
    integer, intent(in) :: n_data

    allocate(aq_rxn_line%data(n_data))

  end subroutine aq_rxn_line_allocate_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Frees all storage.
  subroutine aq_rxn_line_deallocate(aq_rxn_line)

    !> Struct to free.
    type(aq_rxn_line_t), intent(inout) :: aq_rxn_line

    deallocate(aq_rxn_line%data)

  end subroutine aq_rxn_line_deallocate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Copies a aq_rxn_line.
  subroutine aq_rxn_line_copy(from_aq_rxn_line, to_aq_rxn_line)

    !> Original aq_rxn_line.
    type(aq_rxn_line_t), intent(in) :: from_aq_rxn_line
    !> Destination, already alloced.
    type(aq_rxn_line_t), intent(inout) :: to_aq_rxn_line

    if (size(to_aq_rxn_line%data) /= size(from_aq_rxn_line%data)) then
       call aq_rxn_line_deallocate(to_aq_rxn_line)
       call aq_rxn_line_allocate_size(to_aq_rxn_line, &
            size(from_aq_rxn_line%data))
    end if
    to_aq_rxn_line%data = from_aq_rxn_line%data

  end subroutine aq_rxn_line_copy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Strip the comments from a string. Comments are everything after
  !> the "COMMENT" keyword.
  subroutine aq_rxn_line_strip_comment(string)

    !> Complete input string.
    character(len=*), intent(inout) :: string

    integer :: comment_index

    comment_index = index(string, "COMMENT")
    if (comment_index > 0) then
       string = string(1:(comment_index - 1))
    end if

  end subroutine aq_rxn_line_strip_comment

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Expand all tabs in a string into single spaces (one tab makes one
  !> space).
  subroutine aq_rxn_line_tabs_to_spaces(string)

    !> Complete input string.
    character(len=*), intent(inout) :: string

    integer :: i

    do i = 1,len(string)
       if (ichar(string(i:i)) == 9) then
          string(i:i) = ' '
       end if
    end do

  end subroutine aq_rxn_line_tabs_to_spaces

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Strip leading spaces from a string.
  subroutine aq_rxn_line_strip_leading_spaces(string)

    !> Complete input string.
    character(len=*), intent(inout) :: string

    integer :: i

    if (len_trim(string) > 0) then
       i = verify(string, ' ') ! first character not a space
       string = string(i:)
    end if

  end subroutine aq_rxn_line_strip_leading_spaces

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_aq_rxn_line