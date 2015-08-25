! Copyright (C) Matthew Dawson
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_aq_rxn_file module.

!> Reading formatted text input.
module pmc_aq_rxn_file

  use pmc_aq_rxn_line
  use pmc_aq_rxn_data
  use pmc_aq_spec_data
  use pmc_util

  implicit none

  !> Maximum number of lines in an array.
  integer, parameter :: AQ_RXN_FILE_MAX_LIST_LINES = 1000

  !> An input file with extra data for printing messages.
  !!
  !! An aq_rxn_file_t is just a simple wrapper around a Fortran unit
  !! number together with the filename and current line number. The
  !! line number is updated manually by the various \c aq_rxn_file_*()
  !! subroutine. To maintain its validity all file accesses must be
  !! done via the \c aq_rxn_file_*() subroutines, and no data should be
  !! accessed directly via \c aq_rxn_file%%unit.
  type aq_rxn_file_t
     !> Filename.
     character(len=AQ_RXN_LINE_MAX_VAR_LEN) :: name
     !> Attached unit.
     integer :: unit
     !> Current line number.
     integer :: line_num
  end type aq_rxn_file_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Exit with an error message containing filename and line number.
  subroutine aq_rxn_file_die_msg(code, file, msg)

    !> Failure status code.
    integer, intent(in) :: code
    !> Aqueous reaction file.
    type(aq_rxn_file_t), intent(in) :: file
    !> Error message.
    character(len=*), intent(in) :: msg

    call aq_rxn_file_assert_msg(code, file, .false., msg)

  end subroutine aq_rxn_file_die_msg

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Exit with an error message containing filename and line number
  !> if \c condition_ok is \c .false.
  subroutine aq_rxn_file_assert_msg(code, file, condition_ok, msg)

    !> Failure status code.
    integer, intent(in) :: code
    !> Aqueous reaction file.
    type(aq_rxn_file_t), intent(in) :: file
    !> Whether the assertion is ok.
    logical, intent(in) :: condition_ok
    !> Error message.
    character(len=*), intent(in) :: msg

    if (.not. condition_ok) then
       call die_msg(code, "file " // trim(file%name) // " line " &
            // trim(integer_to_string(file%line_num)) // ": " // trim(msg))
    end if

  end subroutine aq_rxn_file_assert_msg

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Open an aqueous reaction file for reading.
  subroutine aq_rxn_file_open(filename, file)

    !> Name of file to open.
    character(len=*), intent(in) :: filename
    !> Aqueous reaction file.
    type(aq_rxn_file_t), intent(out) :: file

    integer :: ios

    file%name = trim(filename)
    file%unit = get_unit()
    open(unit=file%unit, status='old', file=file%name, iostat=ios)
    if (ios /= 0) then
       call die_msg(602785622, "unable to open file " // trim(file%name) &
            // " for reading: IOSTAT = " // trim(integer_to_string(ios)))
    end if
    file%line_num = 0

  end subroutine aq_rxn_file_open

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Close an aqueous reaction file.
  subroutine aq_rxn_file_close(file)

    !> Aqueous reaction file.
    type(aq_rxn_file_t), intent(in) :: file

    close(file%unit)
    call free_unit(file%unit)

  end subroutine aq_rxn_file_close

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read a single line from an aqueous reaction file, signaling if we have hit EOF.
  subroutine aq_rxn_file_read_line_raw(file, line, eof)

    !> Aqueous reaction file.
    type(aq_rxn_file_t), intent(inout) :: file
    !> Complete line read.
    character(len=*), intent(out) :: line
    !> True if at EOF.
    logical, intent(out) :: eof

    integer :: ios

    file%line_num = file%line_num + 1
    eof = .false.
    line = "" ! needed for pgf95 for reading blank lines
    read(unit=file%unit, fmt='(a)', advance='no', end=100, eor=110, &
         iostat=ios) line
    if (ios /= 0) then
       call aq_rxn_file_die_msg(603104955, file, &
            'error reading: IOSTAT = ' // trim(integer_to_string(ios)))
    end if
    ! only reach here if we didn't hit end-of-record (end-of-line) in
    ! the above read, meaning the line was too long
    call aq_rxn_file_die_msg(603222604, file, &
         'line exceeds length: ' // trim(integer_to_string(len(line))))

100 line = "" ! goto here if end-of-file was encountered immediately
    eof = .true.

110 return ! goto here if end-of-record, meaning everything is ok

  end subroutine aq_rxn_file_read_line_raw

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read the next line from the aqueous reaction file that contains useful data
  !> (stripping comments and blank lines).
  subroutine aq_rxn_file_read_next_data_line(file, line, eof)

    !> Aqueous reaction file.
    type(aq_rxn_file_t), intent(inout) :: file
    !> Complete line read.
    character(len=*), intent(out) :: line
    !> True if EOF encountered.
    logical, intent(out) :: eof

    logical :: done

    done = .false.
    do while (.not. done)
       call aq_rxn_file_read_line_raw(file, line, eof)
       if (eof) then
          done = .true.
       else
          call aq_rxn_line_strip_comment(line)
          call aq_rxn_line_tabs_to_spaces(line)
          call aq_rxn_line_strip_leading_spaces(line)
          if (len_trim(line) > 0) then
             done = .true.
          end if
       end if
    end do

  end subroutine aq_rxn_file_read_next_data_line

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read an aq_rxn_line from the aq_rxn_file.
  subroutine aq_rxn_file_read_line(file, line, eof)

    !> Aqueous reaction file.
    type(aq_rxn_file_t), intent(inout) :: file
    !> Aqueous reaction line.
    type(aq_rxn_line_t), intent(inout) :: line
    !> True if EOF encountered.
    logical, intent(out) :: eof

    character(len=aq_rxn_LINE_MAX_LEN) :: line_string, rest
    integer i, n_data
    logical done

    call aq_rxn_file_read_next_data_line(file, line_string, eof)
    if (eof) return

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
          call aq_rxn_line_strip_leading_spaces(rest)
       end if
    end do

    ! allocate the data and read out the data items
    call aq_rxn_line_deallocate(line)
    call aq_rxn_line_allocate_size(line, n_data)
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
             call aq_rxn_file_die_msg(603827656, file, &
                  'internal processing error')
          end if
          if (i >= aq_rxn_LINE_MAX_VAR_LEN) then
             call aq_rxn_file_die_msg(603945305, file, &
                  'data element ' // trim(integer_to_string(n_data)) &
                  // ' longer than: ' &
                  // trim(integer_to_string(AQ_RXN_LINE_MAX_VAR_LEN)))
          end if
          line%data(n_data) = rest(1:(i-1))
          rest = rest(i:)
          call aq_rxn_line_strip_leading_spaces(rest)
       end if
    end do

  end subroutine aq_rxn_file_read_line

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read an aq_rxn_line from the aq_rxn_file. This will always succeed or
  !> error out, so should only be called if we know there should be a
  !> valid line coming.
  subroutine aq_rxn_file_read_line_no_eof(file, line)

    !> Aqueous reaction file.
    type(aq_rxn_file_t), intent(inout) :: file
    !> Aqueous reaction line.
    type(aq_rxn_line_t), intent(inout) :: line

    logical :: eof

    call aq_rxn_file_read_line(file, line, eof)
    if (eof) then
       call aq_rxn_file_die_msg(604113375, file, 'unexpected end of file')
    end if

  end subroutine aq_rxn_file_read_line_no_eof

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Count the number of reactions in an aqueous reaction file
  subroutine aq_rxn_file_count_rxns(file, n_rxn)

    !> Aqueous reaction file.
    type(aq_rxn_file_t), intent(inout) :: file
    !> Number of reactions in file
    integer, intent(out) :: n_rxn

    logical :: eof
    type(aq_rxn_line_t) :: temp_line

    ! Start at the beginning of the file
    rewind(file%unit)
    file%line_num = 0

    ! read file, working out how many reactions we have
    n_rxn = 0
    eof = .false.
    call aq_rxn_line_allocate(temp_line)
    do while (.not. eof)
        call aq_rxn_file_read_line(file, temp_line, eof)
        if (.not.eof) then
            if (temp_line%data(1)=="CLASS:") then
                n_rxn = n_rxn + 1
            endif
        endif
    end do
    call aq_rxn_line_deallocate(temp_line)

    ! Reset file position
    rewind(file%unit)
    file%line_num = 0

  end subroutine aq_rxn_file_count_rxns

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Count the number of reactants and products in a rxn_line
  subroutine aq_rxn_file_count_spec(rxn_line, n_reactant, n_product)

    !> Line containing reaction equation
    type(aq_rxn_line_t), intent(in) :: rxn_line
    !> Number of reactants
    integer, intent(out) :: n_reactant
    !> Number of products
    integer, intent(out) :: n_product
    !> Equation section: 1 for reactant; 2 for product
    integer :: react_or_prod

    integer :: i

    react_or_prod = 1
    n_reactant = 1
    n_product = 1

    do i=1,size(rxn_line%data)

        if (rxn_line%data(i) .eq. "+" .or. rxn_line%data(i) .eq. "-") then
            if (react_or_prod .eq. 1) then
                n_reactant = n_reactant + 1
            else
                n_product = n_product + 1
            endif
        endif

        if (rxn_line%data(i) .eq. "=") then
            react_or_prod = 2
        endif

    enddo

  end subroutine aq_rxn_file_count_spec

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read in a single reaction from the aqueous reaction file
  subroutine aq_rxn_file_read_rxn(file, rxn, aq_spec_data)

    !> Aqueous reaction file.
    type(aq_rxn_file_t), intent(inout) :: file
    !> Reaction data.
    type(aq_rxn_data_t), intent(inout) :: rxn
    !> Mechanism Species
    type(aq_spec_data_t), intent(inout) :: aq_spec_data

    ! Line holding class type
    type(aq_rxn_line_t) :: class_line
    ! Line holding reaction equation
    type(aq_rxn_line_t) :: rxn_line
    ! Line holding rate info
    type(aq_rxn_line_t) :: rate_line

    ! Number of reactants
    integer :: n_reactant
    ! Number of products
    integer :: n_product
    ! Number of rate parameters
    integer :: n_rate_param
    ! Class index
    integer :: class_index
    ! Rate type index
    integer :: rate_index
    ! Species index
    integer :: temp_spec_index
    ! line element index
    integer :: temp_elem_index

    integer :: i, curr_prod
    real :: yield_sign

    ! Prepare the line variables
    call aq_rxn_line_allocate(class_line)
    call aq_rxn_line_allocate(rxn_line)
    call aq_rxn_line_allocate(rate_line)

    ! Read in reaction data from input file
    call aq_rxn_file_read_line_no_eof(file, class_line)
    call aq_rxn_file_read_line_no_eof(file, rxn_line)
    call aq_rxn_file_read_line_no_eof(file, rate_line)

    ! Make sure the CLASS: keyword starts the first line
    if (class_line%data(1).ne."CLASS:") then
       call die_msg(604432708, "expected keyword CLASS: in " // trim(file%name) &
            // " at line " // trim(integer_to_string(file%line_num-2)))
    endif

    ! Get the number of reactants and products
    call aq_rxn_file_count_spec(rxn_line, n_reactant, n_product)

    ! Determine the reaction type
    class_index = aq_rxn_data_get_class_index(class_line%data(2))
    if (class_index .eq. 0) then
        call die_msg(604550357, "found unknown class type: " // trim(class_line%data(2)) &
            // " in " // trim(file%name) // " at line " // &
            trim(integer_to_string(file%line_num-2)))
    endif

    ! Determine the rate type
    rate_index = aq_rxn_data_get_rate_index(rate_line%data(1))
    if (rate_index .eq. 0) then
        call die_msg(604668006, "found unknown rate type: " // trim(rate_line%data(1)) &
            // " in " // trim(file%name) // " at line " // &
            trim(integer_to_string(file%line_num)))
    endif

    ! Determine number of rate parameters
    n_rate_param = aq_rxn_data_get_num_param(rate_index)
    if (n_rate_param .eq. 0) then
        call die_msg(604768848, "Internal error finding number of rate parameters for rate type: " &
            // trim(rate_line%data(1)) // " in " // trim(file%name) // " at line " // &
            trim(integer_to_string(file%line_num)))
    endif

    ! Allocate space for the reaction
    call aq_rxn_data_allocate_size(rxn, n_reactant, n_product, n_rate_param)

    ! Set reaction data
    rxn%class_index = class_index
    rxn%rate_index = rate_index

    ! ... reactants
    do i=1,n_reactant
        temp_elem_index = 1 + 2*(i-1)
        temp_spec_index = aq_spec_data_spec_by_name(aq_spec_data, rxn_line%data(temp_elem_index))
        if (temp_spec_index .eq. 0) then
            temp_spec_index = aq_spec_data_add_spec(aq_spec_data, rxn_line%data(temp_elem_index))
            !> In the capram input files, constant species are surrounded by []
            if (aq_spec_data%name(temp_spec_index)(1:1).eq.'[') then
                aq_spec_data%const_conc(temp_spec_index) = .true.
            endif
        endif
        rxn%reactant(i) = temp_spec_index
    enddo

    ! ... products
    curr_prod = 1
    i = 2*(n_reactant-1)+3
    do while (i.le.size(rxn_line%data))

        ! Check if the first product has a negative sign
        if (curr_prod .eq. 1 .and. rxn_line%data(i) .eq. "-") then
            i = i + 1
        endif

        ! Get the sign of the product yield
        if (rxn_line%data(i-1) .eq. "-") then
            yield_sign = -1.0
        else
            yield_sign = 1.0
        endif

        ! Find out if a yield is included
        if (i+1 .gt. size(rxn_line%data)) then
            ! No yield
            rxn%prod_yield(curr_prod) = 1.0
        else if (rxn_line%data(i+1) .eq. "+" .or. rxn_line%data(i+1) .eq. "-") then
            ! No yield
            rxn%prod_yield(curr_prod) = 1.0
        else if (i+2 .gt. size(rxn_line%data)) then
            ! Yes yield
            if (AQ_RXN_DATA_RATE_YIELD_OK(rate_index).eq.0) then
                read(rxn_line%data(i),*) rxn%prod_yield(curr_prod)
                i = i+1
            else
                call die_msg(605020953, "Yield included in forward/backward rxn in " &
                    // trim(file%name) // " at line " // trim(integer_to_string(file%line_num)))
            endif
        else if (rxn_line%data(i+2) .eq. "+" .or. rxn_line%data(i+2) .eq. "-") then
            ! Yes yield
            if (AQ_RXN_DATA_RATE_YIELD_OK(rate_index).eq.0) then
                read(rxn_line%data(i),*) rxn%prod_yield(curr_prod)
                i = i+1
            else
                call die_msg(605155409, "Yield included in forward/backward rxn in " &
                    // trim(file%name) // " at line " // trim(integer_to_string(file%line_num)))
            endif
        else
            call die_msg(605256251, "Too many elements in product term in " &
                // trim(file%name) // " at line " // trim(integer_to_string(file%line_num)))
        endif

        rxn%prod_yield(curr_prod) = rxn%prod_yield(curr_prod) * yield_sign

        temp_spec_index = aq_spec_data_spec_by_name(aq_spec_data, rxn_line%data(i))
        if (temp_spec_index .eq. 0) then
            temp_spec_index = aq_spec_data_add_spec(aq_spec_data, rxn_line%data(i))
            ! In the capram input files, constant species are surrounded by []
            if (aq_spec_data%name(temp_spec_index)(1:1).eq.'[') then
                aq_spec_data%const_conc(temp_spec_index) = .true.
            endif
        endif
        rxn%product(curr_prod) = temp_spec_index
        i = i+2
        curr_prod = curr_prod + 1
    enddo

    ! ... rate parameters
    do i=1,n_rate_param
        temp_elem_index = 1 + 2*i
        read(rate_line%data(temp_elem_index),*) rxn%rate_param(i)
    enddo

  end subroutine aq_rxn_file_read_rxn

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


end module pmc_aq_rxn_file