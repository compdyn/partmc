! Copyright (C) 2007, 2008 Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_process_spec module.

!> The process_spec_t structure and associated subroutines.
module pmc_process_spec

  use pmc_inout
  use pmc_bin_grid
  use pmc_aero_data
  use pmc_aero_state
  use pmc_gas_data
  use pmc_gas_state
  use pmc_env_state
  use pmc_util

  !> Maximum length of process_spec%%type.
  integer, parameter :: PROCESS_SPEC_TYPE_LEN = 40
  !> Maximum length of process_spec%%name.
  integer, parameter :: PROCESS_SPEC_NAME_LEN = 200

  !> Specification of the processing required to turn internal
  !> particle data into output data.
  !!
  !! Internal per-particle data in an aero_state_t is binned into
  !! arrays before being output. The output routines in
  !! output_processed.f90 use the process_spec_t structure to define
  !! how an aero_state_t is transformed into one output data array.
  !!
  !! There are different types of processing that can be done, as
  !! stored in process_spec%%type. We could have had different types
  !! of process_spec_t for each different \c type, storing only the
  !! data appropriate to that type of processing. Instead we have only
  !! a single process_spec_t that is the union of all the possible
  !! necessary data fields. This means that for some types of
  !! processing some of the entries of the process_spec_t structure
  !! will be ignored.
  type process_spec_t
     !> Processing type.
     character(len=PROCESS_SPEC_TYPE_LEN) :: type
     !> Output variable name.
     character(len=PROCESS_SPEC_NAME_LEN) :: name
     !> Number of steps for histogram.
     integer :: n_step
     !> Minimum histogram value.
     real*8 :: min_val
     !> Maximum histogram value.
     real*8 :: max_val
     !> Use a log-scale for histogram?.
     logical :: log_scale
     !> Composition A species.
     character(len=AERO_NAME_LEN), pointer :: a_species(:)
     !> Composition B species.
     character(len=AERO_NAME_LEN), pointer :: b_species(:)
  end type process_spec_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Allocate and initialize storage.
  subroutine process_spec_alloc(process_spec)

    !> Structure to alloc.
    type(process_spec_t), intent(out) :: process_spec

    process_spec%type = "none"
    process_spec%name = "none"
    process_spec%n_step = 0
    process_spec%min_val = 0d0
    process_spec%max_val = 0d0
    process_spec%log_scale = .false.
    allocate(process_spec%a_species(0))
    allocate(process_spec%b_species(0))

  end subroutine process_spec_alloc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Free storage.
  subroutine process_spec_free(process_spec)

    !> Structure to free.
    type(process_spec_t), intent(inout) :: process_spec

    deallocate(process_spec%a_species)
    deallocate(process_spec%b_species)

  end subroutine process_spec_free

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Free storage.
  subroutine process_spec_list_free(process_spec_list)

    !> Structure to free.
    type(process_spec_t), pointer :: process_spec_list(:)

    integer :: i

    do i = 1,size(process_spec_list)
       call process_spec_free(process_spec_list(i))
    end do
    deallocate(process_spec_list)

  end subroutine process_spec_list_free

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Copy all data.
  subroutine process_spec_copy(process_spec_from, process_spec_to)

    !> Source of data.
    type(process_spec_t), intent(in) :: process_spec_from
    !> Dest of data.
    type(process_spec_t), intent(inout) :: process_spec_to

    process_spec_to%type = process_spec_from%type
    process_spec_to%name = process_spec_from%name
    process_spec_to%n_step = process_spec_from%n_step
    process_spec_to%min_val = process_spec_from%min_val
    process_spec_to%max_val = process_spec_from%max_val
    process_spec_to%log_scale = process_spec_from%log_scale
    allocate(process_spec_to%a_species(size(process_spec_from%a_species)))
    allocate(process_spec_to%b_species(size(process_spec_from%b_species)))
    process_spec_to%a_species = process_spec_from%a_species
    process_spec_to%b_species = process_spec_from%b_species

  end subroutine process_spec_copy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read process_spec_list from filename on line in file.
  subroutine spec_read_process_spec_list_filename(file, name, &
       process_spec_list)

    !> Inout file.
    type(inout_file_t), intent(inout) :: file
    !> Name of data line for filename.
    character(len=*), intent(in) :: name
    !> Data to read.
    type(process_spec_t), pointer :: process_spec_list(:)

    character(len=MAX_CHAR_LEN) :: read_name
    type(inout_file_t) :: read_file

    call inout_read_string(file, name, read_name)
    call inout_open_read(read_name, read_file)
    call inout_read_process_spec_list(read_file, process_spec_list)
    call inout_close(read_file)

  end subroutine spec_read_process_spec_list_filename
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read a list of process_specs from the given file.
  subroutine inout_read_process_spec_list(file, process_spec_list)

    !> File to read.
    type(inout_file_t), intent(inout) :: file
    !> Data to read.
    type(process_spec_t), pointer :: process_spec_list(:)

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
          call process_spec_copy(process_spec_list(i), &
               new_process_spec_list(i))
          call process_spec_free(process_spec_list(i))
       end do
       deallocate(process_spec_list)
       process_spec_list => new_process_spec_list
       nullify(new_process_spec_list)
       call inout_read_process_spec(file, process_spec, eof)
    end do
    
  end subroutine inout_read_process_spec_list

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read from an inout_file.
  subroutine inout_read_process_spec(file, process_spec, eof)

    !> Inout file.
    type(inout_file_t), intent(inout) :: file
    !> Data to read.
    type(process_spec_t), intent(out) :: process_spec
    !> If eof instead of reading data.
    logical :: eof

    type(inout_line_t) :: line

    call inout_read_line(file, line, eof)
    if (.not. eof) then
       call inout_check_line_name(file, line, "process")
       call inout_check_line_length(file, line, 1)
       call process_spec_alloc(process_spec)
       process_spec%type = line%data(1)(1:PROCESS_SPEC_TYPE_LEN)
       call inout_read_string(file, "name", process_spec%name)
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

  !> Read kappa spec from an inout_file.
  subroutine inout_read_process_spec_kappa(file, process_spec)

    !> Inout file.
    type(inout_file_t), intent(inout) :: file
    !> Data to read.
    type(process_spec_t), intent(out) :: process_spec

    call inout_read_integer(file, "n_step", process_spec%n_step)
    call inout_read_real(file, "min", process_spec%min_val)
    call inout_read_real(file, "max", process_spec%max_val)
    process_spec%log_scale = .true.

  end subroutine inout_read_process_spec_kappa

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read comp spec from an inout_file.
  subroutine inout_read_process_spec_comp(file, process_spec)

    !> Inout file.
    type(inout_file_t), intent(inout) :: file
    !> Data to read.
    type(process_spec_t), intent(out) :: process_spec

    type(inout_line_t) :: line
    integer :: i

    call inout_read_integer(file, "n_step", process_spec%n_step)
    call inout_read_real(file, "min", process_spec%min_val)
    call inout_read_real(file, "max", process_spec%max_val)

    call inout_read_line_no_eof(file, line)
    call inout_check_line_name(file, line, "a_species")
    deallocate(process_spec%a_species)
    allocate(process_spec%a_species(size(line%data)))
    do i = 1,size(line%data)
       process_spec%a_species(i) = line%data(i)(1:AERO_NAME_LEN)
    end do
    call inout_line_free(line)

    call inout_read_line_no_eof(file, line)
    call inout_check_line_name(file, line, "b_species")
    deallocate(process_spec%b_species)
    allocate(process_spec%b_species(size(line%data)))
    do i = 1,size(line%data)
       process_spec%b_species(i) = line%data(i)(1:AERO_NAME_LEN)
    end do
    call inout_line_free(line)

  end subroutine inout_read_process_spec_comp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read n_orig_part spec from an inout_file.
  subroutine inout_read_process_spec_n_orig_part(file, process_spec)

    !> Inout file.
    type(inout_file_t), intent(inout) :: file
    !> Data to read.
    type(process_spec_t), intent(out) :: process_spec

    call inout_read_real(file, "min", process_spec%min_val)
    call inout_read_real(file, "max", process_spec%max_val)
    process_spec%n_step = nint(process_spec%max_val - process_spec%min_val)

  end subroutine inout_read_process_spec_n_orig_part

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read optic spec from an inout_file.
  subroutine inout_read_process_spec_optic(file, process_spec)

    !> Inout file.
    type(inout_file_t), intent(inout) :: file
    !> Data to read.
    type(process_spec_t), intent(out) :: process_spec

    call inout_read_integer(file, "n_step", process_spec%n_step)
    call inout_read_real(file, "min", process_spec%min_val)
    call inout_read_real(file, "max", process_spec%max_val)
    process_spec%log_scale = .true.

  end subroutine inout_read_process_spec_optic

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
end module pmc_process_spec
