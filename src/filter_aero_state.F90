! Copyright (C) 2026 Jeffrey Curtis
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The filter_aero_state program.

!> Read a PartMC aero_state NetCDF output, drop a subset of particles
!> based on a user-supplied criterion, and write the result to a new
!> NetCDF file. Two filter criteria are supported (exactly one is
!> required):
!>
!>   --class N        keep only particles whose weight_class is N
!>                    (suitable as a pre-filter for the
!>                    particles_from_file aero_mode type)
!>   --indices FILE   keep only particles whose 1-based index appears
!>                    in FILE (one integer per line; blank lines and
!>                    lines starting with '#' are ignored, with inline
!>                    '#' comments also stripped)
program filter_aero_state

  use pmc_aero_data
  use pmc_aero_state
  use pmc_gas_data
  use pmc_gas_state
  use pmc_env_state
  use pmc_output
  use pmc_util
  use pmc_mpi
  use getopt_m

  integer, parameter :: MODE_CLASS   = 1
  integer, parameter :: MODE_INDICES = 2

  character(len=PMC_MAX_FILENAME_LEN) :: in_filename, out_prefix
  character(len=PMC_MAX_FILENAME_LEN) :: indices_filename
  type(aero_data_t) :: aero_data
  type(aero_state_t) :: aero_state
  type(gas_data_t) :: gas_data
  type(gas_state_t) :: gas_state
  type(env_state_t) :: env_state
  integer :: index, i_repeat, keep_class, n_before, n_after, n_kept
  integer :: filter_mode
  integer, allocatable :: keep_indices(:)
  real(kind=dp) :: time, del_t
  character(len=PMC_UUID_LEN) :: uuid
  type(option_s) :: opts(4)

  call pmc_mpi_init()

  opts(1) = option_s("help",    .false., 'h')
  opts(2) = option_s("class",   .true.,  'c')
  opts(3) = option_s("indices", .true.,  'i')
  opts(4) = option_s("out",     .true.,  'o')

  out_prefix = ""
  indices_filename = ""
  keep_class = -1
  filter_mode = 0

  do
     select case(getopt("hc:i:o:", opts))
     case(char(0))
        exit
     case('h')
        call print_help()
        stop
     case('c')
        call set_filter_mode(filter_mode, MODE_CLASS)
        keep_class = string_to_integer(optarg)
     case('i')
        call set_filter_mode(filter_mode, MODE_INDICES)
        indices_filename = optarg
     case('o')
        out_prefix = optarg
     case( '?' )
        call print_help()
        call die_msg(173820492, 'unknown option: ' // trim(optopt))
     case default
        call print_help()
        call die_msg(641029537, 'unhandled option: ' // trim(optopt))
     end select
  end do

  if (optind /= command_argument_count()) then
     call print_help()
     call die_msg(517284903, &
          'expected exactly one non-option NetCDF filename argument')
  end if
  if (filter_mode == 0) then
     call print_help()
     call die_msg(208461379, &
          'must supply exactly one of --class <N> or --indices <file>')
  end if
  if (filter_mode == MODE_CLASS .and. keep_class < 1) then
     call print_help()
     call die_msg(620038174, '--class <N> must be >= 1')
  end if

  call get_command_argument(optind, in_filename)
  if (out_prefix == "") then
     out_prefix = default_prefix(in_filename, filter_mode, keep_class)
  end if

  call input_state(in_filename, index, time, del_t, i_repeat, uuid, &
       aero_data=aero_data, aero_state=aero_state, &
       gas_data=gas_data, gas_state=gas_state, env_state=env_state)

  n_before = aero_state_n_part(aero_state)

  if (filter_mode == MODE_CLASS) then
     call aero_state_filter_by_weight_class(aero_state, keep_class)
     n_after = aero_state_n_part(aero_state)
     write(*,'(a)') "Input:  " // trim(in_filename)
     write(*,'(a,i0)') "  particles in:  ", n_before
     write(*,'(a,i0,a,i0)') "  kept (class=", keep_class, "): ", n_after
     write(*,'(a,i0)') "  dropped:       ", n_before - n_after
  else
     call read_indices(indices_filename, keep_indices)
     n_kept = size(keep_indices)
     call assert_msg(740583921, n_kept >= 1, &
          'indices file ' // trim(indices_filename) // &
          ' contains no indices')
     call aero_state_filter_by_indices(aero_state, keep_indices)
     n_after = aero_state_n_part(aero_state)
     write(*,'(a)') "Input:   " // trim(in_filename)
     write(*,'(a)') "Indices: " // trim(indices_filename)
     write(*,'(a,i0)') "  particles in:            ", n_before
     write(*,'(a,i0)') "  indices requested:       ", n_kept
     write(*,'(a,i0)') "  particles kept (unique): ", n_after
     write(*,'(a,i0)') "  particles dropped:       ", n_before - n_after
  end if

  if (n_after == 0) then
     call die_msg(391760284, &
          'no particles match the requested filter')
  end if

  write(*,'(a)') "Output prefix: " // trim(out_prefix)
  call output_state(out_prefix, OUTPUT_TYPE_SINGLE, aero_data, aero_state, &
       gas_data, gas_state, env_state, index, time, del_t, i_repeat, &
       record_removals=.false., record_optical=.true., uuid=uuid)

  call pmc_mpi_finalize()

contains

  !> Set the active filter mode, erroring if a different mode was
  !> already set (the user supplied two mutually-exclusive flags).
  subroutine set_filter_mode(current, requested)

    integer, intent(inout) :: current
    integer, intent(in) :: requested

    if (current /= 0 .and. current /= requested) then
       call print_help()
       call die_msg(816274093, &
            '--class and --indices are mutually exclusive')
    end if
    current = requested

  end subroutine set_filter_mode

  !> Read 1-based particle indices from a plain-text file. Blank lines
  !> and lines beginning with '#' are skipped; inline '#' comments on
  !> integer lines are stripped. Two-pass: count valid integers first,
  !> then allocate and parse.
  subroutine read_indices(filename, indices)

    character(len=*), intent(in) :: filename
    integer, allocatable, intent(out) :: indices(:)

    integer :: unit, ios, hash_pos, n, i
    character(len=200) :: line, payload

    do i = 1, 2
       call open_file_read(filename, unit)
       n = 0
       do
          read(unit, '(a)', iostat=ios) line
          if (ios /= 0) exit
          hash_pos = scan(line, '#')
          if (hash_pos > 0) then
             payload = line(1:hash_pos - 1)
          else
             payload = line
          end if
          payload = adjustl(payload)
          if (len_trim(payload) == 0) cycle
          n = n + 1
          if (i == 2) then
             read(payload, *, iostat=ios) indices(n)
             call assert_msg(389417520, ios == 0, &
                  "could not parse integer from line '" &
                  // trim(line) // "' in " // trim(filename))
          end if
       end do
       close(unit)
       if (i == 1) allocate(indices(n))
    end do

  end subroutine read_indices

  !> Build a default output prefix from an input filename of the form
  !> <stem>_NNNN_MMMMMMMM.nc. Strips the trailing _NNNN_MMMMMMMM.nc
  !> and appends a criterion-specific suffix ("_classN" for class
  !> filtering, "_filtered" for index filtering).
  function default_prefix(in_filename, mode, keep_class) result(prefix)

    character(len=*), intent(in) :: in_filename
    integer, intent(in) :: mode
    integer, intent(in) :: keep_class
    character(len=PMC_MAX_FILENAME_LEN) :: prefix

    integer :: dot_idx, last_underscore_idx, second_last_underscore_idx
    character(len=PMC_MAX_FILENAME_LEN) :: stem

    dot_idx = scan(in_filename, '.', back=.true.)
    if (dot_idx > 0) then
       stem = in_filename(1:dot_idx - 1)
    else
       stem = in_filename
    end if

    last_underscore_idx = scan(stem, '_', back=.true.)
    if (last_underscore_idx > 0) then
       stem = stem(1:last_underscore_idx - 1)
       second_last_underscore_idx = scan(stem, '_', back=.true.)
       if (second_last_underscore_idx > 0) then
          stem = stem(1:second_last_underscore_idx - 1)
       end if
    end if

    if (mode == MODE_CLASS) then
       write(prefix, '(a,a,i0)') trim(stem), "_class", keep_class
    else
       prefix = trim(stem) // "_filtered"
    end if

  end function default_prefix

  subroutine print_help()

    write(*,'(a)') 'Usage: filter_aero_state [options] <netcdf_file>'
    write(*,'(a)') ''
    write(*,'(a)') 'Filter a PartMC aero_state NetCDF down to a subset of'
    write(*,'(a)') 'particles. Exactly one filter criterion is required:'
    write(*,'(a)') ''
    write(*,'(a)') '  -c, --class <N>         keep particles with weight_class = N'
    write(*,'(a)') '  -i, --indices <file>    keep particles whose 1-based index'
    write(*,'(a)') '                          appears in <file> (one integer per'
    write(*,'(a)') '                          line; blank / # lines ignored)'
    write(*,'(a)') ''
    write(*,'(a)') 'Other options:'
    write(*,'(a)') '  -h, --help              Print this help message.'
    write(*,'(a)') '  -o, --out <prefix>      Output prefix. Defaults to'
    write(*,'(a)') '                          <stem>_classN (class mode) or'
    write(*,'(a)') '                          <stem>_filtered (indices mode);'
    write(*,'(a)') '                          output_state appends the source'
    write(*,'(a)') '                          file''s i_repeat and index.'
    write(*,'(a)') ''
    write(*,'(a)') 'Examples:'
    write(*,'(a)') '  filter_aero_state --class 2 src_0001_00000005.nc'
    write(*,'(a)') '  filter_aero_state -i keep.txt -o subset src.nc'

  end subroutine print_help

end program filter_aero_state
