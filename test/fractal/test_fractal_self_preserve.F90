! Copyright (C) 2009-2012, 2017 Matthew West
! Copyright (C) 2012 Jian Tian
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The test_fractal_self_preserve program.

!> Read NetCDF output files and write out the self-preserving
!> size distributions in text format.
program test_fractal_self_preserve

  use pmc_aero_state
  use pmc_aero_particle
  use pmc_output
  use pmc_mpi
  use pmc_util
  use pmc_constants
  use getopt_m

  character(len=PMC_MAX_FILENAME_LEN) :: in_prefix, out_filename
  character(len=PMC_MAX_FILENAME_LEN), allocatable :: filename_list(:)
  character(len=1000) :: tmp_str
  type(bin_grid_t) :: dimless_vol_grid
  type(aero_data_t) :: aero_data
  type(aero_state_t) :: aero_state
  integer :: index, i_repeat, i_part, i_spec, out_unit
  integer :: i_file, n_file, i_bin, n_bin
  real(kind=dp) :: time, del_t
  type(aero_particle_t), pointer :: aero_particle
  real(kind=dp) :: dimless_vol_min, dimless_vol_max, log_width
  character(len=PMC_UUID_LEN) :: uuid, run_uuid
  real(kind=dp), allocatable :: diameters(:), num_concs(:)
  real(kind=dp), allocatable :: vol_concs(:), hist(:)
  real(kind=dp), allocatable :: dimless_vol(:), dimless_num_conc(:)
  real(kind=dp), allocatable :: aero_dist(:,:)
  real(kind=dp) :: total_num_conc, total_vol_conc
  type(option_s) :: opts(1)

  call pmc_mpi_init()

  opts(1) = option_s("output", .true., 'o')

  dimless_vol_min = 1d-3
  dimless_vol_max = 10d0
  n_bin = 100
  out_filename = ""

  do
     select case(getopt("o", opts))
     case(char(0))
        exit
     case('o')
        out_filename = optarg
     end select
  end do

  if (optind /= command_argument_count()) then
     call die_msg(927637140, 'expected exactly one non-option prefix argument')
  end if

  call get_command_argument(optind, in_prefix)

  if (out_filename == "") then
     out_filename = trim(in_prefix) // "_self_preserve.txt"
  end if

  allocate(filename_list(0))
  call input_filename_list(in_prefix, filename_list)
  n_file = size(filename_list)
  call assert_msg(958886039, n_file > 0, &
       "no NetCDF files found with prefix: " // trim(in_prefix))

  call bin_grid_make(dimless_vol_grid, BIN_GRID_TYPE_LOG, n_bin, &
       dimless_vol_min, dimless_vol_max)
  allocate(aero_dist(n_bin, n_file))
  allocate(hist(n_bin))
  allocate(diameters(0))
  allocate(vol_concs(0))
  allocate(dimless_vol(0))
  allocate(dimless_num_conc(0))

  call input_state(filename_list(n_file), index, time, del_t, i_repeat, &
       uuid, aero_data=aero_data, aero_state=aero_state)

  run_uuid = uuid
  call assert_msg(551984571, uuid == run_uuid, &
       "UUID mismatch between " // trim(filename_list(1)) // " and " &
       // trim(filename_list(n_file)))

  num_concs = aero_state_num_concs(aero_state, aero_data)
  total_num_conc = aero_state_total_num_conc(aero_state, aero_data)
  call ensure_real_array_size(dimless_vol, aero_state%apa%n_part)
  call ensure_real_array_size(dimless_num_conc, aero_state%apa%n_part)
  call ensure_real_array_size(vol_concs, aero_state%apa%n_part)
  vol_concs = aero_particle_volume( &
       aero_state%apa%particle(1:aero_state%apa%n_part)) * num_concs
  total_vol_conc = sum(vol_concs)
  dimless_vol = aero_particle_volume( &
       aero_state%apa%particle(1:aero_state%apa%n_part)) * total_num_conc &
       / total_vol_conc
  dimless_num_conc = num_concs * total_vol_conc / total_num_conc**2
  hist = bin_grid_histogram_1d_self_preserve(dimless_vol_grid, dimless_vol, &
       dimless_num_conc, total_vol_conc, total_num_conc)
  aero_dist(:, n_file) = hist

  write(*,'(a)') "Output file: " // trim(out_filename)
  write(*,'(a)') "  Each row of output is one size bin."
  write(*,'(a)') "  The columns of output are:"
  write(*,'(a)') "    column 1: dimentionless particle volume"
  write(*,'(a)') "    column 2: dimentionless number concentration"
  write(*,'(a)') "  Dimensionless volume bins have logarithmic width:"
  write(*,'(a,e20.10)') "    log_width = ln(dimless_vol(i+1)) " &
       // "- ln(dimless_vol(i)) =", dimless_vol_grid%widths(1)

  call open_file_write(out_filename, out_unit)
  do i_bin = 1,n_bin
     write(out_unit, '(e30.15e3)', advance='no') &
          dimless_vol_grid%centers(i_bin)
     write(out_unit, '(e30.15e3)', advance='no') aero_dist(i_bin, n_file)
     write(out_unit, '(a)') ''
  end do

  call close_file(out_unit)

  deallocate(filename_list)
  deallocate(aero_dist)
  deallocate(hist)
  deallocate(diameters)
  deallocate(num_concs)
  deallocate(vol_concs)
  deallocate(dimless_vol)
  deallocate(dimless_num_conc)

  call pmc_mpi_finalize()

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Make a histogram with of the given weighted data, scaled by the
  !> non-logarithmic bin sizes, specially for self preserving
  !> size distributions.
  function bin_grid_histogram_1d_self_preserve(x_bin_grid, x_data, &
       weight_data, c_v, c_n)

    !> x-axis bin grid.
    type(bin_grid_t), intent(in) :: x_bin_grid
    !> Data values on the x-axis.
    real(kind=dp), intent(in) :: x_data(:)
    !> Data value weights.
    real(kind=dp), intent(in) :: weight_data(size(x_data))
    !> Total volume and number concentration.
    real(kind=dp), intent(in) :: c_v, c_n

    !> Return histogram.
    real(kind=dp) :: &
         bin_grid_histogram_1d_self_preserve(bin_grid_size(x_bin_grid))

    integer :: i_data, x_bin

    bin_grid_histogram_1d_self_preserve = 0d0
    do i_data = 1,size(x_data)
       x_bin = bin_grid_find(x_bin_grid, x_data(i_data))
       if ((x_bin >= 1) .and. (x_bin <= bin_grid_size(x_bin_grid))) then
          bin_grid_histogram_1d_self_preserve(x_bin) = &
               bin_grid_histogram_1d_self_preserve(x_bin) &
               + weight_data(i_data) * c_n / c_v / x_bin_grid%widths(x_bin)
       end if
    end do

  end function bin_grid_histogram_1d_self_preserve

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program test_fractal_self_preserve
