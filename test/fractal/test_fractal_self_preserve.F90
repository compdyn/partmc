! Copyright (C) 2009-2012 Matthew West
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
  type(env_state_t) :: env_state
  integer :: index, i_repeat, i_part, i_spec, out_unit
  integer :: i_file, n_file, i_bin, n_bin
  real(kind=dp) :: time, del_t
  type(aero_particle_t), pointer :: aero_particle
  real(kind=dp) :: dimless_vol_min, dimless_vol_max, log_width
  character(len=PMC_UUID_LEN) :: uuid, run_uuid
  real(kind=dp), allocatable :: diameters(:), num_concs(:), vol_concs(:), hist(:)
  real(kind=dp), allocatable :: dimless_vol(:), dimless_num_conc(:)
  real(kind=dp), allocatable :: aero_dist(:,:)
  real(kind=dp) :: total_num_conc, total_vol_conc
  type(option_s) :: opts(5)

  call pmc_mpi_init()

  opts(1) = option_s("help", .false., 'h')
  opts(2) = option_s("dimless_vol_min", .true., 'N')
  opts(3) = option_s("dimless_vol_max", .true., 'X')
  opts(4) = option_s("nbin", .true., 'b')
  opts(5) = option_s("output", .true., 'o')

  dimless_vol_min = 1d-3
  dimless_vol_max = 10d0
  n_bin = 100
  out_filename = ""

  do
     select case(getopt("hN:X:b:o:", opts))
     case(char(0))
        exit
     case('h')
        call print_help()
        stop
     case('N')
        dimless_vol_min = string_to_real(optarg)
     case('X')
        dimless_vol_max = string_to_real(optarg)
     case('b')
        n_bin = string_to_integer(optarg)
     case('o')
        out_filename = optarg
     case( '?' )
        call print_help()
        call die_msg(956456220, 'unknown option: ' // trim(optopt))
     case default
        call print_help()
        call die_msg(203991511, 'unhandled option: ' // trim(optopt))
     end select
  end do

  if (optind /= command_argument_count()) then
     call print_help()
     call die_msg(533171694, 'expected exactly one non-option prefix argument')
  end if

  call get_command_argument(optind, in_prefix)

  if (out_filename == "") then
     out_filename = trim(in_prefix) // "_self_preserve.txt"
  end if

  call bin_grid_allocate(dimless_vol_grid)
  call aero_data_allocate(aero_data)
  call aero_state_allocate(aero_state)
  call env_state_allocate(env_state)

  allocate(filename_list(0))
  call input_filename_list(in_prefix, filename_list)
  n_file = size(filename_list)
  call assert_msg(554271458, n_file > 0, &
       "no NetCDF files found with prefix: " // trim(in_prefix))

  call bin_grid_make(dimless_vol_grid, n_bin, dimless_vol_min, dimless_vol_max)
  allocate(aero_dist(n_bin, n_file))
  allocate(hist(n_bin))
  allocate(diameters(0))
  allocate(num_concs(0))
  allocate(vol_concs(0))
  allocate(dimless_vol(0))
  allocate(dimless_num_conc(0))

  call input_state(filename_list(n_file), index, time, del_t, i_repeat, &
       uuid, aero_data=aero_data, aero_state=aero_state, env_state=env_state)

  run_uuid = uuid
  call assert_msg(657993562, uuid == run_uuid, &
       "UUID mismatch between " // trim(filename_list(1)) // " and " &
       // trim(filename_list(n_file)))

  call aero_state_num_concs(aero_state, aero_data, env_state, num_concs)
  total_num_conc = aero_state_total_num_conc(aero_state, aero_data, env_state)
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
  call bin_grid_histogram_1d_self_preserve(dimless_vol_grid, dimless_vol, &
       dimless_num_conc, hist, total_vol_conc, total_num_conc)
  aero_dist(:, n_file) = hist

  write(*,'(a)') "Output file: " // trim(out_filename)
  write(*,'(a)') "  Each row of output is one size bin."
  write(*,'(a)') "  The columns of output are:"
  write(*,'(a)') "    column 1: dimentionless particle volume"
  write(*,'(a)') "    column 2: dimentionless number concentration"
  write(*,'(a)') "  Dimensionless volume bins have logarithmic width:"
  write(*,'(a,e20.10)') "    log_width = ln(dimless_vol(i+1)) - ln(dimless_vol(i)) =", &
       dimless_vol_grid%log_width

  call open_file_write(out_filename, out_unit)
  do i_bin = 1,n_bin
     write(out_unit, '(e30.15e3)', advance='no') &
          dimless_vol_grid%center_radius(i_bin)
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
  call bin_grid_allocate(dimless_vol_grid)
  call aero_data_deallocate(aero_data)
  call aero_state_deallocate(aero_state)

  call pmc_mpi_finalize()

contains

  !> Make a histogram with of the given weighted data, scaled by the
  !> non-logarithmic bin sizes, specially for self preserving
  !> size distributions.
  subroutine bin_grid_histogram_1d_self_preserve(x_bin_grid, x_data, &
       weight_data, hist, c_v, c_n)

    !> x-axis bin grid.
    type(bin_grid_t), intent(in) :: x_bin_grid
    !> Data values on the x-axis.
    real(kind=dp), intent(in) :: x_data(:)
    !> Data value weights.
    real(kind=dp), intent(in) :: weight_data(size(x_data))
    !> Histogram to compute.
    real(kind=dp), intent(out) :: hist(x_bin_grid%n_bin)
    !> Total volume and number concentration.
    real(kind=dp), intent(in) :: c_v, c_n

    integer :: i_data, i_bin

    hist = 0d0
    do i_data = 1,size(x_data)
       i_bin = bin_grid_particle_in_bin(x_bin_grid, x_data(i_data))
       if ((i_bin >= 1) .and. (i_bin <= x_bin_grid%n_bin)) then
          hist(i_bin) = hist(i_bin) &
               + weight_data(i_data) / (x_bin_grid%edge_radius(i_bin + 1) &
               * c_v / c_n - x_bin_grid%edge_radius(i_bin) * c_v / c_n)
       end if
    end do

  end subroutine bin_grid_histogram_1d_self_preserve

  subroutine print_help()

    write(*,'(a)') 'Usage: test_fractal_self_preserve [options] <netcdf_prefix>'
    write(*,'(a)') ''
    write(*,'(a)') 'options are:'
    write(*,'(a)') '  -h, --help        Print this help message.'
    write(*,'(a)') '  -N, --dimless_vol_min <V>    Minimum dimensionless volume.'
    write(*,'(a)') '  -X, --dimless_vol_max <V>    Maximum dimensionless volume.'
    write(*,'(a)') '  -b, --nbin <N>    Number of size bins.'
    write(*,'(a)') '  -o, --out <file>  Output filename.'
    write(*,'(a)') ''
    write(*,'(a)') 'Examples:'
    write(*,'(a)') '  test_fractal_self_preserve ' &
         // '--dimless_vol_min 1e-3 --dimless_vol_max 10 --nbin 100 data_0001'
    write(*,'(a)') ''

  end subroutine print_help

end program test_fractal_self_preserve