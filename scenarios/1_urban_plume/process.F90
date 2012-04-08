! Copyright (C) 2009-2012 Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The process program.

!> Read NetCDF output files and process them.
program process

  use pmc_output

  character(len=PMC_MAX_FILENAME_LEN), parameter :: prefix &
       = "out/urban_plume_wc_0001"

  integer, parameter :: diam_n_bin = 100
  real(kind=dp), parameter :: diam_min = 1d-9
  real(kind=dp), parameter :: diam_max = 1d-3

  integer, parameter :: bc_n_bin = 50
  real(kind=dp), parameter :: bc_min = 0d0
  real(kind=dp), parameter :: bc_max = 1d0

  character(len=PMC_MAX_FILENAME_LEN), allocatable :: filename_list(:)
  character(len=PMC_MAX_FILENAME_LEN) :: out_filename
  type(bin_grid_t) :: diam_grid, bc_grid
  type(aero_data_t) :: aero_data
  type(aero_state_t) :: aero_state
  integer :: ncid, i_file, index, i_repeat
  real(kind=dp) :: time, del_t
  character(len=PMC_UUID_LEN) :: uuid, run_uuid
  real(kind=dp), allocatable :: dry_diameters(:), num_concs(:), dry_masses(:)
  real(kind=dp), allocatable :: bc_masses(:), bc_fracs(:), num_dist(:)
  real(kind=dp), allocatable :: diam_bc_dist(:,:)

  call pmc_mpi_init()

  call bin_grid_allocate(diam_grid)
  call bin_grid_allocate(bc_grid)
  call aero_data_allocate(aero_data)
  call aero_state_allocate(aero_state)

  call input_filename_list(prefix, filename_list)

  call bin_grid_make(diam_grid, BIN_GRID_TYPE_LOG, diam_n_bin, diam_min, &
       diam_max)
  call bin_grid_make(bc_grid, BIN_GRID_TYPE_LINEAR, bc_n_bin, bc_min, bc_max)

  do i_file = 1,size(filename_list)
     write(*,*) "Processing " // trim(filename_list(i_file))
     call input_state(filename_list(i_file), index, time, del_t, i_repeat, &
          uuid, aero_data=aero_data, aero_state=aero_state)

     ! FIXME: add UUID check into input_state(), keyed off of index or
     ! time or something?

     call aero_state_dry_diameters(aero_state, aero_data, dry_diameters)
     call aero_state_num_concs(aero_state, num_concs)
     call aero_state_masses(aero_state, aero_data, dry_masses, &
          exclude=(/"H2O"/))
     call aero_state_masses(aero_state, aero_data, bc_masses, &
          include=(/"BC"/))
     bc_fracs = bc_masses / dry_masses

     call bin_grid_histogram_1d(diam_grid, dry_diameters, num_concs, num_dist)
     call bin_grid_histogram_2d(diam_grid, dry_diameters, bc_grid, bc_fracs, &
          num_concs, diam_bc_dist)

     call make_filename(out_filename, prefix, "_process.nc", index)
     call pmc_nc_open_write(out_filename, ncid)
     call pmc_nc_write_info(ncid, uuid, "1_urban_plume process")
     call bin_grid_output_netcdf(diam_grid, ncid, "diam", unit="m")
     call bin_grid_output_netcdf(bc_grid, ncid, "bc_frac", unit="1")
     call pmc_nc_write_real_1d(ncid, num_dist, "num_dist", &
          dim_name="diam", unit="m^{-3}")
     call pmc_nc_write_real_2d(ncid, diam_bc_dist, "diam_bc_dist", &
          dim_name_1="diam", dim_name_2="bc_frac", unit="m^{-3}")
     call pmc_nc_close(ncid)
  end do

  call bin_grid_allocate(diam_grid)
  call aero_data_deallocate(aero_data)
  call aero_state_deallocate(aero_state)

  call pmc_mpi_finalize()

end program process
