! Copyright (C) 2009-2012 Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The process program.

!> Read NetCDF output files and process them.
program process

  use pmc_output

  character(len=PMC_MAX_FILENAME_LEN), parameter :: prefix &
       = "out/urban_plume_wc"

  integer, parameter :: diam_n_bin = 100
  real(kind=dp), parameter :: diam_min = 1d-9
  real(kind=dp), parameter :: diam_max = 1d-3

  integer, parameter :: bc_n_bin = 50
  real(kind=dp), parameter :: bc_min = 0d0
  real(kind=dp), parameter :: bc_max = 1d0

  character(len=PMC_MAX_FILENAME_LEN) :: in_filename, out_filename
  type(bin_grid_t) :: diam_grid, bc_grid
  type(aero_data_t) :: aero_data
  type(aero_state_t) :: aero_state
  integer :: ncid, index, repeat, i_index, i_repeat, n_index, n_repeat
  real(kind=dp) :: time, del_t, tot_num_conc, tot_mass_conc
  character(len=PMC_UUID_LEN) :: uuid, run_uuid
  real(kind=dp), allocatable :: dry_diameters(:), num_concs(:), dry_masses(:)
  real(kind=dp), allocatable :: masses(:), bc_masses(:), bc_fracs(:)
  real(kind=dp), allocatable :: num_dist(:), num_dist_mean(:)
  real(kind=dp), allocatable :: num_dist_var(:), num_dist_ci_offset(:)
  real(kind=dp), allocatable :: diam_bc_dist(:,:), diam_bc_dist_mean(:,:)
  real(kind=dp), allocatable :: diam_bc_dist_var(:,:)
  real(kind=dp), allocatable :: diam_bc_dist_ci_offset(:,:)
  real(kind=dp), allocatable :: times(:)
  real(kind=dp), allocatable :: tot_num_conc_mean(:), tot_num_conc_var(:)
  real(kind=dp), allocatable :: tot_num_conc_ci_offset(:)
  real(kind=dp), allocatable :: tot_mass_conc_mean(:), tot_mass_conc_var(:)
  real(kind=dp), allocatable :: tot_mass_conc_ci_offset(:)

  call pmc_mpi_init()

  call bin_grid_allocate(diam_grid)
  call bin_grid_allocate(bc_grid)
  call aero_data_allocate(aero_data)
  call aero_state_allocate(aero_state)

  call input_n_files(prefix, n_repeat, n_index)

  call bin_grid_make(diam_grid, BIN_GRID_TYPE_LOG, diam_n_bin, diam_min, &
       diam_max)
  call bin_grid_make(bc_grid, BIN_GRID_TYPE_LINEAR, bc_n_bin, bc_min, bc_max)

  allocate(times(n_index))
  allocate(tot_num_conc_mean(n_index), tot_num_conc_var(n_index))
  allocate(tot_mass_conc_mean(n_index), tot_mass_conc_var(n_index))

  do i_index = 1,n_index
     do i_repeat = 1,n_repeat
        call make_filename(in_filename, prefix, ".nc", i_index, i_repeat)
        write(*,*) "Processing " // trim(in_filename)
        call input_state(in_filename, index, time, del_t, repeat, &
             uuid, aero_data=aero_data, aero_state=aero_state)

        ! FIXME: add UUID check into input_state(), keyed off of index or
        ! time or something?

        call aero_state_dry_diameters(aero_state, aero_data, dry_diameters)
        call aero_state_num_concs(aero_state, num_concs)
        call aero_state_masses(aero_state, aero_data, masses)
        call aero_state_masses(aero_state, aero_data, dry_masses, &
             exclude=(/"H2O"/))
        call aero_state_masses(aero_state, aero_data, bc_masses, &
             include=(/"BC"/))
        bc_fracs = bc_masses / dry_masses
        tot_num_conc = sum(num_concs)
        tot_mass_conc = sum(masses * num_concs)

        times(i_index) = time
        call update_mean_var(tot_num_conc_mean(i_index), &
             tot_num_conc_var(i_index), tot_num_conc, i_repeat)
        call update_mean_var(tot_mass_conc_mean(i_index), &
             tot_mass_conc_var(i_index), tot_mass_conc, i_repeat)
        call bin_grid_histogram_1d(diam_grid, dry_diameters, num_concs, &
             num_dist)
        call update_mean_var_1d(num_dist_mean, num_dist_var, num_dist, &
             i_repeat)
        call bin_grid_histogram_2d(diam_grid, dry_diameters, bc_grid, &
             bc_fracs, num_concs, diam_bc_dist)
        call update_mean_var_2d(diam_bc_dist_mean, diam_bc_dist_var, &
             diam_bc_dist, i_repeat)
     end do

     call conf_95_offset_1d(num_dist_var, n_repeat, num_dist_ci_offset)
     call conf_95_offset_2d(diam_bc_dist_var, n_repeat, diam_bc_dist_ci_offset)

     call make_filename(out_filename, prefix, "_process.nc", index)
     call pmc_nc_open_write(out_filename, ncid)
     call pmc_nc_write_info(ncid, uuid, "1_urban_plume process")
     call bin_grid_output_netcdf(diam_grid, ncid, "diam", unit="m")
     call bin_grid_output_netcdf(bc_grid, ncid, "bc_frac", unit="1")
     call pmc_nc_write_real_1d(ncid, num_dist_mean, "num_dist", &
          dim_name="diam", unit="m^{-3}")
     call pmc_nc_write_real_1d(ncid, num_dist_ci_offset, &
          "num_dist_ci_offset", dim_name="diam", unit="m^{-3}")
     call pmc_nc_write_real_2d(ncid, diam_bc_dist_mean, "diam_bc_dist", &
          dim_name_1="diam", dim_name_2="bc_frac", unit="m^{-3}")
     call pmc_nc_write_real_2d(ncid, diam_bc_dist_ci_offset, &
          "diam_bc_dist_ci_offset", dim_name_1="diam", dim_name_2="bc_frac", &
          unit="m^{-3}")
     call pmc_nc_close(ncid)
  end do

  call conf_95_offset_1d(tot_num_conc_var, n_repeat, tot_num_conc_ci_offset)
  call conf_95_offset_1d(tot_mass_conc_var, n_repeat, tot_mass_conc_ci_offset)

  call make_filename(out_filename, prefix, "_process.nc")
  call pmc_nc_open_write(out_filename, ncid)
  call pmc_nc_write_info(ncid, uuid, "1_urban_plume process")
  call pmc_nc_write_real_1d(ncid, times, "time", dim_name="time", unit="s")
  call pmc_nc_write_real_1d(ncid, tot_num_conc_mean, "tot_num_conc", &
       dim_name="time", unit="m^{-3}")
  call pmc_nc_write_real_1d(ncid, tot_num_conc_ci_offset, &
       "tot_num_conc_ci_offset", dim_name="time", unit="m^{-3}")
  call pmc_nc_write_real_1d(ncid, tot_mass_conc_mean, "tot_mass_conc", &
       dim_name="time", unit="kg m^{-3}")
  call pmc_nc_write_real_1d(ncid, tot_mass_conc_ci_offset, &
       "tot_mass_conc_ci_offset", dim_name="time", unit="kg m^{-3}")
  call pmc_nc_close(ncid)

  call bin_grid_allocate(diam_grid)
  call aero_data_deallocate(aero_data)
  call aero_state_deallocate(aero_state)

  call pmc_mpi_finalize()

end program process
