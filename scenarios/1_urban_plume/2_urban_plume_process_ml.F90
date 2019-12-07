! Copyright (C) 2009-2013 Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The process program.

!> Read NetCDF output files and process them.
program process

  use pmc_output
  use pmc_stats

  character(len=PMC_MAX_FILENAME_LEN), parameter :: prefix &
       = "out/urban_plume"

  character(len=PMC_MAX_FILENAME_LEN) :: in_filename, out_filename
  type(bin_grid_t) :: diam_grid, avg_bin_grid
  type(aero_data_t) :: aero_data
  type(aero_state_t) :: aero_state, aero_state_averaged
  type(env_state_t) :: env_state
  integer :: ncid, index, repeat, i_index, i_repeat, n_index, n_repeat
  real(kind=dp) :: time, del_t, tot_num_conc, tot_mass_conc

  character(len=PMC_UUID_LEN) :: uuid
  real(kind=dp), allocatable :: times(:), dry_diameters(:), wet_diameters(:), &
       num_concs(:), masses(:), dry_masses(:), soa_masses(:), oc_masses(:), bc_masses(:), &
       so4_masses(:), no3_masses(:), &
       num_dist(:), num_concs_averaged(:), masses_averaged(:), dry_masses_averaged(:)
  real(kind=dp) :: total_so4, total_no3
  type(stats_1d_t) :: stats_num_dist, stats_entropy_dist, stats_tot_num_conc, &
       stats_tot_mass_conc, stats_tot_entropy, stats_tot_entropy_averaged
 
  call pmc_mpi_init()

  call bin_grid_allocate(diam_grid)
  call aero_data_allocate(aero_data)
  call aero_state_allocate(aero_state)
  call aero_state_allocate(aero_state_averaged)
  call bin_grid_allocate(avg_bin_grid)
  call env_state_allocate(env_state)

  call input_n_files(prefix, n_repeat, n_index)

  call bin_grid_make(diam_grid, BIN_GRID_TYPE_LOG, 180, 1d-9, 1d-3)
  call bin_grid_make(avg_bin_grid, BIN_GRID_TYPE_LOG, 1, 1d-30, 1d10)

  allocate(times(n_index))

  do i_index = 1,n_index
     do i_repeat = 1,n_repeat
        call make_filename(in_filename, prefix, ".nc", i_index, i_repeat)
        write(*,*) "Processing " // trim(in_filename)
        ! FIXME: add UUID check into input_state(), keyed off of index or
        ! time or something? Or init to "" and check if not this.
        call input_state(in_filename, index, time, del_t, repeat, &
             uuid, aero_data=aero_data, aero_state=aero_state, &
             env_state=env_state)

        times(i_index) = time

        dry_diameters = aero_state_dry_diameters(aero_state, aero_data)
        wet_diameters = aero_state_diameters(aero_state)
        temperature = env_state%temp
                  
        num_concs = aero_state_num_concs(aero_state)
        num_dist = bin_grid_histogram_1d(diam_grid, dry_diameters, num_concs)
        call stats_1d_add(stats_num_dist, num_dist)

        tot_num_conc = sum(num_concs)
        call stats_1d_add_entry(stats_tot_num_conc, tot_num_conc, i_index)

        masses = aero_state_masses(aero_state, aero_data)
        tot_mass_conc = sum(masses * num_concs)
        call stats_1d_add_entry(stats_tot_mass_conc, tot_mass_conc, i_index)

        dry_masses = aero_state_masses(aero_state, aero_data, &
             exclude=(/"H2O"/))
        bc_masses = aero_state_masses(aero_state, aero_data, &
             include=(/"BC"/))
        so4_masses = aero_state_masses(aero_state, aero_data, &
             include=(/"SO4"/))
        no3_masses = aero_state_masses(aero_state, aero_data, &
             include=(/"NO3"/))

        soa_masses = aero_state_masses(aero_state, aero_data, &
             include=(/"ARO1", "ARO2", "ALK1", "OLE1", "API1", "API2", "LIM1", "LIM2"/))
        oc_masses = aero_state_masses(aero_state, aero_data, &
             include=(/"OC"/))

                
        total_so4 = sum(so4_masses * num_concs)
        total_no3 = sum(no3_masses * num_concs)

        call aero_state_copy(aero_state, aero_state_averaged)
        call aero_state_bin_average_comp(aero_state_averaged, avg_bin_grid, &
             aero_data)
        num_concs_averaged = aero_state_num_concs(aero_state_averaged)
        masses_averaged = aero_state_masses(aero_state_averaged, aero_data)

        dry_masses_averaged = aero_state_masses(aero_state_averaged, &
             aero_data, exclude=(/"H2O"/))

     end do

     call make_filename(out_filename, prefix, "_process.nc", index)
     call pmc_nc_open_write(out_filename, ncid)
     call pmc_nc_write_info(ncid, uuid, "1_urban_plume process")
     call bin_grid_output_netcdf(diam_grid, ncid, "diam", unit="m")

     call stats_1d_output_netcdf(stats_num_dist, ncid, "num_dist", &
          dim_name="diam", unit="m^{-3}")
     call stats_1d_clear(stats_num_dist)

     call pmc_nc_close(ncid)
  end do

  call make_filename(out_filename, prefix, "_process.nc")
  call pmc_nc_open_write(out_filename, ncid)
  call pmc_nc_write_info(ncid, uuid, "1_urban_plume process")
  call pmc_nc_write_real_1d(ncid, times, "time", dim_name="time", unit="s")
  call stats_1d_output_netcdf(stats_tot_num_conc, ncid, "tot_num_conc", &
       dim_name="time", unit="m^{-3}")
  call stats_1d_output_netcdf(stats_tot_mass_conc, ncid, "tot_mass_conc", &
       dim_name="time", unit="m^{-3}")
  call pmc_nc_close(ncid)

  call bin_grid_deallocate(diam_grid)
  call aero_data_deallocate(aero_data)
  call aero_state_deallocate(aero_state)
  call env_state_deallocate(env_state)

  call pmc_mpi_finalize()

end program process
