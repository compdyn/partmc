!> Read NetCDF output files from a modal runn and process them.

program process

  use pmc_output
  use pmc_stats
  use pmc_aero_dist
  use pmc_scenario

  character(len=PMC_MAX_FILENAME_LEN), parameter :: prefix & 
= "out/modal_test"

  character(len=PMC_MAX_FILENAME_LEN) :: in_filename, out_filename
  type(bin_grid_t) :: bin_grid
  type(aero_dist_t) :: aero_dist
  type(aero_binned_t) :: aero_binned
  type(aero_data_t) :: aero_data
  type(env_state_t) :: env_state
  type(gas_data_t) :: gas_data
  type(gas_state_t) :: gas_state
  type(scenario_t) :: scenario
  integer :: ncid, index, i_mode, i_index, dum, n_index
  real(kind=dp) :: time, del_t, tot_num_conc, density, tot_mass_conc
  character(len=PMC_UUID_LEN) :: uuid
  real(kind=dp), allocatable :: times(:), num_conc(:)
  type(stats_1d_t) :: stats_tot_num_conc, stats_num_conc, stats_rates_0, &
                      stats_rates_3, stats_tot_mass_conc
  real(kind=dp), allocatable :: rates_0(:), rates_3(:)
  character(len=PMC_MAX_FILENAME_LEN), allocatable :: file_list(:)

  call pmc_mpi_init()

  call input_filename_list(prefix, file_list)
  n_index = size(file_list)
  
  allocate(times(n_index))

  do i_index = 1,n_index
    call make_filename(in_filename, prefix, ".nc", i_index)
    write(*,*) "Processing " // trim(in_filename)
    call input_modal(in_filename, index, time, del_t, uuid, aero_dist=aero_dist, &
        aero_binned=aero_binned, aero_data=aero_data, env_state=env_state, &
        gas_data=gas_data, gas_state=gas_state, bin_grid=bin_grid, scenario=scenario)
    times(i_index) = time

    density = aero_data%density(1)

    num_conc = aero_binned%num_conc * bin_grid%widths
    tot_num_conc = sum(num_conc)

    tot_mass_conc = sum(aero_binned%vol_conc(:,1) * bin_grid%widths) * aero_data%density(1)

    if (.not. allocated(rates_0)) allocate(rates_0(aero_dist_n_mode(aero_dist)))
    if (.not. allocated(rates_3)) allocate(rates_3(aero_dist_n_mode(aero_dist)))

    call scenario_modal_dry_dep_rates(scenario, aero_dist, 0.0d0, density, &
      env_state, rates_0)
    call stats_1d_add(stats_rates_0, rates_0)
    call scenario_modal_dry_dep_rates(scenario, aero_dist, 3.0d0, density, &
      env_state, rates_3)
    call stats_1d_add(stats_rates_3, rates_3)

    call stats_1d_add(stats_num_conc, num_conc)
    call stats_1d_add_entry(stats_tot_num_conc, tot_num_conc, i_index)

    call stats_1d_add_entry(stats_tot_mass_conc, tot_mass_conc, i_index)

    call make_filename(out_filename, prefix, "_process.nc", index)
    write(*,*) "Writing " // trim(out_filename)
    call pmc_nc_open_write(out_filename, ncid)
    call pmc_nc_write_info(ncid, uuid, "1_urban_plume modal process")
    call env_state_output_netcdf(env_state, ncid)
    call aero_data_output_netcdf(aero_data, ncid)
    call aero_dist_output_netcdf(aero_dist, ncid)
    call aero_binned_output_netcdf(aero_binned, ncid, bin_grid, aero_data)

    call stats_1d_output_netcdf(stats_num_conc, ncid, "num concs. per bin", &
      dim_name="diam", unit="m^{-3}")
    call stats_1d_output_netcdf(stats_rates_0, ncid, "loss_rates_0", &
      dim_name="modes", unit="m s^{-1}")
    call stats_1d_output_netcdf(stats_rates_3, ncid, "loss_rates_3", &
      dim_name="modes", unit="m s^{-1}")

    call aero_binned_zero(aero_binned)

    call stats_1d_clear(stats_num_conc)
    call stats_1d_clear(stats_rates_0)
    call stats_1d_clear(stats_rates_3)

    call pmc_nc_close(ncid) 
  end do

  call make_filename(out_filename, prefix, "_process.nc")
  write(*,*) "Writing " // trim(out_filename)
  call pmc_nc_open_write(out_filename, ncid)
  call pmc_nc_write_info(ncid, uuid, "1_urban_plume modal process")
  call pmc_nc_write_real_1d(ncid, times, "time", dim_name="time", unit="s")
  call stats_1d_output_netcdf(stats_tot_num_conc, ncid, "tot_num_conc", &
      dim_name="time", unit="m^{-3}")
  call stats_1d_output_netcdf(stats_tot_mass_conc, ncid, "tot_mass_conc", &
      dim_name="time", unit="ug m^{-3}")
  call pmc_nc_close(ncid)

  call pmc_mpi_finalize()

end program process