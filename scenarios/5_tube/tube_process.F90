! Copyright (C) 2009-2013, 2016, 2017 Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The process program.

!> Read NetCDF output files and process them.
program process

  use pmc_output
  use pmc_stats
  use pmc_aero_state

  character(len=PMC_MAX_FILENAME_LEN), parameter :: prefix &
       = "/Users/nriemer/git/partmc/scenarios/5_tube/out_pfr_suc2as_2/urban_plume"

  character(len=PMC_MAX_FILENAME_LEN) :: in_filename, out_filename
  type(bin_grid_t) :: diam_grid, oc_grid, sc_grid, avg_bin_grid, kappa_grid
  type(aero_data_t) :: aero_data
  type(aero_state_t) :: aero_state, aero_state_averaged
  type(env_state_t) :: env_state
  integer :: ncid, index, repeat, i_index, i_repeat, n_index, n_repeat
  real(kind=dp) :: time, del_t, tot_num_conc, tot_mass_conc
  real(kind=dp) :: d_alpha, d_gamma, chi
  integer       :: n2o5_type 
  character(len=PMC_UUID_LEN) :: uuid
  real(kind=dp), allocatable :: times(:), dry_diameters(:), num_concs(:), &
       dry_masses(:), wet_masses(:), oc_masses(:), oc_fracs(:), &
       so4_masses(:), no3_masses(:), nh4_masses(:), wet_diameters(:), &
       crit_rhs(:), scs(:), num_dist(:), kappa_dist(:), &
       mass_oc_dist(:), mass_so4_dist(:), mass_no3_dist(:), mass_nh4_dist(:), &
       diam_oc_dist(:,:), diam_sc_dist(:,:), &
       surf_area_pr(:), surf_area_comp(:), surf_area_dist_pr(:), &
       surf_area_dist_comp(:), h2o_masses(:), h2o_masses_avg(:), &
       kappas(:), diam_kappa_dist(:,:)
  type(stats_1d_t) :: stats_num_dist, stats_kappa_dist, stats_d_alpha, stats_tot_num_conc, &
       stats_tot_mass_conc, stats_d_gamma, stats_chi, &
       stats_mass_oc_dist, stats_mass_so4_dist, stats_mass_no3_dist, stats_mass_nh4_dist
  type(stats_2d_t) :: stats_diam_oc_dist, stats_diam_sc_dist, stats_diam_kappa_dist

  call pmc_mpi_init()

  call input_n_files(prefix, n_repeat, n_index)

  call bin_grid_make(diam_grid, BIN_GRID_TYPE_LOG, 180, 1d-9, 1d-3)
  call bin_grid_make(avg_bin_grid, BIN_GRID_TYPE_LOG, 1, 1d-30, 1d10)
  call bin_grid_make(oc_grid, BIN_GRID_TYPE_LINEAR, 50, 0d0, 1d0)
  call bin_grid_make(sc_grid, BIN_GRID_TYPE_LOG, 100, 1d-4, 1d0)
  call bin_grid_make(kappa_grid, BIN_GRID_TYPE_LINEAR, 50, 0d0, 1d0)

  allocate(times(n_index))

  scs = [ real(kind=dp) :: ] ! silence compiler warnings
  oc_fracs = [ real(kind=dp) :: ]

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

        !!!!Create an averaged aero_state!!!!
        aero_state_averaged = aero_state
        call aero_state_bin_average_comp(aero_state_averaged, avg_bin_grid, &
            aero_data)

        dry_diameters = aero_state_dry_diameters(aero_state, aero_data)
        wet_diameters = aero_state_diameters(aero_state, aero_data)
        num_concs = aero_state_num_concs(aero_state, aero_data)
        num_dist = bin_grid_histogram_1d(diam_grid, dry_diameters, num_concs)
        call stats_1d_add(stats_num_dist, num_dist)

        tot_num_conc = sum(num_concs)
        call stats_1d_add_entry(stats_tot_num_conc, tot_num_conc, i_index)

        wet_masses = aero_state_masses(aero_state, aero_data)
        tot_mass_conc = sum(wet_masses * num_concs)
        call stats_1d_add_entry(stats_tot_mass_conc, tot_mass_conc, i_index)

        ! Calcluate masses of different species
        oc_masses = aero_state_masses(aero_state, aero_data, &
             include=(/"API1"/))
        so4_masses = aero_state_masses(aero_state, aero_data, &
             include=(/"SO4"/))
        no3_masses = aero_state_masses(aero_state, aero_data, &
             include=(/"NO3"/))
        nh4_masses = aero_state_masses(aero_state, aero_data, &
             include=(/"NH4"/))
        h2o_masses = aero_state_masses(aero_state, aero_data, &
             include=(/"H2O"/))

        h2o_masses_avg = aero_state_masses(aero_state_averaged, aero_data, &
             include=(/"H2O"/))
        
        kappas = aero_state_solute_kappas(aero_state, aero_data)
        kappa_dist = bin_grid_histogram_1d(kappa_grid, kappas, num_concs)
        call stats_1d_add(stats_kappa_dist, kappa_dist)
        
        diam_kappa_dist = bin_grid_histogram_2d(diam_grid, dry_diameters, &
             kappa_grid, kappas, num_concs)
        call stats_2d_add(stats_diam_kappa_dist, diam_kappa_dist)

        ! Make distribution for different species
        mass_oc_dist  = bin_grid_histogram_1d(diam_grid, dry_diameters, oc_masses * num_concs)
        mass_so4_dist = bin_grid_histogram_1d(diam_grid, dry_diameters, so4_masses * num_concs)
        mass_no3_dist = bin_grid_histogram_1d(diam_grid, dry_diameters, no3_masses * num_concs)
        mass_nh4_dist = bin_grid_histogram_1d(diam_grid, dry_diameters, nh4_masses * num_concs)
        
        call stats_1d_add(stats_mass_oc_dist,  mass_oc_dist)
        call stats_1d_add(stats_mass_so4_dist, mass_so4_dist)
        call stats_1d_add(stats_mass_no3_dist, mass_no3_dist)
        call stats_1d_add(stats_mass_nh4_dist, mass_nh4_dist)
        
        dry_masses = aero_state_masses(aero_state, aero_data, &
             exclude=(/"H2O"/))
        oc_fracs = oc_masses / dry_masses
        diam_oc_dist = bin_grid_histogram_2d(diam_grid, dry_diameters, &
             oc_grid, oc_fracs, num_concs)
        call stats_2d_add(stats_diam_oc_dist, diam_oc_dist)

        crit_rhs = aero_state_crit_rel_humids(aero_state, aero_data, &
             env_state)
        scs = crit_rhs - 1d0
        diam_sc_dist = bin_grid_histogram_2d(diam_grid, dry_diameters, &
             sc_grid, scs, num_concs)
        call stats_2d_add(stats_diam_sc_dist, diam_sc_dist)

        call aero_state_mixing_state_metrics(aero_state, aero_data, &
             d_alpha, d_gamma, chi, group=["NH4", "SO4"], exclude=(/"H2O"/))

        call stats_1d_add_entry(stats_d_alpha, d_alpha, i_index)
        call stats_1d_add_entry(stats_d_gamma, d_gamma, i_index)
        call stats_1d_add_entry(stats_chi, chi, i_index)
    
     end do

     call make_filename(out_filename, prefix, "_process.nc", index)
     write(*,*) "Writing " // trim(out_filename)
     call pmc_nc_open_write(out_filename, ncid)
     call pmc_nc_write_info(ncid, uuid, "1_urban_plume process")
     call bin_grid_output_netcdf(diam_grid, ncid, "diam", unit="m")
     call bin_grid_output_netcdf(oc_grid, ncid, "oc_frac", unit="1")
     call bin_grid_output_netcdf(sc_grid, ncid, "sc", unit="1")
     call bin_grid_output_netcdf(kappa_grid, ncid, "kappa", unit="1")

     call stats_1d_output_netcdf(stats_num_dist, ncid, "num_dist", &
          dim_name="diam", unit="m^{-3}")
     call stats_1d_clear(stats_num_dist)

     call stats_1d_output_netcdf(stats_kappa_dist, ncid, "kappa_dist", &
          dim_name="kappa", unit="m^{-3}")
     call stats_1d_clear(stats_kappa_dist)

     call stats_1d_output_netcdf(stats_mass_oc_dist, ncid, "mass_oc_dist", &
          dim_name="diam", unit="m^{-3}")
     call stats_1d_clear(stats_mass_oc_dist)
     call stats_1d_output_netcdf(stats_mass_so4_dist, ncid, "mass_so4_dist", &
          dim_name="diam", unit="m^{-3}")
     call stats_1d_clear(stats_mass_so4_dist)
     call stats_1d_output_netcdf(stats_mass_no3_dist, ncid, "mass_no3_dist", &
          dim_name="diam", unit="m^{-3}")
     call stats_1d_clear(stats_mass_no3_dist)
     call stats_1d_output_netcdf(stats_mass_nh4_dist, ncid, "mass_nh4_dist", &
          dim_name="diam", unit="m^{-3}")
     call stats_1d_clear(stats_mass_nh4_dist)

     call stats_2d_output_netcdf(stats_diam_oc_dist, ncid, "diam_oc_dist", &
          dim_name_1="diam", dim_name_2="oc_frac", unit="m^{-3}")
     call stats_2d_clear(stats_diam_oc_dist)

     call stats_2d_output_netcdf(stats_diam_sc_dist, ncid, "diam_sc_dist", &
          dim_name_1="diam", dim_name_2="sc", unit="m^{-3}")
     call stats_2d_clear(stats_diam_sc_dist)

     call stats_2d_output_netcdf(stats_diam_kappa_dist, ncid, "diam_kappa_dist", &
          dim_name_1="diam", dim_name_2="kappa", unit="m^{-3}")
     call stats_2d_clear(stats_diam_kappa_dist)

     call pmc_nc_close(ncid)
  end do

  call make_filename(out_filename, prefix, "_process.nc")
  write(*,*) "Writing " // trim(out_filename)
  call pmc_nc_open_write(out_filename, ncid)
  call pmc_nc_write_info(ncid, uuid, "1_urban_plume process")
  call pmc_nc_write_real_1d(ncid, times, "time", dim_name="time", unit="s")
  call stats_1d_output_netcdf(stats_tot_num_conc, ncid, "tot_num_conc", &
       dim_name="time", unit="m^{-3}")
  call stats_1d_output_netcdf(stats_tot_mass_conc, ncid, "tot_mass_conc", &
       dim_name="time", unit="kg m^{-3}")
  call stats_1d_output_netcdf(stats_d_alpha, ncid, "d_alpha", &
       dim_name="time", unit="1")
  call stats_1d_output_netcdf(stats_d_gamma, ncid, &
       "d_gamma", dim_name="time", unit="1")
  call stats_1d_output_netcdf(stats_chi, ncid, "chi", &
       dim_name="time", unit="1")
  call pmc_nc_close(ncid)

  call pmc_mpi_finalize()

end program process
