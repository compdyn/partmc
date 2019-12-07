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
  type(bin_grid_t) :: diam_grid, bc_grid, sc_grid, entropy_grid, avg_bin_grid, &
                      gamma_grid, coating_grid
  type(aero_data_t) :: aero_data
  type(aero_state_t) :: aero_state, aero_state_averaged
  type(env_state_t) :: env_state
  integer :: ncid, index, repeat, i_index, i_repeat, n_index, n_repeat
  real(kind=dp) :: time, del_t, tot_num_conc, tot_mass_conc, tot_entropy
  real(kind=dp) :: tot_entropy_averaged, gamma_1, gamma_2, R, H_D_org, pi
  real(kind=dp) :: temperature, m_n2o5, v_n2o5
  character(len=PMC_UUID_LEN) :: uuid
  real(kind=dp), allocatable :: times(:), dry_diameters(:), wet_diameters(:), &
       num_concs(:), soa_masses(:), oc_masses(:), soa_volume(:), oc_volume(:), &
       organic_volume(:), total_volume(:), core_volume(:), core_diameters(:), &
       coating_thickness(:), wet_surfaces(:), surface_frac(:), &
       dry_masses(:), masses(:), bc_masses(:), so4_masses(:), no3_masses(:), &
       f_fraction(:), gamma_core(:), gamma_coat(:), gamma(:), bc_fracs(:), &
       num_concs_averaged(:), dry_masses_averaged(:), masses_averaged(:), &
       entropies(:), entropies_averaged(:), crit_rhs(:), scs(:), num_dist(:), &
       diam_bc_dist(:,:), diam_sc_dist(:,:), diam_gamma_dist(:,:), entropy_dist(:), &
       diam_coating_dist(:,:)
  real(kind=dp) :: total_so4, total_no3, bulk_f_fraction, bulk_gamma, total_surface, &
       gamma_population
  type(stats_1d_t) :: stats_num_dist, stats_entropy_dist, stats_tot_num_conc, &
       stats_tot_mass_conc, stats_tot_entropy, stats_tot_entropy_averaged
  type(stats_2d_t) :: stats_diam_bc_dist, stats_diam_sc_dist, stats_diam_gamma_dist, &
       stats_diam_coating_dist

  call pmc_mpi_init()

  call bin_grid_allocate(diam_grid)
  call bin_grid_allocate(bc_grid)
  call bin_grid_allocate(sc_grid)
  call bin_grid_allocate(entropy_grid)
  call bin_grid_allocate(gamma_grid)
  call bin_grid_allocate(coating_grid)
  call aero_data_allocate(aero_data)
  call aero_state_allocate(aero_state)
  call aero_state_allocate(aero_state_averaged)
  call bin_grid_allocate(avg_bin_grid)
  call env_state_allocate(env_state)

  call input_n_files(prefix, n_repeat, n_index)

  call bin_grid_make(diam_grid, BIN_GRID_TYPE_LOG, 180, 1d-9, 1d-3)
  call bin_grid_make(bc_grid, BIN_GRID_TYPE_LINEAR, 50, 0d0, 1d0)
  call bin_grid_make(sc_grid, BIN_GRID_TYPE_LOG, 50, 1d-4, 1d0)
  call bin_grid_make(gamma_grid, BIN_GRID_TYPE_LINEAR, 50, 0d0, 0.02d0)
  call bin_grid_make(coating_grid, BIN_GRID_TYPE_LINEAR, 50, 0d0, 1d-8)
  call bin_grid_make(entropy_grid, BIN_GRID_TYPE_LINEAR, 50, 0d0, 1d0)
  call bin_grid_make(avg_bin_grid, BIN_GRID_TYPE_LOG, 1, 1d-30, 1d10)

  allocate(times(n_index))

  gamma_1 = 0.02
  gamma_2 = 0.002
  R = 8.3145 ! J/(mol*K)
  pi = 3.14159265359
  H_D_org = 0.03 * 5000 / (1e-3 * 101300) * 1e-9   ! H_aq = 5000 M/atms, D_aq = 1e-9 m^2/s, H_org*D_org=0.03*H_aq*D_aq
  m_n2o5 = 108e-3 ! kg / mol
  
  
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
        v_n2o5 = (8 * R * temperature / (pi * m_n2o5))**0.5
                
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

        soa_volume = soa_masses / 1400d0
        oc_volume = oc_masses / 1000e0

        organic_volume = soa_volume + oc_volume
        total_volume = 3.14159 / 6d0 * wet_diameters**3d0

        core_volume = total_volume - organic_volume
        core_diameters = (6d0 / 3.14159 * core_volume)**(1./3.)

        coating_thickness = (wet_diameters - core_diameters) / 2.
        write(6,*)'coating ', coating_thickness(50:55)

        wet_surfaces = pi * wet_diameters**2.
        total_surface = sum(wet_surfaces * num_concs)
        surface_frac = wet_surfaces / total_surface
                
        total_so4 = sum(so4_masses * num_concs)
        total_no3 = sum(no3_masses * num_concs)
        bulk_f_fraction = total_so4 / (total_so4 + total_no3)
        bulk_gamma = bulk_f_fraction * gamma_1 + (1 - bulk_f_fraction) * gamma_2
                
        bc_fracs = bc_masses / dry_masses

        f_fraction = so4_masses / (so4_masses + no3_masses)

        gamma_core = f_fraction * gamma_1 + (1 - f_fraction) * gamma_2

        gamma_coat = 4 * R * temperature * H_D_org * core_diameters / &
             (v_n2o5 * coating_thickness * wet_diameters) 

        gamma = 1 / (1 / gamma_core + 1 / gamma_coat)
        gamma_population = sum(surface_frac * gamma)
        write(*,*) 'surface_frac ', sum(surface_frac)
        write(*,*) 'gamma population and bulk ', i_index, bulk_gamma, gamma_population
        
        diam_coating_dist = bin_grid_histogram_2d(diam_grid, dry_diameters, &
             coating_grid, coating_thickness, num_concs)
        call stats_2d_add(stats_diam_coating_dist, diam_coating_dist)

        
        diam_gamma_dist = bin_grid_histogram_2d(diam_grid, dry_diameters, &
             gamma_grid, gamma, num_concs)
        call stats_2d_add(stats_diam_gamma_dist, diam_gamma_dist)
        
        diam_bc_dist = bin_grid_histogram_2d(diam_grid, dry_diameters, &
             bc_grid, bc_fracs, num_concs)
        call stats_2d_add(stats_diam_bc_dist, diam_bc_dist)

        crit_rhs = aero_state_crit_rel_humids(aero_state, aero_data, &
             env_state)
        scs = crit_rhs - 1d0
        diam_sc_dist = bin_grid_histogram_2d(diam_grid, dry_diameters, &
             sc_grid, scs, num_concs)
        call stats_2d_add(stats_diam_sc_dist, diam_sc_dist)

        entropies = aero_state_mass_entropies(aero_state, aero_data) !, &
             !exclude=["H2O"]) !, group=["BC"])
        entropy_dist = bin_grid_histogram_1d(entropy_grid, entropies, &
             num_concs)
        call stats_1d_add(stats_entropy_dist, entropy_dist)

        tot_entropy = sum(entropies * masses * num_concs) &
             / sum(masses * num_concs)
        call stats_1d_add_entry(stats_tot_entropy, tot_entropy, i_index)

        call aero_state_copy(aero_state, aero_state_averaged)
        call aero_state_bin_average_comp(aero_state_averaged, avg_bin_grid, &
             aero_data)
        num_concs_averaged = aero_state_num_concs(aero_state_averaged)
        masses_averaged = aero_state_masses(aero_state_averaged, aero_data)
        dry_masses_averaged = aero_state_masses(aero_state_averaged, &
             aero_data, exclude=(/"H2O"/))
        entropies_averaged = aero_state_mass_entropies(aero_state_averaged, &
             aero_data) !, exclude=["H2O"]) !, group=["BC"])
        tot_entropy_averaged &
             = sum(entropies_averaged * masses_averaged &
             * num_concs_averaged) &
             / sum(masses_averaged * num_concs_averaged)
        call stats_1d_add_entry(stats_tot_entropy_averaged, &
             tot_entropy_averaged, i_index)
     end do

     call make_filename(out_filename, prefix, "_process.nc", index)
     call pmc_nc_open_write(out_filename, ncid)
     call pmc_nc_write_info(ncid, uuid, "1_urban_plume process")
     call bin_grid_output_netcdf(diam_grid, ncid, "diam", unit="m")
     call bin_grid_output_netcdf(bc_grid, ncid, "bc_frac", unit="1")
     call bin_grid_output_netcdf(sc_grid, ncid, "sc", unit="1")
     call bin_grid_output_netcdf(entropy_grid, ncid, "entropy", unit="1")
     call bin_grid_output_netcdf(gamma_grid, ncid, "gamma", unit="1")
     call bin_grid_output_netcdf(coating_grid, ncid, "coating", unit="m")

     call stats_1d_output_netcdf(stats_num_dist, ncid, "num_dist", &
          dim_name="diam", unit="m^{-3}")
     call stats_1d_clear(stats_num_dist)

     call stats_2d_output_netcdf(stats_diam_bc_dist, ncid, "diam_bc_dist", &
          dim_name_1="diam", dim_name_2="bc_frac", unit="m^{-3}")
     call stats_2d_clear(stats_diam_bc_dist)

     call stats_2d_output_netcdf(stats_diam_sc_dist, ncid, "diam_sc_dist", &
          dim_name_1="diam", dim_name_2="sc", unit="m^{-3}")
     call stats_2d_clear(stats_diam_sc_dist)

     call stats_2d_output_netcdf(stats_diam_gamma_dist, ncid, "diam_gamma_dist", &
          dim_name_1="diam", dim_name_2="gamma", unit="m^{-3}")
     call stats_2d_clear(stats_diam_gamma_dist)

     call stats_2d_output_netcdf(stats_diam_coating_dist, ncid, "diam_coating_dist", &
          dim_name_1="diam", dim_name_2="coating", unit="m^{-3}")
     call stats_2d_clear(stats_diam_coating_dist)
     
     call stats_1d_output_netcdf(stats_entropy_dist, ncid, "entropy_dist", &
          dim_name="entropy", unit="m^{-3}")
     call stats_1d_clear(stats_entropy_dist)

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
  call stats_1d_output_netcdf(stats_tot_entropy, ncid, "tot_entropy", &
       dim_name="time", unit="m^{-3}")
  call stats_1d_output_netcdf(stats_tot_entropy_averaged, ncid, &
       "tot_entropy_averaged", dim_name="time", unit="m^{-3}")
  call pmc_nc_close(ncid)

  call bin_grid_deallocate(diam_grid)
  call aero_data_deallocate(aero_data)
  call aero_state_deallocate(aero_state)
  call env_state_deallocate(env_state)

  call pmc_mpi_finalize()

end program process
