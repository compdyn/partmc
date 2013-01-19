! Copyright (C) 2009-2012 Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The process program.

!> Read NetCDF output files and process them.
program process

  use pmc_output
  use pmc_stats

  character(len=PMC_MAX_FILENAME_LEN), parameter :: prefix &
       = "out/urban_plume2"

  character(len=PMC_MAX_FILENAME_LEN) :: in_filename, out_filename
  type(bin_grid_t) :: diam_grid, bc_grid, no3_grid, h2o_grid, sc_grid, &
       entropy_grid, avg_bin_grid, time_grid
  type(aero_data_t) :: aero_data
  type(aero_state_t) :: aero_state, aero_state_averaged, aero_state_source
  type(env_state_t) :: env_state
  integer :: ncid, index, repeat, i_index, i_repeat, n_index, n_repeat
  real(kind=dp) :: time, del_t, tot_num_conc, tot_mass_conc, avg_part_entropy, ccn
  real(kind=dp) :: entropy_of_avg_part
  character(len=PMC_UUID_LEN) :: uuid
  real(kind=dp), allocatable :: times(:), dry_diameters(:), num_concs(:), &
       dry_masses(:), masses(:), bc_masses(:), bc_fracs(:), &
       h2o_masses(:), h2o_fracs(:), no3_masses(:), no3_fracs(:), &
       oc_masses(:), oc_fracs(:), so4_masses(:), so4_fracs(:), &
       num_concs_averaged(:), dry_masses_averaged(:), masses_averaged(:), &
       entropies(:), entropies_of_avg_part(:), crit_rhs(:), scs(:), num_dist(:), &
       diam_bc_dist(:,:), diam_no3_dist(:,:), diam_h2o_dist(:,:), &
       diam_sc_dist(:,:), entropy_dist(:), diam_entropy_dist(:,:), time_entropy_dist(:,:)
  real(kind=dp), allocatable :: diesel(:), cooking(:), gasoline(:), init_small(:), &
       init_large(:), back_small(:), back_large(:)
  real(kind=dp) :: tot_entropy_conc, tot_entropy_of_avg_conc, & 
       tot_entropy_ratio
  real(kind=dp), allocatable :: dist_ratio_to_entropy_of_avg_part(:), &
       dist_ratio_to_avg_part_entropy(:)
  real(kind=dp), allocatable :: least_create_times(:), greatest_create_times(:)
  type(stats_1d_t) :: stats_num_dist, stats_entropy_dist, &
       stats_dist_ratio_to_entropy_of_avg_part, stats_dist_ratio_to_avg_part_entropy, &
       stats_tot_num_conc, stats_tot_mass_conc, stats_avg_part_entropy, &
       stats_entropy_of_avg_part, stats_ccn(3), stats_tot_entropy_conc, &
       stats_tot_entropy_of_avg_conc, stats_tot_entropy_ratio
  type(stats_2d_t) :: stats_diam_bc_dist, stats_diam_no3_dist, stats_diam_h2o_dist, &
       stats_diam_sc_dist, stats_diam_entropy_dist, stats_time_entropy_dist

  real(kind=dp) :: s_env(3)
  real(kind=dp), allocatable :: num_concs_active(:)
  logical, allocatable :: is_active(:), is_source(:), is_not_source(:)

  s_env(1) = 0.001
  s_env(2) = 0.003
  s_env(3) = 0.006

  call pmc_mpi_init()

  call bin_grid_allocate(diam_grid)
  call bin_grid_allocate(bc_grid)
  call bin_grid_allocate(no3_grid)
  call bin_grid_allocate(h2o_grid)
  call bin_grid_allocate(sc_grid)
  call bin_grid_allocate(entropy_grid)
  call bin_grid_allocate(time_grid)
  call aero_data_allocate(aero_data)
  call aero_state_allocate(aero_state)
  call aero_state_allocate(aero_state_source)
  call aero_state_allocate(aero_state_averaged)
  call bin_grid_allocate(avg_bin_grid)
  call env_state_allocate(env_state)

  call input_n_files(prefix, n_repeat, n_index)

  call bin_grid_make(diam_grid, BIN_GRID_TYPE_LOG, 180, 1d-9, 1d-3)
  call bin_grid_make(bc_grid, BIN_GRID_TYPE_LINEAR, 50, 0d0, 1d0)
  call bin_grid_make(no3_grid, BIN_GRID_TYPE_LINEAR, 50, 0d0, 1d0)
  call bin_grid_make(h2o_grid, BIN_GRID_TYPE_LINEAR, 50, 0d0, 2d0)
  call bin_grid_make(sc_grid, BIN_GRID_TYPE_LOG, 50, 1d-4, 1d0)
  call bin_grid_make(entropy_grid, BIN_GRID_TYPE_LINEAR, 200, 0d0, 5d0)
  call bin_grid_make(time_grid, BIN_GRID_TYPE_LINEAR, 49, 0d0, 48d0)
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

        dry_masses = aero_state_masses(aero_state, aero_data, &
             exclude=(/"H2O"/))
        bc_masses = aero_state_masses(aero_state, aero_data, &
             include=(/"BC"/))
        bc_fracs = bc_masses / dry_masses
        h2o_masses = aero_state_masses(aero_state, aero_data, &
             include=(/"H2O"/))
        h2o_fracs = h2o_masses / dry_masses
        oc_masses = aero_state_masses(aero_state, aero_data, &
             include=(/"OC"/))
        oc_fracs = oc_masses / dry_masses
        no3_masses = aero_state_masses(aero_state, aero_data, &
             include=(/"NO3"/))
        no3_fracs = no3_masses / dry_masses
        so4_masses = aero_state_masses(aero_state, aero_data, &
             include=(/"SO4"/))
        so4_fracs = so4_masses / dry_masses

        diesel = aero_state_n_orig_part(aero_state, aero_data, "diesel")
        cooking = aero_state_n_orig_part(aero_state, aero_data, "cooking")
        gasoline = aero_state_n_orig_part(aero_state, aero_data, "gasoline")
        init_small = aero_state_n_orig_part(aero_state, aero_data, "init_small")
        init_large = aero_state_n_orig_part(aero_state, aero_data, "init_large")
        back_small = aero_state_n_orig_part(aero_state, aero_data, "back_small")
        back_large = aero_state_n_orig_part(aero_state, aero_data, "back_large")
 
        call aero_state_copy(aero_state, aero_state_source)
        is_not_source = (back_large == 0d0)

        call aero_state_remove_particles(aero_state_source, is_not_source)

        times(i_index) = time

        dry_diameters = aero_state_dry_diameters(aero_state_source, aero_data)
        num_concs = aero_state_num_concs(aero_state_source)
        num_dist = bin_grid_histogram_1d(diam_grid, dry_diameters, num_concs)
        call stats_1d_add(stats_num_dist, num_dist)

        tot_num_conc = sum(num_concs)
        call stats_1d_add_entry(stats_tot_num_conc, tot_num_conc, i_index)

        masses = aero_state_masses(aero_state_source, aero_data)
        tot_mass_conc = sum(masses * num_concs)
        call stats_1d_add_entry(stats_tot_mass_conc, tot_mass_conc, i_index)
        
        least_create_times = aero_state_least_create_times(aero_state_source)
	greatest_create_times = aero_state_greatest_create_times(aero_state_source) 
        
!> per-particle mixing entropy (H_i)
        entropies = aero_state_mass_entropies(aero_state_source, aero_data, &
             exclude=["H2O"])!, group=["BC"])

!> 1d distribution of H_i
        entropy_dist = bin_grid_histogram_1d(entropy_grid, entropies, &
             num_concs)
        call stats_1d_add(stats_entropy_dist, entropy_dist)

!> 2d distribution of number conc. in terms of H_i and dry diameter
        diam_entropy_dist = bin_grid_histogram_2d(diam_grid, dry_diameters, &
             entropy_grid, entropies, num_concs)
        call stats_2d_add(stats_diam_entropy_dist, diam_entropy_dist)

!> 2d distribution of number conc. in terms of H_i and time
        call stats_2d_add_col(stats_time_entropy_dist, entropy_dist, i_index)

!> average per-particle entropy (\bar{H})
        avg_part_entropy = sum(entropies * masses * num_concs) &
             / sum(masses * num_concs)
        call stats_1d_add_entry(stats_avg_part_entropy, avg_part_entropy, i_index)

!> 1d distribution of ratio of per-particle entropies to average per-particle entropy (\bar{R_i})
        dist_ratio_to_avg_part_entropy = bin_grid_histogram_1d(entropy_grid, &
             entropies/avg_part_entropy, num_concs)
        call stats_1d_add(stats_dist_ratio_to_avg_part_entropy, dist_ratio_to_avg_part_entropy)

!> entropy concentration (c_H)
        tot_entropy_conc = sum(entropies * masses * num_concs)
        call stats_1d_add_entry(stats_tot_entropy_conc, tot_entropy_conc, i_index)

!> composition-averaging
        call aero_state_copy(aero_state_source, aero_state_averaged)
        call aero_state_bin_average_comp(aero_state_averaged, avg_bin_grid, &
             aero_data, dry_volume=.false.)
        num_concs_averaged = aero_state_num_concs(aero_state_averaged)
        masses_averaged = aero_state_masses(aero_state_averaged, aero_data)
        dry_masses_averaged = aero_state_masses(aero_state_averaged, &
             aero_data, exclude=(/"H2O"/))

!> per-particle mixing entropy after composition-averaging  (\hat{H_i})
        entropies_of_avg_part = aero_state_mass_entropies(aero_state_averaged, &
             aero_data, exclude=["H2O"])!, group=["BC"])

!> 1d distribution of ratio of per-particle entropies to entropies after composition-averaging (\hat{R_i})
        dist_ratio_to_entropy_of_avg_part = bin_grid_histogram_1d(entropy_grid, &
             entropies/entropies_of_avg_part, num_concs)
        call stats_1d_add(stats_dist_ratio_to_entropy_of_avg_part, dist_ratio_to_entropy_of_avg_part)

!> average per-particle entropy after composition-averaging (\bar{\hat{H}})
        entropy_of_avg_part &
             = sum(entropies_of_avg_part * masses_averaged &
             * num_concs_averaged) &
             / sum(masses_averaged * num_concs_averaged)
        call stats_1d_add_entry(stats_entropy_of_avg_part, &
             entropy_of_avg_part, i_index)

!> entropy concentration after composition-averaging (c_{\hat{H}})
        tot_entropy_of_avg_conc = sum(entropies_of_avg_part * masses * num_concs)
        call stats_1d_add_entry(stats_tot_entropy_of_avg_conc, tot_entropy_of_avg_conc, i_index)

!> ratio of entropy conc. and entropy conc. after composition-averaging (R)
        tot_entropy_ratio = tot_entropy_conc / tot_entropy_of_avg_conc
        call stats_1d_add_entry(stats_tot_entropy_ratio, tot_entropy_ratio, i_index)

        call make_filename(out_filename, prefix, "_process_source.nc", index, repeat)
        call pmc_nc_open_write(out_filename, ncid)
        call pmc_nc_write_info(ncid, uuid, "1_urban_plume process")
        call pmc_nc_write_real_1d(ncid, entropies, "entropies", dim_name="particle", unit="1")
        call pmc_nc_write_real_1d(ncid, dry_diameters, "dry_diam", dim_name="particle", unit="m")
        call pmc_nc_write_real_1d(ncid, least_create_times, "least_ct", dim_name="particle", unit="1")
        call pmc_nc_write_real_1d(ncid, greatest_create_times, "greatest_ct", dim_name="particle", unit="1")
        call pmc_nc_close(ncid)
     end do
     
     call make_filename(out_filename, prefix, "_process_source.nc", index)
     call pmc_nc_open_write(out_filename, ncid)
     call pmc_nc_write_info(ncid, uuid, "1_urban_plume process")
     call bin_grid_output_netcdf(diam_grid, ncid, "diam", unit="m")
     call bin_grid_output_netcdf(entropy_grid, ncid, "entropy", unit="1")

     call stats_1d_output_netcdf(stats_num_dist, ncid, "num_dist", &
          dim_name="diam", unit="m^{-3}")
     call stats_1d_clear(stats_num_dist)

     call stats_1d_output_netcdf(stats_entropy_dist, ncid, "entropy_dist", &
          dim_name="entropy", unit="m^{-3}")
     call stats_1d_clear(stats_entropy_dist)

     call stats_1d_output_netcdf(stats_dist_ratio_to_entropy_of_avg_part, ncid, &
          "dist_ratio_to_entropy_of_avg_part", dim_name="entropy", unit="1")
     call stats_1d_clear(stats_dist_ratio_to_entropy_of_avg_part)

     call stats_1d_output_netcdf(stats_dist_ratio_to_avg_part_entropy, ncid, &
          "dist_ratio_to_avg_part_entropy", dim_name="entropy", unit="1")
     call stats_1d_clear(stats_dist_ratio_to_avg_part_entropy)

     call stats_2d_output_netcdf(stats_diam_entropy_dist, ncid, "diam_entropy_dist", &
          dim_name_1="diam", dim_name_2="entropy", unit="m^{-3}")
     call stats_2d_clear(stats_diam_entropy_dist)

     call pmc_nc_close(ncid)
  end do

  call make_filename(out_filename, prefix, "_process_source.nc")
  write(6,*) out_filename
  call pmc_nc_open_write(out_filename, ncid)
  call pmc_nc_write_info(ncid, uuid, "1_urban_plume process")
  call pmc_nc_write_real_1d(ncid, times, "time", dim_name="time", unit="s")
  call bin_grid_output_netcdf(entropy_grid, ncid, "entropy", unit="1")
  call bin_grid_output_netcdf(time_grid, ncid, "time_grid", unit="h")
  call stats_1d_output_netcdf(stats_tot_num_conc, ncid, "tot_num_conc", &
       dim_name="time", unit="m^{-3}")
  call stats_1d_output_netcdf(stats_tot_mass_conc, ncid, "tot_mass_conc", &
       dim_name="time", unit="m^{-3}")
  call stats_1d_output_netcdf(stats_avg_part_entropy, ncid, "avg_part_entropy", &
       dim_name="time", unit="m^{-3}")
  call stats_1d_output_netcdf(stats_entropy_of_avg_part, ncid, &
       "entropy_of_avg_part", dim_name="time", unit="m^{-3}")
  call stats_1d_output_netcdf(stats_tot_entropy_conc, ncid, "tot_entropy_conc", &
       dim_name="time", unit="m^{-3}")
  call stats_1d_output_netcdf(stats_tot_entropy_of_avg_conc, ncid, "tot_entropy_of_avg_conc", &
       dim_name="time", unit="m^{-3}")
  call stats_1d_output_netcdf(stats_tot_entropy_ratio, ncid, "tot_entropy_ratio", &
       dim_name="time", unit="1")
  call stats_2d_output_netcdf(stats_time_entropy_dist, ncid, "time_entropy_dist", &
       dim_name_1="entropy", dim_name_2="time", unit="m^{-3}")
  call stats_2d_clear(stats_time_entropy_dist)


  call pmc_nc_close(ncid)

  call bin_grid_deallocate(diam_grid)
  call aero_data_deallocate(aero_data)
  call aero_state_deallocate(aero_state_source)
  call env_state_deallocate(env_state)

  call pmc_mpi_finalize()

end program process
