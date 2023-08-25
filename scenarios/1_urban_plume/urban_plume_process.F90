! Copyright (C) 2009-2013, 2016, 2017, 2021 Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The process program.

!> Read NetCDF output files and process them.
program process

  use pmc_output
  use pmc_stats
  use pmc_aero_state
  use pmc_aero_data
  use pmc_env_state
  use pmc_aero_particle
  use pmc_aero_particle_array

  character(len=PMC_MAX_FILENAME_LEN), parameter :: prefix &
       = "out/urban_plume"
  character(len=PMC_MAX_FILENAME_LEN), parameter :: prefix_new &
       = "out/single_particle"
       
  character(len=PMC_MAX_FILENAME_LEN) :: in_filename, out_filename, &
                                         single_particle_output
  type(bin_grid_t) :: diam_grid, bc_grid, oc_grid, sc_grid, &
                      sc_varying_sigma_grid, &
                      so4_grid, no3_grid, nh4_grid, soa_grid, &
                      cl_grid, msa_grid, aro1_grid, aro2_grid, &
                      alk1_grid, ole1_grid, api1_grid, api2_grid, &
                      lim1_grid, lim2_grid, co3_grid, na_grid, &
                      ca_grid, oin_grid, h2o_grid
  type(aero_data_t) :: aero_data
  type(aero_particle_t) :: aero_particle
  type(aero_state_t) :: aero_state
  type(env_state_t) :: env_state
  integer :: ncid, index, repeat, i_index, i_repeat, n_index, n_repeat
  integer :: number
  real(kind=dp) :: crit_diam
  real(kind=dp) :: time, del_t, tot_num_conc, tot_mass_conc, &
               tot_dry_mass_conc, tot_bc_mass_conc, tot_oc_mass_conc, &
               tot_so4_mass_conc, tot_no3_mass_conc, tot_nh4_mass_conc, &
               tot_soa_mass_conc, tot_cl_mass_conc, tot_msa_mass_conc, &
               tot_aro1_mass_conc, tot_aro2_mass_conc, tot_alk1_mass_conc, &
               tot_ole1_mass_conc, tot_api1_mass_conc, tot_api2_mass_conc, &
               tot_lim1_mass_conc, tot_lim2_mass_conc, tot_co3_mass_conc, &
               tot_na_mass_conc, tot_ca_mass_conc, tot_oin_mass_conc, &
               tot_h2o_mass_conc
  real(kind=dp) :: d_alpha, d_gamma, chi
  character(len=PMC_UUID_LEN) :: uuid
  real(kind=dp), allocatable :: times(:), dry_diameters(:), num_concs(:), &
       dry_masses(:), masses(:), bc_masses(:), bc_fracs(:), &
       diam_bc_dist(:,:), oc_masses(:), oc_fracs(:), diam_oc_dist(:,:), &
       crit_rhs(:), scs(:), num_dist(:), diam_sc_dist(:,:), &
       sc_dist(:), sc_varying_sigma_dist(:), &
       crit_rhs_varying_sigma(:), scs_varying_sigma(:), &
       diam_sc_varying_sigma_dist(:,:), &
       so4_masses(:), so4_fracs(:), diam_so4_dist(:,:), &
       no3_masses(:), no3_fracs(:), diam_no3_dist(:,:), &
       nh4_masses(:), nh4_fracs(:), diam_nh4_dist(:,:), &
       soa_masses(:), soa_fracs(:), diam_soa_dist(:,:), &
       cl_masses(:), cl_fracs(:), diam_cl_dist(:,:), &
       msa_masses(:), msa_fracs(:), diam_msa_dist(:,:), &
       aro1_masses(:), aro1_fracs(:), diam_aro1_dist(:,:), &
       aro2_masses(:), aro2_fracs(:), diam_aro2_dist(:,:), &
       alk1_masses(:), alk1_fracs(:), diam_alk1_dist(:,:), &
       ole1_masses(:), ole1_fracs(:), diam_ole1_dist(:,:), &
       api1_masses(:), api1_fracs(:), diam_api1_dist(:,:), &
       api2_masses(:), api2_fracs(:), diam_api2_dist(:,:), &
       lim1_masses(:), lim1_fracs(:), diam_lim1_dist(:,:), &
       lim2_masses(:), lim2_fracs(:), diam_lim2_dist(:,:), &
       co3_masses(:), co3_fracs(:), diam_co3_dist(:,:), &
       na_masses(:), na_fracs(:), diam_na_dist(:,:), &
       ca_masses(:), ca_fracs(:), diam_ca_dist(:,:), &
       oin_masses(:), oin_fracs(:), diam_oin_dist(:,:), &
       h2o_masses(:), h2o_fracs(:), diam_h2o_dist(:,:)

  type(stats_1d_t) :: stats_num_dist, stats_d_alpha, stats_tot_num_conc, &
                      stats_tot_mass_conc, stats_d_gamma, stats_chi, & 
                      stats_tot_bc_mass_conc, stats_tot_oc_mass_conc, &
                      stats_tot_so4_mass_conc, stats_tot_no3_mass_conc, &
                      stats_tot_nh4_mass_conc, stats_tot_dry_mass_conc, &
                      stats_tot_soa_mass_conc, stats_tot_cl_mass_conc, &  
                      stats_tot_msa_mass_conc, stats_tot_aro1_mass_conc, &  
                      stats_tot_aro2_mass_conc, stats_tot_alk1_mass_conc, &  
                      stats_tot_ole1_mass_conc, stats_tot_api1_mass_conc, &  
                      stats_tot_api2_mass_conc, stats_tot_lim1_mass_conc, &  
                      stats_tot_lim2_mass_conc, stats_tot_co3_mass_conc, &  
                      stats_tot_na_mass_conc, stats_tot_ca_mass_conc, &  
                      stats_tot_oin_mass_conc, stats_tot_h2o_mass_conc, &
                      stats_sc_dist, stats_sc_varying_sigma_dist        

  type(stats_2d_t) :: stats_diam_sc_dist, stats_diam_bc_dist, &
                      stats_diam_sc_varying_sigma_dist, &
                      stats_diam_oc_dist, stats_diam_so4_dist, &
                      stats_diam_no3_dist, stats_diam_nh4_dist, &
                      stats_diam_soa_dist, stats_diam_cl_dist, &
                      stats_diam_msa_dist, stats_diam_aro1_dist, &
                      stats_diam_aro2_dist, stats_diam_alk1_dist, &
                      stats_diam_ole1_dist, stats_diam_api1_dist, &
                      stats_diam_api2_dist, stats_diam_lim1_dist, &
                      stats_diam_lim2_dist, stats_diam_co3_dist, &
                      stats_diam_na_dist, stats_diam_oin_dist, &
                      stats_diam_ca_dist, stats_diam_h2o_dist
       
  character(len=AERO_NAME_LEN), allocatable :: mixing_state_groups(:,:)

  call pmc_mpi_init()

  call input_n_files(prefix, n_repeat, n_index)

  call bin_grid_make(diam_grid, BIN_GRID_TYPE_LOG, 180, 1d-9, 1d-3)
  call bin_grid_make(bc_grid, BIN_GRID_TYPE_LINEAR, 50, 0d0, 1d0)
  call bin_grid_make(oc_grid, BIN_GRID_TYPE_LINEAR, 50, 0d0, 1d0)
  call bin_grid_make(so4_grid, BIN_GRID_TYPE_LINEAR, 50, 0d0, 1d0)
  call bin_grid_make(no3_grid, BIN_GRID_TYPE_LINEAR, 50, 0d0, 1d0)
  call bin_grid_make(nh4_grid, BIN_GRID_TYPE_LINEAR, 50, 0d0, 1d0)
  call bin_grid_make(soa_grid, BIN_GRID_TYPE_LINEAR, 50, 0d0, 1d0)
  call bin_grid_make(cl_grid, BIN_GRID_TYPE_LINEAR, 50, 0d0, 1d0)
  call bin_grid_make(msa_grid, BIN_GRID_TYPE_LINEAR, 50, 0d0, 1d0)
  call bin_grid_make(aro1_grid, BIN_GRID_TYPE_LINEAR, 50, 0d0, 1d0)
  call bin_grid_make(aro2_grid, BIN_GRID_TYPE_LINEAR, 50, 0d0, 1d0)
  call bin_grid_make(alk1_grid, BIN_GRID_TYPE_LINEAR, 50, 0d0, 1d0)
  call bin_grid_make(ole1_grid, BIN_GRID_TYPE_LINEAR, 50, 0d0, 1d0)
  call bin_grid_make(api1_grid, BIN_GRID_TYPE_LINEAR, 50, 0d0, 1d0)
  call bin_grid_make(api2_grid, BIN_GRID_TYPE_LINEAR, 50, 0d0, 1d0)
  call bin_grid_make(lim1_grid, BIN_GRID_TYPE_LINEAR, 50, 0d0, 1d0)
  call bin_grid_make(lim2_grid, BIN_GRID_TYPE_LINEAR, 50, 0d0, 1d0)
  call bin_grid_make(co3_grid, BIN_GRID_TYPE_LINEAR, 50, 0d0, 1d0)
  call bin_grid_make(na_grid, BIN_GRID_TYPE_LINEAR, 50, 0d0, 1d0)
  call bin_grid_make(ca_grid, BIN_GRID_TYPE_LINEAR, 50, 0d0, 1d0)
  call bin_grid_make(oin_grid, BIN_GRID_TYPE_LINEAR, 50, 0d0, 1d0)
  call bin_grid_make(h2o_grid, BIN_GRID_TYPE_LINEAR, 50, 0d0, 1d0)
  call bin_grid_make(sc_grid, BIN_GRID_TYPE_LOG, 50, 1d-4, 1d0)
  call bin_grid_make(sc_varying_sigma_grid, BIN_GRID_TYPE_LOG, 50, 1d-4, 1d0)

  allocate(times(n_index))

  allocate(mixing_state_groups(3, 4)) ! 3 groups, max 4 species per group
  mixing_state_groups(1,:) = ["OC    ", "BC    ", "      ", "      "]
  mixing_state_groups(2,:) = ["API1  ", "API2  ", "LIM1  ", "LIM2  "]
  mixing_state_groups(3,:) = ["SO4   ", "NO3   ", "NH4   ", "      "]


  scs = [ real(kind=dp) :: ] ! silence compiler warnings
  scs_varying_sigma = [ real(kind=dp) :: ]
  bc_fracs = [ real(kind=dp) :: ]
  oc_fracs = [ real(kind=dp) :: ]
  so4_fracs = [ real(kind=dp) :: ]
  no3_fracs = [ real(kind=dp) :: ]
  nh4_fracs = [ real(kind=dp) :: ]
  cl_fracs = [ real(kind=dp) :: ]
  msa_fracs = [ real(kind=dp) :: ]
  aro1_fracs = [ real(kind=dp) :: ]
  aro2_fracs = [ real(kind=dp) :: ]
  alk1_fracs = [ real(kind=dp) :: ]
  ole1_fracs = [ real(kind=dp) :: ]
  api1_fracs = [ real(kind=dp) :: ]
  api2_fracs = [ real(kind=dp) :: ]
  lim1_fracs = [ real(kind=dp) :: ]
  lim2_fracs = [ real(kind=dp) :: ]
  co3_fracs = [ real(kind=dp) :: ]
  na_fracs = [ real(kind=dp) :: ]
  ca_fracs = [ real(kind=dp) :: ]
  oin_fracs = [ real(kind=dp) :: ]
  h2o_fracs = [ real(kind=dp) :: ]

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
        num_concs = aero_state_num_concs(aero_state, aero_data)
        num_dist = bin_grid_histogram_1d(diam_grid, dry_diameters, num_concs)
        call stats_1d_add(stats_num_dist, num_dist)

        tot_num_conc = sum(num_concs)
        call stats_1d_add_entry(stats_tot_num_conc, tot_num_conc, i_index)

        masses = aero_state_masses(aero_state, aero_data)
        tot_mass_conc = sum(masses * num_concs)
        call stats_1d_add_entry(stats_tot_mass_conc, tot_mass_conc, i_index)
!for dry mass
        dry_masses = aero_state_masses(aero_state, aero_data, &
             exclude=(/"H2O"/))
        tot_dry_mass_conc = sum(dry_masses * num_concs)
        call stats_1d_add_entry(stats_tot_dry_mass_conc, tot_dry_mass_conc, i_index)
!for bc
        bc_masses = aero_state_masses(aero_state, aero_data, &
             include=(/"BC"/))
        tot_bc_mass_conc = sum(bc_masses * num_concs)
        call stats_1d_add_entry(stats_tot_bc_mass_conc, tot_bc_mass_conc, i_index)
    
        bc_fracs = bc_masses / dry_masses
        diam_bc_dist = bin_grid_histogram_2d(diam_grid, dry_diameters, &
             bc_grid, bc_fracs, num_concs)
        call stats_2d_add(stats_diam_bc_dist, diam_bc_dist)
!for OC
        oc_masses = aero_state_masses(aero_state, aero_data, &
             include=(/"OC"/))
        tot_oc_mass_conc = sum(oc_masses * num_concs)
        call stats_1d_add_entry(stats_tot_oc_mass_conc, tot_oc_mass_conc, i_index)
    
        oc_fracs = oc_masses / dry_masses
        diam_oc_dist = bin_grid_histogram_2d(diam_grid, dry_diameters, &
             oc_grid, oc_fracs, num_concs)
        call stats_2d_add(stats_diam_oc_dist, diam_oc_dist)
!for SO4
        so4_masses = aero_state_masses(aero_state, aero_data, &
             include=(/"SO4"/))
        tot_so4_mass_conc = sum(so4_masses * num_concs)
        call stats_1d_add_entry(stats_tot_so4_mass_conc, tot_so4_mass_conc, i_index)
    
        so4_fracs = so4_masses / dry_masses
        diam_so4_dist = bin_grid_histogram_2d(diam_grid, dry_diameters, &
             so4_grid, so4_fracs, num_concs)
        call stats_2d_add(stats_diam_so4_dist, diam_so4_dist)
!for NO3
        no3_masses = aero_state_masses(aero_state, aero_data, &
             include=(/"NO3"/))
        tot_no3_mass_conc = sum(no3_masses * num_concs)
        call stats_1d_add_entry(stats_tot_no3_mass_conc, tot_no3_mass_conc, i_index)
    
        no3_fracs = no3_masses / dry_masses
        diam_no3_dist = bin_grid_histogram_2d(diam_grid, dry_diameters, &
             no3_grid, no3_fracs, num_concs)
        call stats_2d_add(stats_diam_no3_dist, diam_no3_dist)
!for NH4
        nh4_masses = aero_state_masses(aero_state, aero_data, &
             include=(/"NH4"/))
        tot_nh4_mass_conc = sum(nh4_masses * num_concs)
        call stats_1d_add_entry(stats_tot_nh4_mass_conc, tot_nh4_mass_conc, i_index)
   
        nh4_fracs = nh4_masses / dry_masses
        diam_nh4_dist = bin_grid_histogram_2d(diam_grid, dry_diameters, &
             nh4_grid, nh4_fracs, num_concs)
        call stats_2d_add(stats_diam_nh4_dist, diam_nh4_dist)
!for SOA
        soa_masses = aero_state_masses(aero_state, aero_data, &
             include=(/"ARO1", "ARO2", "ALK1", "OLE1", "API1", "API2", "LIM1", "LIM2"/))
        tot_soa_mass_conc = sum(soa_masses * num_concs)
        call stats_1d_add_entry(stats_tot_soa_mass_conc, tot_soa_mass_conc, i_index)
    
        soa_fracs = soa_masses / dry_masses
        diam_soa_dist = bin_grid_histogram_2d(diam_grid, dry_diameters, &
             soa_grid, soa_fracs, num_concs)
        call stats_2d_add(stats_diam_soa_dist, diam_soa_dist)
!for Cl
        cl_masses = aero_state_masses(aero_state, aero_data, &
             include=(/"Cl"/))
        tot_cl_mass_conc = sum(cl_masses * num_concs)
        call stats_1d_add_entry(stats_tot_cl_mass_conc, tot_cl_mass_conc, i_index)
    
        cl_fracs = cl_masses / dry_masses
        diam_cl_dist = bin_grid_histogram_2d(diam_grid, dry_diameters, &
             cl_grid, cl_fracs, num_concs)
        call stats_2d_add(stats_diam_cl_dist, diam_cl_dist)
!for MSA
        msa_masses = aero_state_masses(aero_state, aero_data, &
             include=(/"MSA"/))
        tot_msa_mass_conc = sum(msa_masses * num_concs)
        call stats_1d_add_entry(stats_tot_msa_mass_conc, tot_msa_mass_conc, i_index)
    
        msa_fracs = msa_masses / dry_masses
        diam_msa_dist = bin_grid_histogram_2d(diam_grid, dry_diameters, &
             msa_grid, msa_fracs, num_concs)
        call stats_2d_add(stats_diam_msa_dist, diam_msa_dist)
!for aro1
        aro1_masses = aero_state_masses(aero_state, aero_data, &
             include=(/"ARO1"/))
        tot_aro1_mass_conc = sum(aro1_masses * num_concs)
        call stats_1d_add_entry(stats_tot_aro1_mass_conc, tot_aro1_mass_conc, i_index)
    
        aro1_fracs = aro1_masses / dry_masses
        diam_aro1_dist = bin_grid_histogram_2d(diam_grid, dry_diameters, &
             aro1_grid, aro1_fracs, num_concs)
        call stats_2d_add(stats_diam_aro1_dist, diam_aro1_dist) 
!for aro2
        aro2_masses = aero_state_masses(aero_state, aero_data, &
             include=(/"ARO2"/))
        tot_aro2_mass_conc = sum(aro2_masses * num_concs)
        call stats_1d_add_entry(stats_tot_aro2_mass_conc, tot_aro2_mass_conc, i_index)
    
        aro2_fracs = aro2_masses / dry_masses
        diam_aro2_dist = bin_grid_histogram_2d(diam_grid, dry_diameters, &
             aro2_grid, aro2_fracs, num_concs)
        call stats_2d_add(stats_diam_aro2_dist, diam_aro2_dist)        
!for alk1
        alk1_masses = aero_state_masses(aero_state, aero_data, &
             include=(/"ALK1"/))
        tot_alk1_mass_conc = sum(alk1_masses * num_concs)
        call stats_1d_add_entry(stats_tot_alk1_mass_conc, tot_alk1_mass_conc, i_index)
    
        alk1_fracs = alk1_masses / dry_masses
        diam_alk1_dist = bin_grid_histogram_2d(diam_grid, dry_diameters, &
             alk1_grid, alk1_fracs, num_concs)
        call stats_2d_add(stats_diam_alk1_dist, diam_alk1_dist)  
!for ole1
        ole1_masses = aero_state_masses(aero_state, aero_data, &
             include=(/"OLE1"/))
        tot_ole1_mass_conc = sum(ole1_masses * num_concs)
        call stats_1d_add_entry(stats_tot_ole1_mass_conc, tot_ole1_mass_conc, i_index)
    
        ole1_fracs = ole1_masses / dry_masses
        diam_ole1_dist = bin_grid_histogram_2d(diam_grid, dry_diameters, &
             ole1_grid, ole1_fracs, num_concs)
        call stats_2d_add(stats_diam_ole1_dist, diam_ole1_dist)  
!for ole1
        ole1_masses = aero_state_masses(aero_state, aero_data, &
             include=(/"OLE1"/))
        tot_ole1_mass_conc = sum(ole1_masses * num_concs)
        call stats_1d_add_entry(stats_tot_ole1_mass_conc, tot_ole1_mass_conc, i_index)
    
        ole1_fracs = ole1_masses / dry_masses
        diam_ole1_dist = bin_grid_histogram_2d(diam_grid, dry_diameters, &
             ole1_grid, ole1_fracs, num_concs)
        call stats_2d_add(stats_diam_ole1_dist, diam_ole1_dist) 
!for api1
        api1_masses = aero_state_masses(aero_state, aero_data, &
             include=(/"API1"/))
        tot_api1_mass_conc = sum(api1_masses * num_concs)
        call stats_1d_add_entry(stats_tot_api1_mass_conc, tot_api1_mass_conc, i_index)
    
        api1_fracs = api1_masses / dry_masses
        diam_api1_dist = bin_grid_histogram_2d(diam_grid, dry_diameters, &
             api1_grid, api1_fracs, num_concs)
        call stats_2d_add(stats_diam_api1_dist, diam_api1_dist) 
!for api2
        api2_masses = aero_state_masses(aero_state, aero_data, &
             include=(/"API2"/))
        tot_api2_mass_conc = sum(api2_masses * num_concs)
        call stats_1d_add_entry(stats_tot_api2_mass_conc, tot_api2_mass_conc, i_index)
    
        api2_fracs = api2_masses / dry_masses
        diam_api2_dist = bin_grid_histogram_2d(diam_grid, dry_diameters, &
             api2_grid, api2_fracs, num_concs)
        call stats_2d_add(stats_diam_api2_dist, diam_api2_dist)      
!for lim1
        lim1_masses = aero_state_masses(aero_state, aero_data, &
             include=(/"LIM1"/))
        tot_lim1_mass_conc = sum(lim1_masses * num_concs)
        call stats_1d_add_entry(stats_tot_lim1_mass_conc, tot_lim1_mass_conc, i_index)
    
        lim1_fracs = lim1_masses / dry_masses
        diam_lim1_dist = bin_grid_histogram_2d(diam_grid, dry_diameters, &
             lim1_grid, lim1_fracs, num_concs)
        call stats_2d_add(stats_diam_lim1_dist, diam_lim1_dist) 
!for lim2
        lim2_masses = aero_state_masses(aero_state, aero_data, &
             include=(/"LIM2"/))
        tot_lim2_mass_conc = sum(lim2_masses * num_concs)
        call stats_1d_add_entry(stats_tot_lim2_mass_conc, tot_lim2_mass_conc, i_index)
    
        lim2_fracs = lim2_masses / dry_masses
        diam_lim2_dist = bin_grid_histogram_2d(diam_grid, dry_diameters, &
             lim2_grid, lim2_fracs, num_concs)
        call stats_2d_add(stats_diam_lim2_dist, diam_lim2_dist) 
!for co3
        co3_masses = aero_state_masses(aero_state, aero_data, &
             include=(/"CO3"/))
        tot_co3_mass_conc = sum(co3_masses * num_concs)
        call stats_1d_add_entry(stats_tot_co3_mass_conc, tot_co3_mass_conc, i_index)
    
        co3_fracs = co3_masses / dry_masses
        diam_co3_dist = bin_grid_histogram_2d(diam_grid, dry_diameters, &
             co3_grid, co3_fracs, num_concs)
        call stats_2d_add(stats_diam_co3_dist, diam_co3_dist) 
!for na
        na_masses = aero_state_masses(aero_state, aero_data, &
             include=(/"Na"/))
        tot_na_mass_conc = sum(na_masses * num_concs)
        call stats_1d_add_entry(stats_tot_na_mass_conc, tot_na_mass_conc, i_index)
    
        na_fracs = na_masses / dry_masses
        diam_na_dist = bin_grid_histogram_2d(diam_grid, dry_diameters, &
             na_grid, na_fracs, num_concs)
        call stats_2d_add(stats_diam_na_dist, diam_na_dist)       
!for ca
        ca_masses = aero_state_masses(aero_state, aero_data, &
             include=(/"Ca"/))
        tot_ca_mass_conc = sum(ca_masses * num_concs)
        call stats_1d_add_entry(stats_tot_ca_mass_conc, tot_ca_mass_conc, i_index)
    
        ca_fracs = ca_masses / dry_masses
        diam_ca_dist = bin_grid_histogram_2d(diam_grid, dry_diameters, &
             ca_grid, ca_fracs, num_concs)
        call stats_2d_add(stats_diam_ca_dist, diam_ca_dist) 
!for oin
        oin_masses = aero_state_masses(aero_state, aero_data, &
             include=(/"OIN"/))
        tot_oin_mass_conc = sum(oin_masses * num_concs)
        call stats_1d_add_entry(stats_tot_oin_mass_conc, tot_oin_mass_conc, i_index)
    
        oin_fracs = oin_masses / dry_masses
        diam_oin_dist = bin_grid_histogram_2d(diam_grid, dry_diameters, &
             oin_grid, oin_fracs, num_concs)
        call stats_2d_add(stats_diam_oin_dist, diam_oin_dist) 
!for h2o
        h2o_masses = aero_state_masses(aero_state, aero_data, &
             include=(/"H2O"/))
        tot_h2o_mass_conc = sum(h2o_masses * num_concs)
        call stats_1d_add_entry(stats_tot_h2o_mass_conc, tot_h2o_mass_conc, i_index)
    
        h2o_fracs = h2o_masses / dry_masses
        diam_h2o_dist = bin_grid_histogram_2d(diam_grid, dry_diameters, &
             h2o_grid, h2o_fracs, num_concs)
        call stats_2d_add(stats_diam_h2o_dist, diam_h2o_dist) 

!for critical RH supersaturation     
        crit_rhs = aero_state_crit_rel_humids(aero_state, aero_data, &
             env_state)
        scs = crit_rhs - 1d0
        diam_sc_dist = bin_grid_histogram_2d(diam_grid, dry_diameters, &
             sc_grid, scs, num_concs)
        sc_dist = bin_grid_histogram_1d(sc_grid, scs, num_concs)
        call stats_1d_add(stats_sc_dist, sc_dist)
        call stats_2d_add(stats_diam_sc_dist, diam_sc_dist)

!for critical RH supersaturation with varying sigma
        crit_rhs_varying_sigma = aero_state_crit_rel_humids_varying_sigma( &
             aero_state, aero_data, env_state)
        scs_varying_sigma = crit_rhs_varying_sigma - 1d0
        diam_sc_varying_sigma_dist = bin_grid_histogram_2d(diam_grid, dry_diameters, &
             sc_varying_sigma_grid, scs_varying_sigma, num_concs)
        sc_varying_sigma_dist = bin_grid_histogram_1d(sc_varying_sigma_grid, & 
             scs_varying_sigma, num_concs)
        call stats_1d_add(stats_sc_varying_sigma_dist, sc_varying_sigma_dist)
        call stats_2d_add(stats_diam_sc_varying_sigma_dist, diam_sc_varying_sigma_dist)

        call for_single_particle(single_particle_output, prefix_new, ".csv", index)
        open(15, file= trim(single_particle_output))
        write(15, *) 'dry_diameters', dry_diameters
        write(15, *) 'masses', masses
        write(15, *) 'dry_masses', dry_masses
        write(15, *) 'bc_masses', bc_masses
        write(15, *) 'oc_masses', oc_masses
        write(15, *) 'so4_masses', so4_masses
        write(15, *) 'no3_masses', no3_masses
        write(15, *) 'nh4_masses', nh4_masses
        write(15, *) 'soa_masses', soa_masses
        write(15, *) 'cl_masses', cl_masses
        write(15, *) 'msa_masses', msa_masses
        write(15, *) 'aro1_masses', aro1_masses
        write(15, *) 'aro2_masses', aro2_masses
        write(15, *) 'alk1_masses', alk1_masses
        write(15, *) 'ole1_masses', ole1_masses
        write(15, *) 'api1_masses', api1_masses
        write(15, *) 'api2_masses', api2_masses
        write(15, *) 'lim1_masses', lim1_masses
        write(15, *) 'lim2_masses', lim2_masses
        write(15, *) 'co3_masses', co3_masses
        write(15, *) 'na_masses', na_masses
        write(15, *) 'ca_masses', ca_masses
        write(15, *) 'oin_masses', oin_masses
        write(15, *) 'h2o_masses', h2o_masses
        write(15, *) 'crit_rhs', crit_rhs
        write(15, *) 'scs', scs
        write(15, *) 'num_con', num_concs
        write(15, *) 'crit_rhs_varying_sigma', crit_rhs_varying_sigma
        write(15, *) 'scs_varying_sigma', scs_varying_sigma
        close(15)
        
         call aero_state_mixing_state_metrics(aero_state, aero_data, &
             d_alpha, d_gamma, chi, groups=mixing_state_groups)

        call stats_1d_add_entry(stats_d_alpha, d_alpha, i_index)
        call stats_1d_add_entry(stats_d_gamma, d_gamma, i_index)
        call stats_1d_add_entry(stats_chi, chi, i_index)

     end do

     call make_filename(out_filename, prefix, "_process.nc", index)
     write(*,*) "Writing " // trim(out_filename)
     call pmc_nc_open_write(out_filename, ncid)
     call pmc_nc_write_info(ncid, uuid, "1_urban_plume process")
     call bin_grid_output_netcdf(diam_grid, ncid, "diam", unit="m")
     call bin_grid_output_netcdf(bc_grid, ncid, "bc_frac", unit="1")
     call bin_grid_output_netcdf(oc_grid, ncid, "oc_frac", unit="1") 
     call bin_grid_output_netcdf(so4_grid, ncid, "so4_frac", unit="1") 
     call bin_grid_output_netcdf(no3_grid, ncid, "no3_frac", unit="1") 
     call bin_grid_output_netcdf(nh4_grid, ncid, "nh4_frac", unit="1")
     call bin_grid_output_netcdf(soa_grid, ncid, "soa_frac", unit="1")
     call bin_grid_output_netcdf(cl_grid, ncid, "cl_frac", unit="1")
     call bin_grid_output_netcdf(msa_grid, ncid, "msa_frac", unit="1")
     call bin_grid_output_netcdf(aro1_grid, ncid, "aro1_frac", unit="1")
     call bin_grid_output_netcdf(aro2_grid, ncid, "aro2_frac", unit="1")
     call bin_grid_output_netcdf(alk1_grid, ncid, "alk1_frac", unit="1")
     call bin_grid_output_netcdf(ole1_grid, ncid, "ole1_frac", unit="1")
     call bin_grid_output_netcdf(api1_grid, ncid, "api1_frac", unit="1")
     call bin_grid_output_netcdf(api2_grid, ncid, "api2_frac", unit="1")
     call bin_grid_output_netcdf(lim1_grid, ncid, "lim1_frac", unit="1")
     call bin_grid_output_netcdf(lim2_grid, ncid, "lim2_frac", unit="1")
     call bin_grid_output_netcdf(co3_grid, ncid, "co3_frac", unit="1")
     call bin_grid_output_netcdf(na_grid, ncid, "na_frac", unit="1")
     call bin_grid_output_netcdf(ca_grid, ncid, "ca_frac", unit="1")
     call bin_grid_output_netcdf(oin_grid, ncid, "oin_frac", unit="1")
     call bin_grid_output_netcdf(h2o_grid, ncid, "h2o_frac", unit="1")
     call bin_grid_output_netcdf(sc_grid, ncid, "sc", unit="1")
     call bin_grid_output_netcdf(sc_varying_sigma_grid, ncid, & 
                                 "sc_varying_sigma", unit="1")

     call stats_1d_output_netcdf(stats_num_dist, ncid, "num_dist", &
          dim_name="diam", unit="m^{-3}")
     call stats_1d_clear(stats_num_dist)

     call stats_2d_output_netcdf(stats_diam_bc_dist, ncid, "diam_bc_dist", &
          dim_name_1="diam", dim_name_2="bc_frac", unit="m^{-3}")
     call stats_2d_clear(stats_diam_bc_dist)

     call stats_2d_output_netcdf(stats_diam_oc_dist, ncid, "diam_oc_dist", &
          dim_name_1="diam", dim_name_2="oc_frac", unit="m^{-3}")
     call stats_2d_clear(stats_diam_oc_dist)

     call stats_2d_output_netcdf(stats_diam_so4_dist, ncid, "diam_so4_dist", &
          dim_name_1="diam", dim_name_2="so4_frac", unit="m^{-3}")
     call stats_2d_clear(stats_diam_so4_dist)

     call stats_2d_output_netcdf(stats_diam_no3_dist, ncid, "diam_no3_dist", &
          dim_name_1="diam", dim_name_2="no3_frac", unit="m^{-3}")
     call stats_2d_clear(stats_diam_no3_dist)

     call stats_2d_output_netcdf(stats_diam_nh4_dist, ncid, "diam_nh4_dist", &
     dim_name_1="diam", dim_name_2="nh4_frac", unit="m^{-3}")
     call stats_2d_clear(stats_diam_nh4_dist)

     call stats_2d_output_netcdf(stats_diam_soa_dist, ncid, "diam_soa_dist", &
     dim_name_1="diam", dim_name_2="soa_frac", unit="m^{-3}")
     call stats_2d_clear(stats_diam_soa_dist)

     call stats_2d_output_netcdf(stats_diam_cl_dist, ncid, "diam_cl_dist", &
     dim_name_1="diam", dim_name_2="cl_frac", unit="m^{-3}")
     call stats_2d_clear(stats_diam_cl_dist)

     call stats_2d_output_netcdf(stats_diam_msa_dist, ncid, "diam_msa_dist", &
     dim_name_1="diam", dim_name_2="msa_frac", unit="m^{-3}")
     call stats_2d_clear(stats_diam_cl_dist)

     call stats_2d_output_netcdf(stats_diam_aro1_dist, ncid, "diam_aro1_dist", &
     dim_name_1="diam", dim_name_2="aro1_frac", unit="m^{-3}")
     call stats_2d_clear(stats_diam_aro1_dist)

     call stats_2d_output_netcdf(stats_diam_aro2_dist, ncid, "diam_aro2_dist", &
     dim_name_1="diam", dim_name_2="aro2_frac", unit="m^{-3}")
     call stats_2d_clear(stats_diam_aro2_dist)
     
     call stats_2d_output_netcdf(stats_diam_alk1_dist, ncid, "diam_alk1_dist", &
     dim_name_1="diam", dim_name_2="alk1_frac", unit="m^{-3}")
     call stats_2d_clear(stats_diam_alk1_dist)

     call stats_2d_output_netcdf(stats_diam_ole1_dist, ncid, "diam_ole1_dist", &
     dim_name_1="diam", dim_name_2="ole1_frac", unit="m^{-3}")
     call stats_2d_clear(stats_diam_ole1_dist)

     call stats_2d_output_netcdf(stats_diam_api1_dist, ncid, "diam_api1_dist", &
     dim_name_1="diam", dim_name_2="api1_frac", unit="m^{-3}")
     call stats_2d_clear(stats_diam_api1_dist)

     call stats_2d_output_netcdf(stats_diam_api2_dist, ncid, "diam_api2_dist", &
     dim_name_1="diam", dim_name_2="api2_frac", unit="m^{-3}")
     call stats_2d_clear(stats_diam_api2_dist)

     call stats_2d_output_netcdf(stats_diam_lim1_dist, ncid, "diam_lim1_dist", &
     dim_name_1="diam", dim_name_2="lim1_frac", unit="m^{-3}")
     call stats_2d_clear(stats_diam_lim1_dist)

     call stats_2d_output_netcdf(stats_diam_lim2_dist, ncid, "diam_lim2_dist", &
     dim_name_1="diam", dim_name_2="lim2_frac", unit="m^{-3}")
     call stats_2d_clear(stats_diam_lim2_dist)
     
     call stats_2d_output_netcdf(stats_diam_co3_dist, ncid, "diam_co3_dist", &
     dim_name_1="diam", dim_name_2="co3_frac", unit="m^{-3}")
     call stats_2d_clear(stats_diam_co3_dist)

     call stats_2d_output_netcdf(stats_diam_na_dist, ncid, "diam_na_dist", &
     dim_name_1="diam", dim_name_2="na_frac", unit="m^{-3}")
     call stats_2d_clear(stats_diam_na_dist)

     call stats_2d_output_netcdf(stats_diam_ca_dist, ncid, "diam_ca_dist", &
     dim_name_1="diam", dim_name_2="ca_frac", unit="m^{-3}")
     call stats_2d_clear(stats_diam_ca_dist)

     call stats_2d_output_netcdf(stats_diam_oin_dist, ncid, "diam_oin_dist", &
     dim_name_1="diam", dim_name_2="oin_frac", unit="m^{-3}")
     call stats_2d_clear(stats_diam_oin_dist)

     call stats_2d_output_netcdf(stats_diam_h2o_dist, ncid, "diam_h2o_dist", &
     dim_name_1="diam", dim_name_2="h2o_frac", unit="m^{-3}")
     call stats_2d_clear(stats_diam_h2o_dist)

     call stats_2d_output_netcdf(stats_diam_sc_dist, ncid, "diam_sc_dist", &
     dim_name_1="diam", dim_name_2="sc", unit="m^{-3}")
     call stats_2d_clear(stats_diam_sc_dist)

     call stats_2d_output_netcdf(stats_diam_sc_varying_sigma_dist, ncid, & 
                                 "diam_sc_varying_sigma_dist", &
     dim_name_1="diam", dim_name_2="sc_varying_sigma", unit="m^{-3}")
     call stats_2d_clear(stats_diam_sc_varying_sigma_dist)

     call stats_1d_output_netcdf(stats_sc_dist, ncid, "sc_dist", &
          dim_name="sc", unit="C")
     call stats_1d_clear(stats_sc_dist)

     call stats_1d_output_netcdf(stats_sc_varying_sigma_dist, ncid, "sc_varying_sigma_dist", &
          dim_name="sc", unit="C")
     call stats_1d_clear(stats_sc_varying_sigma_dist)

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
  call stats_1d_output_netcdf(stats_tot_dry_mass_conc, ncid, "tot_dry_mass_conc",&
       dim_name="time", unit="kg m^{-3}")       
  call stats_1d_output_netcdf(stats_tot_bc_mass_conc, ncid, "tot_bc_mass_conc", &
       dim_name="time", unit="kg m^{-3}")
  call stats_1d_output_netcdf(stats_tot_oc_mass_conc, ncid, "tot_oc_mass_conc", &
       dim_name="time", unit="kg m^{-3}")
  call stats_1d_output_netcdf(stats_tot_so4_mass_conc, ncid, "tot_so4_mass_conc", &
       dim_name="time", unit="kg m^{-3}")
  call stats_1d_output_netcdf(stats_tot_no3_mass_conc, ncid, "tot_no3_mass_conc", &
       dim_name="time", unit="kg m^{-3}")
  call stats_1d_output_netcdf(stats_tot_nh4_mass_conc, ncid, "tot_nh4_mass_conc", &
       dim_name="time", unit="kg m^{-3}")
  call stats_1d_output_netcdf(stats_tot_soa_mass_conc, ncid, "tot_soa_mass_conc", &
       dim_name="time", unit="kg m^{-3}")
  call stats_1d_output_netcdf(stats_tot_cl_mass_conc, ncid, "tot_cl_mass_conc", &
       dim_name="time", unit="kg m^{-3}")
  call stats_1d_output_netcdf(stats_tot_msa_mass_conc, ncid, "tot_msa_mass_conc", &
       dim_name="time", unit="kg m^{-3}")
  call stats_1d_output_netcdf(stats_tot_aro1_mass_conc, ncid, "tot_aro1_mass_conc", &
       dim_name="time", unit="kg m^{-3}")
  call stats_1d_output_netcdf(stats_tot_aro2_mass_conc, ncid, "tot_aro2_mass_conc", &
       dim_name="time", unit="kg m^{-3}")
  call stats_1d_output_netcdf(stats_tot_alk1_mass_conc, ncid, "tot_alk1_mass_conc", &
       dim_name="time", unit="kg m^{-3}")
  call stats_1d_output_netcdf(stats_tot_ole1_mass_conc, ncid, "tot_ole1_mass_conc", &
       dim_name="time", unit="kg m^{-3}")
  call stats_1d_output_netcdf(stats_tot_api1_mass_conc, ncid, "tot_api1_mass_conc", &
       dim_name="time", unit="kg m^{-3}")      
  call stats_1d_output_netcdf(stats_tot_api2_mass_conc, ncid, "tot_api2_mass_conc", &
       dim_name="time", unit="kg m^{-3}")          
  call stats_1d_output_netcdf(stats_tot_lim1_mass_conc, ncid, "tot_lim1_mass_conc", &
       dim_name="time", unit="kg m^{-3}")   
  call stats_1d_output_netcdf(stats_tot_lim2_mass_conc, ncid, "tot_lim2_mass_conc", &
       dim_name="time", unit="kg m^{-3}") 
  call stats_1d_output_netcdf(stats_tot_co3_mass_conc, ncid, "tot_co3_mass_conc", &
       dim_name="time", unit="kg m^{-3}")     
  call stats_1d_output_netcdf(stats_tot_na_mass_conc, ncid, "tot_na_mass_conc", &
       dim_name="time", unit="kg m^{-3}")   
  call stats_1d_output_netcdf(stats_tot_ca_mass_conc, ncid, "tot_ca_mass_conc", &
       dim_name="time", unit="kg m^{-3}")   
  call stats_1d_output_netcdf(stats_tot_oin_mass_conc, ncid, "tot_oin_mass_conc", &
       dim_name="time", unit="kg m^{-3}") 
  call stats_1d_output_netcdf(stats_tot_h2o_mass_conc, ncid, "tot_h2o_mass_conc", &
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
