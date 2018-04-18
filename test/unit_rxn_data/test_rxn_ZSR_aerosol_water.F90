! Copyright (C) 2017 Matt Dawson
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_test_ZSR_aerosol_water program

!> Test of ZSR_aerosol_water reaction module
program pmc_test_ZSR_aerosol_water

  use iso_c_binding

  use pmc_util,                         only: i_kind, dp, assert, &
                                              almost_equal, string_t, &
                                              warn_msg
  use pmc_phlex_core
  use pmc_phlex_state
  use pmc_aero_rep_factory
  use pmc_aero_rep_single_particle
#ifdef PMC_USE_JSON
  use json_module
#endif
  use pmc_mpi

  implicit none

  ! Number of RHs to calculate aerosol water for
  integer(kind=i_kind) :: NUM_RH_STEP = 100

  ! initialize mpi
  call pmc_mpi_init()

  if (run_ZSR_aerosol_water_tests()) then
    write(*,*) "ZSR aerosol water reaction tests - PASS"
  else
    write(*,*) "ZSR aerosol water reaction tests - FAIL"
  end if

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Run all pmc_chem_mech_solver tests
  logical function run_ZSR_aerosol_water_tests() result(passed)

    use pmc_phlex_solver_data

    type(phlex_solver_data_t), pointer :: phlex_solver_data

    phlex_solver_data => phlex_solver_data_t()

    if (phlex_solver_data%is_solver_available()) then
      passed = run_ZSR_aerosol_water_test()
    else
      call warn_msg(713064651, "No solver available")
      passed = .true.
    end if

  end function run_ZSR_aerosol_water_tests

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Solve a mechanism consisting of one aersol water calculation with two ion pairs
  !!
  !! JACOBSON molality parameters are for NaCl from Jacobson et al. 
  !! \cite{Jacobson1996} Table 2. (CaCl2 is used in the test just to test code
  !! when the number of anions and cations differs.) EQSAM molality parameters
  !! are also for NaCl from EQSAM_v03d.
  logical function run_ZSR_aerosol_water_test()

    use pmc_constants

    type(phlex_core_t), pointer :: phlex_core
    type(phlex_state_t), target :: phlex_state
    character(len=:), allocatable :: input_file_path
    type(string_t), allocatable, dimension(:) :: output_file_path

    real(kind=dp), dimension(0:NUM_RH_STEP, 25) :: model_conc, true_conc
    integer(kind=i_kind) :: idx_phase, idx_aero_rep
    integer(kind=i_kind) :: idx_H2O, idx_Na_p, idx_Na_p_act, idx_Cl_m, idx_Cl_m_act, &
           idx_Ca_pp, idx_Ca_pp_act, idx_H2O_aq,  idx_H2O_act
    character(len=:), allocatable :: key
    integer(kind=i_kind) :: i_RH, i_spec
    real(kind=dp) :: RH_step, RH
    real(kind=dp) :: ppm_to_RH
    real(kind=dp) :: molal_NaCl, molal_CaCl2, NaCl_conc, CaCl2_conc, water_NaCl, water_CaCl2

    ! Parameters for calculating true concentrations
    real(kind=dp) :: temp, pressure 

    run_ZSR_aerosol_water_test = .true.

    ! Set the environmental and aerosol test conditions
    temp = 272.5d0              ! temperature (K)
    pressure = 101253.3d0       ! pressure (Pa)

    ! Set RH step (unitless)
    RH_step = 1.0/NUM_RH_STEP

    ! Get the ZSR_aerosol_water reaction mechanism json file
    input_file_path = 'test_ZSR_aerosol_water_config.json'

    ! Construct a phlex_core variable
    phlex_core => phlex_core_t(input_file_path)

    ! Initialize the model
    call phlex_core%initialize()

    ! Initialize the solver
    call phlex_core%solver_initialize()

    ! Get a model state variable
    phlex_state = phlex_core%new_state()

    ! Set the environmental conditions
    phlex_state%env_state%temp = temp
    phlex_state%env_state%pressure = pressure
    call phlex_state%update_env_state()

    ! Find the aerosol representation
    call assert(110830690, size(phlex_core%aero_rep).eq.3)
    idx_aero_rep = 2

    ! Get species indices
    key = "H2O"
    idx_H2O = phlex_core%chem_spec_data%gas_state_id(key);
    key = "aqueous aerosol.H2O_aq"
    idx_H2O_aq = phlex_core%aero_rep(idx_aero_rep)%val%spec_state_id(key);
    idx_H2O_act = phlex_core%aero_rep(idx_aero_rep)%val%activity_coeff_state_id(key);
    key = "aqueous aerosol.Na_p"
    idx_Na_p = phlex_core%aero_rep(idx_aero_rep)%val%spec_state_id(key);
    idx_Na_p_act = phlex_core%aero_rep(idx_aero_rep)%val%activity_coeff_state_id(key);
    key = "aqueous aerosol.Cl_m"
    idx_Cl_m = phlex_core%aero_rep(idx_aero_rep)%val%spec_state_id(key);
    idx_Cl_m_act = phlex_core%aero_rep(idx_aero_rep)%val%activity_coeff_state_id(key);
    key = "aqueous aerosol.Ca_pp"
    idx_Ca_pp = phlex_core%aero_rep(idx_aero_rep)%val%spec_state_id(key);
    idx_Ca_pp_act = phlex_core%aero_rep(idx_aero_rep)%val%activity_coeff_state_id(key);

    ! Make sure the expected species are in the model
    call assert(213525011, idx_H2O.gt.0)
    call assert(943368106, idx_H2O_aq.gt.0)
    call assert(155686452, idx_H2O_act.gt.0)
    call assert(202996397, idx_Na_p.gt.0)
    call assert(315314742, idx_Na_p_act.gt.0)
    call assert(427633087, idx_Cl_m.gt.0)
    call assert(487377180, idx_Cl_m_act.gt.0)
    call assert(317220276, idx_Ca_pp.gt.0)
    call assert(147063372, idx_Ca_pp_act.gt.0)

    ! Save the initial concentrations
    true_conc(:,:) = 0.0
    true_conc(:,idx_H2O) = 0.0
    true_conc(:,idx_H2O_aq) = 0.0
    true_conc(:,idx_H2O_act) = 1.0
    true_conc(:,idx_Na_p) = 2.5
    true_conc(:,idx_Na_p_act) = 1.0
    true_conc(:,idx_Cl_m) = 5.3
    true_conc(:,idx_Cl_m_act) = 0.0
    true_conc(:,idx_Ca_pp) = 1.3
    true_conc(:,idx_Ca_pp_act) = 0.0
    model_conc(:,:) = true_conc(:,:)

    ! Set up the ppm->RH (0-1) conversion
    ! (From MOSAIC code, references Seinfeld and Pandis pg. 181)
    ppm_to_RH = 1.0d0 - 373.15d0/temp
    ppm_to_RH = (((-0.1299d0*ppm_to_RH - 0.6445d0)*ppm_to_RH - 1.976d0)*ppm_to_RH &
            + 13.3185d0)*ppm_to_RH
    ppm_to_RH = exp(ppm_to_RH)  ! VP of water (atm)
    ppm_to_RH = (pressure/101325.0d0) / ppm_to_RH * 1.0d-6 ! ppm -> RH (0-1)

    ! Integrate the mechanism
    do i_RH = 1, NUM_RH_STEP

      ! Calculate the current [H20] (ppm)
      true_conc(i_RH,idx_H2O) = i_RH * RH_step / ppm_to_RH
      model_conc(i_RH,idx_H2O) = true_conc(i_RH, idx_H2O)

      ! Set the initial state in the model
      phlex_state%state_var(:) = model_conc(i_RH,:)

      ! Get the modeled conc
      call phlex_core%solve(phlex_state, real(1.0, kind=dp)) ! time step is arbitrary - equilibrium calculatuions only
      model_conc(i_RH,:) = phlex_state%state_var(:)

      ! Get the analytic conc
      ! Jacobson molality (eq. 29 in \cite{Jacobson1996}) :
      !  m_a^(1/2) = Y0 + Y1*aw + Y2*aw^2 + ... (aw = RH)
      RH = i_RH * RH_step
      RH = MAX(RH, 0.43d0)
      molal_CaCl2 = -1.918004e2 +  2.001540e3 * RH - 8.557205e3 * RH**2 + 1.987670e4 * RH**3 &
              - 2.717192e4 * RH**4 + 2.187103e4 * RH**5 - 9.591577e3 * RH**6 + 1.763672e3 * RH**7
      molal_CaCl2 = molal_CaCl2**2
      ! EQSAM molality (from EQSAM_v03d)
      !  m_i = (NW_i * MW_H2O/MW_i * (1.0/RH-1.0))^ZW_i   where MW_H2O is defined as 55.51*18.01
      RH = i_RH * RH_step
      molal_NaCl = (2.0d0 * 55.51D0 * 18.01d0 / 58.5d0 * (1.0d0/RH - 1.0d0))**0.67d0
      CaCl2_conc = MIN(true_conc(i_RH,idx_Ca_pp)/40.078d0, true_conc(i_RH,idx_Cl_m)/2.0d0/35.453d0) ! (umol/m^3_air = mol/cm^3_air)
      NaCl_conc = true_conc(i_RH,idx_Cl_m)/35.453d0 ! (umol/m^3_air = mol/cm^3_air)
      ! Water content is (eq 28 \cite{Jacobson1996}) :
      ! cw = 1000 / MW_H2O * sum_i (c_i / m_i)   with cw and c_i in (mol_i/cm^3_air) and m_i in (mol_i/kg_H2O)
      water_CaCl2 =  CaCl2_conc / molal_CaCl2 * 1000.0d0 ! (ug_H2O/m^3_air)
      water_NaCl = NaCl_conc / molal_NaCl * 1000.0d0 ! (ug_H2O/m^3_air)
      true_conc(i_RH,idx_H2O_aq) = water_CaCl2 + water_NaCl

    end do

    ! Save the results
    open(unit=7, file="out/ZSR_aerosol_water_results.txt", status="replace", action="write")
    do i_RH = 0, NUM_RH_STEP
      write(7,*) i_RH*RH_step, &
            ' ', true_conc(i_RH, idx_H2O),' ', model_conc(i_RH, idx_H2O), &
            ' ', true_conc(i_RH, idx_H2O_aq),' ', model_conc(i_RH, idx_H2O_aq), &
            ' ', true_conc(i_RH, idx_Na_p),' ', model_conc(i_RH, idx_Na_p), &
            ' ', true_conc(i_RH, idx_Cl_m),' ', model_conc(i_RH, idx_Cl_m), &
            ' ', true_conc(i_RH, idx_Ca_pp),' ', model_conc(i_RH, idx_Ca_pp)
    end do
    close(7)

    ! Analyze the results
    do i_RH = 1, NUM_RH_STEP
      do i_spec = 1, 17
        ! Skip the first aerosol phase
        if (i_spec.ge.2.and.i_spec.le.9) cycle
        call assert_msg(106356995, &
          almost_equal(model_conc(i_RH, i_spec), true_conc(i_RH, i_spec), &
          real(1.0e-2, kind=dp)), "RH: "//to_string(i_RH)//"; species: "// &
          to_string(i_spec)//"; mod: "//to_string(model_conc(i_RH, i_spec))// &
          "; true: "//to_string(true_conc(i_RH, i_spec)))
      end do
    end do

    run_ZSR_aerosol_water_test = .true.

  end function run_ZSR_aerosol_water_test

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program pmc_test_ZSR_aerosol_water
