! Copyright (C) 2017 Matt Dawson
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_test_PDFiTE_activity program

!> Test of PDFiTE_activity reaction module
program pmc_test_PDFiTE_activity

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

  ! Number of timesteps to output in mechanisms
  integer(kind=i_kind) :: NUM_RH_STEP = 101

  ! initialize mpi
  call pmc_mpi_init()

  if (run_PDFiTE_activity_tests()) then
    write(*,*) " PD-FiTE activity reaction tests - PASS"
  else
    write(*,*) " PD-FiTE activity reaction tests - FAIL"
  end if

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Run all pmc_chem_mech_solver tests
  logical function run_PDFiTE_activity_tests() result(passed)

    use pmc_phlex_solver_data

    type(phlex_solver_data_t), pointer :: phlex_solver_data

    phlex_solver_data => phlex_solver_data_t()

    if (phlex_solver_data%is_solver_available()) then
      passed = run_PDFiTE_activity_test()
    else
      call warn_msg(713064651, "No solver available")
      passed = .true.
    end if

  end function run_PDFiTE_activity_tests

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Solve a PD-FiTE activity coefficient calculation 
  !!
  !! The scenario matches that described for the H+ - NH4+ - SO42- - NO3-
  !! system described by equations 16 and 17 in \cite{Topping2009}.
  logical function run_PDFiTE_activity_test()

    use pmc_constants

    type(phlex_core_t), pointer :: phlex_core
    type(phlex_state_t), target :: phlex_state
    character(len=:), allocatable :: input_file_path
    type(string_t), allocatable, dimension(:) :: output_file_path

    real(kind=dp), dimension(0:NUM_RH_STEP, 28) :: model_conc, true_conc
    integer(kind=i_kind) :: idx_H2O, idx_H2O_aq, idx_H_p, idx_NH4_p, idx_SO4_mm, &
            idx_NO3_m, idx_NH42_SO4, idx_NH4_NO3, idx_H_NO3, idx_H2_SO4 
    integer(kind=i_kind) :: idx_phase, idx_aero_rep
    character(len=:), allocatable :: key
    integer(kind=i_kind) :: i_RH, i_spec
    real(kind=dp) :: time_step, time, ppm_to_RH, omega, ln_gamma
    real(kind=dp) :: HNO3_LRH_B0, HNO3_LRH_B1, HNO3_LRH_B2, HNO3_LRH_B3, HNO3_LRH_B4
    real(kind=dp) :: HNO3_HRH_B0, HNO3_HRH_B1, HNO3_HRH_B2, HNO3_HRH_B3
    real(kind=dp) :: H2SO4_B0, H2SO4_B1, H2SO4_B2, H2SO4_B3
    real(kind=dp) :: NH42SO4_B0, NH42SO4_B1, NH42SO4_B2, NH42SO4_B3
    real(kind=dp) :: NH4NO3_B0, NH4NO3_B1, NH4NO3_B2, NH4NO3_B3

    ! Parameters for calculating true concentrations
    real(kind=dp) :: temp, pressure, a_w 

    run_PDFiTE_activity_test = .true.

    ! Set the environmental and aerosol test conditions
    temp = 272.5d0              ! temperature (K)
    pressure = 101253.3d0       ! pressure (Pa)

    ! Set output time step (s)
    time_step = 1.0d0

    ! Get the PDFiTE_activity reaction mechanism json file
    input_file_path = 'test_PDFiTE_activity_config.json'

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
    call assert(917299782, size(phlex_core%aero_rep).eq.3)
    idx_aero_rep = 2

    ! Get species indices
    key = "H2O"
    idx_H2O = phlex_core%chem_spec_data%gas_state_id(key);
    key = "aqueous aerosol.H2O_aq"
    idx_H2O_aq = phlex_core%aero_rep(idx_aero_rep)%val%spec_state_id(key);
    key = "aqueous aerosol.H_p"
    idx_H_p = phlex_core%aero_rep(idx_aero_rep)%val%spec_state_id(key);
    key = "aqueous aerosol.NH4_p"
    idx_NH4_p = phlex_core%aero_rep(idx_aero_rep)%val%spec_state_id(key);
    key = "aqueous aerosol.SO4_mm"
    idx_SO4_mm = phlex_core%aero_rep(idx_aero_rep)%val%spec_state_id(key);
    key = "aqueous aerosol.NO3_m"
    idx_NO3_m = phlex_core%aero_rep(idx_aero_rep)%val%spec_state_id(key);
    key = "aqueous aerosol.(NH4)2-SO4"
    idx_NH42_SO4 = phlex_core%aero_rep(idx_aero_rep)%val%spec_state_id(key);
    key = "aqueous aerosol.NH4-NO3"
    idx_NH4_NO3 = phlex_core%aero_rep(idx_aero_rep)%val%spec_state_id(key);
    key = "aqueous aerosol.H-NO3"
    idx_H_NO3 = phlex_core%aero_rep(idx_aero_rep)%val%spec_state_id(key);
    key = "aqueous aerosol.H2-SO4"
    idx_H2_SO4 = phlex_core%aero_rep(idx_aero_rep)%val%spec_state_id(key);

    ! Make sure the expected species are in the model
    call assert(447873764, idx_H2O.gt.0)
    call assert(225142608, idx_H2O_aq.gt.0)
    call assert(221382734, idx_H_p.gt.0)
    call assert(163544175, idx_NH4_p.gt.0)
    call assert(210854120, idx_SO4_mm.gt.0)
    call assert(388180865, idx_NO3_m.gt.0)
    call assert(218023961, idx_NH42_SO4.gt.0)
    call assert(947867056, idx_NH4_NO3.gt.0)
    call assert(442660651, idx_H_NO3.gt.0)
    call assert(272503747, idx_H2_SO4.gt.0)

    ! Save the initial concentrations
    true_conc(:,:)            = 0.0
    true_conc(:,idx_H2O)      = 1.0
    true_conc(:,idx_H2O_aq)   = 1.0
    true_conc(:,idx_H_p)      = 1.0
    true_conc(:,idx_NH4_p)    = 1.0
    true_conc(:,idx_SO4_mm)   = 1.0
    true_conc(:,idx_NO3_m)    = 1.0
    true_conc(:,idx_NH42_SO4) = 1.0
    true_conc(:,idx_NH4_NO3)  = 1.0
    true_conc(:,idx_H_NO3)    = 1.0
    true_conc(:,idx_H2_SO4)   = 1.0
    model_conc(:,:) = true_conc(:,:)

    ! Set the polynomial parameters
    
    ! H-NO3 0.1 - 0.4 RH
    HNO3_LRH_B0 = 0.12091
    HNO3_LRH_B1 = 13.497
    HNO3_LRH_B2 = -67.771
    HNO3_LRH_B3 = 144.01
    HNO3_LRH_B4 = -117.97

    ! H-NO3 0.4 - 0.9 RH
    HNO3_HRH_B0 = 1.3424
    HNO3_HRH_B1 = -0.8197
    HNO3_HRH_B2 = -0.52983
    HNO3_HRH_B3 = -0.37335

    ! H2-SO4
    H2SO4_B0 = 9.3948
    H2SO4_B1 = -26.808
    H2SO4_B2 = 35.7654
    H2SO4_B3 = -18.5094

    ! (NH4)2-SO4
    NH42SO4_B0 = -40.4136
    NH42SO4_B1 = 108.798
    NH42SO4_B2 = -170.346
    NH42SO4_B3 = 100.926

    ! NH4-NO3
    NH4NO3_B0 = -17.0372
    NH4NO3_B1 = 59.232
    NH4NO3_B2 = -86.312
    NH4NO3_B3 = 44.04

    ! Set up the ppm->RH (0-1) conversion
    ! (From MOSAIC code, references Seinfeld and Pandis pg. 181)
    ppm_to_RH = 1.0d0 - 373.15d0/temp
    ppm_to_RH = (((-0.1299d0*ppm_to_RH - 0.6445d0)*ppm_to_RH - 1.976d0)*ppm_to_RH &
            + 13.3185d0)*ppm_to_RH
    ppm_to_RH = exp(ppm_to_RH)  ! VP of water (atm)
    ppm_to_RH = (pressure/101325.0d0) / ppm_to_RH * 1.0d-6 ! ppm -> RH (0-1)

    ! Set the initial state in the model
    phlex_state%state_var(:) = model_conc(0,:)

    ! Integrate the mechanism
    do i_RH = 1, NUM_RH_STEP

      ! Set the RH
      a_w = real(1.0, kind=dp)/(NUM_RH_STEP-1)*(i_RH-1)
      true_conc(i_RH, idx_H2O) = a_w / ppm_to_RH
      phlex_state%state_var(idx_H2O) = true_conc(i_RH, idx_H2O)

      ! Get the modeled conc
      call phlex_core%solve(phlex_state, time_step)
      model_conc(i_RH,:) = phlex_state%state_var(:)

      ! Calculate the mean binary activity for H-NO3

      ! Calculate omega
      omega = 6.0d0 * true_conc(i_RH,idx_H_p)   * true_conc(i_RH,idx_SO4_mm) + &
              6.0d0 * true_conc(i_RH,idx_NH4_p) * true_conc(i_RH,idx_SO4_mm) + &
              4.0d0 * true_conc(i_RH,idx_NH4_p) * true_conc(i_RH,idx_NO3_m)

      ! Contribution to ln(gamma_HNO3) from ln(gamma_0_HNO3)
      if (a_w.le.0.1d0) then
        ln_gamma = HNO3_LRH_B0 + HNO3_LRH_B1*0.1d0 + HNO3_LRH_B2*0.1d0**2 + &
                    HNO3_LRH_B3*0.1d0**3 + HNO3_LRH_B4*0.1d0**4
      else if (a_w.le.0.4d0) then
        ln_gamma = HNO3_LRH_B0 + HNO3_LRH_B1*a_w + HNO3_LRH_B2*a_w**2 + &
                    HNO3_LRH_B3*a_w**3 + HNO3_LRH_B4*a_w**4
      else if (a_w.le.0.9d0) then
        ln_gamma = HNO3_HRH_B0 + HNO3_HRH_B1*a_w + HNO3_HRH_B2*a_w**2 + &
                    HNO3_HRH_B3*a_w**3
      else
        ln_gamma = HNO3_HRH_B0 + HNO3_HRH_B1*0.9d0 + HNO3_HRH_B2*0.9d0**2 + &
                    HNO3_HRH_B3*0.9d0**3
      end if

      ! ... from (d(ln(gamma_HNO3))/d(N_Hp N_SO4mm) N_Hp N_SO4mm) / omega
      if (a_w.le.0.1d0) then
        ln_gamma = ln_gamma + true_conc(i_RH,idx_H_p) * true_conc(i_RH,idx_SO4_mm) / omega * &
                   ( H2SO4_B0 + H2SO4_B1*0.1d0 + H2SO4_B2*0.1d0**2 + &
                   H2SO4_B3*0.1d0**3 )
      else if (a_w.le.0.99d0) then
        ln_gamma = ln_gamma + true_conc(i_RH,idx_H_p) * true_conc(i_RH,idx_SO4_mm) / omega * &
                   ( H2SO4_B0 + H2SO4_B1*a_w + H2SO4_B2*a_w**2 + &
                   H2SO4_B3*a_w**3 )
      else
        ln_gamma = ln_gamma + true_conc(i_RH,idx_H_p) * true_conc(i_RH,idx_SO4_mm) / omega * &
                   ( H2SO4_B0 + H2SO4_B1*0.99d0 + H2SO4_B2*0.99d0**2 + &
                   H2SO4_B3*0.99d0**3 )
      end if

      ! ... from (d(ln(gamma_HNO3))/d(N_NH4p N_SO4mm) N_NH4p N_SO4mm) / omega
      if (a_w.le.0.1d0) then
        ln_gamma = ln_gamma + true_conc(i_RH,idx_NH4_p) * true_conc(i_RH,idx_SO4_mm) / omega * &
                   ( NH42SO4_B0 + NH42SO4_B1*0.1d0 + NH42SO4_B2*0.1d0**2 + &
                   NH42SO4_B3*0.1d0**3 )
      else if (a_w.le.0.99d0) then
        ln_gamma = ln_gamma + true_conc(i_RH,idx_NH4_p) * true_conc(i_RH,idx_SO4_mm) / omega * &
                   ( NH42SO4_B0 + NH42SO4_B1*a_w + NH42SO4_B2*a_w**2 + &
                   NH42SO4_B3*a_w**3 )
      else
        ln_gamma = ln_gamma + true_conc(i_RH,idx_NH4_p) * true_conc(i_RH,idx_SO4_mm) / omega * &
                   ( NH42SO4_B0 + NH42SO4_B1*0.99d0 + NH42SO4_B2*0.99d0**2 + &
                   NH42SO4_B3*0.99d0**3 )
      end if

      ! ... from (d(ln(gamma_HNO3))/d(N_NH4p N_NO3m) N_NH4p N_NO3m) / omega
      if (a_w.le.0.1d0) then
        ln_gamma = ln_gamma + true_conc(i_RH,idx_NH4_p) * true_conc(i_RH,idx_NO3_m) / omega * &
                   ( NH4NO3_B0 + NH4NO3_B1*0.1d0 + NH4NO3_B2*0.1d0**2 + &
                   NH4NO3_B3*0.1d0**3 )
      else if (a_w.le.0.99d0) then
        ln_gamma = ln_gamma + true_conc(i_RH,idx_NH4_p) * true_conc(i_RH,idx_NO3_m) / omega * &
                   ( NH4NO3_B0 + NH4NO3_B1*a_w + NH4NO3_B2*a_w**2 + &
                   NH4NO3_B3*a_w**3 )
      else
        ln_gamma = ln_gamma + true_conc(i_RH,idx_NH4_p) * true_conc(i_RH,idx_NO3_m) / omega * &
                   ( NH4NO3_B0 + NH4NO3_B1*0.99d0 + NH4NO3_B2*0.99d0**2 + &
                   NH4NO3_B3*0.99d0**3 )
      end if

      ! Update the H-NO3 mean binary activity
      true_conc(i_RH, idx_H_NO3) = exp(ln_gamma)

    end do

    ! Save the results
    open(unit=7, file="out/PDFiTE_activity_results.txt", status="replace", action="write")
    do i_RH = 0, NUM_RH_STEP
      write(7,*) i_RH*time_step, &
            ' ', true_conc(i_RH, idx_H2O),' ',      model_conc(i_RH, idx_H2O), &
            ' ', true_conc(i_RH, idx_H2O_aq),' ',   model_conc(i_RH, idx_H2O_aq), &
            ' ', true_conc(i_RH, idx_H_p),' ',      model_conc(i_RH, idx_H_p), &
            ' ', true_conc(i_RH, idx_NH4_p),' ',    model_conc(i_RH, idx_NH4_p), &
            ' ', true_conc(i_RH, idx_SO4_mm),' ',   model_conc(i_RH, idx_SO4_mm), &
            ' ', true_conc(i_RH, idx_NO3_m),' ',    model_conc(i_RH, idx_NO3_m), &
            ' ', true_conc(i_RH, idx_NH42_SO4),' ', model_conc(i_RH, idx_NH42_SO4), &
            ' ', true_conc(i_RH, idx_NH4_NO3),' ',  model_conc(i_RH, idx_NH4_NO3), &
            ' ', true_conc(i_RH, idx_H_NO3),' ',    model_conc(i_RH, idx_H_NO3), &
            ' ', true_conc(i_RH, idx_H2_SO4),' ',   model_conc(i_RH, idx_H2_SO4)
    end do
    close(7)

    ! Analyze the results
    do i_RH = 1, NUM_RH_STEP
      do i_spec = 1, size(model_conc,2)  
        call assert_msg(510461551, &
          almost_equal(model_conc(i_RH, i_spec), true_conc(i_RH, i_spec), &
          real(1.0e-2, kind=dp)), "time: "//to_string(i_RH)//"; species: "// &
          to_string(i_spec)//"; mod: "//to_string(model_conc(i_RH, i_spec))// &
          "; true: "//to_string(true_conc(i_RH, i_spec)))
      end do
    end do

    run_PDFiTE_activity_test = .true.

  end function run_PDFiTE_activity_test

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program pmc_test_PDFiTE_activity
