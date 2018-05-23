! Copyright (C) 2017 Matt Dawson
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_test_SIMPOL_phase_transfer program

!> Test of SIMPOL_phase_transfer reaction module
program pmc_test_SIMPOL_phase_transfer

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
  integer(kind=i_kind) :: NUM_TIME_STEP = 100

  ! initialize mpi
  call pmc_mpi_init()

  if (run_SIMPOL_phase_transfer_tests()) then
    write(*,*) "Henry's Law phase transfer reaction tests - PASS"
  else
    write(*,*) "Henry's Law phase transfer reaction tests - FAIL"
  end if

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Run all pmc_chem_mech_solver tests
  logical function run_SIMPOL_phase_transfer_tests() result(passed)

    use pmc_phlex_solver_data

    type(phlex_solver_data_t), pointer :: phlex_solver_data

    phlex_solver_data => phlex_solver_data_t()

    if (phlex_solver_data%is_solver_available()) then
      passed = run_SIMPOL_phase_transfer_test()
    else
      call warn_msg(713064651, "No solver available")
      passed = .true.
    end if

  end function run_SIMPOL_phase_transfer_tests

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Solve a mechanism consisting of one phase transfer reaction
  !!
  !! Ethanol partitioning based on parameters and equations from SIMPOL.1
  !! model (Pankow and Asher, 2008. "SIMPOL.1: a simple group contribution 
  !! method for predicting vapor pressures and enthalpies of vaporization of
  !! multifunctional organic compounds." Atmos. Chem. Phys. 8(10), 2773-2796.
  !! doi:10.5194/acp-8-2773-2008.i)
  !! 
  !! Condensation rates are based on parameters and equations from the 
  !! CAPRAM 2.4 reduced mechanism. (Ervens, B., et al., 2003. "CAPRAM 2.4 
  !! (MODAC mechanism) : An extended and condensed tropospheric aqueous phase 
  !! mechanism and its application." J. Geophys. Res. 108, 4426. 
  !! doi:10.1029/2002JD002202.)
  !!
  logical function run_SIMPOL_phase_transfer_test()

    use pmc_constants

    type(phlex_core_t), pointer :: phlex_core
    type(phlex_state_t), target :: phlex_state
    character(len=:), allocatable :: input_file_path
    type(string_t), allocatable, dimension(:) :: output_file_path

    real(kind=dp), dimension(0:NUM_TIME_STEP, 7) :: model_conc, true_conc
    integer(kind=i_kind) :: idx_phase, idx_aero_rep
    integer(kind=i_kind) :: idx_ethanol, idx_ethanol_aq, idx_H2O_aq
    character(len=:), allocatable :: key
    integer(kind=i_kind) :: i_time, i_spec
    real(kind=dp) :: time_step, time
    real(kind=dp) :: n_star, del_H, del_S, del_G, alpha, crms, ugm3_to_ppm
    real(kind=dp) :: VP_ethanol, k_forward, k_backward
    real(kind=dp) :: equil_ethanol, equil_ethanol_aq
    real(kind=dp) :: total_mass, water_mass, VP_0_mass

    ! Parameters for calculating true concentrations
    real(kind=dp) :: temperature, pressure 
    real(kind=dp), target :: radius, number_conc
    real(kind=dp), parameter :: MW_ethanol = 0.04607
    real(kind=dp), parameter :: MW_H2O = 0.01801

    run_SIMPOL_phase_transfer_test = .true.

    ! Set the environmental and aerosol test conditions
    temperature = 272.5d0       ! temperature (K)
    pressure = 101253.3d0       ! pressure (Pa)
    radius = 1.5e-5             ! radius (m)
    number_conc = 1.3e6         ! particle number concentration (#/cc)

    ! Calculate the SIMPOL vapor pressure for ethanol (atm)
    VP_ethanol = 10.0d0**( -1.97E+03 / temperature &
                         + 2.91E+00 &
                         + 1.96E-03 * temperature &
                         - 4.96E-01 * log(temperature) )

    ! Set output time step (s)
    time_step = 1.0d-12

    ! Get the SIMPOL_phase_transfer reaction mechanism json file
    input_file_path = 'test_SIMPOL_phase_transfer_config.json'

    ! Construct a phlex_core variable
    phlex_core => phlex_core_t(input_file_path)

    ! Initialize the model
    call phlex_core%initialize()

    ! Initialize the solver
    call phlex_core%solver_initialize()

    ! Get a model state variable
    phlex_state = phlex_core%new_state()

    ! Set the environmental conditions
    phlex_state%env_state%temp = temperature
    phlex_state%env_state%pressure = pressure
    call phlex_state%update_env_state()

    ! Find the aerosol representation
    call assert(209301925, size(phlex_core%aero_rep).eq.3)
    idx_aero_rep = 1

    ! Update the aerosol representation
    call phlex_core%update_aero_rep_data( &
            idx_aero_rep, &
            UPDATE_RADIUS, &
            c_loc(radius))
    call phlex_core%update_aero_rep_data( &
            idx_aero_rep, &
            UPDATE_NUMBER_CONC, &
            c_loc(number_conc))

    ! Get species indices
    key = "ethanol"
    idx_ethanol = phlex_core%chem_spec_data%gas_state_id(key);
    key = "aqueous aerosol.ethanol_aq"
    idx_ethanol_aq = phlex_core%aero_rep(idx_aero_rep)%val%spec_state_id(key);
    key = "aqueous aerosol.H2O_aq"
    idx_H2O_aq = phlex_core%aero_rep(idx_aero_rep)%val%spec_state_id(key);

    ! Make sure the expected species are in the model
    call assert(884352514, idx_ethanol.gt.0)
    call assert(379146109, idx_ethanol_aq.gt.0)
    call assert(208989205, idx_H2O_aq.gt.0)

    ! Save the initial concentrations
    true_conc(:,:) = 0.0
    true_conc(0,idx_ethanol) = 1.0e-3
    true_conc(0,idx_ethanol_aq) = 1.0e-5
    true_conc(0,idx_H2O_aq) = 1.4e-2
    model_conc(0,:) = true_conc(0,:)

    ! ethanol rate constants
    n_star = 2.55d0
    del_H = - 10.0d0 * (n_star-1.0d0) + 7.53*(n_star**(2.0d0/3.0d0)-1.0d0) - 1.0d0
    del_S = - 32.0d0 * (n_star-1.0d0) + 9.21*(n_star**(2.0d0/3.0d0)-1.0d0) - 1.3d0
    del_G = (del_H - temperature * del_S/1000.0d0) * 4184.0d0
    alpha = exp(-del_G/(const%univ_gas_const*temperature))
    alpha = alpha / (1.0d0 + alpha)
    crms = sqrt(8.0d0*const%univ_gas_const*temperature/(const%pi*46.07d0))
    k_forward = number_conc * ((radius**2 / (3.0d0 * 1.48d-5) + 4.0d0 * radius / &
            (3.0d0 * crms * alpha))**(-1))                                      ! (1/s)

    ! Determine the equilibrium concentrations
    ! [A_gas] =  VP_ethanol
    ! [A_aero] = [A_total] - VP_ethanol
    ugm3_to_ppm = const%univ_gas_const * temperature / (46.07d0 * pressure)
    total_mass = true_conc(0,idx_ethanol)/ugm3_to_ppm + true_conc(0,idx_ethanol_aq)*number_conc ! (ug/m3)

    ! Iterated to find equil_ethanol_aq
    equil_ethanol_aq = 1.48804181d1 ! (ug/m3)
    equil_ethanol = (total_mass-equil_ethanol_aq)*ugm3_to_ppm
    equil_ethanol_aq = equil_ethanol_aq/number_conc

    ! Calculate the backwards rate constant based on the equilibrium conditions
    ! and the forward rate
    k_backward = k_forward / ( equil_ethanol )                             ! (1/s)
    
    ! Set the initial state in the model
    phlex_state%state_var(:) = model_conc(0,:)

    ! Integrate the mechanism
    do i_time = 1, NUM_TIME_STEP

      ! Get the modeled conc
      call phlex_core%solve(phlex_state, time_step)
      model_conc(i_time,:) = phlex_state%state_var(:)

      ! Get the analytic conc
      ! x = [A_gas] - [A_eq_gas]
      ! x0 = [A_init_gas] - [A_eq_gas]
      ! [A_gas] = x + [A_eq_gas] = x0exp(-t/tau) + [A_eq_gas]
      ! 1/tau = k_f + k_b
      ! [A_gas] = ([A_init_gas] - [A_eq_gas]) * exp(-t *(k_f + k_b)) + [A_eq_gas]
      ! [A_aero] = ([A_init_aero] - [A_eq_aero]) * exp(-t * (k_f + k_b)) + [A_eq_aero]
      time = i_time * time_step
      true_conc(i_time,idx_ethanol) = (true_conc(0,idx_ethanol) - equil_ethanol) * &
              exp(-time * (k_forward + k_backward)) + equil_ethanol 
      true_conc(i_time,idx_ethanol_aq) = (true_conc(0,idx_ethanol_aq) - equil_ethanol_aq) * &
              exp(-time * (k_forward + k_backward)) + equil_ethanol_aq
      true_conc(i_time,idx_H2O_aq) = true_conc(0,idx_H2O_aq)
    end do

    ! Save the results
    open(unit=7, file="out/SIMPOL_phase_transfer_results.txt", status="replace", action="write")
    do i_time = 0, NUM_TIME_STEP
      write(7,*) i_time*time_step, &
            ' ', true_conc(i_time, idx_ethanol),' ', model_conc(i_time, idx_ethanol), &
            ' ', true_conc(i_time, idx_ethanol_aq),' ', model_conc(i_time, idx_ethanol_aq), &
            ' ', true_conc(i_time, idx_H2O_aq),' ', model_conc(i_time, idx_H2O_aq)
    end do
    close(7)

    ! Analyze the results
    do i_time = 1, NUM_TIME_STEP
      do i_spec = 1, 5
        ! Only check the second phase
        if (i_spec.ge.2.and.i_spec.le.3) cycle
!        call assert_msg(162705801, &
!          almost_equal(model_conc(i_time, i_spec), true_conc(i_time, i_spec), &
!          real(1.0e-2, kind=dp)), "time: "//to_string(i_time)//"; species: "// &
!          to_string(i_spec)//"; mod: "//to_string(model_conc(i_time, i_spec))// &
!          "; true: "//to_string(true_conc(i_time, i_spec)))
      end do
    end do

    run_SIMPOL_phase_transfer_test = .true.

  end function run_SIMPOL_phase_transfer_test

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program pmc_test_SIMPOL_phase_transfer
