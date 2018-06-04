! Copyright (C) 2017 Matt Dawson
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_test_HL_phase_transfer program

!> Test of HL_phase_transfer reaction module
program pmc_test_HL_phase_transfer

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

  if (run_HL_phase_transfer_tests()) then
    write(*,*) "Henry's Law phase transfer reaction tests - PASS"
  else
    write(*,*) "Henry's Law phase transfer reaction tests - FAIL"
  end if

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Run all pmc_chem_mech_solver tests
  logical function run_HL_phase_transfer_tests() result(passed)

    use pmc_phlex_solver_data

    type(phlex_solver_data_t), pointer :: phlex_solver_data

    phlex_solver_data => phlex_solver_data_t()

    if (phlex_solver_data%is_solver_available()) then
      passed = run_HL_phase_transfer_test()
    else
      call warn_msg(713064651, "No solver available")
      passed = .true.
    end if

    deallocate(phlex_solver_data)

  end function run_HL_phase_transfer_tests

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Solve a mechanism consisting of two phase transfer reactions
  !!
  !! Ozone (O3) and hydrogen peroxide (H2O2) partitioning based on parameters
  !! and equations from CAPRAM 2.4 reduced mechanism.
  !! (Ervens, B., et al., 2003. "CAPRAM 2.4 (MODAC mechanism) : An extended
  !! and condensed tropospheric aqueous phase mechanism and its
  !! application." J. Geophys. Res. 108, 4426. doi:10.1029/2002JD002202
  logical function run_HL_phase_transfer_test()

    use pmc_constants

    type(phlex_core_t), pointer :: phlex_core
    type(phlex_state_t), pointer :: phlex_state
    character(len=:), allocatable :: input_file_path
    type(string_t), allocatable, dimension(:) :: output_file_path

    real(kind=dp), dimension(0:NUM_TIME_STEP, 11) :: model_conc, true_conc
    integer(kind=i_kind) :: idx_phase, idx_aero_rep, idx_O3, idx_O3_aq, &
            idx_H2O2, idx_H2O2_aq, idx_H2O_aq, i_time, i_spec
    character(len=:), allocatable :: key
    real(kind=dp) :: time_step, time, n_star, del_H, del_S, del_G, alpha, &
            crms, M_to_ppm, ugm3_to_ppm, K_eq_O3, K_eq_H2O2, k_O3_forward, &
            k_O3_backward, k_H2O2_forward, k_H2O2_backward, equil_O3, &
            equil_O3_aq, equil_H2O2, equil_H2O2_aq, temp, pressure
    real(kind=dp), target :: radius, number_conc

    run_HL_phase_transfer_test = .true.

    ! Set the environmental and aerosol test conditions
    temp = 272.5d0              ! temperature (K)
    pressure = 101253.3d0       ! pressure (Pa)
    radius = 1.5e-5             ! radius (m)
    number_conc = 1.3e6         ! particle number concentration (#/cc)

    ! Henry's Law equilibrium constants (M/ppm)
    ! O3 HLC Equil Const (M/ppm)
    K_eq_O3 = 1.14d-2 * exp(2300.0d0 * (1.0d0/temp - 1.0d0/298.0d0)) / 1.0d6
    ! H2O2 HLC Equil Const (M/ppm)
    K_eq_H2O2 = 1.025d5 * exp(6340.0d0 * (1.0d0/temp - 1.0d0/298.0d0)) / 1.0d6 
    
    ! Set output time step (s)
    time_step = 1.0d-13

    ! Get the HL_phase_transfer reaction mechanism json file
    input_file_path = 'test_HL_phase_transfer_config.json'

    ! Construct a phlex_core variable
    phlex_core => phlex_core_t(input_file_path)

    ! Initialize the model
    call phlex_core%initialize()

    ! Initialize the solver
    call phlex_core%solver_initialize()

    ! Get a model state variable
    phlex_state => phlex_core%new_state()

    ! Set the environmental conditions
    phlex_state%env_state%temp = temp
    phlex_state%env_state%pressure = pressure
    call phlex_state%update_env_state()

    ! Find the aerosol representation
    call assert(116793129, size(phlex_core%aero_rep).eq.3)
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
    key = "O3"
    idx_O3 = phlex_core%chem_spec_data%gas_state_id(key);
    key = "aqueous aerosol.O3_aq"
    idx_O3_aq = phlex_core%aero_rep(idx_aero_rep)%val%spec_state_id(key);
    key = "H2O2"
    idx_H2O2 = phlex_core%chem_spec_data%gas_state_id(key);
    key = "aqueous aerosol.H2O2_aq"
    idx_H2O2_aq = phlex_core%aero_rep(idx_aero_rep)%val%spec_state_id(key);
    key = "aqueous aerosol.H2O_aq"
    idx_H2O_aq = phlex_core%aero_rep(idx_aero_rep)%val%spec_state_id(key);

    ! Make sure the expected species are in the model
    call assert(202593939, idx_O3.gt.0)
    call assert(262338032, idx_O3_aq.gt.0)
    call assert(374656377, idx_H2O2.gt.0)
    call assert(704441571, idx_H2O2_aq.gt.0)
    call assert(758921357, idx_H2O_aq.gt.0)

    ! Save the initial concentrations
    true_conc(:,:) = 0.0
    true_conc(0,idx_O3) = 0.0
    true_conc(0,idx_O3_aq) = 1.0e-3
    true_conc(0,idx_H2O2) = 1.0
    true_conc(0,idx_H2O2_aq) = 0.0
    true_conc(0,idx_H2O_aq) = 1.4e-2
    model_conc(0,:) = true_conc(0,:)

    ! Determine the M -> ppm conversion using the total aerosol water
    M_to_ppm = number_conc * 1.0d-3 * true_conc(0,idx_H2O_aq) * &
            const%univ_gas_const * temp / pressure

    ! O3 rate constants
    n_star = 1.89d0
    del_H = - 10.0d0 * (n_star-1.0d0) + &
          7.53*(n_star**(2.0d0/3.0d0)-1.0d0) - 1.0d0
    del_S = - 32.0d0 * (n_star-1.0d0) + &
          9.21*(n_star**(2.0d0/3.0d0)-1.0d0) - 1.3d0
    del_G = (del_H - temp * del_S/1000.0d0) * 4184.0d0
    alpha = exp(-del_G/(const%univ_gas_const*temp))
    alpha = alpha / (1.0d0 + alpha)
    crms = sqrt(8.0d0*const%univ_gas_const*temp/(const%pi*48.0d0))
    k_O3_forward = number_conc * ((radius**2 / (3.0d0 * 1.48d-5) + &
          4.0d0 * radius / (3.0d0 * crms * alpha))**(-1))           ! (1/s)
    k_O3_backward = k_O3_forward / (K_eq_O3 * M_to_ppm)             ! (1/s)

    ! H2O2 rate constants
    n_star = 1.74d0
    del_H = - 10.0d0 * (n_star-1.0d0) + &
          7.53*(n_star**(2.0d0/3.0d0)-1.0d0) - 1.0d0
    del_S = - 32.0d0 * (n_star-1.0d0) + &
          9.21*(n_star**(2.0d0/3.0d0)-1.0d0) - 1.3d0
    del_G = (del_H - temp * del_S/1000.0d0) * 4184.0d0
    alpha = exp(-del_G/(const%univ_gas_const*temp))
    alpha = alpha / (1.0d0 + alpha)
    crms = sqrt(8.0d0*const%univ_gas_const*temp/(const%pi*34.0d0))
    k_H2O2_forward = number_conc * ((radius**2 / (3.0d0 * 1.46d-5) + &
          4.0d0 * radius / (3.0d0 * crms * alpha))**(-1))           ! (1/s)
    k_H2O2_backward = k_H2O2_forward / (K_eq_H2O2 * M_to_ppm)       ! (1/s)

    ! Determine the equilibrium concentrations
    ! [A_gas] = [A_total] / (1 + 1/K_HL)
    ! [A_aero] = [A_total] / (K_HL + 1)
    ugm3_to_ppm = const%univ_gas_const * temp / (48.0d0 * pressure)
    equil_O3 = (true_conc(0,idx_O3) + &
            true_conc(0,idx_O3_aq)*number_conc*ugm3_to_ppm) / &
            (K_eq_O3*M_to_ppm + 1.0d0)
    equil_O3_aq = (true_conc(0,idx_O3)/ugm3_to_ppm/number_conc + &
            true_conc(0,idx_O3_aq)) / &
            (1.0d0 + 1.0d0/(K_eq_O3*M_to_ppm))
    
    ugm3_to_ppm = const%univ_gas_const * temp / (34.0d0 * pressure)
    equil_H2O2 = (true_conc(0,idx_H2O2) + &
            true_conc(0,idx_H2O2_aq)*number_conc*ugm3_to_ppm) / &
            (K_eq_H2O2*M_to_ppm + 1.0d0)
    equil_H2O2_aq = (true_conc(0,idx_H2O2)/ugm3_to_ppm/number_conc + &
            true_conc(0,idx_H2O2_aq)) / &
            (1.0d0 + 1.0d0/(K_eq_H2O2*M_to_ppm))

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
      ! [A_gas] = ([A_init_gas] - [A_eq_gas]) 
      !     * exp(-t *(k_f + k_b)) + [A_eq_gas]
      ! [A_aero] = ([A_init_aero] - [A_eq_aero]) 
      !     * exp(-t * (k_f + k_b)) + [A_eq_aero]
      time = i_time * time_step
      true_conc(i_time,idx_O3) = (true_conc(0,idx_O3) - equil_O3) * &
              exp(-time * (k_O3_forward + k_O3_backward)) + equil_O3 
      true_conc(i_time,idx_O3_aq) = (true_conc(0,idx_O3_aq) - equil_O3_aq) * &
              exp(-time * (k_O3_forward + k_O3_backward)) + equil_O3_aq
      true_conc(i_time,idx_H2O2) = (true_conc(0,idx_H2O2) - equil_H2O2) * &
              exp(-time * (k_H2O2_forward + k_H2O2_backward)) + equil_H2O2 
      true_conc(i_time,idx_H2O2_aq) = &
              (true_conc(0,idx_H2O2_aq) - equil_H2O2_aq) * &
              exp(-time * (k_H2O2_forward + k_H2O2_backward)) + equil_H2O2_aq
      true_conc(i_time,idx_H2O_aq) = true_conc(0,idx_H2O_aq)
    end do

    ! Save the results
    open(unit=7, file="out/HL_phase_transfer_results.txt", status="replace", &
            action="write")
    do i_time = 0, NUM_TIME_STEP
      write(7,*) i_time*time_step, &
            ' ', true_conc(i_time, idx_O3), &
            ' ', model_conc(i_time, idx_O3), &
            ' ', true_conc(i_time, idx_O3_aq),&
            ' ', model_conc(i_time, idx_O3_aq), &
            ' ', true_conc(i_time, idx_H2O2), &
            ' ', model_conc(i_time, idx_H2O2), &
            ' ', true_conc(i_time, idx_H2O2_aq), &
            ' ', model_conc(i_time, idx_H2O2_aq), &
            ' ', true_conc(i_time, idx_H2O_aq), &
            ' ', model_conc(i_time, idx_H2O_aq)
    end do
    close(7)

    ! Analyze the results
    do i_time = 1, NUM_TIME_STEP
      do i_spec = 1, 11
        ! Only check the second phase
        if (i_spec.ge.2.and.i_spec.le.8) cycle
        call assert_msg(114526423, &
          almost_equal(model_conc(i_time, i_spec), &
          true_conc(i_time, i_spec), real(1.0e-2, kind=dp)), "time: "// &
          to_string(i_time)//"; species: "//to_string(i_spec)//"; mod: "// &
          to_string(model_conc(i_time, i_spec))//"; true: "// &
          to_string(true_conc(i_time, i_spec)))
      end do
    end do

    deallocate(phlex_state)
    deallocate(phlex_core)

    run_HL_phase_transfer_test = .true.

  end function run_HL_phase_transfer_test

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program pmc_test_HL_phase_transfer
