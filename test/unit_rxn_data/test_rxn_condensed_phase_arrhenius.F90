! Copyright (C) 2017 Matt Dawson
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_test_condensed_phase_arrhenius program

!> Test of condensed_phase_arrhenius reaction module
program pmc_test_condensed_phase_arrhenius

  use iso_c_binding

  use pmc_util,                         only: i_kind, dp, assert, &
                                              almost_equal, string_t, &
                                              warn_msg
  use pmc_phlex_core
  use pmc_phlex_state
  use pmc_aero_rep_data
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

  if (run_condensed_phase_arrhenius_tests()) then
    write(*,*) "Condensed-phase Arrhenius reaction tests - PASS"
  else
    write(*,*) "Condensed-phase Arrhenius reaction tests - FAIL"
  end if

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Run all pmc_chem_mech_solver tests
  logical function run_condensed_phase_arrhenius_tests() result(passed)

    use pmc_phlex_solver_data

    type(phlex_solver_data_t), pointer :: phlex_solver_data

    phlex_solver_data => phlex_solver_data_t()

    if (phlex_solver_data%is_solver_available()) then
      passed = run_condensed_phase_arrhenius_test()
    else
      call warn_msg(187891224, "No solver available")
      passed = .true.
    end if

    deallocate(phlex_solver_data)

  end function run_condensed_phase_arrhenius_tests

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Solve a mechanism consisting of two sets of condensed-phase reactions
  !!
  !! The mechanism is of the form:
  !!
  !!   A + D -k1-> B -k2-> C + D
  !!
  !! where k1 and k2 are condensed-phase Arrhenius reaction rate constants,
  !! and D is a constant species. One set of reactions is done with units
  !! of 'M' and one is with 'ug/m3'.
  logical function run_condensed_phase_arrhenius_test()

    use pmc_constants

    type(phlex_core_t), pointer :: phlex_core
    type(phlex_state_t), pointer :: phlex_state
    character(len=:), allocatable :: input_file_path, key
    type(string_t), allocatable, dimension(:) :: output_file_path

    class(aero_rep_data_t), pointer :: aero_rep_ptr
    integer(kind=i_kind), parameter :: NUM_STATE_VAR = 27
    real(kind=dp), dimension(0:NUM_TIME_STEP, NUM_STATE_VAR) :: model_conc, &
            true_conc
    integer(kind=i_kind) :: idx_A_aq, idx_B_aq, idx_C_aq, idx_D_aq, idx_H2O, &
            idx_A_org, idx_B_org, idx_C_org, idx_D_org, idx_aq_phase, &
            idx_org_phase, i_time, i_spec
    real(kind=dp) :: time_step, time, conc_D, conc_water, MW_A, MW_B, MW_C, &
            MW_D, k1_aq, k2_aq, k1_org, k2_org, temp, pressure

    run_condensed_phase_arrhenius_test = .true.

    ! Set the environmental and aerosol test conditions
    temp = 272.5d0              ! temperature (K)
    pressure = 101253.3d0       ! pressure (Pa)

    ! Set the rate constants (for calculating the true values)
    conc_D = 1.2d-2
    conc_water = 2.3d0
    MW_A = 0.1572
    MW_B = 0.0219
    MW_C = 0.2049
    MW_D = 0.0345
    k1_aq = conc_D / MW_D / conc_water + 1476.0d0 * &
            exp( -5.5d-21 / (const%boltzmann * temp) ) * &
            (temp/300.0d0)**(150.0d0) * (1.0d0 + 0.15d0 * pressure) / 60.0d0
    k2_aq = 21.0d0 * exp( -4000.0d0/temp ) * (temp/315d0)**(11.0d0) * &
            (1.0d0 + 0.05d0 * pressure)
    k1_org = conc_D / (MW_D*1.0e9) + 1476.0d0 * &
             exp( -5.5d-21 / (const%boltzmann * temp) ) * &
             (temp/300.0d0)**(150.0d0) * (1.0d0 + 0.15d0 * pressure) / 60.0d0
    k2_org = 21.0d0 * exp( -4000.0d0/temp ) * (temp/315d0)**(11.0d0) * &
             (1.0d0 + 0.05d0 * pressure)

    ! Set output time step (s)
    time_step = 1.0d0

    ! Get the condensed_phase_arrhenius reaction mechanism json file
    input_file_path = 'test_condensed_phase_arrhenius_config.json'

    ! Construct a phlex_core variable
    phlex_core => phlex_core_t(input_file_path)

    ! Initialize the model
    call phlex_core%initialize()

    ! Initialize the solver
    call phlex_core%solver_initialize()

    ! Get a model state variable
    phlex_state => phlex_core%new_state()

    ! Check the size of the state array
    call assert(235226766, size(phlex_state%state_var).eq.NUM_STATE_VAR)

    ! Set the environmental conditions
    phlex_state%env_state%temp = temp
    phlex_state%env_state%pressure = pressure
    call phlex_state%update_env_state()

    ! Find the aerosol representation
    key = "my aero rep 2"
    call assert(421062613, phlex_core%get_aero_rep(key, aero_rep_ptr))

    ! Get species indices
    key = "aqueous aerosol.A"
    idx_A_aq = aero_rep_ptr%spec_state_id(key);
    key = "aqueous aerosol.B"
    idx_B_aq = aero_rep_ptr%spec_state_id(key);
    key = "aqueous aerosol.C"
    idx_C_aq = aero_rep_ptr%spec_state_id(key);
    key = "aqueous aerosol.D"
    idx_D_aq = aero_rep_ptr%spec_state_id(key);
    key = "aqueous aerosol.H2O_aq"
    idx_H2O = aero_rep_ptr%spec_state_id(key);
    key = "organic aerosol.A"
    idx_A_org = aero_rep_ptr%spec_state_id(key);
    key = "organic aerosol.B"
    idx_B_org = aero_rep_ptr%spec_state_id(key);
    key = "organic aerosol.C"
    idx_C_org = aero_rep_ptr%spec_state_id(key);
    key = "organic aerosol.D"
    idx_D_org = aero_rep_ptr%spec_state_id(key);

    ! Make sure the expected species are in the model
    call assert(643455452, idx_A_aq.gt.0)
    call assert(917307621, idx_B_aq.gt.0)
    call assert(747150717, idx_C_aq.gt.0)
    call assert(241944312, idx_D_aq.gt.0)
    call assert(971787407, idx_H2O.gt.0)
    call assert(801630503, idx_A_org.gt.0)
    call assert(631473599, idx_B_org.gt.0)
    call assert(743791944, idx_C_org.gt.0)
    call assert(573635040, idx_D_org.gt.0)

    ! Save the initial concentrations
    true_conc(:,:) = 0.0
    true_conc(0,idx_A_aq) = 1.0
    true_conc(0,idx_B_aq) = 0.0
    true_conc(0,idx_C_aq) = 0.0
    true_conc(0,idx_D_aq) = conc_D
    true_conc(:,idx_H2O) = conc_water
    true_conc(0,idx_A_org) = 1.0
    true_conc(0,idx_B_org) = 0.0
    true_conc(0,idx_C_org) = 0.0
    true_conc(0,idx_D_org) = conc_D
    model_conc(0,:) = true_conc(0,:)

    ! Set the initial state in the model
    phlex_state%state_var(:) = model_conc(0,:)

    ! Integrate the mechanism
    do i_time = 1, NUM_TIME_STEP

      ! Get the modeled conc
      call phlex_core%solve(phlex_state, time_step)
      model_conc(i_time,:) = phlex_state%state_var(:)

      ! Get the analytic concentrations
      time = i_time * time_step
      true_conc(i_time,idx_A_aq) = true_conc(0,idx_A_aq) * exp(-k1_aq*time)
      true_conc(i_time,idx_B_aq) = true_conc(0,idx_A_aq) * &
              (k1_aq/(k2_aq-k1_aq)) * &
              (exp(-k1_aq*time) - exp(-k2_aq*time)) * MW_B / MW_A
      true_conc(i_time,idx_C_aq) = true_conc(0,idx_A_aq) * MW_C / MW_A * &
              (1.0d0 + (k1_aq*exp(-k2_aq*time) - &
              k2_aq*exp(-k1_aq*time))/(k2_aq-k1_aq))
      true_conc(i_time,idx_D_aq) = true_conc(0,idx_D_aq)
      true_conc(i_time,idx_H2O) = true_conc(0,idx_H2O)
      true_conc(i_time,idx_A_org) = true_conc(0,idx_A_org) * exp(-k1_org*time)
      true_conc(i_time,idx_B_org) = true_conc(0,idx_A_org) * &
              (k1_org/(k2_org-k1_org)) * &
              (exp(-k1_org*time) - exp(-k2_org*time)) * MW_B / MW_A
      true_conc(i_time,idx_C_org) = true_conc(0,idx_A_org) * MW_C / MW_A * &
              (1.0d0 + (k1_org*exp(-k2_org*time) - &
              k2_org*exp(-k1_org*time))/(k2_org-k1_org))
      true_conc(i_time,idx_D_org) = true_conc(0,idx_D_org)

    end do

    ! Save the results
    open(unit=7, file="out/condensed_phase_arrhenius_results.txt", &
            status="replace", action="write")
    do i_time = 0, NUM_TIME_STEP
      write(7,*) i_time*time_step, &
            ' ', true_conc(i_time, idx_A_aq), &
            ' ', model_conc(i_time, idx_A_aq), &
            ' ', true_conc(i_time, idx_B_aq), &
            ' ', model_conc(i_time, idx_B_aq), &
            ' ', true_conc(i_time, idx_C_aq), &
            ' ', model_conc(i_time, idx_C_aq), &
            ' ', true_conc(i_time, idx_D_aq), &
            ' ', model_conc(i_time, idx_D_aq), &
            ' ', true_conc(i_time, idx_H2O), &
            ' ', model_conc(i_time, idx_H2O), &
            ' ', true_conc(i_time, idx_A_org), &
            ' ', model_conc(i_time, idx_A_org), &
            ' ', true_conc(i_time, idx_B_org), &
            ' ', model_conc(i_time, idx_B_org), &
            ' ', true_conc(i_time, idx_C_org), &
            ' ', model_conc(i_time, idx_C_org), &
            ' ', true_conc(i_time, idx_D_org), &
            ' ', model_conc(i_time, idx_D_org)
    end do
    close(7)

    ! Analyze the results
    do i_time = 1, NUM_TIME_STEP
      do i_spec = 1, size(model_conc, 2)
        call assert_msg(923311346, &
          almost_equal(model_conc(i_time, i_spec),  &
          true_conc(i_time, i_spec), real(1.0e-2, kind=dp)), "time: "// &
          to_string(i_time)//"; species: "//to_string(i_spec)//"; mod: "// &
          to_string(model_conc(i_time, i_spec))//"; true: "// &
          to_string(true_conc(i_time, i_spec)))
      end do
    end do

    deallocate(phlex_state)
    deallocate(phlex_core)

    run_condensed_phase_arrhenius_test = .true.

  end function run_condensed_phase_arrhenius_test

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program pmc_test_condensed_phase_arrhenius
