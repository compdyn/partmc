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
  use pmc_solver_stats
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
    if (pmc_mpi_rank().eq.0) write(*,*) "Condensed-phase Arrhenius reaction tests - PASS"
  else
    if (pmc_mpi_rank().eq.0) write(*,*) "Condensed-phase Arrhenius reaction tests - FAIL"
    stop 3
  end if

  ! finalize mpi
  call pmc_mpi_finalize()

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Run all pmc_chem_mech_solver tests
  logical function run_condensed_phase_arrhenius_tests() result(passed)

    use pmc_phlex_solver_data

    type(phlex_solver_data_t), pointer :: phlex_solver_data

    phlex_solver_data => phlex_solver_data_t()

    if (phlex_solver_data%is_solver_available()) then
      passed = run_condensed_phase_arrhenius_test(1)
      passed = passed .and. run_condensed_phase_arrhenius_test(2)
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
  !!
  !! One of two scenarios is tested, depending on the passed integer:
  !! (1) single-particle aerosol representation and fixed water concentration
  !! (2) modal aerosol representation and ZSR-calculated water concentration
  logical function run_condensed_phase_arrhenius_test(scenario)

    use pmc_constants

    !> Scenario flag
    integer, intent(in) :: scenario

    type(phlex_core_t), pointer :: phlex_core
    type(phlex_state_t), pointer :: phlex_state
    character(len=:), allocatable :: input_file_path, key, idx_prefix
    type(string_t), allocatable, dimension(:) :: output_file_path

    class(aero_rep_data_t), pointer :: aero_rep_ptr
    integer(kind=i_kind) :: num_state_var, state_size
    real(kind=dp), allocatable, dimension(:,:) :: model_conc, true_conc
    integer(kind=i_kind) :: idx_A_aq, idx_B_aq, idx_C_aq, idx_D_aq, idx_H2O, &
            idx_A_org, idx_B_org, idx_C_org, idx_D_org, idx_aq_phase, &
            idx_org_phase, i_time, i_spec
    real(kind=dp) :: time_step, time, conc_D, conc_water, MW_A, MW_B, MW_C, &
            MW_D, k1_aq, k2_aq, k1_org, k2_org, temp, pressure
#ifdef PMC_USE_MPI
    character, allocatable :: buffer(:), buffer_copy(:)
    integer(kind=i_kind) :: pack_size, pos, i_elem, results
#endif

    type(solver_stats_t), target :: solver_stats

    call assert_msg(619416982, scenario.ge.1 .and. scenario.le.2, &
                    "Invalid scenario specified: "//to_string( scenario ))

    run_condensed_phase_arrhenius_test = .true.

    ! Allocate space for the results
    if (scenario.eq.1) then
      num_state_var = 27
    else if (scenario.eq.2) then
      num_state_var = 14
    end if
    allocate(model_conc(0:NUM_TIME_STEP, num_state_var))
    allocate(true_conc(0:NUM_TIME_STEP, num_state_var))

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

#ifdef PMC_USE_MPI
    ! Load the model data on root process and pass it to process 1 for solving
    if (pmc_mpi_rank().eq.0) then
#endif

      ! Get the condensed_phase_arrhenius reaction mechanism json file
      if (scenario.eq.1) then
        input_file_path = 'test_condensed_phase_arrhenius_config.json'
      else if (scenario.eq.2) then
        input_file_path = 'test_condensed_phase_arrhenius_config_2.json'
      end if

      ! Construct a phlex_core variable
      phlex_core => phlex_core_t(input_file_path)

      deallocate(input_file_path)

      ! Initialize the model
      call phlex_core%initialize()

      ! Find the aerosol representation
      key = "my aero rep 2"
      call assert(421062613, phlex_core%get_aero_rep(key, aero_rep_ptr))

      ! Get species indices
      if (scenario.eq.1) then
        idx_prefix = ""
      else if (scenario.eq.2) then
        idx_prefix = "the mode."
      end if
      key = idx_prefix//"aqueous aerosol.A"
      idx_A_aq = aero_rep_ptr%spec_state_id(key);
      key = idx_prefix//"aqueous aerosol.B"
      idx_B_aq = aero_rep_ptr%spec_state_id(key);
      key = idx_prefix//"aqueous aerosol.C"
      idx_C_aq = aero_rep_ptr%spec_state_id(key);
      key = idx_prefix//"aqueous aerosol.D"
      idx_D_aq = aero_rep_ptr%spec_state_id(key);
      key = idx_prefix//"aqueous aerosol.H2O_aq"
      idx_H2O = aero_rep_ptr%spec_state_id(key);
      key = idx_prefix//"organic aerosol.A"
      idx_A_org = aero_rep_ptr%spec_state_id(key);
      key = idx_prefix//"organic aerosol.B"
      idx_B_org = aero_rep_ptr%spec_state_id(key);
      key = idx_prefix//"organic aerosol.C"
      idx_C_org = aero_rep_ptr%spec_state_id(key);
      key = idx_prefix//"organic aerosol.D"
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

#ifdef PMC_USE_MPI
      ! pack the phlex core
      pack_size = phlex_core%pack_size()
      allocate(buffer(pack_size))
      pos = 0
      call phlex_core%bin_pack(buffer, pos)
      call assert(636035849, pos.eq.pack_size)
    end if

    ! broadcast the species ids
    call pmc_mpi_bcast_integer(idx_A_aq)
    call pmc_mpi_bcast_integer(idx_B_aq)
    call pmc_mpi_bcast_integer(idx_C_aq)
    call pmc_mpi_bcast_integer(idx_D_aq)
    call pmc_mpi_bcast_integer(idx_H2O)
    call pmc_mpi_bcast_integer(idx_A_org)
    call pmc_mpi_bcast_integer(idx_B_org)
    call pmc_mpi_bcast_integer(idx_C_org)
    call pmc_mpi_bcast_integer(idx_D_org)

    ! broadcast the buffer size
    call pmc_mpi_bcast_integer(pack_size)

    if (pmc_mpi_rank().eq.1) then
      ! allocate the buffer to receive data
      allocate(buffer(pack_size))
    end if

    ! broadcast the data
    call pmc_mpi_bcast_packed(buffer)

    if (pmc_mpi_rank().eq.1) then
      ! unpack the data
      phlex_core => phlex_core_t()
      pos = 0
      call phlex_core%bin_unpack(buffer, pos)
      call assert(913246791, pos.eq.pack_size)
      allocate(buffer_copy(pack_size))
      pos = 0
      call phlex_core%bin_pack(buffer_copy, pos)
      call assert(408040386, pos.eq.pack_size)
      do i_elem = 1, pack_size
        call assert_msg(185309230, buffer(i_elem).eq.buffer_copy(i_elem), &
                "Mismatch in element :"//trim(to_string(i_elem)))
      end do

      ! solve and evaluate results on process 1
#endif

      ! Initialize the solver
      call phlex_core%solver_initialize()

      ! Get a model state variable
      phlex_state => phlex_core%new_state()

      ! Check the size of the state array
      state_size = size(phlex_state%state_var)
      call assert_msg(235226766, state_size.eq.NUM_STATE_VAR, &
                      "Wrong state size: "//to_string( state_size ))

      ! Set the environmental conditions
      phlex_state%env_state%temp = temp
      phlex_state%env_state%pressure = pressure
      call phlex_state%update_env_state()

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

#ifdef PMC_DEBUG
      ! Evaluate the Jacobian during solving
      solver_stats%eval_Jac = .true.
#endif

      ! Integrate the mechanism
      do i_time = 1, NUM_TIME_STEP

        ! Get the modeled conc
        call phlex_core%solve(phlex_state, time_step, &
                              solver_stats = solver_stats)
        model_conc(i_time,:) = phlex_state%state_var(:)

#ifdef PMC_DEBUG
        ! Check the Jacobian evaluations
        call assert_msg(772386254, solver_stats%Jac_eval_fails.eq.0, &
                        trim( to_string( solver_stats%Jac_eval_fails ) )// &
                        " Jacobian evaluation failures at time step "// &
                        trim( to_string( i_time ) ) )
#endif

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
      if (scenario.eq.1) then
        open(unit=7, file="out/condensed_phase_arrhenius_results.txt", &
              status="replace", action="write")
      else if (scenario.eq.2) then
        open(unit=7, file="out/condensed_phase_arrhenius_results_2.txt", &
              status="replace", action="write")
      end if
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

      ! Analyze the results (single-particle only)
      if (scenario.eq.1) then
        do i_time = 1, NUM_TIME_STEP
          do i_spec = 1, size(model_conc, 2)
            call assert_msg(248109045, &
              almost_equal(model_conc(i_time, i_spec), &
              true_conc(i_time, i_spec), real(1.0e-2, kind=dp)).or. &
              (model_conc(i_time, i_spec).lt.1e-5*model_conc(1, i_spec).and. &
              true_conc(i_time, i_spec).lt.1e-5*true_conc(1, i_spec)).or. &
              (model_conc(i_time, i_spec).lt.1.0e-30.and. &
              true_conc(i_time, i_spec).lt.1.0e-30), &
              "time: "//trim(to_string(i_time))//"; species: "// &
              trim(to_string(i_spec))//"; mod: "// &
              trim(to_string(model_conc(i_time, i_spec)))//"; true: "// &
              trim(to_string(true_conc(i_time, i_spec))))
          end do
        end do
      end if

      deallocate(phlex_state)

#ifdef PMC_USE_MPI
      ! convert the results to an integer
      if (run_condensed_phase_arrhenius_test) then
        results = 0
      else
        results = 1
      end if
    end if

    ! Send the results back to the primary process
    call pmc_mpi_transfer_integer(results, results, 1, 0)

    ! convert the results back to a logical value
    if (pmc_mpi_rank().eq.0) then
      if (results.eq.0) then
        run_condensed_phase_arrhenius_test = .true.
      else
        run_condensed_phase_arrhenius_test = .false.
      end if
    end if

    deallocate(buffer)
#endif

    deallocate(phlex_core)

  end function run_condensed_phase_arrhenius_test

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program pmc_test_condensed_phase_arrhenius
