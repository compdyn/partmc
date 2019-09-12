! Copyright (C) 2017 Matt Dawson
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_test_aqueous_equilibrium program

!> Test of aqueous_equilibrium reaction module
program pmc_test_aqueous_equilibrium

  use iso_c_binding

  use pmc_util,                         only: i_kind, dp, assert, &
                                              almost_equal, string_t, &
                                              warn_msg
  use pmc_camp_core
  use pmc_camp_state
  use pmc_aero_rep_data
  use pmc_aero_rep_factory
  use pmc_aero_rep_single_particle
  use pmc_aero_rep_modal_binned_mass
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

  if (run_aqueous_equilibrium_tests()) then
    if (pmc_mpi_rank().eq.0) write(*,*) &
            "Aqueous equilibrium reaction tests - PASS"
  else
    if (pmc_mpi_rank().eq.0) write(*,*) &
            "Aqueous equilibrium reaction tests - FAIL"
    stop 3
  end if

  ! finalize mpi
  call pmc_mpi_finalize()

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Run all pmc_chem_mech_solver tests
  logical function run_aqueous_equilibrium_tests() result(passed)

    use pmc_camp_solver_data

    type(camp_solver_data_t), pointer :: camp_solver_data

    camp_solver_data => camp_solver_data_t()

    if (camp_solver_data%is_solver_available()) then
      passed = run_aqueous_equilibrium_test(1)
      passed = passed .and. run_aqueous_equilibrium_test(2)
    else
      call warn_msg(713064651, "No solver available")
      passed = .true.
    end if

    deallocate(camp_solver_data)

  end function run_aqueous_equilibrium_tests

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Solve a mechanism consisting of two aqueous equilibrium reactions
  !!
  !! One of two scenarios is tested, depending on the passed integer:
  !! (1) single-particle aerosol representation and fixed water concentration
  !! (2) modal aerosol representation and ZSR-calculated water concentration
  logical function run_aqueous_equilibrium_test(scenario)

    use pmc_constants

    !> Scenario flag
    integer, intent(in) :: scenario

    type(camp_core_t), pointer :: camp_core
    type(camp_state_t), pointer :: camp_state
    character(len=:), allocatable :: input_file_path, key, idx_prefix
    type(string_t), allocatable, dimension(:) :: output_file_path

    class(aero_rep_data_t), pointer :: aero_rep_ptr
    real(kind=dp), allocatable, dimension(:,:) :: model_conc, true_conc
    integer(kind=i_kind) :: idx_A, idx_B, idx_C, idx_D, idx_E, idx_F, idx_G, &
            idx_H, idx_BC_act, idx_H2O, idx_phase, i_time, i_spec
    real(kind=dp) :: time_step, time, Keq_1, Keq_2, Keq_3, k1_forward, &
            k2_forward, k3_forward, k1_reverse, k2_reverse, k3_reverse, &
            total_init, equil_A, equil_B, equil_C, equil_D, equil_E, &
            equil_F, equil_G, equil_H, x, x0, temp, pressure
    real(kind=dp), target :: radius, number_conc
#ifdef PMC_USE_MPI
    character, allocatable :: buffer(:), buffer_copy(:)
    integer(kind=i_kind) :: pack_size, pos, i_elem, results
#endif

    type(solver_stats_t), target :: solver_stats

    integer(kind=i_kind) :: i_sect_unused, i_sect_the_mode
    type(aero_rep_factory_t) :: aero_rep_factory
    type(aero_rep_update_data_modal_binned_mass_GMD_t) :: update_data_GMD
    type(aero_rep_update_data_modal_binned_mass_GSD_t) :: update_data_GSD

    call assert_msg(533630504, scenario.ge.1 .and. scenario.le.2, &
                    "Invalid scenario specified: "//to_string( scenario ))

    run_aqueous_equilibrium_test = .true.


    ! Allocate space for the results
    if (scenario.eq.1) then
      allocate(model_conc(0:NUM_TIME_STEP, 30))
      allocate(true_conc(0:NUM_TIME_STEP, 30))
    else if (scenario.eq.2) then
      allocate(model_conc(0:NUM_TIME_STEP, 15))
      allocate(true_conc(0:NUM_TIME_STEP, 15))
    end if

    ! Set the environmental and aerosol test conditions
    temp = 272.5d0              ! temperature (K)
    pressure = 101253.3d0       ! pressure (Pa)

    ! Set output time step (s)
    time_step = 1.0d0

#ifdef PMC_USE_MPI
    ! Load the model data on the root process and pass it to process 1 for solving
    if (pmc_mpi_rank().eq.0) then
#endif

      ! Get the aqueous_equilibrium reaction mechanism json file
      if (scenario.eq.1) then
        input_file_path = 'test_aqueous_equilibrium_config.json'
      else if (scenario.eq.2) then
        input_file_path = 'test_aqueous_equilibrium_config_2.json'
      end if

      ! Construct a camp_core variable
      camp_core => camp_core_t(input_file_path)

      deallocate(input_file_path)

      ! Initialize the model
      call camp_core%initialize()

      ! Find the aerosol representation
      key = "my aero rep 2"
      call assert(750324390, camp_core%get_aero_rep(key, aero_rep_ptr))

      if (scenario.eq.2) then

        ! Set the aerosol representation id
        select type (aero_rep_ptr)
          type is (aero_rep_modal_binned_mass_t)
            call aero_rep_factory%initialize_update_data( aero_rep_ptr, &
                                                          update_data_GMD)
            call aero_rep_factory%initialize_update_data( aero_rep_ptr, &
                                                          update_data_GSD)
            call assert_msg(225977671, &
                  aero_rep_ptr%get_section_id("unused mode", i_sect_unused), &
                  "Could not get section id for the unused mode")
            call assert_msg(738353917, &
                  aero_rep_ptr%get_section_id("the mode", i_sect_the_mode), &
                  "Could not get section id for the unused mode")

          class default
            call die_msg(804363261, "Wrong aero rep type")
        end select

      end if

      ! Get species indices
      if (scenario.eq.1) then
        idx_prefix = ""
      else if (scenario.eq.2) then
        idx_prefix = "the mode."
      end if
      key = idx_prefix//"aqueous aerosol.A"
      idx_A = aero_rep_ptr%spec_state_id(key);
      key = idx_prefix//"aqueous aerosol.B"
      idx_B = aero_rep_ptr%spec_state_id(key);
      key = idx_prefix//"aqueous aerosol.C"
      idx_C = aero_rep_ptr%spec_state_id(key);
      key = idx_prefix//"aqueous aerosol.B-C"
      idx_BC_act = aero_rep_ptr%spec_state_id(key);
      key = idx_prefix//"aqueous aerosol.D"
      idx_D = aero_rep_ptr%spec_state_id(key);
      key = idx_prefix//"aqueous aerosol.E"
      idx_E = aero_rep_ptr%spec_state_id(key);
      key = idx_prefix//"aqueous aerosol.F"
      idx_F = aero_rep_ptr%spec_state_id(key);
      key = idx_prefix//"aqueous aerosol.G"
      idx_G = aero_rep_ptr%spec_state_id(key);
      key = idx_prefix//"aqueous aerosol.H"
      idx_H = aero_rep_ptr%spec_state_id(key);
      key = idx_prefix//"aqueous aerosol.H2O_aq"
      idx_H2O = aero_rep_ptr%spec_state_id(key);

      ! Make sure the expected species are in the model
      call assert(503629528, idx_A.gt.0)
      call assert(445790969, idx_B.gt.0)
      call assert(217795506, idx_C.gt.0)
      call assert(612589100, idx_BC_act.gt.0)
      call assert(159956947, idx_D.gt.0)
      call assert(784651538, idx_E.gt.0)
      call assert(444337730, idx_F.gt.0)
      call assert(779060521, idx_G.gt.0)
      call assert(438746713, idx_H.gt.0)
      call assert(446243264, idx_H2O.gt.0)

#ifdef PMC_USE_MPI
      ! pack the camp core
      pack_size = camp_core%pack_size() &
                  + update_data_GMD%pack_size() &
                  + update_data_GSD%pack_size()
      allocate(buffer(pack_size))
      pos = 0
      call camp_core%bin_pack(buffer, pos)
      call update_data_GMD%bin_pack(buffer, pos)
      call update_data_GSD%bin_pack(buffer, pos)
      call assert(672194140, pos.eq.pack_size)
    end if

    ! broadcast the species ids
    call pmc_mpi_bcast_integer(idx_A)
    call pmc_mpi_bcast_integer(idx_B)
    call pmc_mpi_bcast_integer(idx_BC_act)
    call pmc_mpi_bcast_integer(idx_C)
    call pmc_mpi_bcast_integer(idx_D)
    call pmc_mpi_bcast_integer(idx_E)
    call pmc_mpi_bcast_integer(idx_F)
    call pmc_mpi_bcast_integer(idx_G)
    call pmc_mpi_bcast_integer(idx_H)
    call pmc_mpi_bcast_integer(idx_H2O)
    call pmc_mpi_bcast_integer(i_sect_unused)
    call pmc_mpi_bcast_integer(i_sect_the_mode)

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
      camp_core => camp_core_t()
      pos = 0
      call camp_core%bin_unpack(buffer, pos)
      call update_data_GMD%bin_unpack(buffer, pos)
      call update_data_GSD%bin_unpack(buffer, pos)
      call assert(101979335, pos.eq.pack_size)
      allocate(buffer_copy(pack_size))
      pos = 0
      call camp_core%bin_pack(buffer_copy, pos)
      call update_data_GMD%bin_pack(buffer_copy, pos)
      call update_data_GSD%bin_pack(buffer_copy, pos)
      call assert(161723428, pos.eq.pack_size)
      do i_elem = 1, pack_size
        call assert_msg(556517022, buffer(i_elem).eq.buffer_copy(i_elem), &
                "Mismatch in element: "//trim(to_string(i_elem)))
      end do

      ! solve and evaluate results on process 1
#endif

      ! Initialize the solver
      call camp_core%solver_initialize()

      ! Get a model state variable
      camp_state => camp_core%new_state()

      ! Update the GMD and GSD for the aerosol modes
      if (scenario.eq.2) then
        ! unused mode
        call update_data_GMD%set_GMD(i_sect_unused, 1.2d-6)
        call update_data_GSD%set_GSD(i_sect_unused, 1.2d0)
        call camp_core%update_aero_rep_data(update_data_GMD)
        call camp_core%update_aero_rep_data(update_data_GSD)
        ! the mode
        call update_data_GMD%set_GMD(i_sect_the_mode, 9.3d-7)
        call update_data_GSD%set_GSD(i_sect_the_mode, 0.9d0)
        call camp_core%update_aero_rep_data(update_data_GMD)
        call camp_core%update_aero_rep_data(update_data_GSD)
      end if

      ! Set the environmental conditions
      camp_state%env_state%temp = temp
      camp_state%env_state%pressure = pressure
      call camp_state%update_env_state()

      ! Save the initial concentrations
      true_conc(:,:) = 0.0
      true_conc(0,idx_A) = 13.5
      true_conc(0,idx_B) = 0.0
      true_conc(0,idx_C) = 0.0
      true_conc(:,idx_BC_act) = 0.5
      true_conc(0,idx_D) = 0.0
      true_conc(0,idx_E) = 0.0
      true_conc(0,idx_F) = 8.0
      true_conc(0,idx_G) = 12.0
      true_conc(0,idx_H) = 0.0
      true_conc(:,idx_H2O) = 1.4e-2
      model_conc(0,:) = true_conc(0,:)

      ! Henry's Law equilibrium constants (M/ppm)
      Keq_1 = 1.14d-2 * exp(2300.0d0 * (1.0d0/temp - 1.0d0/298.0d0)) ! (M^2/M^2)
      Keq_2 = 12.3                                                   ! (M^3/M^2)
      Keq_3 = 2.35 * exp(1245.7d0 * (1.0d0/temp - 1.0d0/298.0d0))    ! (M/M)

      ! Calculate the forward and reverse rate constants for reaction 1
      k1_reverse = 0.32d0 * true_conc(0,idx_BC_act)                  ! (1/M/s)
      k1_forward = Keq_1 * k1_reverse                                ! (1/M/s)

      ! Calculate the forward and reverse rate constants for reaction 2
      k2_reverse = 3.25e-3                                           ! (1/M/M/s)
      k2_forward = Keq_2 * k2_reverse                                ! (1/M/s)

      ! Calculate the forward and reverse rate constants for reaction 3
      k3_reverse = 1.56e-4                                           ! (1/s)
      k3_forward = Keq_3 * k3_reverse                                ! (1/s)

      ! Determine the equilibrium concentrations (ug/m3)
      !
      ! Reaction 1 (equil values in M)
      !
      ! K_eq = ([B][C])/([A]^2)
      ! [Total] = [A]i
      ! [B] = [C] = [Total] * (sqrt(1/Keq)-2)/(1/Keq-4)
      ! [A] = [Total] - [B] - [C]
      total_init = true_conc(0,idx_A)/true_conc(0,idx_H2O) * 1000.0d0/48.0d0
      equil_B = (total_init * (sqrt(1.0d0/Keq_1)-2.0d0) / (1.0d0/Keq_1-4.0d0))
      equil_C = (total_init * (sqrt(1.0d0/Keq_1)-2.0d0) / (1.0d0/Keq_1-4.0d0))
      equil_A = (total_init * (1.0d0 - 2.0d0*(sqrt(1.0d0/Keq_1)-2.0d0) / &
              (1.0d0/Keq_1-4.0d0)))

      ! Reaction 2
      !
      ! K_eq = [F]^2/([D][E])
      ! [Total] = [F]i
      ! [D] = [E] = [Total] * (sqrt(Keq)-2)/(Keq-4)
      ! [F] = [Total] - [D] - [E]
      total_init = true_conc(0,idx_F)/true_conc(0,idx_H2O) * 1000.0d0/28.0d0
      equil_D = (total_init * (sqrt(Keq_2)-2.0d0) / (Keq_2-4.0d0)) * &
              true_conc(0,idx_H2O) * 27.6d0 / 1000.0d0
      equil_E = (total_init * (sqrt(Keq_2)-2.0d0) / (Keq_2-4.0d0)) * &
            true_conc(0,idx_H2O) * 202.4d0 / 1000.0d0
      equil_F = (total_init * (1.0d0 - 2.0d0*(sqrt(Keq_2)-2.0d0) / &
              (Keq_2-4.0d0))) * true_conc(0,idx_H2O) * 28.0d0 / 1000.0d0

      ! Reaction 3
      !
      ! K_eq = [G]/[H]
      ! [Total] = [G]i
      ! [G] = [Total] / (1 + 1/Keq)
      ! [G] = [Total] / (Keq + 1)
      total_init = true_conc(0,idx_G)/true_conc(0,idx_H2O) * 1000.0d0/35.67d0
      equil_H = (total_init / (1.0d0 + 1.0d0/Keq_3)) * &
              true_conc(0,idx_H2O) * 284.2 / 1000.0d0
      equil_G = (total_init / (Keq_3 + 1.0d0)) * &
              true_conc(0,idx_H2O) * 35.67d0 / 1000.0d0

      ! Set the initial state in the model
      camp_state%state_var(:) = model_conc(0,:)

#ifdef PMC_DEBUG
      ! Evaluate the Jacobian during solving
      solver_stats%eval_Jac = .true.
#endif

      ! Integrate the mechanism
      do i_time = 1, NUM_TIME_STEP

        ! Get the modeled conc
        call camp_core%solve(camp_state, time_step, &
                              solver_stats = solver_stats)
        model_conc(i_time,:) = camp_state%state_var(:)

#ifdef PMC_DEBUG
        ! Check the Jacobian evaluations
        call assert_msg(164937823, solver_stats%Jac_eval_fails.eq.0, &
                        trim( to_string( solver_stats%Jac_eval_fails ) )// &
                        " Jacobian evaluation failures at time step "// &
                        trim( to_string( i_time ) ) )
#endif

        ! Get the analytic concentrations
        !
        ! Two-reactant, two-product reactions
        ! FIXME This is not the analytic answer (waiting to use MatLab to solve)
        ! x = [A] - [A_eq]
        ! x0 = [A_init] - [A_eq]
        ! [A] = x + [A_eq] = x0exp(-t/tau) + [A_eq]
        ! 1/tau = k_f + k_b
        ! [A] = ([A_init] - [A_eq]) * exp(-t *(k_f + k_b)) + [A_eq]
        time = i_time * time_step
        x0 = true_conc(0,idx_A) * 1000.0d0 / 48.0d0 / true_conc(0,idx_H2O)
        x = 2.0 / ((1.0d0+2.0d0/x0)* &
                exp(-2.0d0*(k1_forward*equil_A+k1_reverse*equil_B)*time)-1.0d0)
        true_conc(i_time,idx_A) = &
                (equil_A + x) * true_conc(0,idx_H2O) * 48.0d0 / 1000.0d0
        true_conc(i_time,idx_B) = &
                (equil_B - x) * true_conc(0,idx_H2O) * 32.67d0 / 1000.0d0
        true_conc(i_time,idx_C) = &
                (equil_C - x) * true_conc(0,idx_H2O) * 114.3d0 / 1000.0d0
        true_conc(i_time,idx_D) = (true_conc(0,idx_D) - equil_D) * &
                exp(-time * (k2_forward + k2_reverse)) + equil_D
        true_conc(i_time,idx_E) = (true_conc(0,idx_E) - equil_E) * &
                exp(-time * (k2_forward + k2_reverse)) + equil_E
        true_conc(i_time,idx_F) = (true_conc(0,idx_F) - equil_F) * &
                exp(-time * (k2_forward + k2_reverse)) + equil_F
        !
        ! One-reactant, one-product reaction
        ! x = [G] - [G_eq]
        ! x0 = [G_init] - [G_eq]
        ! [G] = x + [G_eq] = x0 * exp(-t/tau) + [G_eq]
        ! 1/tau = k_f + f_r
        ! [G] = ([G_init] - [G_eq]) * exp(-t * (k_f + k_r)) + [G_eq]
        ! [H] = ([H_init] - [H_eq]) * exp(-t * (k_f + k_r)) + [H_eq]
        true_conc(i_time,idx_G) = (true_conc(0,idx_G) - equil_G) * &
                exp(-time * (k3_forward + k3_reverse)) + equil_G
        true_conc(i_time,idx_H) = (true_conc(0,idx_H) - equil_H) * &
                exp(-time * (k3_forward + k3_reverse)) + equil_H

      end do

      ! Save the results
      if (scenario.eq.1) then
        open(unit=7, file="out/aqueous_equilibrium_results.txt", &
                status="replace", action="write")
      else if (scenario.eq.2) then
        open(unit=7, file="out/aqueous_equilibrium_results_2.txt", &
                status="replace", action="write")
      end if
      do i_time = 0, NUM_TIME_STEP
        write(7,*) i_time*time_step, &
              ' ', true_conc(i_time, idx_A),' ', model_conc(i_time, idx_A), &
              ' ', true_conc(i_time, idx_B),' ', model_conc(i_time, idx_B), &
              ' ', true_conc(i_time, idx_C),' ', model_conc(i_time, idx_C), &
              ' ', true_conc(i_time, idx_D),' ', model_conc(i_time, idx_D), &
              ' ', true_conc(i_time, idx_E),' ', model_conc(i_time, idx_E), &
              ' ', true_conc(i_time, idx_F),' ', model_conc(i_time, idx_F), &
              ' ', true_conc(i_time, idx_G),' ', model_conc(i_time, idx_G), &
              ' ', true_conc(i_time, idx_H),' ', model_conc(i_time, idx_H), &
              ' ', true_conc(i_time, idx_H2O),' ', model_conc(i_time, idx_H2O)
      end do
      close(7)

      ! Analyze the results (single-particle only)
      if (scenario.eq.1) then
        do i_time = 1, NUM_TIME_STEP
          do i_spec = 1, size(model_conc, 2)
            ! FIXME Check all species once a true value is found for the 4
            ! component reactions
            if (i_spec.ne.idx_G.and.i_spec.ne.idx_H) cycle
            call assert_msg(703162515, &
              almost_equal(model_conc(i_time, i_spec), &
              true_conc(i_time, i_spec), real(1.0e-2, kind=dp)).or. &
              (model_conc(i_time, i_spec).lt.1e-5*model_conc(1, i_spec).and. &
              true_conc(i_time, i_spec).lt.1e-5*true_conc(1, i_spec)), &
              "time: "//trim(to_string(i_time))//"; species: "// &
              trim(to_string(i_spec))//"; mod: "// &
              trim(to_string(model_conc(i_time, i_spec)))//"; true: "// &
              trim(to_string(true_conc(i_time, i_spec))))
          end do
        end do
      end if

      deallocate(camp_state)

#ifdef PMC_USE_MPI
      ! convert the results to an integer
      if (run_aqueous_equilibrium_test) then
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
        run_aqueous_equilibrium_test = .true.
      else
        run_aqueous_equilibrium_test = .false.
      end if
    end if

    deallocate(buffer)
#endif

    deallocate(camp_core)
  end function run_aqueous_equilibrium_test

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program pmc_test_aqueous_equilibrium
