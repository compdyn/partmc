! Copyright (C) 2017 Matt Dawson
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_test_aqueous_equilibrium program

!> Test of aqueous_equilibrium reaction module
program pmc_test_aqueous_equilibrium

  use iso_c_binding

  use pmc_util,                         only: phlex_real, phlex_int, &
                                              assert, almost_equal, &
                                              string_t, warn_msg
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
  integer(kind=phlex_int) :: NUM_TIME_STEP = 100

  ! initialize mpi
  call pmc_mpi_init()

  if (run_aqueous_equilibrium_tests()) then
    if (pmc_mpi_rank().eq.0) write(*,*) &
            "Aqueous equilibrium reaction tests - PASS"
  else
    if (pmc_mpi_rank().eq.0) write(*,*) &
            "Aqueous equilibrium reaction tests - FAIL"
  end if

  ! finalize mpi
  call pmc_mpi_finalize()

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Run all pmc_chem_mech_solver tests
  logical function run_aqueous_equilibrium_tests() result(passed)

    use pmc_phlex_solver_data

    type(phlex_solver_data_t), pointer :: phlex_solver_data

    phlex_solver_data => phlex_solver_data_t()

    if (phlex_solver_data%is_solver_available()) then
      passed = run_aqueous_equilibrium_test()
    else
      call warn_msg(713064651, "No solver available")
      passed = .true.
    end if

    deallocate(phlex_solver_data)

  end function run_aqueous_equilibrium_tests

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Solve a mechanism consisting of two aqueous equilibrium reactions
  logical function run_aqueous_equilibrium_test()

    use pmc_constants

    type(phlex_core_t), pointer :: phlex_core
    type(phlex_state_t), pointer :: phlex_state
    character(len=:), allocatable :: input_file_path, key
    type(string_t), allocatable, dimension(:) :: output_file_path

    class(aero_rep_data_t), pointer :: aero_rep_ptr
    real(kind=phlex_real), dimension(0:NUM_TIME_STEP, 30) :: model_conc, true_conc
    integer(kind=phlex_int) :: idx_A, idx_B, idx_C, idx_D, idx_E, idx_F, idx_G, &
            idx_H, idx_BC_act, idx_H2O, idx_phase, i_time, i_spec
    real(kind=phlex_real) :: time_step, time, Keq_1, Keq_2, Keq_3, k1_forward, &
            k2_forward, k3_forward, k1_reverse, k2_reverse, k3_reverse, &
            total_init, equil_A, equil_B, equil_C, equil_D, equil_E, &
            equil_F, equil_G, equil_H, x, x0, temp, pressure
    real(kind=phlex_real), target :: radius, number_conc
#ifdef PMC_USE_MPI
    character, allocatable :: buffer(:), buffer_copy(:)
    integer(kind=phlex_int) :: pack_size, pos, i_elem, results
#endif

    run_aqueous_equilibrium_test = .true.

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
      input_file_path = 'test_aqueous_equilibrium_config.json'

      ! Construct a phlex_core variable
      phlex_core => phlex_core_t(input_file_path)

      deallocate(input_file_path)

      ! Initialize the model
      call phlex_core%initialize()

      ! Find the aerosol representation
      key = "my aero rep 2"
      call assert(750324390, phlex_core%get_aero_rep(key, aero_rep_ptr))

      ! Get species indices
      key = "aqueous aerosol.A"
      idx_A = aero_rep_ptr%spec_state_id(key);
      key = "aqueous aerosol.B"
      idx_B = aero_rep_ptr%spec_state_id(key);
      key = "aqueous aerosol.C"
      idx_C = aero_rep_ptr%spec_state_id(key);
      key = "aqueous aerosol.B-C"
      idx_BC_act = aero_rep_ptr%spec_state_id(key);
      key = "aqueous aerosol.D"
      idx_D = aero_rep_ptr%spec_state_id(key);
      key = "aqueous aerosol.E"
      idx_E = aero_rep_ptr%spec_state_id(key);
      key = "aqueous aerosol.F"
      idx_F = aero_rep_ptr%spec_state_id(key);
      key = "aqueous aerosol.G"
      idx_G = aero_rep_ptr%spec_state_id(key);
      key = "aqueous aerosol.H"
      idx_H = aero_rep_ptr%spec_state_id(key);
      key = "aqueous aerosol.H2O_aq"
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
      ! pack the phlex core
      pack_size = phlex_core%pack_size()
      allocate(buffer(pack_size))
      pos = 0
      call phlex_core%bin_pack(buffer, pos)
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
      call assert(101979335, pos.eq.pack_size)
      allocate(buffer_copy(pack_size))
      pos = 0
      call phlex_core%bin_pack(buffer_copy, pos)
      call assert(161723428, pos.eq.pack_size)
      do i_elem = 1, pack_size
        call assert_msg(556517022, buffer(i_elem).eq.buffer_copy(i_elem), &
                "Mismatch in element: "//trim(to_string(i_elem)))
      end do

      ! solve and evaluate results on process 1
#endif

      ! Initialize the solver
      call phlex_core%solver_initialize()

      ! Get a model state variable
      phlex_state => phlex_core%new_state()

      ! Set the environmental conditions
      phlex_state%env_state%temp = temp
      phlex_state%env_state%pressure = pressure
      call phlex_state%update_env_state()

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
      phlex_state%state_var(:) = model_conc(0,:)

      ! Integrate the mechanism
      do i_time = 1, NUM_TIME_STEP

        ! Get the modeled conc
        call phlex_core%solve(phlex_state, time_step)
        model_conc(i_time,:) = phlex_state%state_var(:)

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
      open(unit=7, file="out/aqueous_equilibrium_results.txt", &
              status="replace", action="write")
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

      ! Analyze the results
      do i_time = 1, NUM_TIME_STEP
        do i_spec = 1, size(model_conc, 2)
          ! FIXME Check all species once a true value is found for the 4 
          ! component reactions
          if (i_spec.ne.idx_G.and.i_spec.ne.idx_H) cycle
          call assert_msg(703162515, &
            almost_equal(model_conc(i_time, i_spec), &
            true_conc(i_time, i_spec), real(1.0e-2, kind=phlex_real)).or. &
            (model_conc(i_time, i_spec).lt.1e-5*model_conc(1, i_spec).and. &
            true_conc(i_time, i_spec).lt.1e-5*true_conc(1, i_spec)), &
            "time: "//trim(to_string(i_time))//"; species: "// &
            trim(to_string(i_spec))//"; mod: "// &
            trim(to_string(model_conc(i_time, i_spec)))//"; true: "// &
            trim(to_string(true_conc(i_time, i_spec))))
        end do
      end do

      deallocate(phlex_state)

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

    deallocate(phlex_core)
  end function run_aqueous_equilibrium_test

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program pmc_test_aqueous_equilibrium
