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
  use pmc_chem_spec_data
  use pmc_aero_rep_data
  use pmc_aero_rep_factory
  use pmc_aero_rep_modal_binned_mass
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

  if (run_HL_phase_transfer_tests()) then
    if (pmc_mpi_rank().eq.0) write(*,*) "Henry's Law phase transfer reaction tests - PASS"
  else
    if (pmc_mpi_rank().eq.0) write(*,*) "Henry's Law phase transfer reaction tests - FAIL"
    stop 3
  end if

  ! finalize mpi
  call pmc_mpi_finalize()

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Run all pmc_chem_mech_solver tests
  logical function run_HL_phase_transfer_tests() result(passed)

    use pmc_phlex_solver_data

    type(phlex_solver_data_t), pointer :: phlex_solver_data

    phlex_solver_data => phlex_solver_data_t()

    if (phlex_solver_data%is_solver_available()) then
      passed = run_HL_phase_transfer_test(1)
      passed = passed .and. run_HL_phase_transfer_test(2)
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
  !!
  !! One of two scenarios is tested, depending on the passed integer:
  !! one with a single-particle aerosol representation (1)
  !! and one with a modal aerosol representation (2)
  logical function run_HL_phase_transfer_test(scenario)

    use pmc_constants

    !> Scenario flag
    integer, intent(in) :: scenario

    type(phlex_core_t), pointer :: phlex_core
    type(phlex_state_t), pointer :: phlex_state
    character(len=:), allocatable :: input_file_path
    type(string_t), allocatable, dimension(:) :: output_file_path

    type(chem_spec_data_t), pointer :: chem_spec_data
    class(aero_rep_data_t), pointer :: aero_rep_ptr
    real(kind=dp), allocatable, dimension(:,:) :: model_conc, true_conc
    integer(kind=i_kind) :: idx_phase, idx_aero_rep, idx_O3, idx_O3_aq, &
            idx_H2O2, idx_H2O2_aq, idx_H2O_aq, i_time, i_spec
    character(len=:), allocatable :: key, idx_prefix
    real(kind=dp) :: time_step, time, n_star, del_H, del_S, del_G, alpha, &
            crms, M_to_ppm, ugm3_to_ppm, K_eq_O3, K_eq_H2O2, k_O3_forward, &
            k_O3_backward, k_H2O2_forward, k_H2O2_backward, equil_O3, &
            equil_O3_aq, equil_H2O2, equil_H2O2_aq, temp, pressure
    real(kind=dp), target :: radius, number_conc
#ifdef PMC_USE_MPI
    character, allocatable :: buffer(:), buffer_copy(:)
    integer(kind=i_kind) :: pack_size, pos, i_elem, results
#endif

    type(solver_stats_t), target :: solver_stats

    ! For setting particle radius and number concentration
    type(aero_rep_factory_t) :: aero_rep_factory
    type(aero_rep_update_data_single_particle_radius_t) :: radius_update
    type(aero_rep_update_data_single_particle_number_t) :: number_update
    integer(kind=i_kind), parameter :: aero_rep_external_id = 12

    call assert_msg(144071521, scenario.ge.1 .and. scenario.le.2, &
                    "Invalid scenario specified: "//to_string( scenario ) )

    run_HL_phase_transfer_test = .true.

    ! Allocate space for the results
    if (scenario.eq.1) then
      allocate(model_conc(0:NUM_TIME_STEP, 11))
      allocate(true_conc(0:NUM_TIME_STEP, 11))
    else if (scenario.eq.2) then
      allocate(model_conc(0:NUM_TIME_STEP, 6))
      allocate(true_conc(0:NUM_TIME_STEP, 6))
    endif

    ! Set the environmental conditions
    temp = 272.5d0              ! temperature (K)
    pressure = 101253.3d0       ! pressure (Pa)

    ! Henry's Law equilibrium constants (M/ppm)
    ! O3 HLC Equil Const (M/ppm)
    K_eq_O3 = 1.14d-2 * exp(2300.0d0 * (1.0d0/temp - 1.0d0/298.0d0)) / 1.0d6
    ! H2O2 HLC Equil Const (M/ppm)
    K_eq_H2O2 = 1.025d5 * exp(6340.0d0 * (1.0d0/temp - 1.0d0/298.0d0)) / 1.0d6

    ! Set output time step (s)
    time_step = 1.0d-13

#ifdef PMC_USE_MPI
    ! Load the model data on root process and pass it to process 1 for solving
    if (pmc_mpi_rank().eq.0) then
#endif

      ! Get the HL_phase_transfer reaction mechanism json file
      if (scenario.eq.1) then
        input_file_path = 'test_HL_phase_transfer_config.json'
      else if (scenario.eq.2) then
        input_file_path = 'test_HL_phase_transfer_config_2.json'
      end if

      ! Construct a phlex_core variable
      phlex_core => phlex_core_t(input_file_path)

      deallocate(input_file_path)

      ! Initialize the model
      call phlex_core%initialize()

      ! Find the aerosol representation
      key = "my aero rep 2"
      call assert(116793129, phlex_core%get_aero_rep(key, aero_rep_ptr))
      if (scenario.eq.1) then
        select type (aero_rep_ptr)
          type is (aero_rep_single_particle_t)
            call aero_rep_ptr%set_id(aero_rep_external_id)
          class default
            call die_msg(866102326, "Incorrect aerosol representation type")
        end select
      else if (scenario.eq.2) then
        select type (aero_rep_ptr)
          type is (aero_rep_modal_binned_mass_t)
            call aero_rep_ptr%set_id(aero_rep_external_id)
          class default
            call die_msg(290304323, "Incorrect aerosol representation type")
        end select
      end if

      ! Get the chemical species data
      call assert(191714381, phlex_core%get_chem_spec_data(chem_spec_data))

      ! Get species indices
      if (scenario.eq.1) then
        idx_prefix = ""
      else if (scenario.eq.2) then
        idx_prefix = "the mode."
      end if
      key = "O3"
      idx_O3 = chem_spec_data%gas_state_id(key);
      key = idx_prefix//"aqueous aerosol.O3_aq"
      idx_O3_aq = aero_rep_ptr%spec_state_id(key);
      key = "H2O2"
      idx_H2O2 = chem_spec_data%gas_state_id(key);
      key = idx_prefix//"aqueous aerosol.H2O2_aq"
      idx_H2O2_aq = aero_rep_ptr%spec_state_id(key);
      key = idx_prefix//"aqueous aerosol.H2O_aq"
      idx_H2O_aq = aero_rep_ptr%spec_state_id(key);

      ! Make sure the expected species are in the model
      call assert(202593939, idx_O3.gt.0)
      call assert(262338032, idx_O3_aq.gt.0)
      call assert(374656377, idx_H2O2.gt.0)
      call assert(704441571, idx_H2O2_aq.gt.0)
      call assert(758921357, idx_H2O_aq.gt.0)

#ifdef PMC_USE_MPI
      ! pack the phlex core
      pack_size = phlex_core%pack_size()
      allocate(buffer(pack_size))
      pos = 0
      call phlex_core%bin_pack(buffer, pos)
      call assert(947937853, pos.eq.pack_size)
    end if

    ! broadcast the species ids
    call pmc_mpi_bcast_integer(idx_O3)
    call pmc_mpi_bcast_integer(idx_O3_aq)
    call pmc_mpi_bcast_integer(idx_H2O2)
    call pmc_mpi_bcast_integer(idx_H2O2_aq)
    call pmc_mpi_bcast_integer(idx_H2O_aq)

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
      call assert(711319310, pos.eq.pack_size)
      allocate(buffer_copy(pack_size))
      pos = 0
      call phlex_core%bin_pack(buffer_copy, pos)
      call assert(206112905, pos.eq.pack_size)
      do i_elem = 1, pack_size
        call assert_msg(265856998, buffer(i_elem).eq.buffer_copy(i_elem), &
                "Mismatch in element :"//trim(to_string(i_elem)))
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
      true_conc(0,idx_O3) = 0.0
      true_conc(0,idx_O3_aq) = 1.0e-3
      true_conc(0,idx_H2O2) = 1.0
      true_conc(0,idx_H2O2_aq) = 0.0
      true_conc(0,idx_H2O_aq) = 1.4e-2
      model_conc(0,:) = true_conc(0,:)

      ! Calculate the radius and number concentration to use
      ! ( the real values for the modal representation cannot be calculated
      !   because the number concentrations change sligthly during the run
      !   but the Jacobian checker can be run as a check. )
      if (scenario.eq.1) then
        radius = 1.5e-5             ! radius (m)
        number_conc = 1.3e6         ! particle number concentration (#/cc)
      else if (scenario.eq.2) then
        ! radius (m)
        radius = 9.37e-7 / 2.0 * exp(9.0/2.0 * 0.9 * 0.9)
        ! number conc
        number_conc = 1.0 / (const%pi/6.0 * (9.37e-7)**3.0 * &
                             exp(9.0/2.0 * 0.9 * 0.9))
        number_conc = number_conc * 1.0e-9 * (1.0e-3 + 1.4e-2)
      end if

      ! Update the aerosol representation (single-particle only)
      if (scenario.eq.1) then
        call aero_rep_factory%initialize_update_data(radius_update)
        call aero_rep_factory%initialize_update_data(number_update)
        call radius_update%set_radius(aero_rep_external_id, radius)
        call number_update%set_number(aero_rep_external_id, number_conc)
        call phlex_core%update_aero_rep_data(radius_update)
        call phlex_core%update_aero_rep_data(number_update)
      end if

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
        call assert_msg(470924483, solver_stats%Jac_eval_fails.eq.0, &
                        trim( to_string( solver_stats%Jac_eval_fails ) )// &
                        " Jacobian evaluation failures at time step "// &
                        trim( to_string( i_time ) ) )
#endif

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
      if (scenario.eq.1) then
        open(unit=7, file="out/HL_phase_transfer_results.txt", status="replace", &
              action="write")
      else if (scenario.eq.2) then
        open(unit=7, file="out/HL_phase_transfer_results_2.txt", status="replace", &
              action="write")
      end if
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

      ! Analyze the results (single-particle only)
      if (scenario.eq.1) then
        do i_time = 1, NUM_TIME_STEP
          do i_spec = 1, size(model_conc, 2)
            if (i_spec.ge.2.and.i_spec.le.8) cycle
            call assert_msg(411096108, &
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
      endif

      deallocate(phlex_state)

#ifdef PMC_USE_MPI
      ! convert the results to an integer
      if (run_HL_phase_transfer_test) then
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
        run_HL_phase_transfer_test = .true.
      else
        run_HL_phase_transfer_test = .false.
      end if
    end if

    deallocate(buffer)
#endif

    deallocate(phlex_core)
    deallocate(model_conc)
    deallocate(true_conc)

  end function run_HL_phase_transfer_test

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program pmc_test_HL_phase_transfer
