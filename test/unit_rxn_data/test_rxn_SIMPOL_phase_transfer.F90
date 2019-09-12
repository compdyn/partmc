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
  use pmc_camp_core
  use pmc_camp_state
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

  if (run_SIMPOL_phase_transfer_tests()) then
    if (pmc_mpi_rank().eq.0) write(*,*) &
            "SIMPOL phase transfer reaction tests - PASS"
  else
    if (pmc_mpi_rank().eq.0) write(*,*) &
            "SIMPOL phase transfer reaction tests - FAIL"
    stop 3
  end if

  ! finalize mpi
  call pmc_mpi_finalize()

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Run all pmc_chem_mech_solver tests
  logical function run_SIMPOL_phase_transfer_tests() result(passed)

    use pmc_camp_solver_data

    type(camp_solver_data_t), pointer :: camp_solver_data

    camp_solver_data => camp_solver_data_t()

    if (camp_solver_data%is_solver_available()) then
      passed = run_SIMPOL_phase_transfer_test(1)
      passed = passed .and. run_SIMPOL_phase_transfer_test(2)
    else
      call warn_msg(713064651, "No solver available")
      passed = .true.
    end if

    deallocate(camp_solver_data)

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
  !! One of two scenarios is tested, depending on the passed integer:
  !! one with a single-particle aerosol representation (1)
  !! and one with a modal aerosol representation (2).
  !! Scenario (2) includes UNIFAC activity calculations to test the
  !! Jacobian calculations of the UNIFAC sub model with the Jacobian
  !! checker.
  logical function run_SIMPOL_phase_transfer_test(scenario)

    use pmc_constants

    !> Scenario flag
    integer, intent(in) :: scenario

    type(camp_core_t), pointer :: camp_core
    type(camp_state_t), pointer :: camp_state
    character(len=:), allocatable :: input_file_path, key, idx_prefix
    type(string_t), allocatable, dimension(:) :: output_file_path

    type(chem_spec_data_t), pointer :: chem_spec_data
    class(aero_rep_data_t), pointer :: aero_rep_ptr
    real(kind=dp), allocatable, dimension(:,:) :: model_conc, true_conc
    integer(kind=i_kind) :: idx_phase, idx_aero_rep
    integer(kind=i_kind) :: idx_ethanol, idx_ethanol_aq, idx_H2O_aq
    integer(kind=i_kind) :: i_time, i_spec
    real(kind=dp) :: time_step, time
    real(kind=dp) :: n_star, del_H, del_S, del_G, alpha, crms, ugm3_to_ppm
    real(kind=dp) :: VP_ethanol, k_forward, k_backward
    real(kind=dp) :: equil_ethanol, equil_ethanol_aq
    real(kind=dp) :: total_mass, water_mass, VP_0_mass
#ifdef PMC_USE_MPI
    character, allocatable :: buffer(:), buffer_copy(:)
    integer(kind=i_kind) :: pack_size, pos, i_elem, results
#endif

    ! Parameters for calculating true concentrations
    real(kind=dp) :: temperature, pressure
    real(kind=dp), target :: radius, number_conc
    real(kind=dp), parameter :: MW_ethanol = 0.04607
    real(kind=dp), parameter :: MW_H2O = 0.01801
    real(kind=dp), parameter :: Dg_ethanol = 9.50d-6

    ! For setting particle radius and number concentration
    type(aero_rep_factory_t) :: aero_rep_factory
    type(aero_rep_update_data_single_particle_radius_t) :: radius_update
    type(aero_rep_update_data_single_particle_number_t) :: number_update

    ! For setting the GSD and GMD for modes
    integer(kind=i_kind) :: i_sect_unused, i_sect_the_mode
    type(aero_rep_update_data_modal_binned_mass_GMD_t) :: update_data_GMD
    type(aero_rep_update_data_modal_binned_mass_GSD_t) :: update_data_GSD

    type(solver_stats_t), target :: solver_stats

    call assert_msg(378317062, scenario.ge.1 .and. scenario.le.2, &
                    "Invalid scenario specified: "//to_string( scenario ) )

    run_SIMPOL_phase_transfer_test = .true.

    ! Allocate space for the results
    if (scenario.eq.1) then
      allocate(model_conc(0:NUM_TIME_STEP, 7))
      allocate(true_conc(0:NUM_TIME_STEP, 7))
    else if (scenario.eq.2) then
      allocate(model_conc(0:NUM_TIME_STEP, 6))
      allocate(true_conc(0:NUM_TIME_STEP, 6))
    endif

    ! Set the environmental conditions
    temperature = 272.5d0       ! temperature (K)
    pressure = 101253.3d0       ! pressure (Pa)

    ! Calculate the SIMPOL vapor pressure for ethanol (atm)
    VP_ethanol = 10.0d0**( -1.97E+03 / temperature &
                         + 2.91E+00 &
                         + 1.96E-03 * temperature &
                         - 4.96E-01 * log(temperature) )

    ! Set output time step (s)
    time_step = 1.0d-11

#ifdef PMC_USE_MPI
    ! Load the model data on root process and pass it to process 1 for solving
    if (pmc_mpi_rank().eq.0) then
#endif

      ! Get the SIMPOL_phase_transfer reaction mechanism json file
      if (scenario.eq.1) then
        input_file_path = 'test_SIMPOL_phase_transfer_config.json'
      else if (scenario.eq.2) then
        input_file_path = 'test_SIMPOL_phase_transfer_config_2.json'
      endif

      ! Construct a camp_core variable
      camp_core => camp_core_t(input_file_path)

      deallocate(input_file_path)

      ! Initialize the model
      call camp_core%initialize()

      ! Find the aerosol representation
      key ="my aero rep 2"
      call assert(209301925, camp_core%get_aero_rep(key, aero_rep_ptr))
      if (scenario.eq.1) then
        select type (aero_rep_ptr)
          type is (aero_rep_single_particle_t)
            call aero_rep_factory%initialize_update_data( aero_rep_ptr, &
                                                          radius_update)
            call aero_rep_factory%initialize_update_data( aero_rep_ptr, &
                                                          number_update)
          class default
            call die_msg(261298847, "Incorrect aerosol representation type")
        end select
      else if (scenario.eq.2) then
        select type (aero_rep_ptr)
          type is (aero_rep_modal_binned_mass_t)
            call aero_rep_factory%initialize_update_data( aero_rep_ptr, &
                                                          update_data_GMD)
            call aero_rep_factory%initialize_update_data( aero_rep_ptr, &
                                                          update_data_GSD)
            call assert_msg(883833294, &
                  aero_rep_ptr%get_section_id("unused mode", i_sect_unused), &
                  "Could not get section id for the unused mode")
            call assert_msg(431201141, &
                  aero_rep_ptr%get_section_id("the mode", i_sect_the_mode), &
                  "Could not get section id for the unused mode")
          class default
            call die_msg(304136933, "Incorrect aerosol representation type")
        end select
      endif

      ! Get chemical species data
      call assert(250292358, camp_core%get_chem_spec_data(chem_spec_data))

      ! Get species indices
      if (scenario.eq.1) then
        idx_prefix = ""
      else if (scenario.eq.2) then
        idx_prefix = "the mode."
      endif
      key = "ethanol"
      idx_ethanol = chem_spec_data%gas_state_id(key);
      key = idx_prefix//"aqueous aerosol.ethanol_aq"
      idx_ethanol_aq = aero_rep_ptr%spec_state_id(key);
      key = idx_prefix//"aqueous aerosol.H2O_aq"
      idx_H2O_aq = aero_rep_ptr%spec_state_id(key);

      ! Make sure the expected species are in the model
      call assert(884352514, idx_ethanol.gt.0)
      call assert(379146109, idx_ethanol_aq.gt.0)
      call assert(208989205, idx_H2O_aq.gt.0)

#ifdef PMC_USE_MPI
      ! pack the camp core
      pack_size = camp_core%pack_size()
      if (scenario.eq.1) then
        pack_size = pack_size &
                  + radius_update%pack_size() &
                  + number_update%pack_size()
      else if (scenario.eq.2) then
        pack_size = pack_size &
                  + update_data_GMD%pack_size() &
                  + update_data_GSD%pack_size()
      end if
      allocate(buffer(pack_size))
      pos = 0
      call camp_core%bin_pack(buffer, pos)
      if (scenario.eq.1) then
        call radius_update%bin_pack(buffer, pos)
        call number_update%bin_pack(buffer, pos)
      else if (scenario.eq.2) then
        call update_data_GMD%bin_pack(buffer, pos)
        call update_data_GSD%bin_pack(buffer, pos)
      end if
      call assert(325888214, pos.eq.pack_size)
    end if

    ! broadcast the species ids
    call pmc_mpi_bcast_integer(idx_ethanol)
    call pmc_mpi_bcast_integer(idx_ethanol_aq)
    call pmc_mpi_bcast_integer(idx_H2O_aq)
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
      if (scenario.eq.1) then
        call radius_update%bin_unpack(buffer, pos)
        call number_update%bin_unpack(buffer, pos)
      else if (scenario.eq.2) then
        call update_data_GMD%bin_unpack(buffer, pos)
        call update_data_GSD%bin_unpack(buffer, pos)
      end if
      call assert(320623907, pos.eq.pack_size)
      allocate(buffer_copy(pack_size))
      pos = 0
      call camp_core%bin_pack(buffer_copy, pos)
      if (scenario.eq.1) then
        call radius_update%bin_pack(buffer_copy, pos)
        call number_update%bin_pack(buffer_copy, pos)
      else if (scenario.eq.2) then
        call update_data_GMD%bin_pack(buffer_copy, pos)
        call update_data_GSD%bin_pack(buffer_copy, pos)
      end if
      call assert(432942252, pos.eq.pack_size)
      do i_elem = 1, pack_size
        call assert_msg(827735846, buffer(i_elem).eq.buffer_copy(i_elem), &
                "Mismatch in element :"//trim(to_string(i_elem)))
      end do

      ! solve and evaluate results on process 1
#endif

      ! Initialize the solver
      call camp_core%solver_initialize()

      ! Get a model state variable
      camp_state => camp_core%new_state()

      ! Set the environmental conditions
      camp_state%env_state%temp = temperature
      camp_state%env_state%pressure = pressure
      call camp_state%update_env_state()

      ! Save the initial concentrations
      true_conc(:,:) = 0.0
      true_conc(0,idx_ethanol) = 1.0e-3
      true_conc(0,idx_ethanol_aq) = 1.0e-5
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

      ! Update the aerosol representation (single partile only)
      if (scenario.eq.1) then
        call radius_update%set_radius(radius)
        call number_update%set_number(number_conc)
        call camp_core%update_aero_rep_data(radius_update)
        call camp_core%update_aero_rep_data(number_update)
      end if

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

      ! ethanol rate constants
      n_star = 2.55d0
      del_H = - 10.0d0 * (n_star-1.0d0) + 7.53*(n_star**(2.0d0/3.0d0)-1.0d0) &
              - 1.0d0
      del_S = - 32.0d0 * (n_star-1.0d0) + 9.21*(n_star**(2.0d0/3.0d0)-1.0d0) &
              - 1.3d0
      del_G = (del_H - temperature * del_S/1000.0d0) * 4184.0d0
      alpha = exp(-del_G/(const%univ_gas_const*temperature))
      alpha = alpha / (1.0d0 + alpha)
      crms = sqrt(8.0d0*const%univ_gas_const*temperature/(const%pi*46.07d0))
      k_forward = number_conc * ((radius**2 / (3.0d0 * Dg_ethanol) + &
              4.0d0 * radius / (3.0d0 * crms * alpha))**(-1))     ! (1/s)

      ! Determine the equilibrium concentrations
      ! [A_gas] =  VP_ethanol
      ! [A_aero] = [A_total] - VP_ethanol
      ugm3_to_ppm = const%univ_gas_const * temperature / (46.07d0 * pressure)
      total_mass = true_conc(0,idx_ethanol)/ugm3_to_ppm + &
              true_conc(0,idx_ethanol_aq)*number_conc ! (ug/m3)

      ! Iterated to find equil_ethanol_aq
      equil_ethanol_aq = 1.48804181d1 ! (ug/m3)
      equil_ethanol = (total_mass-equil_ethanol_aq)*ugm3_to_ppm
      equil_ethanol_aq = equil_ethanol_aq/number_conc

      ! Calculate the backwards rate constant based on the equilibrium
      ! conditions and the forward rate (1/s)
      k_backward = k_forward / ( equil_ethanol/ugm3_to_ppm/equil_ethanol_aq )

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
        call assert_msg(173108608, solver_stats%Jac_eval_fails.eq.0, &
                        trim( to_string( solver_stats%Jac_eval_fails ) )// &
                        " Jacobian evaluation failures at time step "// &
                        trim( to_string( i_time ) ) )
#endif

        ! Get the analytic conc
        ! x = [A_gas] - [A_eq_gas]
        ! x0 = [A_init_gas] - [A_eq_gas]
        ! [A_gas] = x + [A_eq_gas] = x0exp(-t/tau) + [A_eq_gas]
        ! 1/tau = k_f + k_b
        ! [A_gas] = ([A_init_gas] - [A_eq_gas]) * exp(-t *(k_f + k_b)) +
        !                 [A_eq_gas]
        ! [A_aero] = ([A_init_aero] - [A_eq_aero]) * exp(-t * (k_f + k_b)) +
        !                 [A_eq_aero]
        time = i_time * time_step
        true_conc(i_time,idx_ethanol) = &
                (true_conc(0,idx_ethanol) - equil_ethanol) * &
                exp(-time * (k_forward + k_backward)) + equil_ethanol
        true_conc(i_time,idx_ethanol_aq) = &
                (true_conc(0,idx_ethanol_aq) - equil_ethanol_aq) * &
                exp(-time * (k_forward + k_backward)) + equil_ethanol_aq
        true_conc(i_time,idx_H2O_aq) = true_conc(0,idx_H2O_aq)
      end do

      ! Save the results
      if (scenario.eq.1) then
        open(unit=7, file="out/SIMPOL_phase_transfer_results.txt", &
                status="replace", action="write")
      else if (scenario.eq.2) then
        open(unit=7, file="out/SIMPOL_phase_transfer_results_2.txt", &
                status="replace", action="write")
      endif
      do i_time = 0, NUM_TIME_STEP
        write(7,*) i_time*time_step, &
              ' ', true_conc(i_time, idx_ethanol),        &
              ' ', model_conc(i_time, idx_ethanol),       &
              ' ', true_conc(i_time, idx_ethanol_aq),     &
              ' ', model_conc(i_time, idx_ethanol_aq),    &
              ' ', true_conc(i_time, idx_H2O_aq),         &
              ' ', model_conc(i_time, idx_H2O_aq)
      end do
      close(7)

      ! Analyze the results (single particle only)
      if (scenario.eq.1) then
        do i_time = 1, NUM_TIME_STEP
          do i_spec = 1, 5
            ! Only check the second phase
            if (i_spec.ge.2.and.i_spec.le.3) cycle
            call assert_msg(237580431, &
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
      if (run_SIMPOL_phase_transfer_test) then
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
        run_SIMPOL_phase_transfer_test = .true.
      else
        run_SIMPOL_phase_transfer_test = .false.
      end if
    end if

    deallocate(buffer)
#endif

    deallocate(model_conc)
    deallocate(true_conc)
    deallocate(camp_core)

  end function run_SIMPOL_phase_transfer_test

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program pmc_test_SIMPOL_phase_transfer
