! Copyright (C) 2017 Matt Dawson
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_test_ZSR_aerosol_water program

!> Test of ZSR_aerosol_water sub model
program pmc_test_ZSR_aerosol_water

  use iso_c_binding

  use pmc_util,                         only: i_kind, dp, assert, &
                                              almost_equal, string_t, &
                                              warn_msg
  use pmc_camp_core
  use pmc_camp_state
  use pmc_chem_spec_data
  use pmc_aero_rep_data
  use pmc_aero_rep_factory
  use pmc_aero_rep_single_particle
  use pmc_solver_stats
#ifdef PMC_USE_JSON
  use json_module
#endif
  use pmc_mpi

  use iso_c_binding
  implicit none

  ! Number of RHs to calculate aerosol water for
  integer(kind=i_kind) :: NUM_RH_STEP = 100

  !> Interface to the c ODE solver and test functions
  interface
    !> Run the c functions tests
    integer(kind=c_int) function run_sub_model_zsr_c_tests(solver_data, &
        state, env) bind (c)
      use iso_c_binding
      !> Pointer to the initialized solver data
      type(c_ptr), value :: solver_data
      !> Pointer to the state array
      type(c_ptr), value :: state
      !> Pointer to the environmental state array
      type(c_ptr), value :: env
    end function run_sub_model_zsr_c_tests
  end interface

  ! initialize mpi
  call pmc_mpi_init()

  if (run_ZSR_aerosol_water_tests()) then
    if (pmc_mpi_rank().eq.0) write(*,*) &
            "ZSR aerosol water sub model tests - PASS"
  else
    if (pmc_mpi_rank().eq.0) write(*,*) &
            "ZSR aerosol water sub model tests - FAIL"
          stop 3
  end if

  ! finalize mpi
  call pmc_mpi_finalize()

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Run all pmc_chem_mech_solver tests
  logical function run_ZSR_aerosol_water_tests() result(passed)

    use pmc_camp_solver_data

    type(camp_solver_data_t), pointer :: camp_solver_data

    camp_solver_data => camp_solver_data_t()

    if (camp_solver_data%is_solver_available()) then
      passed = run_ZSR_aerosol_water_test()
    else
      call warn_msg(713064651, "No solver available")
      passed = .true.
    end if

    deallocate(camp_solver_data)

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

    type(camp_core_t), pointer :: camp_core
    type(camp_state_t), pointer :: camp_state
    character(len=:), allocatable :: input_file_path, key
    type(string_t), allocatable, dimension(:) :: output_file_path

    type(chem_spec_data_t), pointer :: chem_spec_data
    class(aero_rep_data_t), pointer :: aero_rep_ptr
    real(kind=dp), dimension(0:NUM_RH_STEP, 13) :: model_conc, true_conc
    integer(kind=i_kind) :: idx_H2O, idx_Na_p, idx_Na_p_act, idx_Cl_m, &
            idx_Cl_m_act, idx_Ca_pp, idx_Ca_pp_act, idx_H2O_aq, idx_H2O_act, &
            i_RH, i_spec, idx_phase
    real(kind=dp) :: RH_step, RH, ppm_to_RH, molal_NaCl, molal_CaCl2, &
            NaCl_conc, CaCl2_conc, water_NaCl, water_CaCl2, temp, pressure
#ifdef PMC_USE_MPI
    character, allocatable :: buffer(:), buffer_copy(:)
    integer(kind=i_kind) :: pack_size, pos, i_elem, results
#endif

    type(solver_stats_t), target :: solver_stats

    run_ZSR_aerosol_water_test = .true.

    ! Set the environmental and aerosol test conditions
    temp = 272.5d0              ! temperature (K)
    pressure = 101253.3d0       ! pressure (Pa)

    ! Set RH step (unitless)
    RH_step = 1.0/NUM_RH_STEP

#ifdef PMC_USE_MPI
    ! Load the model data on root process and pass it to process 1 for solving
    if (pmc_mpi_rank().eq.0) then
#endif

      ! Get the ZSR_aerosol_water sub model mechanism json file
      input_file_path = 'test_ZSR_aerosol_water_config.json'

      ! Construct a camp_core variable
      camp_core => camp_core_t(input_file_path)

      deallocate(input_file_path)

      ! Initialize the model
      call camp_core%initialize()

      ! Get the chemical species data
      call assert(604467956, camp_core%get_chem_spec_data(chem_spec_data))

      ! Find the aerosol representation
      key = "my aero rep 2"
      call assert(110830690, camp_core%get_aero_rep(key, aero_rep_ptr))

      ! Get species indices
      key = "H2O"
      idx_H2O = chem_spec_data%gas_state_id(key);
      key = "aqueous aerosol.H2O_aq"
      idx_H2O_aq = aero_rep_ptr%spec_state_id(key);
      key = "aqueous aerosol.Na_p"
      idx_Na_p = aero_rep_ptr%spec_state_id(key);
      key = "aqueous aerosol.Cl_m"
      idx_Cl_m = aero_rep_ptr%spec_state_id(key);
      key = "aqueous aerosol.Ca_pp"
      idx_Ca_pp = aero_rep_ptr%spec_state_id(key);

      ! Make sure the expected species are in the model
      call assert(213525011, idx_H2O.gt.0)
      call assert(943368106, idx_H2O_aq.gt.0)
      call assert(202996397, idx_Na_p.gt.0)
      call assert(427633087, idx_Cl_m.gt.0)
      call assert(317220276, idx_Ca_pp.gt.0)

#ifdef PMC_USE_MPI
      ! pack the camp core
      pack_size = camp_core%pack_size()
      allocate(buffer(pack_size))
      pos = 0
      call camp_core%bin_pack(buffer, pos)
      call assert(346780270, pos.eq.pack_size)
    end if

    ! broadcast the species ids
    call pmc_mpi_bcast_integer(idx_H2O)
    call pmc_mpi_bcast_integer(idx_H2O_aq)
    call pmc_mpi_bcast_integer(idx_Na_p)
    call pmc_mpi_bcast_integer(idx_Cl_m)
    call pmc_mpi_bcast_integer(idx_Ca_pp)

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
      call assert(401260056, pos.eq.pack_size)
      allocate(buffer_copy(pack_size))
      pos = 0
      call camp_core%bin_pack(buffer_copy, pos)
      call assert(796053650, pos.eq.pack_size)
      do i_elem = 1, pack_size
        call assert_msg(625896746, buffer(i_elem).eq.buffer_copy(i_elem), &
                "Mismatch in element :"//trim(to_string(i_elem)))
      end do

      ! solve and evaluate results on process 1
#endif

      ! Initialize the solver
      call camp_core%solver_initialize()

      ! Get a model state variable
      camp_state => camp_core%new_state()

      ! Set the environmental conditions
      call camp_state%env_states(1)%set_temperature_K(   temp )
      call camp_state%env_states(1)%set_pressure_Pa( pressure )

      ! Save the initial concentrations
      true_conc(:,:) = 0.0
      true_conc(:,idx_H2O) = 0.0
      true_conc(:,idx_H2O_aq) = 0.0
      true_conc(:,idx_Na_p) = 2.5
      true_conc(:,idx_Cl_m) = 5.3
      true_conc(:,idx_Ca_pp) = 1.3
      model_conc(:,:) = true_conc(:,:)

      ! Set up the ppm->RH (0-1) conversion
      ! (From MOSAIC code, references Seinfeld and Pandis pg. 181)
      ppm_to_RH = 1.0d0 - 373.15d0/temp
      ppm_to_RH = (((-0.1299d0*ppm_to_RH - 0.6445d0)*ppm_to_RH - 1.976d0)* &
              ppm_to_RH + 13.3185d0)*ppm_to_RH
      ppm_to_RH = exp(ppm_to_RH)  ! VP of water (atm)
      ppm_to_RH = (pressure/101325.0d0) / ppm_to_RH * 1.0d-6 ! ppm -> RH (0-1)

#ifdef PMC_DEBUG
      ! Evaluate the Jacobian during solving
      solver_stats%eval_Jac = .true.
#endif

      ! Integrate the mechanism
      do i_RH = 1, NUM_RH_STEP

        ! Calculate the current [H20] (ppm)
        true_conc(i_RH,idx_H2O) = i_RH * RH_step / ppm_to_RH
        model_conc(i_RH,idx_H2O) = true_conc(i_RH, idx_H2O)

        ! Set the initial state in the model
        camp_state%state_var(:) = model_conc(i_RH,:)

        ! Get the modeled conc
        ! time step is arbitrary - equilibrium calculatuions only
        call camp_core%solve(camp_state, real(1.0, kind=dp), &
                              solver_stats = solver_stats)
        model_conc(i_RH,:) = camp_state%state_var(:)

#ifdef PMC_DEBUG
        ! Check the Jacobian evaluations
        call assert_msg(961242557, solver_stats%Jac_eval_fails.eq.0, &
                        trim( to_string( solver_stats%Jac_eval_fails ) )// &
                        " Jacobian evaluation failures for i_RH "// &
                        trim( to_string( i_RH ) ) )
#endif

        ! Get the analytic conc
        ! Jacobson molality (eq. 29 in \cite{Jacobson1996}) :
        !  m_a^(1/2) = Y0 + Y1*aw + Y2*aw^2 + ... (aw = RH)
        RH = i_RH * RH_step
        RH = MAX(RH, 0.43d0)
        molal_CaCl2 = -1.918004e2 +  2.001540e3 * RH - 8.557205e3 * RH**2 + &
              1.987670e4 * RH**3 - 2.717192e4 * RH**4 + 2.187103e4 * RH**5 - &
              9.591577e3 * RH**6 + 1.763672e3 * RH**7
        molal_CaCl2 = molal_CaCl2**2
        ! EQSAM molality (from EQSAM_v03d)
        !  m_i = (NW_i * MW_H2O/MW_i * (1.0/RH-1.0))^ZW_i
        !  where MW_H2O is defined as 55.51*18.01
        RH = i_RH * RH_step
        molal_NaCl = (2.0d0 * 55.51D0 * 18.01d0 / 58.5d0 * &
                (1.0d0/RH - 1.0d0))**0.67d0
        CaCl2_conc = MIN(true_conc(i_RH,idx_Ca_pp)/40.078d0, &
                true_conc(i_RH,idx_Cl_m)/2.0d0/35.453d0)
                ! (umol/m^3_air = mol/cm^3_air)
        NaCl_conc = true_conc(i_RH,idx_Cl_m)/35.453d0
                ! (umol/m^3_air = mol/cm^3_air)
        ! Water content is (eq 28 \cite{Jacobson1996}) :
        ! cw = 1000 / MW_H2O * sum_i (c_i / m_i)
        !   with cw and c_i in (mol_i/cm^3_air) and m_i in (mol_i/kg_H2O)
        water_CaCl2 =  CaCl2_conc / molal_CaCl2 * 1000.0d0 ! (ug_H2O/m^3_air)
        water_NaCl = NaCl_conc / molal_NaCl * 1000.0d0 ! (ug_H2O/m^3_air)
        true_conc(i_RH,idx_H2O_aq) = water_CaCl2 + water_NaCl

      end do

      ! Save the results
      open(unit=7, file="out/ZSR_aerosol_water_results.txt", status="replace", &
              action="write")
      do i_RH = 0, NUM_RH_STEP
        write(7,*) i_RH*RH_step, &
              ' ',true_conc(i_RH, idx_H2O),' ',   model_conc(i_RH, idx_H2O), &
              ' ',true_conc(i_RH, idx_H2O_aq),' ',model_conc(i_RH, idx_H2O_aq),&
              ' ',true_conc(i_RH, idx_Na_p),' ',  model_conc(i_RH, idx_Na_p), &
              ' ',true_conc(i_RH, idx_Cl_m),' ',  model_conc(i_RH, idx_Cl_m), &
              ' ',true_conc(i_RH, idx_Ca_pp),' ', model_conc(i_RH, idx_Ca_pp)
      end do
      close(7)

      ! Analyze the results
      do i_RH = 1, NUM_RH_STEP
        do i_spec = 1, 13
          ! Skip the first aerosol phase
          if (i_spec.ge.2.and.i_spec.le.9) cycle
          call assert_msg(221026833, &
            almost_equal(model_conc(i_RH, i_spec), &
            true_conc(i_RH, i_spec), real(1.0e-2, kind=dp)).or. &
            (model_conc(i_RH, i_spec).lt.1e-5*model_conc(1, i_spec).and. &
            true_conc(i_RH, i_spec).lt.1e-5*true_conc(1, i_spec)), &
            "time: "//trim(to_string(i_RH))//"; species: "// &
            trim(to_string(i_spec))//"; mod: "// &
            trim(to_string(model_conc(i_RH, i_spec)))//"; true: "// &
            trim(to_string(true_conc(i_RH, i_spec))))
        end do
      end do

      deallocate(camp_state)

#ifdef PMC_USE_MPI
      ! convert the results to an integer
      if (run_ZSR_aerosol_water_test) then
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
        run_ZSR_aerosol_water_test = .true.
      else
        run_ZSR_aerosol_water_test = .false.
      end if
    end if

    deallocate(buffer)
#endif

    deallocate(camp_core)

    ! Evaluate the sub model c functions
    run_ZSR_aerosol_water_test = &
      run_ZSR_aerosol_water_test .and. &
      eval_c_func()

  end function run_ZSR_aerosol_water_test

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Evaluate the sub model c functions
  logical function eval_c_func() result(passed)

    type(camp_core_t), pointer :: camp_core
    type(camp_state_t), pointer :: camp_state
    character(len=:), allocatable :: input_file_path

    ! Get the ZSR_aerosol_water sub model mechanism json file
    input_file_path = 'test_ZSR_aerosol_water_config.json'

    ! Construct a camp_core variable
    camp_core => camp_core_t(input_file_path)

    deallocate(input_file_path)

    ! Initialize the model
    call camp_core%initialize()

    ! Initialize the solver
    call camp_core%solver_initialize()

    ! Get a new state variable
    camp_state => camp_core%new_state()

    ! Set the initial conditions
    camp_state%state_var(:) = 0.0
    call camp_state%env_states(1)%set_temperature_K(  298.0d0 )
    call camp_state%env_states(1)%set_pressure_Pa( 101325.0d0 )
    call camp_state%update_env_state( )

    passed = run_sub_model_zsr_c_tests(                                      &
                 camp_core%solver_data_gas_aero%solver_c_ptr,                &
                 c_loc(camp_state%state_var),                                &
                 c_loc(camp_state%env_var)                                   &
                 ) .eq. 0

   deallocate(camp_state)
   deallocate(camp_core)

  end function eval_c_func

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program pmc_test_ZSR_aerosol_water
