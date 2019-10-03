! Copyright (C) 2017 Matt Dawson
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_test_PDFiTE program

!> Test of PDFiTE sub model module
program pmc_test_PDFiTE

  use iso_c_binding

  use pmc_util,                         only: i_kind, dp, assert, &
                                              almost_equal, string_t, &
                                              warn_msg
  use pmc_camp_core
  use pmc_camp_state
  use pmc_aero_rep_data
  use pmc_aero_rep_factory
  use pmc_aero_rep_single_particle
  use pmc_chem_spec_data
  use pmc_solver_stats
#ifdef PMC_USE_JSON
  use json_module
#endif
  use pmc_mpi

  implicit none

  ! Number of timesteps to output in mechanisms
  integer(kind=i_kind) :: NUM_RH_STEP = 101

  ! initialize mpi
  call pmc_mpi_init()

  if (run_PDFiTE_tests()) then
    if (pmc_mpi_rank().eq.0) write(*,*) " PD-FiTE activity sub model tests - PASS"
  else
    if (pmc_mpi_rank().eq.0) write(*,*) " PD-FiTE activity sub model tests - FAIL"
    stop 3
  end if

  ! finalize mpi
  call pmc_mpi_finalize()

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Run all pmc_chem_mech_solver tests
  logical function run_PDFiTE_tests() result(passed)

    use pmc_camp_solver_data

    type(camp_solver_data_t), pointer :: camp_solver_data

    camp_solver_data => camp_solver_data_t()

    if (camp_solver_data%is_solver_available()) then
      passed = run_PDFiTE_test()
    else
      call warn_msg(713064651, "No solver available")
      passed = .true.
    end if

    deallocate(camp_solver_data)

  end function run_PDFiTE_tests

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Solve a PD-FiTE activity coefficient calculation
  !!
  !! The scenario matches that described for the H+ - NH4+ - SO42- - NO3-
  !! system described by equations 16 and 17 in \cite{Topping2009}.
  logical function run_PDFiTE_test()

    use pmc_constants

    type(camp_core_t), pointer :: camp_core
    type(camp_state_t), pointer :: camp_state
    character(len=:), allocatable :: input_file_path, key
    type(string_t), allocatable, dimension(:) :: output_file_path

    type(chem_spec_data_t), pointer :: chem_spec_data
    class(aero_rep_data_t), pointer :: aero_rep_ptr
    real(kind=dp), dimension(0:NUM_RH_STEP, 28) :: model_conc, true_conc
    integer(kind=i_kind) :: idx_H2O, idx_H2O_aq, idx_H_p, idx_NH4_p, &
            idx_SO4_mm, idx_NO3_m, idx_NH42_SO4, idx_NH4_NO3, idx_H_NO3, &
            idx_H2_SO4,  idx_phase, i_RH, i_spec
    real(kind=dp) :: time_step, time, ppm_to_RH, omega, ln_gamma, &
            HNO3_LRH_B0, HNO3_LRH_B1, HNO3_LRH_B2, HNO3_LRH_B3, HNO3_LRH_B4, &
            HNO3_HRH_B0, HNO3_HRH_B1, HNO3_HRH_B2, HNO3_HRH_B3, &
            HNO3_SHRH_B0, HNO3_SHRH_B1, HNO3_SHRH_B2, HNO3_SHRH_B3, &
            H2SO4_B0, H2SO4_B1, H2SO4_B2, H2SO4_B3, NH42SO4_B0, NH42SO4_B1, &
            NH42SO4_B2, NH42SO4_B3, NH4NO3_B0, NH4NO3_B1, NH4NO3_B2, &
            NH4NO3_B3, H_p_mol_m3, NH4_p_mol_m3, SO4_mm_mol_m3, &
            NO3_m_mol_m3, temp, pressure, a_w
#ifdef PMC_USE_MPI
    character, allocatable :: buffer(:), buffer_copy(:)
    integer(kind=i_kind) :: pack_size, pos, i_elem, results
#endif

    type(solver_stats_t), target :: solver_stats

    run_PDFiTE_test = .true.

    ! Set the environmental and aerosol test conditions
    temp = 272.5d0              ! temperature (K)
    pressure = 101253.3d0       ! pressure (Pa)

    ! Set output time step (s)
    time_step = 1.0d0

#ifdef PMC_USE_MPI
    ! Load the model data on root process and pass it to process 1 for solving
    if (pmc_mpi_rank().eq.0) then
#endif

      ! Get the PDFiTE sub model mechanism json file
      input_file_path = 'test_PDFiTE_config.json'

      ! Construct a camp_core variable
      camp_core => camp_core_t(input_file_path)

      deallocate(input_file_path)

      ! Initialize the model
      call camp_core%initialize()

      ! Get the chemical species data
      call assert(585768557, camp_core%get_chem_spec_data(chem_spec_data))

      ! Find the aerosol representation
      key = "my aero rep 2"
      call assert(917299782, camp_core%get_aero_rep(key, aero_rep_ptr))

      ! Get species indices
      key = "H2O"
      idx_H2O = chem_spec_data%gas_state_id(key);
      key = "aqueous aerosol.H2O_aq"
      idx_H2O_aq = aero_rep_ptr%spec_state_id(key);
      key = "aqueous aerosol.H_p"
      idx_H_p = aero_rep_ptr%spec_state_id(key);
      key = "aqueous aerosol.NH4_p"
      idx_NH4_p = aero_rep_ptr%spec_state_id(key);
      key = "aqueous aerosol.SO4_mm"
      idx_SO4_mm = aero_rep_ptr%spec_state_id(key);
      key = "aqueous aerosol.NO3_m"
      idx_NO3_m = aero_rep_ptr%spec_state_id(key);
      key = "aqueous aerosol.(NH4)2-SO4"
      idx_NH42_SO4 = aero_rep_ptr%spec_state_id(key);
      key = "aqueous aerosol.NH4-NO3"
      idx_NH4_NO3 = aero_rep_ptr%spec_state_id(key);
      key = "aqueous aerosol.H-NO3"
      idx_H_NO3 = aero_rep_ptr%spec_state_id(key);
      key = "aqueous aerosol.H2-SO4"
      idx_H2_SO4 = aero_rep_ptr%spec_state_id(key);

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

#ifdef PMC_USE_MPI
      ! pack the camp core
      pack_size = camp_core%pack_size()
      allocate(buffer(pack_size))
      pos = 0
      call camp_core%bin_pack(buffer, pos)
      call assert(304632244, pos.eq.pack_size)
    end if

    ! broadcast the species ids
    call pmc_mpi_bcast_integer(idx_H2O)
    call pmc_mpi_bcast_integer(idx_H2O_aq)
    call pmc_mpi_bcast_integer(idx_H_p)
    call pmc_mpi_bcast_integer(idx_NH4_p)
    call pmc_mpi_bcast_integer(idx_SO4_mm)
    call pmc_mpi_bcast_integer(idx_NO3_m)
    call pmc_mpi_bcast_integer(idx_NH42_SO4)
    call pmc_mpi_bcast_integer(idx_NH4_NO3)
    call pmc_mpi_bcast_integer(idx_H_NO3)
    call pmc_mpi_bcast_integer(idx_H2_SO4)

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
      call assert(359112030, pos.eq.pack_size)
      allocate(buffer_copy(pack_size))
      pos = 0
      call camp_core%bin_pack(buffer_copy, pos)
      call assert(471430375, pos.eq.pack_size)
      do i_elem = 1, pack_size
        call assert_msg(414489565, buffer(i_elem).eq.buffer_copy(i_elem), &
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
      true_conc(:,:)            = 0.0d0
      true_conc(:,idx_H2O)      = 0.0d0
      true_conc(:,idx_H2O_aq)   = 1.0d0 ! TODO remove from sub model parameters
      true_conc(:,idx_H_p)      = 1.0d-9
      true_conc(:,idx_NH4_p)    = 0.62d0
      true_conc(:,idx_SO4_mm)   = 0.25d0
      true_conc(:,idx_NO3_m)    = 0.13d0
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

      ! H-NO3 0.9 - 0.99 RH
      HNO3_SHRH_B0 = -1420.5
      HNO3_SHRH_B1 = 4467.9
      HNO3_SHRH_B2 = -4682.7
      HNO3_SHRH_B3 = 1635.1

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
      ppm_to_RH = (((-0.1299d0*ppm_to_RH - 0.6445d0)*ppm_to_RH - 1.976d0)* &
              ppm_to_RH + 13.3185d0)*ppm_to_RH
      ppm_to_RH = exp(ppm_to_RH)  ! VP of water (atm)
      ppm_to_RH = (pressure/101325.0d0) / ppm_to_RH * 1.0d-6 ! ppm -> RH (0-1)

      ! Set the initial state in the model
      camp_state%state_var(:) = model_conc(0,:)

#ifdef PMC_DEBUG
      ! Evaluate the Jacobian during solving
      solver_stats%eval_Jac = .true.
#endif

      ! Integrate the mechanism
      do i_RH = 1, NUM_RH_STEP

        ! Set the RH (Stay slightly away from transitions to avoid
        ! discrepancies in RH calc from [H2O]_g (ppm)
        a_w = real(1.0, kind=dp)/(NUM_RH_STEP-1)*(i_RH-1) + 1.0d-10
        true_conc(i_RH, idx_H2O) = a_w / ppm_to_RH
        camp_state%state_var(idx_H2O) = true_conc(i_RH, idx_H2O)

        ! Get the modeled conc
        call camp_core%solve(camp_state, time_step, &
                              solver_stats = solver_stats)
        model_conc(i_RH,:) = camp_state%state_var(:)

#ifdef PMC_DEBUG
        ! Check the Jacobian evaluations
        call assert_msg(404462844, solver_stats%Jac_eval_fails.eq.0, &
                        trim( to_string( solver_stats%Jac_eval_fails ) )// &
                        " Jacobian evaluation failures for i_RH "// &
                        trim( to_string( i_RH ) ) )
#endif

        ! Calculate the mean binary activity for H-NO3

        ! Convert concentrations to ug/m^3 -> mol/m^3
        ! (ug/m^3) * (ug/umol) * (umol/mol)
        H_p_mol_m3    = true_conc(i_RH, idx_H_p)    / 1.008d0   * 1.0d-6
        NH4_p_mol_m3  = true_conc(i_RH, idx_NH4_p)  / 18.04d0   * 1.0d-6
        SO4_mm_mol_m3 = true_conc(i_RH, idx_SO4_mm) / 96.06d0   * 1.0d-6
        NO3_m_mol_m3  = true_conc(i_RH, idx_NO3_m)  / 62.0049d0 * 1.0d-6

        ! Calculate omega
        omega = 6.0d0 * H_p_mol_m3   * SO4_mm_mol_m3 + &
                6.0d0 * NH4_p_mol_m3 * SO4_mm_mol_m3 + &
                4.0d0 * NH4_p_mol_m3 * NO3_m_mol_m3

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
        else if (a_w.le.0.99d0) then
          ln_gamma = HNO3_SHRH_B0 + HNO3_SHRH_B1*a_w + HNO3_SHRH_B2*a_w**2 + &
                      HNO3_SHRH_B3*a_w**3
        else
          ln_gamma = HNO3_SHRH_B0 + HNO3_SHRH_B1*0.99d0 + &
                      HNO3_SHRH_B2*0.99d0**2 + HNO3_SHRH_B3*0.99d0**3
        end if

        ! ... from (d(ln(gamma_HNO3))/d(N_Hp N_SO4mm) N_Hp N_SO4mm) / omega
        if (a_w.le.0.1d0) then
          ln_gamma = ln_gamma + H_p_mol_m3 * SO4_mm_mol_m3 / omega * &
                     ( H2SO4_B0 + H2SO4_B1*0.1d0 + H2SO4_B2*0.1d0**2 + &
                     H2SO4_B3*0.1d0**3 )
        else if (a_w.le.0.99d0) then
          ln_gamma = ln_gamma + H_p_mol_m3 * SO4_mm_mol_m3 / omega * &
                     ( H2SO4_B0 + H2SO4_B1*a_w + H2SO4_B2*a_w**2 + &
                     H2SO4_B3*a_w**3 )
        else
        ln_gamma = ln_gamma + H_p_mol_m3 * SO4_mm_mol_m3 / omega * &
                   ( H2SO4_B0 + H2SO4_B1*0.99d0 + H2SO4_B2*0.99d0**2 + &
                   H2SO4_B3*0.99d0**3 )
        end if

        ! ... from (d(ln(gamma_HNO3))/d(N_NH4p N_SO4mm) N_NH4p N_SO4mm) / omega
        if (a_w.le.0.1d0) then
          ln_gamma = ln_gamma + NH4_p_mol_m3 * SO4_mm_mol_m3 / omega * &
                     ( NH42SO4_B0 + NH42SO4_B1*0.1d0 + NH42SO4_B2*0.1d0**2 + &
                     NH42SO4_B3*0.1d0**3 )
        else if (a_w.le.0.99d0) then
          ln_gamma = ln_gamma + NH4_p_mol_m3 * SO4_mm_mol_m3 / omega * &
                     ( NH42SO4_B0 + NH42SO4_B1*a_w + NH42SO4_B2*a_w**2 + &
                     NH42SO4_B3*a_w**3 )
        else
          ln_gamma = ln_gamma + NH4_p_mol_m3 * SO4_mm_mol_m3 / omega * &
                     ( NH42SO4_B0 + NH42SO4_B1*0.99d0 + NH42SO4_B2*0.99d0**2 + &
                     NH42SO4_B3*0.99d0**3 )
        end if

        ! ... from (d(ln(gamma_HNO3))/d(N_NH4p N_NO3m) N_NH4p N_NO3m) / omega
        if (a_w.le.0.1d0) then
          ln_gamma = ln_gamma + NH4_p_mol_m3 * NO3_m_mol_m3 / omega * &
                     ( NH4NO3_B0 + NH4NO3_B1*0.1d0 + NH4NO3_B2*0.1d0**2 + &
                     NH4NO3_B3*0.1d0**3 )
        else if (a_w.le.0.99d0) then
          ln_gamma = ln_gamma + NH4_p_mol_m3 * NO3_m_mol_m3 / omega * &
                     ( NH4NO3_B0 + NH4NO3_B1*a_w + NH4NO3_B2*a_w**2 + &
                     NH4NO3_B3*a_w**3 )
        else
          ln_gamma = ln_gamma + NH4_p_mol_m3 * NO3_m_mol_m3 / omega * &
                     ( NH4NO3_B0 + NH4NO3_B1*0.99d0 + NH4NO3_B2*0.99d0**2 + &
                     NH4NO3_B3*0.99d0**3 )
        end if

        ! Update the H-NO3 mean binary activity
        true_conc(i_RH, idx_H_NO3) = exp(ln_gamma)

      end do

      ! Save the results
      open(unit=7, file="out/PDFiTE_results.txt", status="replace", &
              action="write")
      do i_RH = 0, NUM_RH_STEP
        a_w = real(1.0, kind=dp)/(NUM_RH_STEP-1)*(i_RH-1) + 1.0d-10
        write(7,*) a_w, &
              ' ', true_conc(i_RH, idx_H2O),      &
              ' ', model_conc(i_RH, idx_H2O),     &
              ' ', true_conc(i_RH, idx_H2O_aq),   &
              ' ', model_conc(i_RH, idx_H2O_aq),  &
              ' ', true_conc(i_RH, idx_H_p),      &
              ' ', model_conc(i_RH, idx_H_p),     &
              ' ', true_conc(i_RH, idx_NH4_p),    &
              ' ', model_conc(i_RH, idx_NH4_p),   &
              ' ', true_conc(i_RH, idx_SO4_mm),   &
              ' ', model_conc(i_RH, idx_SO4_mm),  &
              ' ', true_conc(i_RH, idx_NO3_m),    &
              ' ', model_conc(i_RH, idx_NO3_m),   &
              ' ', true_conc(i_RH, idx_NH42_SO4), &
              ' ', model_conc(i_RH, idx_NH42_SO4),&
              ' ', true_conc(i_RH, idx_NH4_NO3),  &
              ' ', model_conc(i_RH, idx_NH4_NO3), &
              ' ', true_conc(i_RH, idx_H_NO3),    &
              ' ', model_conc(i_RH, idx_H_NO3),   &
              ' ', true_conc(i_RH, idx_H2_SO4),   &
              ' ', model_conc(i_RH, idx_H2_SO4)
      end do
      close(7)

      ! Analyze the results
      do i_RH = 1, NUM_RH_STEP
        do i_spec = 1, size(model_conc, 2)
          if (i_spec.ne.1.and.(i_spec.lt.11.or.i_spec.gt.19)) cycle
          call assert_msg(357068617, &
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
      if (run_PDFiTE_test) then
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
        run_PDFiTE_test = .true.
      else
        run_PDFiTE_test = .false.
      end if
    end if

    deallocate(buffer)
#endif

    deallocate(camp_core)

  end function run_PDFiTE_test

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program pmc_test_PDFiTE
