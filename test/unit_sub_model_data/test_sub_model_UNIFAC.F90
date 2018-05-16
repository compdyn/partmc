! Copyright (C) 2017 Matt Dawson
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_test_sub_module_UNIFAC program

!> Test of UNIFAC activity sub module
program pmc_test_sub_module_UNIFAC

  use pmc_util,                         only: i_kind, dp, assert, &
                                              almost_equal, string_t, &
                                              warn_msg
  use pmc_phlex_core
  use pmc_phlex_state
#ifdef PMC_USE_JSON
  use json_module
#endif
  use pmc_mpi

  implicit none

  ! Number of mole_fracsteps to output in mechanisms
  integer(kind=i_kind) :: NUM_MASS_FRAC_STEP = 100

  ! initialize mpi
  call pmc_mpi_init()

  if (run_UNIFAC_tests()) then
    write(*,*) "UNIFAC activity sub module tests - PASS"
  else
    write(*,*) "UNIFAC activity sub module tests - FAIL"
  end if

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Run all pmc_sub_module_UNIFAC tests
  logical function run_UNIFAC_tests() result(passed)

    use pmc_phlex_solver_data

    type(phlex_solver_data_t), pointer :: phlex_solver_data

    phlex_solver_data => phlex_solver_data_t()

    if (phlex_solver_data%is_solver_available()) then
      passed = run_UNIFAC_test()
    else
      call warn_msg(185492126, "No solver available")
      passed = .true.
    end if

  end function run_UNIFAC_tests

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Calculate UNIFAC acitivity coefficients for the n-butanol/water mixture
  !!
  !! Example is from Marcolli and Peter, ACP 5(2), 1501-1527, 2005. (fig 3a)
  !! 
  !! Equations in the test comments are from the same reference.
  !!
  logical function run_UNIFAC_test()

    use pmc_constants

    type(phlex_core_t), pointer :: phlex_core
    type(phlex_state_t), target :: phlex_state
    character(len=:), allocatable :: input_file_path
    type(string_t), allocatable, dimension(:) :: output_file_path

    real(kind=dp), dimension(0:NUM_MASS_FRAC_STEP, 6) :: model_conc, calc_conc
    real(kind=dp), dimension(0:NUM_MASS_FRAC_STEP, 6) :: model_activity, calc_activity
    integer(kind=i_kind) :: idx_butanol, idx_water
    character(len=:), allocatable :: key
    integer(kind=i_kind) :: i_mass_frac, i_spec, i_aero_rep
    real(kind=dp) :: mass_frac_step, mole_frac, mass_frac

    ! Parameters for calculating true concentrations
    real(kind=dp) :: temperature, pressure

    ! Molecular weights
    real(kind=dp), parameter :: mw_butanol = 74.12
    real(kind=dp), parameter :: mw_water = 18.01
    
    ! Number of functional groups
    integer(kind=i_kind), parameter :: num_group = 5

    ! Volume parameters (R_k in eq. 6)                       CH2(-OH) CH2(phb) CH3(phb)  OH     H2O
    real(kind=dp), parameter, dimension(num_group) :: R_k = [ 0.6744, 0.6744, 0.9011, 1.0000, 0.9200 ]
    ! Surface parameters (Q_k in eq. 6)
    real(kind=dp), parameter, dimension(num_group) :: Q_k = [ 0.540,  0.540,  0.848,  1.200,  1.400  ]

    ! Group interactions (a_mn in eq. 9)
    real(kind=dp), parameter, dimension(num_group, num_group) :: a_mn = reshape( [&
            ! CH2(-OH) CH2(phb) CH3(phb)      OH      H2O
                  0.0,     0.0,     0.0,   156.4,  -89.71,   & ! CH2(-OH)
                  0.0,     0.0,     0.0,   156.4,   362.1,   & ! CH2(phb)
                  0.0,     0.0,     0.0,   156.4,   362.1,   & ! CH3(phb)
                986.5,   986.5,   986.5,     0.0,  -153.0,   & ! OH
               2314.0,  1325.0,  1325.0,   276.4,     0.0 ], & ! H2O
            [num_group, num_group] )

    ! Number of groups in each species
    integer(kind=i_kind), parameter :: num_grps_butanol = 5
    integer(kind=i_kind), parameter :: num_grps_water = 1
    integer(kind=i_kind), parameter, dimension(num_group) :: butanol_grps =  [ 1, 2, 1, 1, 0 ]
    integer(kind=i_kind), parameter, dimension(num_group) :: water_grps   =  [ 0, 0, 0, 0, 1 ]

    ! Variables used in UNIFAC calculations
    real(kind=dp) :: r_butanol, r_water
    real(kind=dp) :: q_butanol, q_water
    real(kind=dp) :: l_butanol, l_water
    real(kind=dp), dimension(num_group, num_group) :: PSI_mn
    real(kind=dp), dimension(num_group) :: ln_GAMMA_butanol_pure, ln_GAMMA_water_pure
    real(kind=dp), dimension(num_group) :: ln_GAMMA_mixture
    real(kind=dp) :: sum_Qn_Xn_butanol, sum_Qn_Xn_water, sum_Qn_Xn_mixture
    real(kind=dp) :: PHI_butanol, PHI_water
    real(kind=dp) :: THETA_butanol, THETA_water
    real(kind=dp) :: ln_gamma_C_butanol, ln_gamma_C_water
    real(kind=dp) :: ln_gamma_R_butanol, ln_gamma_R_water

    real(kind=dp), dimension(num_group) :: THETA_m
    real(kind=dp) :: sum_m_A, sum_m_B, sum_n 

    integer(kind=i_kind) :: i, k, m, n

    run_UNIFAC_test = .true.

    ! Set the rate constants (for calculating the true value)
    temperature = 298.15d0
    pressure = 101253.3d0

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! Calculate mole frac-invariant parameters !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Calculate total molecule surface area and volume (eq. 6)
    r_butanol = 0.0d0
    r_water   = 0.0d0
    q_butanol = 0.0d0
    q_water   = 0.0d0
    do k = 1, num_group
      r_butanol = r_butanol + real(butanol_grps(k), kind=dp) * R_k(k)
      r_water   = r_water   + real(water_grps(k),   kind=dp) * R_k(k)
      q_butanol = q_butanol + real(butanol_grps(k), kind=dp) * Q_k(k)
      q_water   = q_water   + real(water_grps(k),   kind=dp) * Q_k(k)
    end do

    ! Calculate l_i (eq. 5)
    l_butanol = 10.0d0 / 2.0d0 * (r_butanol - q_butanol) - (r_butanol - 1.0d0)
    l_water = 10.0d0 / 2.0d0 * (r_water - q_water) - (r_water - 1.0d0)

    ! Calculate the PSI interaction parameters (eq. 9)
    PSI_mn(:,:) = exp( -a_mn(:,:) / temperature )

    ! Calculate the sum (Q_n * X_n) in the denominator of eq 9 for pure liquids
    sum_Qn_Xn_butanol = 0.0d0
    sum_Qn_Xn_water   = 0.0d0
    do n = 1, num_group
      sum_Qn_Xn_butanol = sum_Qn_Xn_butanol + Q_k(n) * real(butanol_grps(n), kind=dp) 
      sum_Qn_Xn_water   = sum_Qn_Xn_water   + Q_k(n) * real(water_grps(n),   kind=dp)
    end do

    ! Calculate group residual acitivty coefficient for pure liquids (ln(GAMMA_k^(i)) in eq 7 & 8)
    ! ... for butanol
    do m = 1, num_group
      THETA_m(m) = Q_k(m) &
                   * real(butanol_grps(m), kind=dp) &
                   / sum_Qn_Xn_butanol
    end do
    do k = 1, num_group
      sum_m_A = 0.0d0 ! ln( sum_m_A  ) term in eq 8
      sum_m_B = 0.0d0 ! last term in eq 8
      do m = 1, num_group
        sum_m_A = sum_m_A + THETA_m(m) * PSI_mn(m,k)
        sum_n = 0.0d0
        do n = 1, num_group
          sum_n = sum_n + THETA_m(n) * PSI_mn(n,m)
        end do
        sum_m_B = sum_m_B + THETA_m(m) * PSI_mn(k,m) / sum_n
      end do
      ln_GAMMA_butanol_pure(k) = Q_k(k) * (1.0d0 - log(sum_m_A) - sum_m_B)
    end do

    ! ... for water
    do m = 1, num_group
      THETA_m(m) = Q_k(m) &
                   * real(water_grps(m), kind=dp) & 
                   / sum_Qn_Xn_water
    end do
    do k = 1, num_group
      sum_m_A = 0.0d0 ! ln( sum_m_A  ) term in eq 8
      sum_m_B = 0.0d0 ! last term in eq 8
      do m = 1, num_group
        sum_m_A = sum_m_A + THETA_m(m) * PSI_mn(m,k)
        sum_n = 0.0d0
        do n = 1, num_group
          sum_n = sum_n + THETA_m(n) * PSI_mn(n,m)
        end do
        sum_m_B = sum_m_B + THETA_m(m) * PSI_mn(k,m) / sum_n
      end do
      ln_GAMMA_water_pure(k) = Q_k(k) * (1.0d0 - log(sum_m_A) - sum_m_B)
    end do

    ! Set output mole faction step (unitless)
    mass_frac_step = 1.0d0 / real(NUM_MASS_FRAC_STEP, kind=dp)

    ! Get the UNIFAC sub module json file
    input_file_path = 'test_UNIFAC_config.json'

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
    call assert(217936937, size(phlex_core%aero_rep).eq.3)
    i_aero_rep = 2

    ! Get species indices
    key = "n-butanol/water mixture.n-butanol"
    idx_butanol = phlex_core%aero_rep(i_aero_rep)%val%spec_state_id(key);
    key = "n-butanol/water mixture.water"
    idx_water = phlex_core%aero_rep(i_aero_rep)%val%spec_state_id(key);

    ! Make sure the expected species are in the model
    call assert(486977094, idx_butanol.gt.0)
    call assert(881770688, idx_water.gt.0)

    ! Integrate the mechanism
    do i_mass_frac = 0, NUM_MASS_FRAC_STEP

      ! Update the concentrations
      mass_frac = i_mass_frac * mass_frac_step
      mole_frac = (1.0d0-mass_frac)/mw_water / &
        ((1.0d0-mass_frac)/mw_water + mass_frac/mw_butanol)
      calc_conc(i_mass_frac, idx_butanol) = (1.0d0 - mole_frac) * mw_butanol
      calc_conc(i_mass_frac, idx_water)   = mole_frac * mw_water
      model_conc(i_mass_frac,:) = calc_conc(i_mass_frac,:)
      
      ! Set the concentrations in the model
      phlex_state%state_var(:) = model_conc(i_mass_frac,:)

      ! Get the modeled conc
      ! TODO finish
      model_activity(i_mass_frac,:) = 0.0d0

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!! Get the UNIFAC activities !!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      ! Equation 4
      PHI_butanol = r_butanol * (1.0d0 - mole_frac) / &
              (r_butanol * (1.0d0 - mole_frac) + r_water * mole_frac)
      PHI_water   = r_water * mole_frac / &
              (r_butanol * (1.0d0 - mole_frac) + r_water * mole_frac)
      THETA_butanol = q_butanol * (1.0d0 - mole_frac) / &
              (q_butanol * (1.0d0 - mole_frac) + q_water * mole_frac)
      THETA_water   = q_water * mole_frac / &
              (q_butanol * (1.0d0 - mole_frac) + q_water * mole_frac)
     
      ! Combinatorial portion (ln(gamma_i^C)) Eq. 3
      if (i_mass_frac.eq.NUM_MASS_FRAC_STEP) then
        ln_gamma_C_butanol = 0.0d0
      else
        ln_gamma_C_butanol = log(PHI_butanol / (1.0d0 - mole_frac)) &
                           + 5.0d0 * q_butanol * log(THETA_butanol/PHI_butanol) &
                           + l_butanol &
                           - PHI_butanol / (1.0d0 - mole_frac) &
                              * ((1.0d0 - mole_frac) * l_butanol + mole_frac * l_water)
      end if
      if (i_mass_frac.eq.0) then
        ln_gamma_C_water = 0.0d0
      else
        ln_gamma_C_water   = log(PHI_water / mole_frac) &
                           + 5.0d0 * q_water * log(THETA_water/PHI_water) &
                           + l_water &
                           - PHI_water / mole_frac &
                              * ((1.0d0 - mole_frac) * l_butanol + mole_frac * l_water) 
      end if

      ! Calculate the sum (Q_n * X_n) in the denominator of eq 9 for the mixture
      sum_Qn_Xn_mixture = 0.0d0
      do n = 1, num_group
        sum_Qn_Xn_mixture = sum_Qn_Xn_mixture &
                            + Q_k(n) &
                            * ( (1.0d0 - mole_frac) * real(butanol_grps(n), kind=dp) &
                                 + mole_frac * real(water_grps(n), kind=dp) )
      end do

      ! Group surface area fraction (Eq. 9)
      do m = 1, num_group
        THETA_m(m) = Q_k(m) &
                     * ( (1.0d0 - mole_frac) * real(butanol_grps(m), kind=dp) &
                          + mole_frac * real(water_grps(m), kind=dp) ) &
                     / sum_Qn_Xn_mixture
      end do
      
      ! Group residual activity coefficients (ln(GAMMA_k)) Eq. 8
      do k = 1, num_group
        sum_m_A = 0.0d0 ! ln( sum_m_A  ) term in eq 8
        sum_m_B = 0.0d0 ! last term in eq 8
        do m = 1, num_group
          sum_m_A = sum_m_A + THETA_m(m) * PSI_mn(m,k)
          sum_n = 0.0d0
          do n = 1, num_group
            sum_n = sum_n + THETA_m(n) * PSI_mn(n,m)
          end do
          sum_m_B = sum_m_B + THETA_m(m) * PSI_mn(k,m) / sum_n
        end do
        ln_GAMMA_mixture(k) = Q_k(k) * (1.0d0 - log(sum_m_A) - sum_m_B)
      end do

      ! Residual term (ln(gamma_i^C)) in Eq. 7
      ln_gamma_R_butanol = 0.0d0
      ln_gamma_R_water = 0.0d0
      do k = 1, num_group
        ln_gamma_R_butanol = ln_gamma_R_butanol &
                             + real(butanol_grps(k), kind=dp) &
                             * ( ln_GAMMA_mixture(k) - ln_GAMMA_butanol_pure(k) )
        ln_gamma_R_water   = ln_gamma_R_water &
                             + real(water_grps(k), kind=dp) &
                             * ( ln_GAMMA_mixture(k) - ln_GAMMA_water_pure(k) )
      end do

      ! Save activities
      calc_activity(i_mass_frac,idx_butanol) = (1.0d0 - mole_frac) &
                                             * exp( ln_gamma_C_butanol + ln_gamma_R_butanol )
      calc_activity(i_mass_frac,idx_water)   = mole_frac &
                                             * exp( ln_gamma_C_water + ln_gamma_R_water )

    end do

    ! Save the results
    open(unit=7, file="out/UNIFAC_activity_results.txt", status="replace", action="write")
    do i_mass_frac = 0, NUM_MASS_FRAC_STEP
      mass_frac = i_mass_frac * mass_frac_step
      write(7,*) mass_frac, &
            ' ', calc_conc(i_mass_frac, idx_butanol),    ' ', model_conc(i_mass_frac, idx_butanol), &
            ' ', calc_conc(i_mass_frac, idx_water),      ' ', model_conc(i_mass_frac, idx_water), &
            ' ', calc_activity(i_mass_frac, idx_butanol),' ', model_activity(i_mass_frac, idx_butanol), &
            ' ', calc_activity(i_mass_frac, idx_water),  ' ', model_activity(i_mass_frac, idx_water)
    end do
    close(7)

    ! Analyze the results
    do i_mass_frac = 1, NUM_MASS_FRAC_STEP
      mass_frac = i_mass_frac * mass_frac_step
     
      ! Check these calculations against the digitized plot for n-butanol/water @ 25C 
      ! (The digitized plot is noisy close to butanol mass fraction = 1.0)
      if (mass_frac < 0.97d0) then
        call assert_msg(666095395, &
          almost_equal(get_water_activity_fig3a(mass_frac), calc_activity(i_mass_frac, idx_water), &
          real(1.0e-1, kind=dp)), "mole_frac: "//trim(to_string(i_mass_frac))//"; species: water"// &
          "; mod: "//trim(to_string(get_water_activity_fig3a(mass_frac)))// &
          "; calc: "//trim(to_string(calc_activity(i_mass_frac, idx_water))))
      end if
    end do

    run_UNIFAC_test = .true.

  end function run_UNIFAC_test

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get the water activity for the n-butanol/water mixture of Marcolli & Peters 2005
  !!
  !! Values for water activity were digitized from Fig. 3(a) solid red lines
  !! for UNIFAC (25 C) calculations
  function get_water_activity_fig3a(mass_frac_butanol) result (water_activity)

    !> Water activity (unitless)
    real(kind=dp) :: water_activity
    !> Mass fraction of butanol (unitless)
    real(kind=dp), intent(in) :: mass_frac_butanol
    !> Number of digitized points
    integer(kind=i_kind), parameter :: NUM_POINTS = 101
    !> Digitized water activities (unitless)
    real(kind=dp), parameter, dimension(NUM_POINTS) :: a_w = [ &
      0.995733, 0.994375, 0.993023, 0.991684, 0.990363, 0.989068, 0.987805, &
      0.986579, 0.985398, 0.984267, 0.983193, 0.982182, 0.981240, 0.980375, &
      0.979592, 0.978898, 0.978295, 0.977784, 0.977369, 0.977051, 0.976832, &
      0.976716, 0.976703, 0.976796, 0.976993, 0.977292, 0.977693, 0.978192, &
      0.978790, 0.979481, 0.980265, 0.981136, 0.982093, 0.983131, 0.984248, &
      0.985440, 0.986704, 0.988037, 0.989436, 0.990897, 0.992414, 0.993984, &
      0.995600, 0.997259, 0.998954, 1.000680, 1.002430, 1.004210, 1.005990, &
      1.007770, 1.009550, 1.011300, 1.013030, 1.014720, 1.016360, 1.017930, &
      1.019420, 1.020810, 1.022090, 1.023240, 1.024240, 1.025090, 1.025750, &
      1.026190, 1.026370, 1.026280, 1.025860, 1.025080, 1.023880, 1.022180, &
      1.019920, 1.017000, 1.013340, 1.008840, 1.003390, 0.996856, 0.989114, &
      0.980037, 0.969500, 0.957374, 0.943524, 0.927812, 0.910098, 0.890250, &
      0.868153, 0.843702, 0.816797, 0.787379, 0.755413, 0.720865, 0.683730, &
      0.644058, 0.601906, 0.557361, 0.510599, 0.461815, 0.411277, 0.359302, &
      0.306229, 0.252423, 0.198251 ]

    integer(kind=i_kind) :: i_a_w

    i_a_w = int(mass_frac_butanol*100.0d0, kind=i_kind) + 1

    water_activity = a_w(i_a_w)

  end function get_water_activity_fig3a

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program pmc_test_sub_module_UNIFAC
