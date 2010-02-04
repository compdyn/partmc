! Copyright (C) 2005-2010 Nicole Riemer and Matthew West
! Copyright (C) 2009 Joseph Ching
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_condense module.

!> Water condense onto aerosol particles.
module pmc_condense

  use pmc_aero_state
  use pmc_bin_grid
  use pmc_env_data
  use pmc_env_state
  use pmc_aero_data
  use pmc_util
  use pmc_aero_particle
  use pmc_constants
  !>DEBUG
  !use dvode_f90_m
  !<DEBUG
  use iso_c_binding

  !> Whether to numerically test the Jacobian-times-vector function
  !> during execution (for debugging only).
  logical, parameter :: CONDENSE_DO_TEST_JTIMES = .false.

  !> Internal-use structure for storing the inputs for the
  !> rate-calculation function.
  type condense_rates_inputs_t
     !> Temperature (K).
     real(kind=dp) :: T
     !> Rate of change of temperature (K s^{-1}).
     real(kind=dp) :: Tdot
     !> Relative humidity (1).
     real(kind=dp) :: H
     !> Pressure (Pa).
     real(kind=dp) :: p
     !> Computational volume (m^3).
     real(kind=dp) :: V_comp
     !> Particle diameter (m).
     real(kind=dp) :: D
     !> Particle dry diameter (m).
     real(kind=dp) :: D_dry
     !> Kappa parameter (1).
     real(kind=dp) :: kappa
  end type condense_rates_inputs_t

  !> Internal-use structure for storing the outputs from the
  !> rate-calculation function.
  type condense_rates_outputs_t
     !> Change rate of diameter (m s^{-1}).
     real(kind=dp) :: Ddot
     !> Change rate of relative humidity due to this particle (s^{-1}).
     real(kind=dp) :: Hdot_i
     !> Change rate of relative humidity due to environment changes (s^{-1}).
     real(kind=dp) :: Hdot_env
     !> Sensitivity of \c Ddot to input \c D (m s^{-1} m^{-1}).
     real(kind=dp) :: dDdot_dD
     !> Sensitivity of \c Ddot to input \c H (m s^{-1}).
     real(kind=dp) :: dDdot_dH
     !> Sensitivity of \c Hdot_i to input \c D (s^{-1} m^{-1}).
     real(kind=dp) :: dHdoti_dD
     !> Sensitivity of \c Hdot_i to input \c D (s^{-1}).
     real(kind=dp) :: dHdoti_dH
     !> Sensitivity of \c Hdot_env to input \c D (s^{-1} m^{-1}).
     real(kind=dp) :: dHdotenv_dD
     !> Sensitivity of \c Hdot_env to input \c D (s^{-1}).
     real(kind=dp) :: dHdotenv_dH
  end type condense_rates_outputs_t

  !> Internal-use variable for storing the aerosol data during calls
  !> to the ODE solver.
  type(aero_data_t) :: condense_saved_aero_data
  !> Internal-use variable for storing the environment data during
  !> calls to the ODE solver.
  type(env_data_t) :: condense_saved_env_data
  !> Internal-use variable for storing the initial environment state
  !> during calls to the ODE solver.
  type(env_state_t) :: condense_saved_env_state_initial
  !> Internal-use variable for storing the inital computational volume
  !> during calls to the ODE solver.
  real(kind=dp) :: condense_saved_V_comp_initial
  !> Internal-use variable for storing the rate of change of the
  !> temperature during calls to the ODE solver.
  real(kind=dp) :: condense_saved_Tdot
  !> Internal-use variable for storing the per-particle kappa values
  !> during calls to the ODE solver.
  real(kind=dp), allocatable :: condense_saved_kappa(:)
  !> Internal-use variable for storing the per-particle dry diameters
  !> during calls to the ODE solver.
  real(kind=dp), allocatable :: condense_saved_D_dry(:)
  !>DEBUG
  integer, save :: condense_n_vf
  integer, save :: condense_n_jtimes
  integer, save :: condense_n_prec
  !<DEBUG

contains
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Do condensation to all the particles for a given time interval,
  !> including updating the environment to account for the lost
  !> vapor.
  subroutine condense_particles(bin_grid, env_state, env_data, &
       aero_data, aero_state, del_t)

    !> Bin grid.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Environment state.
    type(env_state_t), intent(inout) :: env_state
    !> Environment data.
    type(env_data_t), intent(in) :: env_data
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Aerosol state.
    type(aero_state_t), intent(inout) :: aero_state
    !> Total time to integrate.
    real(kind=dp), intent(in) :: del_t
    
    integer :: i_bin, j, new_bin, k, i_part, n_eqn, i_eqn, solver_stat
    type(aero_particle_t), pointer :: aero_particle
    real(kind=dp) :: state(aero_state%n_part + 1), init_time, final_time
    real(kind=dp) :: abs_tol_vector(aero_state%n_part + 1)
    !type(vode_opts) :: options
    type(env_state_t) :: env_state_final
    real(kind=dp) :: water_vol_initial, water_vol_final, d_water_vol
    real(kind=dp) :: vapor_vol_initial, vapor_vol_final, d_vapor_vol
    real(kind=dp) :: V_comp_final, water_rel_error
    real(kind=c_double), target :: state_f(aero_state%n_part + 1)
    real(kind=c_double), target :: abstol_f(aero_state%n_part + 1)
    type(c_ptr) :: state_f_p, abstol_f_p
    !type(c_funptr) :: condense_vf_f_p, condense_jtimes_f_p
    integer(kind=c_int) :: n_eqn_f
    real(kind=c_double) :: reltol_f, t_initial_f, t_final_f

    interface
       integer(kind=c_int) function condense_solver(neq, x_f, abstol_f, reltol_f, &
            t_initial_f, t_final_f) bind(c)
         use iso_c_binding
         integer(kind=c_int), value :: neq
         type(c_ptr), value :: x_f
         type(c_ptr), value :: abstol_f
         real(kind=c_double), value :: reltol_f
         real(kind=c_double), value :: t_initial_f
         real(kind=c_double), value :: t_final_f
         !type(c_funptr) :: condense_vf_f_p
         !type(c_funptr) :: condense_jtimes_f_p
       end function condense_solver
    end interface

    ! initial water volume in the aerosol particles
    water_vol_initial = 0d0
    do i_bin = 1,bin_grid%n_bin
       do j = 1,aero_state%bin(i_bin)%n_part
          aero_particle => aero_state%bin(i_bin)%particle(j)
          water_vol_initial = water_vol_initial &
               + aero_particle%vol(aero_data%i_water)
       end do
    end do

    ! save data for use within the timestepper
    call aero_data_allocate(condense_saved_aero_data)
    call env_data_allocate(condense_saved_env_data)
    call env_state_allocate(condense_saved_env_state_initial)

    call aero_data_copy(aero_data, condense_saved_aero_data)
    call env_data_copy(env_data, condense_saved_env_data)
    call env_state_copy(env_state, condense_saved_env_state_initial)

    condense_saved_V_comp_initial = aero_state%comp_vol
    
    call env_state_allocate(env_state_final)
    call env_state_copy(env_state, env_state_final)
    call env_data_update_state(env_data, env_state_final, &
         env_state_final%elapsed_time + del_t, update_rel_humid = .false.)
    condense_saved_Tdot = (env_state_final%temp - env_state%temp) / del_t

    ! construct initial state vector from aero_state and env_state
    allocate(condense_saved_kappa(aero_state%n_part))
    allocate(condense_saved_D_dry(aero_state%n_part))
    i_part = 0
    do i_bin = 1,bin_grid%n_bin
       do j = 1,aero_state%bin(i_bin)%n_part
          i_part = i_part + 1
          aero_particle => aero_state%bin(i_bin)%particle(j)
          condense_saved_kappa(i_part) &
               = aero_particle_solute_kappa(aero_particle, aero_data)
          condense_saved_D_dry(i_part) = vol2diam(&
               aero_particle_solute_volume(aero_particle, aero_data))
          state(i_part) = aero_particle_diameter(aero_particle)
          abs_tol_vector(i_part) = max(1d-30, &
                1d-8 * (state(i_part) - condense_saved_D_dry(i_part)))
       end do
    end do
    state(aero_state%n_part + 1) = env_state%rel_humid
    abs_tol_vector(aero_state%n_part + 1) = 1d-10

    ! call VODE solver
    !n_eqn = aero_state%n_part + 1
    !init_time = 0d0
    !final_time = del_t
    !itask = 1 ! just output val at final_time
    !istate = 1 ! first call for this ODE
    !options = set_opts(sparse_j = .true., &
    !     user_supplied_jacobian = .true., &
    !     abserr_vector = abs_tol_vector, relerr = 1d-8)
    !call dvode_f90(condense_vode_f, n_eqn, state, init_time, final_time, &
    !     itask, istate, options, j_fcn = condense_vode_jac)

    ! call SUNDIALS solver
    n_eqn = aero_state%n_part + 1
    n_eqn_f = int(n_eqn, kind=c_int)
    reltol_f = real(1d-8, kind=c_double)
    t_initial_f = real(0, kind=c_double)
    t_final_f = real(del_t, kind=c_double)
    do i_eqn = 1,n_eqn
       state_f(i_eqn) = real(state(i_eqn), kind=c_double)
       abstol_f(i_eqn) = real(abs_tol_vector(i_eqn), kind=c_double)
    end do
    state_f_p = c_loc(state_f)
    abstol_f_p = c_loc(abstol_f)
    !condense_vf_f_p = c_funloc(condense_vf_f)
    !condense_jtimes_f_p = c_funloc(condense_jtimes_f)
    if (CONDENSE_DO_TEST_JTIMES) then
       call condense_test_jtimes(n_eqn, 0d0, state)
    end if
    !>DEBUG
    condense_n_vf = 0
    condense_n_jtimes = 0
    condense_n_prec = 0
    !<DEBUG
    solver_stat = condense_solver(n_eqn_f, state_f_p, abstol_f_p, reltol_f, &
         t_initial_f, t_final_f) !, condense_vf_f_p, condense_jtimes_f_p)
    !>DEBUG
    write(*,*) 'condense_n_vf ', condense_n_vf
    write(*,*) 'condense_n_jtimes ', condense_n_jtimes
    write(*,*) 'condense_n_prec ', condense_n_prec
    !<DEBUG
    !>DEBUG
    !stop
    !<DEBUG
    do i_eqn = 1,n_eqn
       state(i_eqn) = real(state_f(i_eqn), kind=dp)
    end do

    ! unpack result state vector into aero_state and env_state
    i_part = 0
    do i_bin = 1,bin_grid%n_bin
       do j = 1,aero_state%bin(i_bin)%n_part
          i_part = i_part + 1
          aero_particle => aero_state%bin(i_bin)%particle(j)

          ! translate output back to particle
          aero_particle%vol(aero_data%i_water) = diam2vol(state(i_part)) &
               - aero_particle_solute_volume(aero_particle, aero_data)

          ! ensure volumes stay positive
          aero_particle%vol(aero_data%i_water) = max(0d0, &
               aero_particle%vol(aero_data%i_water))
       end do
    end do
    ! We've modified particle diameters, so we need to update which
    ! bins they are in.
    call aero_state_resort(bin_grid, aero_state)
    env_state%rel_humid = state(aero_state%n_part + 1)

    ! final water volume in the aerosol particles and in the vapor
    water_vol_final = 0d0
    do i_bin = 1,bin_grid%n_bin
       do j = 1,aero_state%bin(i_bin)%n_part
          aero_particle => aero_state%bin(i_bin)%particle(j)
          water_vol_final = water_vol_final &
               + aero_particle%vol(aero_data%i_water)
       end do
    end do

    ! check that water removed from particles equals water added to vapor
    V_comp_final = condense_saved_V_comp_initial &
         * env_state_final%temp / condense_saved_env_state_initial%temp
    vapor_vol_initial = aero_data%molec_weight(aero_data%i_water) &
         / (const%univ_gas_const * condense_saved_env_state_initial%temp) &
         * env_state_sat_vapor_pressure(condense_saved_env_state_initial) &
         * condense_saved_env_state_initial%rel_humid &
         * condense_saved_V_comp_initial &
         / aero_particle_water_density(aero_data)
    vapor_vol_final = aero_data%molec_weight(aero_data%i_water) &
         / (const%univ_gas_const * env_state_final%temp) &
         * env_state_sat_vapor_pressure(env_state_final) &
         * env_state%rel_humid &
         * V_comp_final / aero_particle_water_density(aero_data)
    d_vapor_vol = vapor_vol_final - vapor_vol_initial
    d_water_vol = water_vol_final - water_vol_initial
    water_rel_error = (d_vapor_vol + d_water_vol) &
         / (vapor_vol_final + water_vol_final)
    call warn_assert_msg(477865387, abs(water_rel_error) < 1d-6, &
         "condensation water imbalance too high: " &
         // real_to_string(water_rel_error))

    deallocate(condense_saved_kappa)
    deallocate(condense_saved_D_dry)
    call env_state_deallocate(env_state_final)
    call aero_data_deallocate(condense_saved_aero_data)
    call env_data_deallocate(condense_saved_env_data)
    call env_state_deallocate(condense_saved_env_state_initial)

  end subroutine condense_particles

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Fills in the \c env_state with the current environment state,
  !> taken from the \c state vector and from global variables.
  subroutine condense_current_env_state(n_eqn, time, state, env_state)

    !> Length of state vector.
    integer, intent(in) :: n_eqn
    !> Current time (s).
    real(kind=dp), intent(in) :: time
    !> Current state vector.
    real(kind=dp), intent(in) :: state(n_eqn)
    !> Current environment state.
    type(env_state_t), intent(inout) :: env_state

    call env_state_copy(condense_saved_env_state_initial, env_state)
    call env_data_update_state(condense_saved_env_data, &
         env_state, env_state%elapsed_time + time, &
         update_rel_humid = .false.)
    env_state%rel_humid = state(n_eqn)

  end subroutine condense_current_env_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Compute the rate of change of particle diameter and relative
  !> humidity for a single particle.
  subroutine condense_rates(inputs, outputs)

    !> Inputs to rates.
    type(condense_rates_inputs_t), intent(in) :: inputs
    !> Outputs rates.
    type(condense_rates_outputs_t), intent(out) :: outputs

    real(kind=dp) :: rho_w, M_w, P_0, dP0_dT, rho_air, k_a, D_v, U, V
    real(kind=dp) :: dV_dT, W, X, Y, Z, k_ap, dkap_dD, D_vp, dDvp_dD
    real(kind=dp) :: a_w, daw_dD, delta_star, h, dh_ddelta, dh_dD
    real(kind=dp) :: dh_dH, ddeltastar_dD, ddeltastar_dH, dVcomp_dT
    integer :: newton_step

    rho_w = const%water_density
    M_w = const%water_molec_weight
    P_0 = const%water_eq_vap_press &
         * 10d0**(7.45d0 * (inputs%T - const%water_freeze_temp) &
         / (inputs%T - 38d0))
    dP0_dT = P_0 * 7.45d0 * log(10d0) * (const%water_freeze_temp - 38d0) &
         / (inputs%T - 38d0)**2
    rho_air = const%air_molec_weight * inputs%p &
         / (const%univ_gas_const * inputs%T)
    dVcomp_dT = inputs%V_comp / inputs%T

    k_a = 1d-3 * (4.39d0 + 0.071d0 * inputs%T)
    D_v = 0.211d-4 / (inputs%p / const%air_std_press) &
         * (inputs%T / 273d0)**1.94d0
    U = const%water_latent_heat * rho_w / (4d0 * inputs%T)
    V = 4d0 * M_w * P_0 / (rho_w * const%univ_gas_const * inputs%T)
    dV_dT = (dP0_dT / P_0 - 1d0 / inputs%T) * V
    W = const%water_latent_heat * M_w / (const%univ_gas_const * inputs%T)
    X = 4d0 * M_w * const%water_surf_eng &
         / (const%univ_gas_const * inputs%T * rho_w) 
    Y = 2d0 * k_a / (const%accom_coeff * rho_air &
         * const%air_spec_heat) &
         * sqrt(2d0 * const%pi * const%air_molec_weight &
         / (const%univ_gas_const * inputs%T))
    Z = 2d0 * D_v / const%accom_coeff * sqrt(2d0 * const%pi * M_w &
         / (const%univ_gas_const * inputs%T))

    outputs%Hdot_env = - dV_dT * inputs%Tdot * inputs%H / V &
         - dVcomp_dT * inputs%Tdot * inputs%H / inputs%V_comp
    outputs%dHdotenv_dD = 0d0
    outputs%dHdotenv_dH = - dV_dT * inputs%Tdot / V &
         - dVcomp_dT * inputs%Tdot / inputs%V_comp

    if (inputs%D <= inputs%D_dry) then
       k_ap = k_a / (1d0 + Y / inputs%D_dry)
       dkap_dD = 0d0
       D_vp = D_v / (1d0 + Z / inputs%D_dry)
       dDvp_dD = 0d0
       a_w = 0d0
       daw_dD = 0d0

       delta_star = U * V * D_vp * inputs%H / k_ap
       
       outputs%Ddot = k_ap * delta_star / (U * inputs%D_dry)
       outputs%Hdot_i = - 2d0 * const%pi / (V * inputs%V_comp) &
            * inputs%D_dry**2 * outputs%Ddot
       
       dh_ddelta = k_ap
       dh_dD = 0d0
       dh_dH = - U * V * D_vp

       ddeltastar_dD = - dh_dD / dh_ddelta
       ddeltastar_dH = - dh_dH / dh_ddelta
       
       outputs%dDdot_dD = 0d0
       outputs%dDdot_dH = k_ap / (U * inputs%D_dry) * ddeltastar_dH
       outputs%dHdoti_dD = - 2d0 * const%pi / (V * inputs%V_comp) &
            * inputs%D_dry**2 * outputs%dDdot_dD
       outputs%dHdoti_dH = - 2d0 * const%pi / (V * inputs%V_comp) &
            * inputs%D_dry**2 * outputs%dDdot_dH

       return
    end if

    k_ap = k_a / (1d0 + Y / inputs%D)
    dkap_dD = k_a * Y / (inputs%D + Y)**2
    D_vp = D_v / (1d0 + Z / inputs%D)
    dDvp_dD = D_v * Z / (inputs%D + Z)**2
    a_w = (inputs%D**3 - inputs%D_dry**3) &
         / (inputs%D**3 + (inputs%kappa - 1d0) * inputs%D_dry**3)
    daw_dD = 3d0 * inputs%D**2 * inputs%kappa * inputs%D_dry**3 &
         / (inputs%D**3 + (inputs%kappa - 1d0) * inputs%D_dry**3)**2

    delta_star = 0d0
    h = 0d0
    dh_ddelta = 1d0
    do newton_step = 1,5
       ! update delta_star first so when the newton loop ends we have
       ! h and dh_ddelta evaluated at the final delta_star value
       delta_star = delta_star - h / dh_ddelta
       h = k_ap * delta_star - U * V * D_vp &
            * (inputs%H - a_w / (1d0 + delta_star) &
            * exp(W * delta_star / (1d0 + delta_star) &
            + (X / inputs%D) / (1d0 + delta_star)))
       dh_ddelta = &
            k_ap - U * V * D_vp * a_w / (1d0 + delta_star)**2 &
            * (1d0 - W / (1d0 + delta_star) &
            + (X / inputs%D) / (1d0 + delta_star)) &
            * exp(W * delta_star / (1d0 + delta_star) &
            + (X / inputs%D) / (1d0 + delta_star))
    end do
    call warn_assert_msg(387362320, &
         abs(h) < 1d3 * epsilon(1d0) * abs(U * V * D_vp * inputs%H), &
         "condensation newton loop did not satisfy convergence tolerance")

    outputs%Ddot = k_ap * delta_star / (U * inputs%D)
    outputs%Hdot_i = - 2d0 * const%pi / (V * inputs%V_comp) &
         * inputs%D**2 * outputs%Ddot

    dh_dD = dkap_dD * delta_star &
         - U * V * dDvp_dD * inputs%H + U * V &
         * (a_w * dDvp_dD + D_vp * daw_dD &
         - D_vp * a_w * (X / inputs%D**2) / (1d0 + delta_star)) &
         * (1d0 / (1d0 + delta_star)) &
         * exp((W * delta_star) / (1d0 + delta_star) &
         + (X / inputs%D) / (1d0 + delta_star))
    dh_dH = - U * V * D_vp

    ddeltastar_dD = - dh_dD / dh_ddelta
    ddeltastar_dH = - dh_dH / dh_ddelta

    outputs%dDdot_dD = dkap_dD * delta_star / (U * inputs%D) &
         + k_ap * ddeltastar_dD / (U * inputs%D) &
         - k_ap * delta_star / (U * inputs%D**2)
    outputs%dDdot_dH = k_ap / (U * inputs%D) * ddeltastar_dH
    outputs%dHdoti_dD = - 2d0 * const%pi / (V * inputs%V_comp) &
         * (2d0 * inputs%D * outputs%Ddot + inputs%D**2 * outputs%dDdot_dD)
    outputs%dHdoti_dH = - 2d0 * const%pi / (V * inputs%V_comp) &
         * inputs%D**2 * outputs%dDdot_dH

  end subroutine condense_rates

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine condense_vode_f(n_eqn, time, state, state_dot)

    !> Length of state vector.
    integer, intent(in) :: n_eqn
    !> Current time (s).
    real(kind=dp), intent(in) :: time
    !> Current state vector.
    real(kind=dp), intent(in) :: state(n_eqn)
    !> Time derivative of state vector.
    real(kind=dp), intent(out) :: state_dot(n_eqn)

    real(kind=dp) :: Hdot
    integer :: i_part
    type(env_state_t) :: env_state
    type(condense_rates_inputs_t) :: inputs
    type(condense_rates_outputs_t) :: outputs

    call env_state_allocate(env_state)
    call condense_current_env_state(n_eqn, time, state, env_state)

    inputs%T = env_state%temp
    inputs%Tdot = condense_saved_Tdot
    inputs%H = env_state%rel_humid
    inputs%p = env_state%pressure
    inputs%V_comp = condense_saved_V_comp_initial &
         * env_state%temp / condense_saved_env_state_initial%temp

    Hdot = 0d0
    do i_part = 1,(n_eqn - 1)
       inputs%D = state(i_part)
       inputs%D_dry = condense_saved_D_dry(i_part)
       inputs%kappa = condense_saved_kappa(i_part)
       call condense_rates(inputs, outputs)
       state_dot(i_part) = outputs%Ddot
       Hdot = Hdot + outputs%Hdot_i
    end do
    Hdot = Hdot + outputs%Hdot_env

    state_dot(n_eqn) = Hdot

    call env_state_deallocate(env_state)

  end subroutine condense_vode_f

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine condense_vode_jac(n_eqn, time, state, ia, ja, &
       n_non_zero, state_jac)

    !> Length of state vector.
    integer, intent(in) :: n_eqn
    !> Current time (s).
    real(kind=dp), intent(in) :: time
    !> Current state vector.
    real(kind=dp), intent(in) :: state(n_eqn)
    !> Non-zero column offsets.
    integer, intent(out) :: ia(*)
    !> Non-zero row locations.
    integer, intent(out) :: ja(*)
    !> Number of non-zero elements in the Jacobian.
    integer, intent(inout) :: n_non_zero
    !> Sparse Jacobian of time derivative of state vector.
    real(kind=dp), intent(out) :: state_jac(*)

    real(kind=dp) :: dDdot_dD(n_eqn - 1), dDdot_dH(n_eqn - 1)
    real(kind=dp) :: dHdot_dD(n_eqn - 1), dHdot_dH
    integer :: i_nz, i_part
    type(env_state_t) :: env_state
    type(condense_rates_inputs_t) :: inputs
    type(condense_rates_outputs_t) :: outputs

    ! signal from vode to initialize number of non-zeros
    if (n_non_zero == 0) then
       n_non_zero = 3 * n_eqn - 2
       return
    end if

    ! if not initializing, this should be correct
    call assert(395158320, n_non_zero == 3 * n_eqn - 2)

    call env_state_allocate(env_state)
    call condense_current_env_state(n_eqn, time, state, env_state)

    inputs%T = env_state%temp
    inputs%Tdot = condense_saved_Tdot
    inputs%H = env_state%rel_humid
    inputs%p = env_state%pressure
    inputs%V_comp = condense_saved_V_comp_initial &
         * env_state%temp / condense_saved_env_state_initial%temp

    dHdot_dH = 0d0
    do i_part = 1,(n_eqn - 1)
       inputs%D = state(i_part)
       inputs%D_dry = condense_saved_D_dry(i_part)
       inputs%kappa = condense_saved_kappa(i_part)
       call condense_rates(inputs, outputs)
       dDdot_dD(i_part) = outputs%dDdot_dD
       dDdot_dH(i_part) = outputs%dDdot_dH
       dHdot_dD(i_part) = outputs%dHdoti_dD + outputs%dHdotenv_dD
       dHdot_dH = dHdot_dH + outputs%dHdoti_dH
    end do
    dHdot_dH = dHdot_dH + outputs%dHdotenv_dH

    ! Copied from dvode_f90_m documentation:
    !
    ! IA defines the number of nonzeros including the
    ! diagonal in each column of the Jacobian. Define
    ! IA(1) = 1 and for J = 1,..., N,
    ! IA(J+1) = IA(J) + number of nonzeros in column J.
    ! Diagonal elements must be included even if they are
    ! zero. You should check to ensure that IA(N+1)-1 = NZ.
    ! JA defines the rows in which the nonzeros occur.
    ! For I = 1,...,NZ, JA(I) is the row in which the Ith
    ! element of the Jacobian occurs. JA must also include
    ! the diagonal elements of the Jacobian.
    ! PD defines the numerical value of the Jacobian
    ! elements. For I = 1,...,NZ, PD(I) is the numerical
    ! value of the Ith element in the Jacobian. PD must
    ! also include the diagonal elements of the Jacobian.

    ia(1) = 1
    do i_part = 1,(n_eqn - 1)
       ia(1 + i_part) = ia(i_part) + 2
    end do
    ia(n_eqn + 1) = ia(n_eqn) + n_eqn
    call assert(824345254, ia(n_eqn + 1) - 1 == n_non_zero)
    
    i_nz = 0
    do i_part = 1,(n_eqn - 1)
       i_nz = i_nz + 1
       ja(i_nz) = i_part
       state_jac(i_nz) = dDdot_dD(i_part)
       i_nz = i_nz + 1
       ja(i_nz) = n_eqn
       state_jac(i_nz) = dHdot_dD(i_part)
    end do
    do i_part = 1,(n_eqn - 1)
       i_nz = i_nz + 1
       ja(i_nz) = i_part
       state_jac(i_nz) = dDdot_dH(i_part)
    end do
    i_nz = i_nz + 1
    ja(i_nz) = n_eqn
    state_jac(i_nz) = dHdot_dH
    call assert(901219708, i_nz == n_non_zero)

    call env_state_deallocate(env_state)

  end subroutine condense_vode_jac

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Determine the equilibrium state of a single particle.
  subroutine condense_equilib_particle(env_state, aero_data, &
       aero_particle)

    !> Environment state.
    type(env_state_t), intent(in) :: env_state
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Particle.
    type(aero_particle_t), intent(inout) :: aero_particle

    real(kind=dp) :: X, kappa, D_dry, D, g, dg_dD, a_w, daw_dD
    integer :: newton_step

    X = 4d0 * const%water_molec_weight * const%water_surf_eng &
         / (const%univ_gas_const * env_state%temp &
         * const%water_density)
    kappa = aero_particle_solute_kappa(aero_particle, aero_data)
    D_dry = vol2diam(aero_particle_solute_volume(aero_particle, aero_data))

    D = D_dry
    g = 0d0
    dg_dD = 1d0
    do newton_step = 1,20
       D = D - g / dg_dD
       a_w = (D**3 - D_dry**3) / (D**3 + (kappa - 1d0) * D_dry**3)
       daw_dD = 3d0 * D**2 * kappa * D_dry**3 &
            / (D**3 + (kappa - 1d0) * D_dry**3)**2
       g = env_state%rel_humid - a_w * exp(X / D)
       dg_dD = - daw_dD * exp(X / D) + a_w * exp(X / D) * (X / D**2)
    end do
    call warn_assert_msg(426620001, abs(g) < 1d3 * epsilon(1d0), &
         "convergence problem in equilibriation")

    aero_particle%vol(aero_data%i_water) = diam2vol(D) - diam2vol(D_dry)

  end subroutine condense_equilib_particle

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Call condense_equilib_particle() on each particle in the aerosol.
  subroutine condense_equilib_particles(bin_grid, env_state, aero_data, &
       aero_state)

    !> Bin grid.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Environment state.
    type(env_state_t), intent(inout) :: env_state
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Aerosol state.
    type(aero_state_t), intent(inout) :: aero_state
 
    integer :: i_bin, i
    
    do i_bin = 1,bin_grid%n_bin
       do i = 1,aero_state%bin(i_bin)%n_part
          call condense_equilib_particle(env_state, aero_data, &
               aero_state%bin(i_bin)%particle(i))
       end do
    end do

  end subroutine condense_equilib_particles

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine condense_vf_f(n_eqn, time, state_p, state_dot_p) bind(c)
    
    !> Length of state vector.
    integer(kind=c_int), value, intent(in) :: n_eqn
    !> Current time (s).
    real(kind=c_double), value, intent(in) :: time
    !> Pointer to state data.
    type(c_ptr), value, intent(in) :: state_p
    !> Pointer to state_dot data.
    type(c_ptr), value, intent(in) :: state_dot_p

    real(kind=c_double), pointer :: state(:)
    real(kind=c_double), pointer :: state_dot(:)
    real(kind=dp) :: Hdot
    integer :: i_part
    type(env_state_t) :: env_state
    type(condense_rates_inputs_t) :: inputs
    type(condense_rates_outputs_t) :: outputs

    !>DEBUG
    condense_n_vf = condense_n_vf + 1
    !<DEBUG
    call c_f_pointer(state_p, state, (/ n_eqn /))
    call c_f_pointer(state_dot_p, state_dot, (/ n_eqn /))
    
    call env_state_allocate(env_state)
    call condense_current_env_state(n_eqn, time, state, env_state)

    inputs%T = env_state%temp
    inputs%Tdot = condense_saved_Tdot
    inputs%H = env_state%rel_humid
    inputs%p = env_state%pressure
    inputs%V_comp = condense_saved_V_comp_initial &
         * env_state%temp / condense_saved_env_state_initial%temp
    
    Hdot = 0d0
    do i_part = 1,(n_eqn - 1)
       inputs%D = state(i_part)
       inputs%D_dry = condense_saved_D_dry(i_part)
       inputs%kappa = condense_saved_kappa(i_part)
       call condense_rates(inputs, outputs)
       state_dot(i_part) = outputs%Ddot
       Hdot = Hdot + outputs%Hdot_i
    end do
    Hdot = Hdot + outputs%Hdot_env
    
    state_dot(n_eqn) = Hdot
    
    call env_state_deallocate(env_state)
    
    !>DEBUG
    !write(*,*) 'condense_jtimes_f: time = ', time
    !write(*,*) 'condense_vf_f: state(1), state(n_eqn) = ', state(1), state(n_eqn)
    !write(*,*) 'condense_vf_f: state_dot(1), state_dot(n_eqn) = ', state_dot(1), state_dot(n_eqn)
    !<DEBUG
    
  end subroutine condense_vf_f
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine condense_jac(n_eqn, time, state_p, dDdot_dD, dDdot_dH, &
       dHdot_dD, dHdot_dH)

    !> Length of state vector.
    integer(kind=c_int), intent(in) :: n_eqn
    !> Current time (s).
    real(kind=c_double), intent(in) :: time
    !> Pointer to current state vector.
    type(c_ptr), intent(in) :: state_p
    !> Derivative of Ddot with respect to D.
    real(kind=dp), intent(out) :: dDdot_dD(n_eqn - 1)
    !> Derivative of Ddot with respect to H.
    real(kind=dp), intent(out) :: dDdot_dH(n_eqn - 1)
    !> Derivative of Hdot with respect to D.
    real(kind=dp), intent(out) :: dHdot_dD(n_eqn - 1)
    !> Derivative of Hdot with respect to H.
    real(kind=dp), intent(out) :: dHdot_dH

    real(kind=c_double), pointer :: state(:)
    integer :: i_part
    type(env_state_t) :: env_state
    type(condense_rates_inputs_t) :: inputs
    type(condense_rates_outputs_t) :: outputs

    call c_f_pointer(state_p, state, (/ n_eqn /))

    call env_state_allocate(env_state)
    call condense_current_env_state(n_eqn, time, state, env_state)
    
    inputs%T = env_state%temp
    inputs%Tdot = condense_saved_Tdot
    inputs%H = env_state%rel_humid
    inputs%p = env_state%pressure
    inputs%V_comp = condense_saved_V_comp_initial &
         * env_state%temp / condense_saved_env_state_initial%temp
    
    dHdot_dH = 0d0
    do i_part = 1,(n_eqn - 1)
       inputs%D = state(i_part)
       inputs%D_dry = condense_saved_D_dry(i_part)
       inputs%kappa = condense_saved_kappa(i_part)
       call condense_rates(inputs, outputs)
       dDdot_dD(i_part) = outputs%dDdot_dD
       dDdot_dH(i_part) = outputs%dDdot_dH
       dHdot_dD(i_part) = outputs%dHdoti_dD + outputs%dHdotenv_dD
       dHdot_dH = dHdot_dH + outputs%dHdoti_dH
    end do
    dHdot_dH = dHdot_dH + outputs%dHdotenv_dH
    
    call env_state_deallocate(env_state)
    
  end subroutine condense_jac

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine condense_jtimes_f(n_eqn, time, state_p, vec_p, &
       jac_times_vec_p) bind(c)
    
    !> Length of state vector.
    integer(kind=c_int), value, intent(in) :: n_eqn
    !> Current time (s).
    real(kind=c_double), value, intent(in) :: time
    !> Pointer to current state vector.
    type(c_ptr), value, intent(in) :: state_p
    !> Pointer to vector to multiply by the Jacobian.
    type(c_ptr), value, intent(in) :: vec_p
    !> Pointer to jacobian multiplied by the input vector \c vec.
    type(c_ptr), value, intent(in) :: jac_times_vec_p

    real(kind=c_double), pointer :: vec(:)
    real(kind=c_double), pointer :: jac_times_vec(:)
    real(kind=dp) :: dDdot_dD(n_eqn - 1), dDdot_dH(n_eqn - 1)
    real(kind=dp) :: dHdot_dD(n_eqn - 1), dHdot_dH
    integer :: i_part

    !>DEBUG
    condense_n_jtimes = condense_n_jtimes + 1
    !<DEBUG
    call condense_jac(n_eqn, time, state_p, dDdot_dD, dDdot_dH, &
         dHdot_dD, dHdot_dH)

    call c_f_pointer(vec_p, vec, (/ n_eqn /))
    call c_f_pointer(jac_times_vec_p, jac_times_vec, (/ n_eqn /))
    
    jac_times_vec(n_eqn) = dHdot_dH * vec(n_eqn)
    do i_part = 1,(n_eqn - 1)
       jac_times_vec(i_part) = dDdot_dD(i_part) * vec(i_part) &
            + dDdot_dH(i_part) * vec(n_eqn)
       jac_times_vec(n_eqn) = jac_times_vec(n_eqn) &
            + dHdot_dD(i_part) * vec(i_part)
    end do

  end subroutine condense_jtimes_f

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Precondition by solving the system \f$ Pz = r \f$ where \f$ P
  !> \approx M = I - \gamma J \f$ and \f$ J = \partial f / \partial y
  !> \f$.
  subroutine condense_prec_diag_f(n_eqn, time, state_p, state_dot_p, rhs_p, &
       soln_p, gamma, delta, left_or_right) bind(c)

    !> Length of state vector.
    integer(kind=c_int), value, intent(in) :: n_eqn
    !> Current time (s).
    real(kind=c_double), value, intent(in) :: time
    !> Pointer to current state vector.
    type(c_ptr), value, intent(in) :: state_p
    !> Pointer to current state derivative vector.
    type(c_ptr), value, intent(in) :: state_dot_p
    !> Pointer to right-hand-side vector.
    type(c_ptr), value, intent(in) :: rhs_p
    !> Pointer to solution vector.
    type(c_ptr), value, intent(in) :: soln_p
    !> Value of \f$ \gamma \f$ parameter.
    real(kind=c_double), value, intent(in) :: gamma
    !> Input tolerance.
    real(kind=c_double), value, intent(in) :: delta
    !> Whether to use the left (1) or right (2) preconditioner.
    integer(kind=c_int), value, intent(in) :: left_or_right

    real(kind=c_double), pointer :: rhs(:), soln(:)
    real(kind=dp) :: dDdot_dD(n_eqn - 1), dDdot_dH(n_eqn - 1)
    real(kind=dp) :: dHdot_dD(n_eqn - 1), dHdot_dH
    integer :: i_part

    !>DEBUG
    condense_n_prec = condense_n_prec + 1
    !<DEBUG
    call condense_jac(n_eqn, time, state_p, dDdot_dD, dDdot_dH, &
         dHdot_dD, dHdot_dH)

    call c_f_pointer(rhs_p, rhs, (/ n_eqn /))
    call c_f_pointer(soln_p, soln, (/ n_eqn /))
    
    do i_part = 1,(n_eqn - 1)
       soln(i_part) = rhs(i_part) / (real(1d0, kind=c_double) &
            - gamma * real(dDdot_dD(i_part), kind=c_double))
    end do
    soln(n_eqn) = rhs(n_eqn) / (real(1d0, kind=c_double) &
         - gamma * real(dHdot_dH, kind=c_double))

  end subroutine condense_prec_diag_f

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Precondition by solving the system \f$ Pz = r \f$ where \f$ P
  !> \approx M = I - \gamma J \f$ and \f$ J = \partial f / \partial y
  !> \f$.
  subroutine condense_prec_exact_f(n_eqn, time, state_p, state_dot_p, rhs_p, &
       soln_p, gamma, delta, left_or_right) bind(c)

    !> Length of state vector.
    integer(kind=c_int), value, intent(in) :: n_eqn
    !> Current time (s).
    real(kind=c_double), value, intent(in) :: time
    !> Pointer to current state vector.
    type(c_ptr), value, intent(in) :: state_p
    !> Pointer to current state derivative vector.
    type(c_ptr), value, intent(in) :: state_dot_p
    !> Pointer to right-hand-side vector.
    type(c_ptr), value, intent(in) :: rhs_p
    !> Pointer to solution vector.
    type(c_ptr), value, intent(in) :: soln_p
    !> Value of \f$ \gamma \f$ parameter.
    real(kind=c_double), value, intent(in) :: gamma
    !> Input tolerance.
    real(kind=c_double), value, intent(in) :: delta
    !> Whether to use the left (1) or right (2) preconditioner.
    integer(kind=c_int), value, intent(in) :: left_or_right

    real(kind=c_double), pointer :: rhs(:), soln(:)
    real(kind=dp) :: dDdot_dD(n_eqn - 1), dDdot_dH(n_eqn - 1)
    real(kind=dp) :: dHdot_dD(n_eqn - 1), dHdot_dH
    real(kind=dp) :: lhs_n, rhs_n
    integer :: i_part

    !>DEBUG
    condense_n_prec = condense_n_prec + 1
    !<DEBUG
    call condense_jac(n_eqn, time, state_p, dDdot_dD, dDdot_dH, &
         dHdot_dD, dHdot_dH)

    call c_f_pointer(rhs_p, rhs, (/ n_eqn /))
    call c_f_pointer(soln_p, soln, (/ n_eqn /))

    lhs_n = 1d0 - gamma * dHdot_dH
    rhs_n = rhs(n_eqn)
    do i_part = 1,(n_eqn - 1)
       lhs_n = lhs_n - (- gamma * dDdot_dH(i_part)) &
            * (- gamma * dHdot_dD(i_part)) / (1d0 - gamma * dDdot_dD(i_part))
       rhs_n = rhs_n - (- gamma * dHdot_dD(i_part)) * rhs(i_part) &
            / (1d0 - gamma * dDdot_dD(i_part))
    end do
    soln(n_eqn) = rhs_n / lhs_n

    do i_part = 1,(n_eqn - 1)
       soln(i_part) = (rhs(i_part) &
            - (- gamma * dDdot_dH(i_part)) * soln(n_eqn)) &
            / (1d0 - gamma * dDdot_dD(i_part))
    end do

  end subroutine condense_prec_exact_f

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine condense_compare_jtimes(n_eqn, time, state, i_col, eps, &
       jac_col, fd_jac_col)

    !> Length of state vector.
    integer, value, intent(in) :: n_eqn
    !> Current time (s).
    real(kind=dp), value, intent(in) :: time
    !> Current state vector.
    real(kind=dp), intent(in) :: state(n_eqn)
    !> Column of the jacobian to compute.
    integer, intent(in) :: i_col
    !> Scaling for finite difference approximation.
    real(kind=dp), intent(in) :: eps
    !> Exact column of the jacobian.
    real(kind=dp), intent(out) :: jac_col(n_eqn)
    !> Finite difference column of the jacobian.
    real(kind=dp), intent(out) :: fd_jac_col(n_eqn)

    real(kind=c_double), target :: state_f(n_eqn)
    real(kind=c_double), target :: state_plus_eps_f(n_eqn)
    real(kind=c_double), target :: state_dot_f(n_eqn)
    real(kind=c_double), target :: state_plus_eps_dot_f(n_eqn)
    real(kind=c_double), target :: vec_f(n_eqn)
    real(kind=c_double), target :: jac_times_vec_f(n_eqn)
    integer :: i

    do i = 1,n_eqn
       state_f(i) = real(state(i), kind=c_double)
    end do
    state_plus_eps_f = state_f
    state_plus_eps_f(i_col) = state_plus_eps_f(i_col) &
         + real(eps, kind=c_double)
    vec_f = real(0d0, kind=c_double)
    vec_f(i_col) = real(1d0, kind=c_double)

    call condense_vf_f(int(n_eqn, kind=c_int), real(time, kind=c_double), &
         c_loc(state_f), c_loc(state_dot_f))
    call condense_vf_f(int(n_eqn, kind=c_int), real(time, kind=c_double), &
         c_loc(state_plus_eps_f), c_loc(state_plus_eps_dot_f))
    call condense_jtimes_f(int(n_eqn, kind=c_int), &
         real(time, kind=c_double), &
         c_loc(state_f), c_loc(vec_f), c_loc(jac_times_vec_f))

    do i = 1,n_eqn
       jac_col(i) = real(jac_times_vec_f(i), kind=dp)
       fd_jac_col(i) = real(state_plus_eps_dot_f(i) - state_dot_f(i), &
            kind=dp) / eps
    end do

  end subroutine condense_compare_jtimes

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine condense_test_jtimes(n_eqn, time, state)

    !> Length of state vector.
    integer, value, intent(in) :: n_eqn
    !> Current time (s).
    real(kind=dp), value, intent(in) :: time
    !> Current state vector.
    real(kind=dp), intent(in) :: state(n_eqn)

    real(kind=dp), parameter :: eps_big = 1d-6
    real(kind=dp), parameter :: eps_small = 1d-12
    real(kind=dp), parameter :: error_factor = 10d0

    real(kind=dp) :: jac_col(n_eqn), fd_jac_col(n_eqn)
    real(kind=dp) :: fd_jac_col_small(n_eqn), fd_jac_col_big(n_eqn)
    real(kind=dp) :: error_big, error_small, error, alpha, eps
    integer :: i_col, i_row, i_eps
    real(kind=dp) :: error_ratio, expected_error_ratio

    do i_col = 1,n_eqn
       call condense_compare_jtimes(n_eqn, time, state, i_col, eps_big, &
            jac_col, fd_jac_col_big)
       call condense_compare_jtimes(n_eqn, time, state, i_col, eps_small, &
            jac_col, fd_jac_col_small)
       do i_row = 1,n_eqn
          error_big = abs(fd_jac_col_big(i_row) - jac_col(i_row))
          error_small = abs(fd_jac_col_small(i_row) - jac_col(i_row))
          error_ratio = error_small / error_big
          expected_error_ratio = eps_small / eps_big * error_factor
          if (error_ratio > expected_error_ratio) then
             write(0,*) '***********************************************'
             write(0,*) 'row,col = ', i_row, i_col
             write(0,*) 'jac, fd_jac_small, fd_jac_big = ', &
                  jac_col(i_row), fd_jac_col_small(i_row), &
                  fd_jac_col_big(i_row)
             write(0,*) 'error_small, error_big = ', error_small, error_big
             write(0,*) 'error_ratio, expected = ', &
                  error_ratio, expected_error_ratio

             do i_eps = 1,101
                alpha = real(i_eps - 1, kind=dp) / 100d0
                eps = exp((1d0 - alpha) * log(eps_big) &
                     + alpha * log(eps_small))
                call condense_compare_jtimes(n_eqn, time, state, i_col, eps, &
                     jac_col, fd_jac_col)
                error = abs(fd_jac_col(i_row) - jac_col(i_row))
                write(0,*) 'eps, error, rel_error = ', &
                     eps, error, error / jac_col(i_row)
             end do
          end if
       end do
    end do
    stop
  
  end subroutine condense_test_jtimes
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
end module pmc_condense
