! Copyright (C) 2005-2009 Nicole Riemer and Matthew West
! Copyright (C) 2009 Joseph Ching
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_condensation module.

!> Water condensation onto aerosol particles.
module pmc_condensation

  use pmc_aero_state
  use pmc_bin_grid
  use pmc_env_data
  use pmc_env_state
  use pmc_aero_data
  use pmc_util
  use pmc_aero_particle
  use pmc_constants
  use dvode_f90_m
  
  !> Relative error tolerance (scalar) for all solution components.
  real(kind=dp), parameter :: CONDENSE_REL_TOL = 1d-8
  !> Lower bound on the absolute error tolerance for particle diameters (m).
  real(kind=dp), parameter :: CONDENSE_D_ABS_TOL_MIN = 1d-30
  !> Scale factor \$ \gamma \$ to determine per-particle tolerance \$ \Delta D = \gamma (D - D_{\rm dry}) \$.
  real(kind=dp), parameter :: CONDENSE_D_ABS_TOL_FACTOR = 1d-8
  !> Absolute error tolerance for relative humidity.
  real(kind=dp), parameter :: CONDENSE_RH_ABS_TOL = 1d-10
  !> Characteristic time to maintain error control (s).
  real(kind=dp), parameter :: CONDENSE_ERROR_CHAR_TIME = 1000d0

  logical, save :: rh_first = .true.
  integer, save :: d_offset = 1 ! 1 if rh_first, otherwise 0

  type(aero_data_t) :: condense_saved_aero_data
  type(env_data_t) :: condense_saved_env_data
  type(env_state_t) :: condense_saved_env_state_initial
  real(kind=dp) :: condense_saved_V_comp_initial
  real(kind=dp) :: condense_saved_Tdot
  real(kind=dp), pointer, save :: condense_saved_kappa(:)
  real(kind=dp), pointer, save :: condense_saved_D_dry(:)

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

  type condense_equilib_inputs_t
     !> Temperature (K).
     real(kind=dp) :: T
     !> Relative humidity (1).
     real(kind=dp) :: H
     !> Pressure (Pa).
     real(kind=dp) :: p
     !> Particle dry diameter (m).
     real(kind=dp) :: D_dry
     !> Kappa parameter (1).
     real(kind=dp) :: kappa
  end type condense_equilib_inputs_t

  type condense_equilib_outputs_t
     !> Partial diameter (m).
     real(kind=dp) :: D
  end type condense_equilib_outputs_t

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
    
    integer :: i_bin, j, new_bin, k, i_part
    real(kind=dp) :: pre_water, post_water
    type(aero_particle_t), pointer :: aero_particle
    integer :: n_eqn, itask, istate
    real(kind=dp) :: state(aero_state%n_part + 1), init_time, final_time
    real(kind=dp) :: abs_tol_vector(aero_state%n_part + 1)
    type(vode_opts) :: options
    type(env_state_t) :: env_state_final
    real(kind=dp) :: V_comp_final, water_vol_initial, water_vol_final
    real(kind=dp) :: d_water_vol_rh, d_water_vol_particles, water_rel_err

    pre_water = 0d0
    do i_bin = 1,bin_grid%n_bin
       do j = 1,aero_state%bin(i_bin)%n_part
          aero_particle => aero_state%bin(i_bin)%particle(j)
          pre_water = pre_water + aero_particle%vol(aero_data%i_water)
       end do
    end do

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
         env_state_final%elapsed_time + del_t)
    condense_saved_Tdot = (env_state_final%temp - env_state%temp) / del_t

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
          state(i_part + d_offset) = aero_particle_diameter(aero_particle)
          abs_tol_vector(i_part + d_offset) = max(CONDENSE_D_ABS_TOL_MIN, &
                CONDENSE_D_ABS_TOL_FACTOR &
                * (state(i_part + d_offset) - condense_saved_D_dry(i_part)))
       end do
    end do
    if (rh_first) then
       state(1) = env_state%rel_humid
       abs_tol_vector(1) = CONDENSE_RH_ABS_TOL
    else
       state(aero_state%n_part + 1) = env_state%rel_humid
       abs_tol_vector(aero_state%n_part + 1) = CONDENSE_RH_ABS_TOL
    end if

    ! set VODE inputs
    n_eqn = aero_state%n_part + 1
    init_time = 0d0
    final_time = del_t
    itask = 1 ! just output val at final_time
    istate = 1 ! first call for this ODE

    options = set_opts(sparse_j = .true., &
         user_supplied_jacobian = .true., &
         abserr_vector = abs_tol_vector, relerr = CONDENSE_REL_TOL)
    call dvode_f90(condense_vode_f, n_eqn, state, init_time, final_time, &
         itask, istate, options, j_fcn = condense_vode_jac)

    if (rh_first) then
       env_state%rel_humid = state(1)
    else
       env_state%rel_humid = state(aero_state%n_part + 1)
    end if

    i_part = 0
    do i_bin = 1,bin_grid%n_bin
       do j = 1,aero_state%bin(i_bin)%n_part
          i_part = i_part + 1
          aero_particle => aero_state%bin(i_bin)%particle(j)

          ! translate output back to particle
          aero_particle%vol(aero_data%i_water) = diam2vol(state(i_part + d_offset)) &
               - aero_particle_solute_volume(aero_particle, aero_data)

          ! ensure volumes stay positive
          aero_particle%vol(aero_data%i_water) = max(0d0, &
               aero_particle%vol(aero_data%i_water))
       end do
    end do

    ! We resort the particles in the bins after only all particles
    ! have condensation done, otherwise we will lose track of which
    ! ones have had condensation and which have not.
    call aero_state_resort(bin_grid, aero_state)

    post_water = 0d0
    do i_bin = 1,bin_grid%n_bin
       do j = 1,aero_state%bin(i_bin)%n_part
          aero_particle => aero_state%bin(i_bin)%particle(j)
          post_water = post_water + aero_particle%vol(aero_data%i_water)
       end do
    end do

    !>DEBUG
    V_comp_final = condense_saved_V_comp_initial * env_state_final%temp / condense_saved_env_state_initial%temp
    water_vol_initial = aero_data%molec_weight(aero_data%i_water) &
         / (const%univ_gas_const * condense_saved_env_state_initial%temp) &
         * env_state_sat_vapor_pressure(condense_saved_env_state_initial) &
         * condense_saved_env_state_initial%rel_humid &
         * condense_saved_V_comp_initial / aero_particle_water_density(aero_data)
    water_vol_final = aero_data%molec_weight(aero_data%i_water) &
         / (const%univ_gas_const * env_state_final%temp) &
         * env_state_sat_vapor_pressure(env_state_final) &
         * env_state%rel_humid &
         * V_comp_final / aero_particle_water_density(aero_data)
    d_water_vol_rh = water_vol_final - water_vol_initial
    d_water_vol_particles = post_water - pre_water
    water_rel_err = (d_water_vol_rh + d_water_vol_particles) / water_vol_final
    call warn_assert_msg(477865387, water_rel_err < 1d-6, &
         "condensation water balance error exceeded bounds")
    !<DEBUG

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
         env_state, env_state%elapsed_time + time)

    if (rh_first) then
       env_state%rel_humid = state(1)
    else
       env_state%rel_humid = state(n_eqn)
    end if

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
       inputs%D = state(i_part + d_offset)
       inputs%D_dry = condense_saved_D_dry(i_part)
       inputs%kappa = condense_saved_kappa(i_part)
       call condense_rates(inputs, outputs)
       state_dot(i_part + d_offset) = outputs%Ddot
       Hdot = Hdot + outputs%Hdot_i
    end do
    Hdot = Hdot + outputs%Hdot_env

    if (rh_first) then
       state_dot(1) = Hdot
    else
       state_dot(n_eqn) = Hdot
    end if

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
       inputs%D = state(i_part + d_offset)
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
    if (rh_first) then
       ia(1) = 1
       ia(2) = ia(1) + n_eqn
       do i_part = 1,(n_eqn - 1)
          ia(2 + i_part) = ia(1 + i_part) + 2
       end do
       call assert(351522940, ia(n_eqn + 1) - 1 == n_non_zero)

       i_nz = 0
       i_nz = i_nz + 1
       ja(i_nz) = 1
       state_jac(i_nz) = dHdot_dH
       do i_part = 1,(n_eqn - 1)
          i_nz = i_nz + 1
          ja(i_nz) = 1 + i_part
          state_jac(i_nz) = dDdot_dH(i_part)
       end do
       do i_part = 1,(n_eqn - 1)
          i_nz = i_nz + 1
          ja(i_nz) = 1
          state_jac(i_nz) = dHdot_dD(i_part)
          i_nz = i_nz + 1
          ja(i_nz) = i_part + 1
          state_jac(i_nz) = dDdot_dD(i_part)
       end do
       call assert(260941580, i_nz == n_non_zero)
    else
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
    end if

    call env_state_deallocate(env_state)

  end subroutine condense_vode_jac

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Determine the equilibrium state of a single particle.
  subroutine condense_equilib_particle(inputs, outputs)

    !> Inputs to equilibriation.
    type(condense_equilib_inputs_t), intent(in) :: inputs
    !> Outputs equilibrium.
    type(condense_equilib_outputs_t), intent(out) :: outputs

    real(kind=dp) :: X, D, g, dg_dD, a_w, daw_dD
    integer :: newton_step

    X = 4d0 * const%water_molec_weight * const%water_surf_eng &
         / (const%univ_gas_const * inputs%T * const%water_density)

    D = inputs%D_dry
    newton_step = 0
    g = 0d0
    dg_dD = 1d0
    do while (.true.)
       newton_step = newton_step + 1
       D = D - g / dg_dD
       a_w = (D**3 - inputs%D_dry**3) &
            / (D**3 + (inputs%kappa - 1d0) * inputs%D_dry**3)
       daw_dD = 3d0 * D**2 * inputs%kappa * inputs%D_dry**3 &
            / (D**3 + (inputs%kappa - 1d0) * inputs%D_dry**3)**2
       g = inputs%H - a_w * exp(X / D)
       dg_dD = - daw_dD * exp(X / D) + a_w * exp(X / D) * (X / D**2)
       if (abs(g) < 1d3 * epsilon(1d0)) then
          exit
       end if
       if (newton_step > 100) then
          exit
       end if
    end do
    write(*,*) 'newton_step ', newton_step

    outputs%D = D

  end subroutine condense_equilib_particle

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Scalar Newton's method for solving the implicit condensation
  !> functions.
  subroutine cond_newt(x, env_state, aero_data, func, x_tol, f_tol, &
       iter_max, aero_particle)

    !> Variable (set to init value on call).
    real(kind=dp), intent(inout) :: x
    !> Environment state.
    type(env_state_t), intent(in) :: env_state
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> X convergence tolerance.
    real(kind=dp), intent(in) :: x_tol
    !> F convergence tolerance.
    real(kind=dp), intent(in) :: f_tol
    !> Maximum number of iterations.
    integer, intent(in) :: iter_max
    !> Particle.
    type(aero_particle_t), intent(in) :: aero_particle

#ifndef DOXYGEN_SKIP_DOC
    interface
       subroutine func(env_state, aero_data, init, x, f, df, aero_particle)
         use pmc_env_state
         use pmc_aero_data
         use pmc_aero_particle
         !> Environment state.
         type(env_state_t), intent(in) :: env_state
         !> Aerosol data.
         type(aero_data_t), intent(in) :: aero_data
         !> True if first Newton loop.
         logical, intent(in) :: init
         !> Independent variable to solve for.
         real(kind=dp), intent(in) :: x
         !> Function to solve.
         real(kind=dp), intent(out) :: f
         !> Derivative df/dx.
         real(kind=dp), intent(out) :: df
         !> Particle.
         type(aero_particle_t), intent(in) :: aero_particle
       end subroutine func
    end interface
#endif
    
    integer iter, k
    real(kind=dp) delta_f, delta_x, f, old_f, df
    real(kind=dp) :: check_x, check_f, check_df, finite_diff_df, rel_error

    call func(env_state, aero_data, .true., x, f, df, aero_particle)
    old_f = f

    iter = 0
    do
       !write(*,*) 'no. of iteration', iter
       iter = iter + 1
       delta_x = f / df
       !write(*,*) 'old x =' , x , 'old f=' ,old_f
       
       x = x - delta_x
       call func(env_state, aero_data, .false., x, f, df, aero_particle)
       delta_f = f - old_f
       old_f = f

       !write(*,*) 'new x =' , x , 'new f=' , f
       !write(*,*) 'delta x =', delta_x , 'delta_f=', delta_f
       if (iter .ge. iter_max) then
          call die_msg(449343787, 'Newton iteration failed to terminate')
       end if
       
       ! FIXME: gfortran 4.1.1 requires the "then" in the following
       ! statement, rather than using a single-line "if" statement.
       if ((abs(delta_x) .lt. x_tol) &
            .and. (abs(delta_f) .lt. f_tol)) then
          exit
       end if
    end do
 
    !write(*,*) 'id =', aero_particle%id, 'finish using', iter ,'iteration(s)'
  end subroutine cond_newt

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Add water to the particle until it is in equilibrium.
  subroutine equilibriate_particle(env_state, aero_data, aero_particle)

    !> Environment state.
    type(env_state_t), intent(in) :: env_state
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Particle.
    type(aero_particle_t), intent(inout) :: aero_particle

    ! parameters
    !> Maximum iterations.
    integer, parameter :: it_max = 400
    !> Volume relative convergence tolerance.
!    real(kind=dp), parameter :: pv_rel_tol = 1d-6
    real(kind=dp), parameter :: pv_rel_tol = 1d-3
    !> Function convergence tolerance.
    real(kind=dp), parameter :: f_tol = 1d-15
    !> Maximum number of iterations.
    integer, parameter :: iter_max = 100
    !> Initial value.
    !> want to get variable initial guess as the dry diameter
!    real(kind=dp), parameter :: dw_init = 1d-6
    
    real(kind=dp) dw ! wet diameter of particle
    real(kind=dp) dw_tol, pv, di
    real(kind=dp) dw_init     
    real(kind=dp) :: D_init, D_final, D_final_new, kappa
    type(condense_equilib_inputs_t) :: inputs
    type(condense_equilib_outputs_t) :: outputs

    !> total volume of the particle before equilibriation
    pv = aero_particle_volume(aero_particle)
    di = vol2diam(pv)
    !write(*,*) 'id =', aero_particle%id, 'wet diam (before equilibriation)', di  
    !> dry diameter of the particle
    dw_init = vol2diam(aero_particle_solute_volume(aero_particle, aero_data))
    D_init = dw_init
    !write(*,*) 'id =', aero_particle%id, 'dry diam (before iteration) =', dw_init
    dw = dw_init

    ! convert volume relative tolerance to diameter absolute tolerance
    dw_tol = vol2diam(pv * (1d0 + pv_rel_tol)) - vol2diam(pv)
    !write(*,*) 'sat (in equilibriate_particle) BF newton =' , env_state%rel_humid
!    write(*,*) 'dw_tol = ', dw_tol
!    write(*,*) 'id =', aero_particle%id 

    call cond_newt(dw, env_state, aero_data, equilibriate_func, &
         dw_tol, f_tol, iter_max, aero_particle)

    !write(*,*) 'sat (in equilibriate_particle) AF newton =' , env_state%rel_humid
!    aero_particle%vol(aero_data%i_water) = diam2vol(dw) - pv
    !> the amount of water after equilibriation
    !> dw is the wet diameter after equilibriation 
    !> dw_init is the dry diameter before equilibriation
    aero_particle%vol(aero_data%i_water) = diam2vol(dw) - diam2vol(dw_init)
    di=vol2diam(pv)

    D_final = vol2diam(aero_particle_volume(aero_particle))

    inputs%T = env_state%temp
    inputs%H = env_state%rel_humid
    inputs%p = env_state%pressure
    inputs%D_dry = D_init
    inputs%kappa = aero_particle_solute_kappa(aero_particle, aero_data)

    call condense_equilib_particle(inputs, outputs)
    D_final_new = outputs%D

    call compute_kappa(aero_data, aero_particle, kappa)
    write(*,*) 'D_dry,D,D_new,kappa,kappa_new ', D_init, D_final, D_final_new, &
         kappa, inputs%kappa

    !write(*,*)  'id = ', aero_particle%id, 'wet diam (after iteration) = ', dw
  end subroutine equilibriate_particle

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Return the error function value and its derivative for the
  !> implicit function that determines the equilibrium state of a
  !> particle.
  subroutine equilibriate_func(env_state, aero_data, init, dw, f, df, &
       aero_particle)

    !> Environment state.
    type(env_state_t), intent(in) :: env_state
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> True if first Newton loop.
    logical, intent(in) :: init
    !> Wet diameter (m).
    real(kind=dp), intent(in) :: dw
    !> Function value.
    real(kind=dp), intent(out) :: f
    !> Function derivative df/dx.
    real(kind=dp), intent(out) :: df
    !> Particle.
    type(aero_particle_t), intent(inout) :: aero_particle

!    call equilibriate_func_old(env_state, aero_data, init, dw, f, df, &
!       aero_particle)

!    call equilibriate_func_new(env_state, aero_data, init, dw, f, df, &
!       aero_particle)

!     write(*,*) 'sat (in equilibriate_func) 1 =', env_state%rel_humid
    call equilibriate_func_aw(env_state, aero_data, init, dw, f, df, &
       aero_particle)
    
!     write(*,*) 'sat (in equiliraiate_func) 2 =', env_state%rel_humid
!    write(*,*) 'dw = ', dw
!    write(*,*) 'f = ', f
!    write(*,*) 'df = ', df
!    write(*,*) 'id = ', aero_particle%id

    end subroutine equilibriate_func

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Return the error function value and its derivative for the
  !> implicit function that determines the equilibrium state of a
  !> particle.
  subroutine equilibriate_func_old(env_state, aero_data, init, dw, f, df, &
       aero_particle)

    !> Environment state.
    type(env_state_t), intent(in) :: env_state
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> True if first Newton loop.
    logical, intent(in) :: init
    !> Wet diameter (m).
    real(kind=dp), intent(in) :: dw
    !> Function value.
    real(kind=dp), intent(out) :: f
    !> Function derivative df/dx.
    real(kind=dp), intent(out) :: df
    !> Particle.
    type(aero_particle_t), intent(inout) :: aero_particle

    real(kind=dp), save :: c0, c1, c3, c4, dc0, dc2, dc3
    real(kind=dp), save :: A, B
    real(kind=dp), save ::  pv
    real(kind=dp), save :: M_water, M_solute, rho_water, rho_solute
    real(kind=dp), save :: eps, nu, g_water, g_solute

    if (init) then
       ! Start of new Newton loop, compute all constants

       M_water = aero_particle_water_molec_weight(aero_data) ! (kg/mole)
       M_solute = aero_particle_solute_molec_weight(aero_particle, aero_data)
       ! (kg/mole)
       nu = aero_particle_solute_num_ions(aero_particle, aero_data) ! (1)
       eps = aero_particle_solute_solubility(aero_particle, aero_data) ! (1)
       rho_water = aero_particle_water_density(aero_data) ! (kg/m^3)
       rho_solute = aero_particle_solute_density(aero_particle, aero_data)
       ! (kg/m^3)
       g_water = aero_particle_water_mass(aero_particle, aero_data) ! (kg)
       g_solute = aero_particle_solute_mass(aero_particle, aero_data) ! (kg)

       pv = aero_particle_volume(aero_particle) ! (m^3)
    
       A = 4d0 * M_water * const%water_surf_eng &
            / (const%univ_gas_const * env_state%temp * rho_water)
       
       B = nu * eps * M_water * rho_solute &
            * vol2rad(pv)**3.d0 / (M_solute * rho_water)
       
       c4 = log(env_state%rel_humid) / 8d0
       c3 = A / 8d0
       
       dc3 = log(env_state%rel_humid) / 2d0
       dc2 = 3d0 * A / 8d0
       
       c1 = B - log(env_state%rel_humid) * vol2rad(pv)**3d0
       c0 = A * vol2rad(pv)**3d0
       dc0 = c1

    endif
    
    f = c4 * dw**4d0 - c3 * dw**3d0 + c1 * dw + c0
    df = dc3 * dw**3d0 -dc2 * dw**2d0 + dc0

  end subroutine equilibriate_func_old

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Return the error function value and its derivative for the
  !> implicit function that determines the equilibrium state of a
  !> particle.
  subroutine equilibriate_func_new(env_state, aero_data, init, dw, f, df, &
       aero_particle)

    !> Environment state.
    type(env_state_t), intent(in) :: env_state
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> True if first Newton loop.
    logical, intent(in) :: init
    !> Wet diameter (m).
    real(kind=dp), intent(in) :: dw
    !> Function value.
    real(kind=dp), intent(out) :: f
    !> Function derivative df/dx.
    real(kind=dp), intent(out) :: df
    !> Particle.
    type(aero_particle_t), intent(inout) :: aero_particle

    real(kind=dp), save :: c1, c2, c3, c4, c5, c6
    real(kind=dp), save :: A, B
    real(kind=dp), save ::  pv
    real(kind=dp), save :: M_water, M_solute, rho_water, rho_solute
    real(kind=dp), save :: eps, nu, g_water, g_solute

    if (init) then
       ! Start of new Newton loop, compute all constants

       M_water = aero_particle_water_molec_weight(aero_data) ! (kg/mole)
       M_solute = aero_particle_solute_molec_weight(aero_particle, aero_data)
       ! (kg/mole)
       nu = aero_particle_solute_num_ions(aero_particle, aero_data) ! (1)
       eps = aero_particle_solute_solubility(aero_particle, aero_data) ! (1)
       rho_water = aero_particle_water_density(aero_data) ! (kg/m^3)
       rho_solute = aero_particle_solute_density(aero_particle, aero_data)
       ! (kg/m^3)
       g_water = aero_particle_water_mass(aero_particle, aero_data) ! (kg)
       g_solute = aero_particle_solute_mass(aero_particle, aero_data) ! (kg)

       pv = aero_particle_volume(aero_particle) ! (m^3)

       c1=log(env_state%rel_humid)
       c2 = 4d0 * M_water &
            * const%water_surf_eng / (env_state%temp * const%univ_gas_const * rho_water)
       c3=(6d0*nu*M_water) / (const%pi*rho_water*M_solute) * g_solute 
       
!       c4=1d-9                   ! hardcode insoluble diameter
        c4=0d0
       c5=c2*(c4**3d0)
       c6=c1*(c4**3d0) 
      
      endif
    
    f = c1*dw**4d0-c2*dw**3d0+(c3-c6)*dw+c5     
    df= 4d0*c1*dw**3d0 - 3d0*c2*dw**2d0 + (c3-c6) 

    ! NOTE: we could return T_a (the droplet temperature) if we have
    ! any need for it.


  end subroutine equilibriate_func_new

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Return the error function value and its derivative for the
  !> implicit function that determines the equilibrium state of a
  !> particle.
  subroutine equilibriate_func_aw(env_state, aero_data, init, dw, f, df, &
       aero_particle)

    !> Environment state.
    type(env_state_t), intent(in) :: env_state
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> True if first Newton loop.
    logical, intent(in) :: init
    !> Wet diameter (m).
    real(kind=dp), intent(in) :: dw
    !> Function value.
    real(kind=dp), intent(out) :: f
    !> Function derivative df/dx.
    real(kind=dp), intent(out) :: df
    !> Particle.
    type(aero_particle_t), intent(inout) :: aero_particle

    real(kind=dp) :: sat, tot_kappa, kappa
    
    integer :: k
 
!       write(*,*) 'sat (in equilibriate_func_aw) 0 = ', env_state%rel_humid
!       write(*,*) 'init =', init
    if (init) then
          
       ! Start of new Newton loop, compute all constants

       call compute_kappa(aero_data, aero_particle, tot_kappa)   
       kappa=tot_kappa
!       write(*,*) ' id = ', aero_particle%id, 'kappa = ' , kappa       

       sat = env_state%rel_humid
!       write(*,*) 'sat (in equilibriate_func_aw) 1 = ', sat
   
    endif
         
!       write(*,*) 'now start to calculate f and df in equilibration'
!       write(*,*) ' id = ', aero_particle%id, 'kappa = ', kappa

       sat = env_state%rel_humid
       
        f= sat - aw(aero_data, aero_particle, dw) & 
              * exp(kelvin_argu(aero_particle, aero_data, env_state, dw))
       df= -  aw_deri(aero_data, aero_particle, dw) & 
              * exp(kelvin_argu(aero_particle, aero_data, env_state, dw)) &
              + kelvin_deri(aero_particle, aero_data, env_state, dw) &
              * aw(aero_data, aero_particle, dw)

!       write(*,*) 'f = ', f 
!       write(*,*) 'df = ', df
!       write(*,*)  'wet diam (after iteration) = ', dw
!       write(*,*) 'sat (in equilibriate_func_aw) 2 = ', sat
    
   end subroutine equilibriate_func_aw

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> calculate the argument of the Kelvin term.
  real(kind=dp) function kelvin_argu(aero_particle, aero_data, env_state, dw)

    !> Particle.
    type(aero_particle_t), intent(in) :: aero_particle
    !> Aerosol Particle.
    type(aero_data_t), intent(in) :: aero_data
    !> Environment.
    type(env_state_t), intent(in) :: env_state
    !> Wet diameter
    real(kind=dp), intent(in) :: dw

    real(kind=dp), save :: M_water, rho_water

!      write(*,*) 'now in function kelvin_term'
!      write(*,*) 'input dw to the funtion kelvin_argu = ', dw

       M_water = aero_particle_water_molec_weight(aero_data) ! (kg/mole)
       rho_water = aero_particle_water_density(aero_data) ! (kg/m^3)

    kelvin_argu= 4d0 * M_water * const%water_surf_eng &
         / (env_state%temp * const%univ_gas_const * rho_water * dw) 
    
!    write(*,*) 'input dw to the funtion aw = ', dw
!    write(*,*) 'mole weight of h20 =', M_water
!    write(*,*) 'temp =', env_state%temp
!    write(*,*) 'gas const =', const%univ_gas_const
!    write(*,*) 'rho water =', rho_water 
!    write(*,*) 'surface tension =', const%water_surf_eng
!    write(*,*) aero_particle%id, 'kelvin argu= ', kelvin_argu

  end function kelvin_argu

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> calculate the derivative of the Kelvin term.
  real(kind=dp) function kelvin_deri(aero_particle, aero_data, env_state, dw)

    !> Particle.
    type(aero_particle_t), intent(in) :: aero_particle
    !> Aerosol Particle.
    type(aero_data_t), intent(in) :: aero_data
    !> Environment.
    type(env_state_t), intent(in) :: env_state
    !> Wet diameter.
    real(kind=dp), intent(in) :: dw

!      write(*,*) 'now in function kelvin_deri'      
!      write(*,*) 'input dw to the funtion aw = ', dw

    kelvin_deri=exp(kelvin_argu(aero_particle, aero_Data, env_state, dw)) & 
             *kelvin_argu(aero_particle, aero_data, env_state, dw) / dw

!    write(*,*) aero_particle%id, 'kelvin deri= ', kelvin_deri
  end function kelvin_deri

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> calculate the water activity.
  real(kind=dp) function aw(aero_data, aero_particle, dw)

    !> Particle.
    type(aero_data_t), intent(in) :: aero_data
    !> Aerosol data
    type(aero_particle_t), intent(in) :: aero_particle

    real(kind=dp), intent(in) :: dw
    
    real(kind=dp) :: vs, vt, kappa, tot_kappa

    integer :: k

!       write(*,*) 'now in function aw' 

       vs=0d0
       vt=0d0
!   vs is the dry particle volume (non-water components)
       do k=1,aero_data%n_spec 
          if (k /= aero_data%i_water) then
          vs=vs+aero_particle%vol(k)
          endif
          vt=vt+aero_particle%vol(k)
       enddo
    
       call compute_kappa(aero_data, aero_particle, tot_kappa)   
       kappa=tot_kappa

!      write(*,*) aero_particle%id, 'kappa = ', kappa
!      write(*,*) 'input dw to the funtion aw = ', dw
!      write(*,*) 'dry salt volume in the function aw = ', vs
!      write(*,*) 'wet salt volume (dw^3) ', const%pi/6d0*(dw**3d0)     
       
      aw=(const%pi/6d0*(dw**3d0)-vs)/(const%pi/6d0*(dw**3d0)+(kappa-1d0)*vs)
  
!      write(*,*) 'kappa =', kappa   
!      write(*,*) aero_particle%id, 'aw= ', aw
  end function aw

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> calculate the derivative of the water activity.
  real(kind=dp) function aw_deri(aero_data, aero_particle, dw)

    !> Particle.
    type(aero_particle_t), intent(in) :: aero_particle
    !> Aerosol data
    type(aero_data_t), intent(in) :: aero_data

    real(kind=dp), intent(in) :: dw
    
    real(kind=dp) :: vs, vt, kappa, tot_kappa
    
    integer :: k

!       write(*,*) 'now in function aw_deri'

       vs=0d0
       vt=0d0
!   vs is the dry particle volume (non-water components)
       do k=1,aero_data%n_spec 
          if (k /= aero_data%i_water) then
          vs=vs+aero_particle%vol(k)
          endif
          vt=vt+aero_particle%vol(k)
       enddo
    
       call compute_kappa(aero_data, aero_particle, tot_kappa)   
       kappa=tot_kappa
      
!      write(*,*) 'input dw to the funtion aw = ', dw
!      write(*,*) 'dry salt volume in the function aw = ', vs

       aw_deri= (kappa*vs*const%pi/2d0*(dw**2d0))   &
       /((const%pi/6d0*(dw**3d0)+(kappa-1d0)*vs)**2d0)

!    write(*,*) aero_particle%id, 'aw_deri= ', aw_deri
  end function aw_deri

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Call equilibriate_particle() on each particle in the aerosol.
  subroutine aero_state_equilibriate(bin_grid, env_state, aero_data, &
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
         !write(*,*) 'hello id =', aero_particle%id
         env_state%rel_humid=9.50d-1
         !write(*,*) 'RH (in aero_state_equi) =', env_state%rel_humid
          call equilibriate_particle(env_state, aero_data, &
               aero_state%bin(i_bin)%particle(i))
!     temporary stop statement
          !write(*,*) 'stop at particle i', i
!        STOP
       end do
    end do
    stop

  end subroutine aero_state_equilibriate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Compute the Kappa value of individual particle based on volume
  !> fractions
  subroutine compute_kappa(aero_data, aero_particle, tot_kappa)

    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Aerosol particle.
    type(aero_particle_t), intent(in) :: aero_particle
    !> Kappa values for individual particle. 
    real(kind=dp), intent(out) :: tot_kappa
    !> Quantity to total. 
!    real(kind=dp), intent(in) :: quantity(:)

    ! local variables
    real(kind=dp) :: M_water, M_solute, rho_water, rho_solute
    real(kind=dp) :: g_water, g_solute
    real(kind=dp) :: indiv_kappa
    real(kind=dp) :: tot_vol, vol_frac, tot_vol_frac
    real(kind=dp) :: water_vol,other_vol
    integer :: j, k
    
       M_water = aero_particle_water_molec_weight(aero_data) ! (kg/mole)
       M_solute = aero_particle_solute_molec_weight(aero_particle, aero_data)
       ! (kg/mole)
       rho_water = aero_particle_water_density(aero_data) ! (kg/m^3)
       rho_solute = aero_particle_solute_density(aero_particle, aero_data)
       ! (kg/m^3)
       g_water = aero_particle_water_mass(aero_particle, aero_data) ! (kg)
       g_solute = aero_particle_solute_mass(aero_particle, aero_data) ! (kg)

       tot_kappa=0d0
       tot_vol=0d0
       tot_vol_frac=0d0

       do k=1,aero_data%n_spec 
          if (k /= aero_data%i_water) then
          ! total volume of solute (non-water components)
          tot_vol=tot_vol+aero_particle%vol(k)
!          write(*,*) 'k=',k ,'vol= ', aero_particle%vol(k)
          endif
       enddo
       do k=1,aero_data%n_spec 
          if (k /= aero_data%i_water) then
              if(aero_data%kappa(k) > 500d0) then
          indiv_kappa=(M_water/rho_water)*(aero_data%density(k)/aero_data%molec_weight(k))    !volume fraction
!          write(*,*) 'k=', k, 'mole weight =', M_water, 'den =', rho_water
!          write(*,*) 'k=', k, 'aero den =', aero_data%density(k), 'aero weight=', aero_data%molec_weight(k)   
!          write(*,*) 'inorganic', k, indiv_kappa
              else
          indiv_kappa=aero_data%kappa(k)
!          write(*,*) 'organic', k, indiv_kappa
              endif

          vol_frac=aero_particle%vol(k)/tot_vol
!          write(*,*) 'vol= ', vol_frac, 'kappa= ', indiv_kappa
          tot_kappa=tot_kappa+vol_frac*indiv_kappa
          tot_vol_frac=tot_vol_frac+vol_frac
          write(*,*) 'k,vol_frac,indiv_kappa ', k, vol_frac, indiv_kappa
          write(*,*) 'aero_data%kappa(k) ', aero_data%kappa(k)
          endif
       enddo
       stop
!       write(*,*) 'total vol frac=', tot_vol_frac
 
  end subroutine compute_kappa

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
end module pmc_condensation
