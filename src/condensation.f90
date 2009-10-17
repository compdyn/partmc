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

  type(aero_data_t) :: condense_saved_aero_data
  type(env_data_t) :: condense_saved_env_data
  type(env_state_t) :: condense_saved_env_state_initial
  real(kind=dp) :: condense_saved_V_comp_initial
  real(kind=dp) :: condense_saved_Tdot
  real(kind=dp), allocatable :: condense_saved_kappa(:)
  real(kind=dp), allocatable :: condense_saved_D_dry(:)

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
    
    integer :: i_bin, j, new_bin, k, i_part, n_eqn, itask, istate
    type(aero_particle_t), pointer :: aero_particle
    real(kind=dp) :: state(aero_state%n_part + 1), init_time, final_time
    real(kind=dp) :: abs_tol_vector(aero_state%n_part + 1)
    type(vode_opts) :: options
    type(env_state_t) :: env_state_final
    real(kind=dp) :: water_vol_initial, water_vol_final, d_water_vol
    real(kind=dp) :: vapor_vol_initial, vapor_vol_final, d_vapor_vol
    real(kind=dp) :: V_comp_final, water_rel_error

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
         env_state_final%elapsed_time + del_t)
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
    n_eqn = aero_state%n_part + 1
    init_time = 0d0
    final_time = del_t
    itask = 1 ! just output val at final_time
    istate = 1 ! first call for this ODE
    options = set_opts(sparse_j = .true., &
         user_supplied_jacobian = .true., &
         abserr_vector = abs_tol_vector, relerr = 1d-8)
    call dvode_f90(condense_vode_f, n_eqn, state, init_time, final_time, &
         itask, istate, options, j_fcn = condense_vode_jac)

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
         env_state, env_state%elapsed_time + time)
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
    
    !>DEBUG
    env_state%rel_humid=9.50d-1
    call warn_msg(842966545, "resetting rel_humid")
    !>DEBUG
    do i_bin = 1,bin_grid%n_bin
       do i = 1,aero_state%bin(i_bin)%n_part
          call condense_equilib_particle(env_state, aero_data, &
               aero_state%bin(i_bin)%particle(i))
       end do
    end do

  end subroutine condense_equilib_particles

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
end module pmc_condensation
