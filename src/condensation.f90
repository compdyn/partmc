! Copyright (C) 2005-2008 Nicole Riemer and Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_condensation module.

!> Water condensation onto aerosol particles.
module pmc_condensation

  use pmc_aero_state
  use pmc_bin_grid
  use pmc_aero_binned
  use pmc_env_state
  use pmc_aero_data
  use pmc_util
  use pmc_aero_particle
  use pmc_constants
  
contains
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Do condensation to all the particles for a given time interval,
  !> including updating the environment to account for the lost
  !> vapor.
  subroutine condense_particles(bin_grid, aero_binned, env_state, aero_data, &
       aero_state, del_t)

    !> Bin grid.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Binned distributions.
    type(aero_binned_t), intent(inout) :: aero_binned
    !> Environment state.
    type(env_state_t), intent(inout) :: env_state
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Aerosol state.
    type(aero_state_t), intent(inout) :: aero_state
    !> Total time to integrate.
    real*8, intent(in) :: del_t
    
    integer :: i_bin, j, new_bin, k
    real*8 :: pv, pre_water_vol, post_water_vol

    ! FIXME: don't rely on binned data, but rather compute the total
    ! water transfered in condense_particle() and return it
    pre_water_vol = sum(aero_binned%vol_den(:,aero_data%i_water)) &
         * aero_state%comp_vol * bin_grid%dlnr

    do i_bin = 1,bin_grid%n_bin
       do j = 1,aero_state%bin(i_bin)%n_part
          call condense_particle(del_t, env_state, aero_data, &
               aero_state%bin(i_bin)%particle(j))
       end do
    end do

    ! We resort the particles in the bins after only all particles
    ! have condensation done, otherwise we will lose track of which
    ! ones have had condensation and which have not.
    call aero_state_resort(bin_grid, aero_state)

    ! update the bin arrays
    call aero_state_to_binned(bin_grid, aero_data, aero_state, aero_binned)

    ! update the environment due to condensation of water
    post_water_vol = sum(aero_binned%vol_den(:,aero_data%i_water)) &
         * aero_state%comp_vol * bin_grid%dlnr
    call env_state_change_water_volume(env_state, aero_data, &
         (post_water_vol - pre_water_vol) / aero_state%comp_vol)

  end subroutine condense_particles

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Integrate the condensation growth or decay ODE for total time
  !> del_t for a single particle.
  subroutine condense_particle(del_t, env_state, aero_data, aero_particle)

    !> Total time to integrate.
    real*8, intent(in) :: del_t
    !> Environment state.
    type(env_state_t), intent(in) :: env_state
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Particle.
    type(aero_particle_t), intent(inout) :: aero_particle

    real*8 time_step, time
    logical done

    integer i
    real*8 dvdt

    time = 0d0
    done = .false.
    do while (.not. done)
       call condense_step_euler(del_t - time, time_step, done, env_state, &
            aero_data, aero_particle)
       time = time + time_step
    end do
    
  end subroutine condense_particle

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Does one timestep (determined by this subroutine) of the
  !> condensation ODE.
  !!
  !! The timestep will not exceed max_dt, but might be less. If we in
  !! fact step all the way to max_dt then done will be true. This uses
  !! the explicit (forward) Euler integrator.
  subroutine condense_step_euler(max_dt, dt, done, env_state, aero_data, &
       aero_particle)
    
    !> Maximum timestep to integrate.
    real*8, intent(in) :: max_dt
    !> Actual timestep used.
    real*8, intent(out) :: dt
    !> Did we reach the maximum timestep?.
    logical, intent(out) :: done
    !> Environment state.
    type(env_state_t), intent(in) :: env_state
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Particle.
    type(aero_particle_t), intent(inout) :: aero_particle

    real*8 dvdt

    ! get timestep
    done = .false.
    call find_condense_timestep_variable(dt, env_state, aero_data, &
         aero_particle)
    if (dt .ge. max_dt) then
       dt = max_dt
       done = .true.
    end if

    ! do condensation
    call cond_growth_rate(dvdt, env_state, aero_data, aero_particle)
    aero_particle%vol(aero_data%i_water) = &
         aero_particle%vol(aero_data%i_water) + dt * dvdt

    ! ensure volumes stay positive
    aero_particle%vol(aero_data%i_water) = max(0d0, &
         aero_particle%vol(aero_data%i_water))
   
  end subroutine condense_step_euler

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Does one timestep (determined by this subroutine) of the
  !> condensation ODE.
  !!
  !! The timestep will not exceed max_dt, but might be less. If we in
  !! fact step all the way to max_dt then done will be true. This uses
  !! the explicit 4th-order Runge-Kutta integrator.
  subroutine condense_step_rk_fixed(max_dt, &
       dt, done, env_state, aero_data, aero_particle)
    
    !> Maximum timestep to integrate.
    real*8, intent(in) :: max_dt
    !> Actual timestep used.
    real*8, intent(out) :: dt
    !> Did we reach the maximum timestep?.
    logical, intent(out) :: done
    !> Environment state.
    type(env_state_t), intent(in) :: env_state
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Particle.
    type(aero_particle_t), intent(inout) :: aero_particle

    ! get timestep
    done = .false.
    call find_condense_timestep_variable(dt, env_state, aero_data, &
         aero_particle)
    if (dt .ge. max_dt) then
       dt = max_dt
       done = .true.
    end if

    ! do condensation
    call condense_step_rk(dt, env_state, aero_data, aero_particle)
   
  end subroutine condense_step_rk_fixed

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Does one fixed timestep of Runge-Kutta-4.
  subroutine condense_step_rk(dt, env_state, aero_data, aero_particle)

    !> Timestep.
    real*8, intent(out) :: dt
    !> Environment state.
    type(env_state_t), intent(in) :: env_state
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Particle.
    type(aero_particle_t), intent(inout) :: aero_particle

    ! local variables
    real*8 k1, k2, k3, k4
    type(aero_particle_t) :: aero_particle_tmp

    call aero_particle_alloc(aero_particle_tmp, aero_data%n_spec)
    call aero_particle_copy(aero_particle, aero_particle_tmp)

    ! step 1
    call cond_growth_rate(k1, env_state, aero_data, aero_particle_tmp)

    ! step 2
    aero_particle_tmp%vol(aero_data%i_water) = &
         aero_particle%vol(aero_data%i_water) + dt * k1 / 2d0
    aero_particle_tmp%vol(aero_data%i_water) = &
         max(0d0, aero_particle_tmp%vol(aero_data%i_water))
    call cond_growth_rate(k2, env_state, aero_data, aero_particle_tmp)

    ! step 3
    aero_particle_tmp%vol(aero_data%i_water) = &
         aero_particle%vol(aero_data%i_water) + dt * k2 / 2d0
    aero_particle_tmp%vol(aero_data%i_water) = &
         max(0d0, aero_particle_tmp%vol(aero_data%i_water))
    call cond_growth_rate(k3, env_state, aero_data, aero_particle_tmp)

    ! step 4
    aero_particle_tmp%vol(aero_data%i_water) = &
         aero_particle%vol(aero_data%i_water) + dt * k3
    aero_particle_tmp%vol(aero_data%i_water) = &
         max(0d0, aero_particle_tmp%vol(aero_data%i_water))
    call cond_growth_rate(k4, env_state, aero_data, aero_particle_tmp)

    aero_particle%vol(aero_data%i_water) = &
         aero_particle%vol(aero_data%i_water) &
         + dt * (k1 / 6d0 + k2 / 3d0 + k3 / 3d0 + k4 / 6d0)
    aero_particle%vol(aero_data%i_water) = &
         max(0d0, aero_particle%vol(aero_data%i_water))

    call aero_particle_free(aero_particle_tmp)
   
  end subroutine condense_step_rk

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Just returns a constant timestep.
  subroutine find_condense_timestep_constant(dt, env_state, aero_data, &
       aero_particle)

    !> Timestep to use.
    real*8, intent(out) :: dt
    !> Environment state.
    type(env_state_t), intent(in) :: env_state
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Particle.
    type(aero_particle_t), intent(in) :: aero_particle

    dt = 5d-3

  end subroutine find_condense_timestep_constant

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Computes a timestep proportional to V / (dV/dt).
  subroutine find_condense_timestep_variable(dt, env_state, aero_data, &
       aero_particle)

    !> Timestep to use.
    real*8, intent(out) :: dt
    !> Environment state.
    type(env_state_t), intent(in) :: env_state
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Particle.
    type(aero_particle_t), intent(in) :: aero_particle

    !> Scale factor for timestep.
    real*8, parameter :: scale = 0.1d0

    real*8 pv, dvdt

    pv = aero_particle_volume(aero_particle)
    call cond_growth_rate(dvdt, env_state, aero_data, aero_particle)
    dt = abs(scale * pv / dvdt)

  end subroutine find_condense_timestep_variable
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Find the water volume growth rate due to condensation.
  subroutine cond_growth_rate(dvdt, env_state, aero_data, aero_particle)

    !> Dv/dt (m^3 s^{-1}).
    real*8, intent(out) :: dvdt
    !> Environment state.
    type(env_state_t), intent(in) :: env_state
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Particle.
    type(aero_particle_t), intent(in) :: aero_particle

    !> Relative dm/dt convergence tol.
    real*8, parameter :: dmdt_rel_tol = 1d-8
    !> Function convergence tolerance.
    real*8, parameter :: f_tol = 1d-15
    !> Maximum number of iterations.
    integer, parameter :: iter_max = 100

    real*8 dmdt, pm, dmdt_tol

    pm = aero_particle_mass(aero_particle, aero_data)
    dmdt_tol = pm * dmdt_rel_tol

    dmdt = 0d0
    call cond_newt(dmdt, env_state, aero_data, cond_growth_rate_func, &
         dmdt_tol, f_tol, iter_max, aero_particle)
    
    dvdt = dmdt / aero_data%density(aero_data%i_water)

  end subroutine cond_growth_rate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Scalar Newton's method for solving the implicit condensation
  !> functions.
  subroutine cond_newt(x, env_state, aero_data, func, x_tol, f_tol, &
       iter_max, aero_particle)

    !> Variable (set to init value on call).
    real*8, intent(inout) :: x
    !> Environment state.
    type(env_state_t), intent(in) :: env_state
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> X convergence tolerance.
    real*8, intent(in) :: x_tol
    !> F convergence tolerance.
    real*8, intent(in) :: f_tol
    !> Maximum number of iterations.
    integer, intent(in) :: iter_max
    !> Particle.
    type(aero_particle_t), intent(in) :: aero_particle

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
         real*8, intent(in) :: x
         !> Function to solve.
         real*8, intent(out) :: f
         !> Derivative df/dx.
         real*8, intent(out) :: df
         !> Particle.
         type(aero_particle_t), intent(in) :: aero_particle
       end subroutine func
    end interface
    
    integer iter, k
    real*8 delta_f, delta_x, f, old_f, df

    call func(env_state, aero_data, .true., x, f, df, aero_particle)
    old_f = f

    iter = 0
    do
       iter = iter + 1

       delta_x = f / df
       x = x - delta_x
       call func(env_state, aero_data, .false., x, f, df, aero_particle)
       delta_f = f - old_f
       old_f = f
       
       if (iter .ge. iter_max) then
          write(0,*) 'ERROR: Newton iteration failed to terminate'
          write(0,*) 'iter_max = ', iter_max
          write(0,*) 'x = ', x
          write(0,*) 'aero_particle%vol = ', aero_particle%vol
          call exit(1)
       end if
       
       ! FIXME: gfortran 4.1.1 requires the "then" in the following
       ! statement, rather than using a single-line "if" statement.
       if ((abs(delta_x) .lt. x_tol) &
            .and. (abs(delta_f) .lt. f_tol)) then
          exit
       end if
    end do

  end subroutine cond_newt

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Return the error function value and its derivative for the
  !> implicit growth rate function.
  subroutine cond_growth_rate_func(env_state, aero_data, init, dmdt, &
       f, df, aero_particle)

    !> Environment state.
    type(env_state_t), intent(in) :: env_state
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> True if first Newton loop.
    logical, intent(in) :: init
    !> Mass growth rate dm/dt (kg s^{-1}).
    real*8, intent(in) :: dmdt
    !> Error.
    real*8, intent(out) :: f
    !> Derivative of error with respect to x.
    real*8, intent(out) :: df
    !> Particle.
    type(aero_particle_t), intent(in) :: aero_particle
    
    ! local variables
    real*8, save :: k_a, k_ap, k_ap_div, D_v, D_v_div, D_vp, d_p, pv
    real*8, save :: rat, fact1, fact2, c1, c2, c3, c4, c5
    real*8, save :: M_water, M_solute, rho_water, rho_solute
    real*8, save :: eps, nu, g_water, g_solute

    real*8 T_a ! droplet temperature (K), determined as part of solve

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
       d_p = vol2diam(pv) ! (m)
       
       ! molecular diffusion coefficient uncorrected
       D_v = 0.211d-4 / (env_state%pressure / const%air_std_press) &
            * (env_state%temp / 273d0)**1.94d0 ! (m^2 s^{-1})
       ! molecular diffusion coefficient corrected for non-continuum effects
       D_v_div = 1d0 + (2d0 * D_v * 1d-4 / (const%accom_coeff * d_p)) &
            * (2d0 * const%pi * M_water &
            / (const%univ_gas_const * env_state%temp))**0.5d0
       D_vp = D_v / D_v_div
       
       ! TEST: use the basic expression for D_vp
       ! D_vp = D_v                ! (m^2 s^{-1})
       ! FIXME: check whether we can reinstate the correction
       
       ! thermal conductivity uncorrected
       k_a = 1d-3 * (4.39d0 + 0.071d0 &
            * env_state%temp) ! (J m^{-1} s^{-1} K^{-1})
       k_ap_div = 1d0 + 2d0 &  ! dimensionless
            * k_a / (const%accom_coeff * d_p * env_state_air_den(env_state) &
            * const%water_spec_heat) &
            * (2d0 * const%pi * const%air_molec_weight &
            / (const%univ_gas_const * env_state%temp))**0.5d0
       ! thermal conductivity corrected
       k_ap = k_a / k_ap_div     ! (J m^{-1} s^{-1} K^{-1})

       rat = env_state_sat_vapor_pressure(env_state) &
            / (const%univ_gas_const * env_state%temp)
       fact1 = const%water_latent_heat * M_water &
            / (const%univ_gas_const * env_state%temp)
       fact2 = const%water_latent_heat &
            / (2d0 * const%pi * d_p * k_ap * env_state%temp)

       c1 = 2d0 * const%pi * d_p * D_vp * M_water * rat
       c2 = 4d0 * M_water &
            * const%water_surf_eng / (const%univ_gas_const * rho_water * d_p)
       c3 = c1 * fact1 * fact2
       c4 = const%water_latent_heat / (2d0 * const%pi * d_p * k_ap)
       ! incorrect expression from Majeed and Wexler:
       !     c5 = nu * eps * M_water * rho_solute * r_n**3d0 &
       !         / (M_solute * rho_water * ((d_p / 2)**3d0 - r_n**3))
       ! c5 = nu * eps * M_water / M_solute * g_solute / g_water
       ! corrected according to Jim's note:
           c5 = nu * eps * M_water / M_solute * g_solute / &
                (g_water + (rho_water / rho_solute) * eps * g_solute)
    end if

    T_a = env_state%temp + c4 * dmdt ! (K)
    
    f = dmdt - c1 * (env_state%rel_humid - exp(c2 / T_a - c5)) &
         / (1d0 + c3 * exp(c2 / T_a - c5))
    
    df = 1d0 + c1 * env_state%rel_humid &
         * (1d0 + c3 * exp(c2 / T_a -c5))**(-2d0) * c3 * &
         exp(c2 / T_a - c5) * (-1d0) * c2 * c4 / T_a**2d0 + c1 * &
         (exp(c2 / T_a - c5) * (-1d0) * c2 * c4 / T_a**2d0 * (1d0 + c3 &
         * exp(c2 / T_a -c5))**(-1d0) + exp(c2 / T_a - c5) * (-1d0) * &
         (1d0 + c3 * exp(c2 / T_a -c5))**(-2d0) * c3 * exp(c2 / T_a - &
         c5) * (-1d0) * c2 * c4 / T_a**2d0)

    ! NOTE: we could return T_a (the droplet temperature) if we have
    ! any need for it.

  end subroutine cond_growth_rate_func
  
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
    !> Dw relative convergence tolerance.
    real*8, parameter :: pv_rel_tol = 1d-6
    !> Function convergence tolerance.
    real*8, parameter :: f_tol = 1d-15
    !> Maximum number of iterations.
    integer, parameter :: iter_max = 100
    !> Initial value.
    real*8, parameter :: dw_init = 1d0
    
    real*8 dw ! wet diameter of particle
    real*8 dw_tol, pv

    pv = aero_particle_volume(aero_particle)
    dw = dw_init
    ! convert volume relative tolerance to diameter absolute tolerance
    dw_tol = vol2diam(pv * (1d0 + pv_rel_tol)) - vol2diam(pv)
    call cond_newt(dw, env_state, aero_data, equilibriate_func, &
         dw_tol, f_tol, iter_max, aero_particle)

    aero_particle%vol(aero_data%i_water) = diam2vol(dw) - pv

  end subroutine equilibriate_particle

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
    real*8, intent(in) :: dw
    !> Function value.
    real*8, intent(out) :: f
    !> Function derivative df/dx.
    real*8, intent(out) :: df
    !> Particle.
    type(aero_particle_t), intent(inout) :: aero_particle

    real*8, save :: c0, c1, c3, c4, dc0, dc2, dc3
    real*8, save :: A, B
    real*8, save ::  pv
    real*8, save :: M_water, M_solute, rho_water, rho_solute
    real*8, save :: eps, nu, g_water, g_solute

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
    end if
    
    f = c4 * dw**4d0 - c3 * dw**3d0 + c1 * dw + c0
    df = dc3 * dw**3d0 -dc2 * dw**2d0 + dc0

  end subroutine equilibriate_func

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
          call equilibriate_particle(env_state, aero_data, &
               aero_state%bin(i_bin)%particle(i))
       end do
    end do

  end subroutine aero_state_equilibriate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
end module pmc_condensation
