! Copyright (C) 2005-2007 Nicole Riemer and Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.
!
! Condensation routines for water condensing onto particles.

module mod_condensation
contains
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine condense_particles(bin_grid, aero_binned, env, aero_data, &
       aero_state, del_t)

    ! Do condensation to all the particles for a given time interval,
    ! including updating the environment to account for the lost
    ! vapor.

    use mod_aero_state
    use mod_bin_grid
    use mod_aero_binned
    use mod_environ
    use mod_aero_data

    type(bin_grid_t), intent(in) :: bin_grid ! bin grid
    type(aero_binned_t), intent(inout) :: aero_binned ! binned distributions
    type(environ), intent(inout) :: env ! environment state
    type(aero_data_t), intent(in) :: aero_data ! aerosol data
    type(aero_state_t), intent(inout) :: aero_state ! aerosol state
    real*8, intent(in) :: del_t         ! total time to integrate
    
    ! local variables
    integer bin, j, new_bin, k
    real*8 pv, pre_water_vol, post_water_vol

    ! FIXME: don't rely on binned data, but rather compute the total
    ! water transfered in condense_particle() and return it
    pre_water_vol = sum(aero_binned%vol_den(:,aero_data%i_water)) &
         * aero_state%comp_vol * bin_grid%dlnr

    do bin = 1,bin_grid%n_bin
       do j = 1,aero_state%bins(bin)%n_part
          call condense_particle(del_t, env, aero_data, &
               aero_state%bins(bin)%particles(j))
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
    call change_water_volume(env, aero_data, &
         (post_water_vol - pre_water_vol) / aero_state%comp_vol)

  end subroutine condense_particles

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine condense_particle(del_t, env, aero_data, aero_particle)

    ! Integrate the condensation growth or decay ODE for total time
    ! del_t for a single particle.

    use mod_util
    use mod_environ
    use mod_aero_data
    use mod_aero_particle

    real*8, intent(in) :: del_t         ! total time to integrate
    type(environ), intent(in) :: env    ! environment state
    type(aero_data_t), intent(in) :: aero_data ! aerosol data
    type(aero_particle_t), intent(inout) :: aero_particle ! particle

    real*8 time_step, time
    logical done

    integer i
    real*8 dvdt

    time = 0d0
    done = .false.
    do while (.not. done)
       call condense_step_euler(del_t - time, time_step, done, env, &
            aero_data, aero_particle)
       time = time + time_step
    end do
    
  end subroutine condense_particle

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine condense_step_euler(max_dt, dt, done, env, aero_data, &
       aero_particle)

    ! Does one timestep (determined by this subroutine) of the
    ! condensation ODE. The timestep will not exceed max_dt, but might
    ! be less. If we in fact step all the way to max_dt then done will
    ! be true. This uses the explicit (forward) Euler integrator.
    
    use mod_environ
    use mod_aero_data
    use mod_aero_particle

    real*8, intent(in) :: max_dt        ! maximum timestep to integrate
    real*8, intent(out) :: dt           ! actual timestep used
    logical, intent(out) :: done        ! did we reach the maximum timestep?
    type(environ), intent(in) :: env    ! environment state
    type(aero_data_t), intent(in) :: aero_data   ! aerosol data
    type(aero_particle_t), intent(inout) :: aero_particle ! particle

    real*8 dvdt

    ! get timestep
    done = .false.
    call find_condense_timestep_variable(dt, env, aero_data, aero_particle)
    if (dt .ge. max_dt) then
       dt = max_dt
       done = .true.
    end if

    ! do condensation
    call cond_growth_rate(dvdt, env, aero_data, aero_particle)
    aero_particle%vols(aero_data%i_water) = &
         aero_particle%vols(aero_data%i_water) + dt * dvdt

    ! ensure volumes stay positive
    aero_particle%vols(aero_data%i_water) = max(0d0, &
         aero_particle%vols(aero_data%i_water))
   
  end subroutine condense_step_euler

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine condense_step_rk_fixed(max_dt, &
       dt, done, env, aero_data, aero_particle)

    ! Does one timestep (determined by this subroutine) of the
    ! condensation ODE. The timestep will not exceed max_dt, but might
    ! be less. If we in fact step all the way to max_dt then done will
    ! be true. This uses the explicit 4th-order Runge-Kutta integrator.
    
    use mod_environ
    use mod_aero_data
    use mod_aero_particle

    real*8, intent(in) :: max_dt        ! maximum timestep to integrate
    real*8, intent(out) :: dt           ! actual timestep used
    logical, intent(out) :: done        ! did we reach the maximum timestep?
    type(environ), intent(in) :: env    ! environment state
    type(aero_data_t), intent(in) :: aero_data   ! aerosol data
    type(aero_particle_t), intent(inout) :: aero_particle ! particle

    ! get timestep
    done = .false.
    call find_condense_timestep_variable(dt, env, aero_data, aero_particle)
    if (dt .ge. max_dt) then
       dt = max_dt
       done = .true.
    end if

    ! do condensation
    call condense_step_rk(dt, env, aero_data, aero_particle)
   
  end subroutine condense_step_rk_fixed

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine condense_step_rk(dt, env, aero_data, aero_particle)

    ! Does one fixed timestep of Runge-Kutta-4.

    use mod_environ
    use mod_aero_data
    use mod_aero_particle

    real*8, intent(out) :: dt           ! timestep
    type(environ), intent(in) :: env    ! environment state
    type(aero_data_t), intent(in) :: aero_data   ! aerosol data
    type(aero_particle_t), intent(inout) :: aero_particle ! particle

    ! local variables
    real*8 k1, k2, k3, k4
    type(aero_particle_t) :: aero_particle_tmp

    call aero_particle_alloc(aero_data%n_spec, aero_particle_tmp)
    call aero_particle_copy(aero_particle, aero_particle_tmp)

    ! step 1
    call cond_growth_rate(k1, env, aero_data, aero_particle_tmp)

    ! step 2
    aero_particle_tmp%vols(aero_data%i_water) = &
         aero_particle%vols(aero_data%i_water) + dt * k1 / 2d0
    aero_particle_tmp%vols(aero_data%i_water) = &
         max(0d0, aero_particle_tmp%vols(aero_data%i_water))
    call cond_growth_rate(k2, env, aero_data, aero_particle_tmp)

    ! step 3
    aero_particle_tmp%vols(aero_data%i_water) = &
         aero_particle%vols(aero_data%i_water) + dt * k2 / 2d0
    aero_particle_tmp%vols(aero_data%i_water) = &
         max(0d0, aero_particle_tmp%vols(aero_data%i_water))
    call cond_growth_rate(k3, env, aero_data, aero_particle_tmp)

    ! step 4
    aero_particle_tmp%vols(aero_data%i_water) = &
         aero_particle%vols(aero_data%i_water) + dt * k3
    aero_particle_tmp%vols(aero_data%i_water) = &
         max(0d0, aero_particle_tmp%vols(aero_data%i_water))
    call cond_growth_rate(k4, env, aero_data, aero_particle_tmp)

    aero_particle%vols(aero_data%i_water) = &
         aero_particle%vols(aero_data%i_water) &
         + dt * (k1 / 6d0 + k2 / 3d0 + k3 / 3d0 + k4 / 6d0)
    aero_particle%vols(aero_data%i_water) = &
         max(0d0, aero_particle%vols(aero_data%i_water))

    call aero_particle_free(aero_particle_tmp)
   
  end subroutine condense_step_rk

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine find_condense_timestep_constant(dt, env, aero_data, &
       aero_particle)

    ! Just returns a constant timestep.

    use mod_aero_state
    use mod_environ
    use mod_aero_data
    use mod_aero_particle

    real*8, intent(out) :: dt           ! timestep to use
    type(environ), intent(in) :: env    ! environment state
    type(aero_data_t), intent(in) :: aero_data   ! aerosol data
    type(aero_particle_t), intent(in) :: aero_particle ! particle

    dt = 5d-3

  end subroutine find_condense_timestep_constant

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine find_condense_timestep_variable(dt, env, aero_data, &
       aero_particle)

    ! Computes a timestep proportional to V / (dV/dt).

    use mod_aero_state
    use mod_environ
    use mod_aero_data
    use mod_aero_particle

    real*8, intent(out) :: dt           ! timestep to use
    type(environ), intent(in) :: env    ! environment state
    type(aero_data_t), intent(in) :: aero_data   ! aerosol data
    type(aero_particle_t), intent(in) :: aero_particle ! particle

    real*8, parameter :: scale = 0.1d0  ! scale factor for timestep

    real*8 pv, dvdt

    pv = aero_particle_volume(aero_particle)
    call cond_growth_rate(dvdt, env, aero_data, aero_particle)
    dt = abs(scale * pv / dvdt)

  end subroutine find_condense_timestep_variable
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
  subroutine cond_growth_rate(dvdt, env, aero_data, aero_particle)
    
    ! Find the water volume growth rate due to condensation.

    use mod_environ
    use mod_aero_data
    use mod_aero_particle

    real*8, intent(out) :: dvdt         ! dv/dt (m^3 s^{-1})
    type(environ), intent(in) :: env    ! environment state
    type(aero_data_t), intent(in) :: aero_data   ! aerosol data
    type(aero_particle_t), intent(in) :: aero_particle ! particle

    real*8, parameter :: dmdt_rel_tol = 1d-8 ! relative dm/dt convergence tol
    real*8, parameter :: f_tol = 1d-15  ! function convergence tolerance
    integer, parameter :: iter_max = 100 ! maximum number of iterations

    real*8 dmdt, pm, dmdt_tol

    pm = aero_particle_mass(aero_particle, aero_data)
    dmdt_tol = pm * dmdt_rel_tol

    dmdt = 0d0
    call cond_newt(dmdt, env, aero_data, cond_growth_rate_func, &
         dmdt_tol, f_tol, iter_max, aero_particle)
    
    dvdt = dmdt / aero_data%rho(aero_data%i_water)

  end subroutine cond_growth_rate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
  subroutine cond_newt(x, env, aero_data, func, x_tol, f_tol, &
       iter_max, aero_particle)
    
    ! Scalar Newton's method for solving the implicit condensation
    ! functions.

    use mod_environ
    use mod_aero_data
    use mod_aero_particle

    real*8, intent(inout) :: x          ! variable (set to inital value on call)
    type(environ), intent(in) :: env    ! environment state
    type(aero_data_t), intent(in) :: aero_data   ! aerosol data
    real*8, intent(in) :: x_tol         ! x convergence tolerance
    real*8, intent(in) :: f_tol         ! f convergence tolerance
    integer, intent(in) :: iter_max     ! maximum number of iterations
    type(aero_particle_t), intent(in) :: aero_particle ! particle

    interface
       subroutine func(env, aero_data, init, x, f, df, aero_particle)
         use mod_environ
         use mod_aero_data
         use mod_aero_particle
         type(environ), intent(in) :: env  ! environment state
         type(aero_data_t), intent(in) :: aero_data ! aerosol data
         logical, intent(in) :: init       ! true if first Newton loop
         real*8, intent(in) :: x           ! independent variable to solve for
         real*8, intent(out) :: f          ! function to solve
         real*8, intent(out) :: df         ! derivative df/dx
         type(aero_particle_t), intent(in) :: aero_particle ! particle
       end subroutine func
    end interface
    
    integer iter, k
    real*8 delta_f, delta_x, f, old_f, df

    call func(env, aero_data, .true., x, f, df, aero_particle)
    old_f = f

    iter = 0
    do
       iter = iter + 1

       delta_x = f / df
       x = x - delta_x
       call func(env, aero_data, .false., x, f, df, aero_particle)
       delta_f = f - old_f
       old_f = f
       
       if (iter .ge. iter_max) then
          write(0,*) 'ERROR: Newton iteration failed to terminate'
          write(0,*) 'iter_max = ', iter_max
          write(0,*) 'x = ', x
          write(0,*) 'aero_particle%vols = ', aero_particle%vols
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

  subroutine cond_growth_rate_func(env, aero_data, init, dmdt, &
       f, df, aero_particle)

    ! Return the error function value and its derivative for the
    ! implicit growth rate function.

    use mod_aero_state
    use mod_util
    use mod_environ
    use mod_aero_data
    use mod_constants
    use mod_aero_particle

    type(environ), intent(in) :: env    ! environment state
    type(aero_data_t), intent(in) :: aero_data   ! aerosol data
    logical, intent(in) :: init         ! true if first Newton loop
    real*8, intent(in) :: dmdt          ! mass growth rate dm/dt (kg s^{-1})
    real*8, intent(out) :: f            ! error
    real*8, intent(out) :: df           ! derivative of error with respect to x
    type(aero_particle_t), intent(in) :: aero_particle ! particle
    
    ! local variables
    real*8, save :: k_a, k_ap, k_ap_div, D_v, D_v_div, D_vp, d_p, pv
    real*8, save :: rat, fact1, fact2, c1, c2, c3, c4, c5
    real*8, save :: M_water, M_solute, rho_water, rho_solute
    real*8, save :: eps, nu, g_water, g_solute

    real*8 T_a ! droplet temperature (K), determined as part of solve

    if (init) then
       ! Start of new Newton loop, compute all constants

       M_water = aero_particle_water_M_w(aero_data)                  ! (kg/mole)
       M_solute = aero_particle_solute_M_w(aero_data, aero_particle) ! (kg/mole)
       nu = aero_particle_solute_nu(aero_data, aero_particle)        ! (1)
       eps = aero_particle_solute_eps(aero_data, aero_particle)      ! (1)
       rho_water = aero_particle_water_rho(aero_data)                ! (kg/m^3)
       rho_solute = aero_particle_solute_rho(aero_data, aero_particle)! (kg/m^3)
       g_water = aero_particle_water_mass(aero_data, aero_particle)  ! (kg)
       g_solute = aero_particle_solute_mass(aero_data, aero_particle) ! (kg)

       pv = aero_particle_volume(aero_particle) ! m^3
       d_p = vol2diam(pv) ! m
       
       ! molecular diffusion coefficient uncorrected
       D_v = 0.211d-4 / (env%p / const%atm) &
            * (env%T / 273d0)**1.94d0 ! m^2 s^{-1}
       ! molecular diffusion coefficient corrected for non-continuum effects
       D_v_div = 1d0 + (2d0 * D_v * 1d-4 / (const%alpha * d_p)) &
            * (2d0 * const%pi * M_water / (const%R * env%T))**0.5d0
       D_vp = D_v / D_v_div
       
       ! TEST: use the basic expression for D_vp
       ! D_vp = D_v                ! m^2 s^{-1}
       ! FIXME: check whether we can reinstate the correction
       
       ! thermal conductivity uncorrected
       k_a = 1d-3 * (4.39d0 + 0.071d0 * env%T) ! J m^{-1} s^{-1} K^{-1}
       k_ap_div = 1d0 + 2d0 &
            * k_a / (const%alpha * d_p * env%rho_a * const%cp) &
            * (2d0 * const%pi * const%M_a / (const%R * env%T))**0.5d0 ! dim-less
       ! thermal conductivity corrected
       k_ap = k_a / k_ap_div     ! J m^{-1} s^{-1} K^{-1}

       rat = sat_vapor_pressure(env) / (const%R * env%T)
       fact1 = const%L_v * M_water / (const%R * env%T)
       fact2 = const%L_v / (2d0 * const%pi * d_p * k_ap * env%T)

       c1 = 2d0 * const%pi * d_p * D_vp * M_water * rat
       c2 = 4d0 * M_water &
            * const%sig / (const%R * rho_water * d_p)
       c3 = c1 * fact1 * fact2
       c4 = const%L_v / (2d0 * const%pi * d_p * k_ap)
       ! incorrect expression from Majeed and Wexler:
       !     c5 = nu * eps * M_water * rho_solute * r_n**3d0 &
       !         / (M_solute * rho_water * ((d_p / 2)**3d0 - r_n**3))
       ! c5 = nu * eps * M_water / M_solute * g_solute / g_water
       ! corrected according to Jim's note:
           c5 = nu * eps * M_water / M_solute * g_solute / &
                (g_water + (rho_water / rho_solute) * eps * g_solute)
    end if

    T_a = env%T + c4 * dmdt ! K
    
    f = dmdt - c1 * (env%RH - exp(c2 / T_a - c5)) &
         / (1d0 + c3 * exp(c2 / T_a - c5))
    
    df = 1d0 + c1 * env%RH * (1d0 + c3 * exp(c2 / T_a -c5))**(-2d0) * c3 * &
         exp(c2 / T_a - c5) * (-1d0) * c2 * c4 / T_a**2d0 + c1 * &
         (exp(c2 / T_a - c5) * (-1d0) * c2 * c4 / T_a**2d0 * (1d0 + c3 &
         * exp(c2 / T_a -c5))**(-1d0) + exp(c2 / T_a - c5) * (-1d0) * &
         (1d0 + c3 * exp(c2 / T_a -c5))**(-2d0) * c3 * exp(c2 / T_a - &
         c5) * (-1d0) * c2 * c4 / T_a**2d0)

    ! NOTE: we could return T_a (the droplet temperature) if we have
    ! any need for it.

  end subroutine cond_growth_rate_func
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine equilibriate_particle(env, aero_data, aero_particle)

    ! Add water to the particle until it is in equilibrium.

    use mod_util
    use mod_aero_state
    use mod_environ
    use mod_aero_data
    use mod_constants
    use mod_aero_particle

    type(environ), intent(in) :: env    ! environment state
    type(aero_data_t), intent(in) :: aero_data   ! aerosol data
    type(aero_particle_t), intent(inout) :: aero_particle ! particle

    ! parameters
    integer, parameter :: it_max = 400     ! maximum iterations
    real*8, parameter :: pv_rel_tol = 1d-6 ! dw relative convergence tolerance
    real*8, parameter :: f_tol = 1d-15     ! function convergence tolerance
    integer, parameter :: iter_max = 100   ! maximum number of iterations
    real*8, parameter :: dw_init = 1d0     ! initial value
    
    real*8 dw ! wet diameter of particle
    real*8 dw_tol, pv

    pv = aero_particle_volume(aero_particle)
    dw = dw_init
    ! convert volume relative tolerance to diameter absolute tolerance
    dw_tol = vol2diam(pv * (1d0 + pv_rel_tol)) - vol2diam(pv)
    call cond_newt(dw, env, aero_data, equilibriate_func, &
         dw_tol, f_tol, iter_max, aero_particle)

    aero_particle%vols(aero_data%i_water) = diam2vol(dw) - pv

  end subroutine equilibriate_particle

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine equilibriate_func(env, aero_data, init, dw, f, df, &
       aero_particle)

    ! Return the error function value and its derivative for the
    ! implicit function that determines the equilibrium state of a
    ! particle.

    use mod_util
    use mod_aero_state
    use mod_environ
    use mod_aero_data
    use mod_constants
    use mod_aero_particle

    type(environ), intent(in) :: env    ! environment state
    type(aero_data_t), intent(in) :: aero_data   ! aerosol data
    logical, intent(in) :: init         ! true if first Newton loop
    real*8, intent(in) :: dw            ! wet diameter (m)
    real*8, intent(out) :: f            ! function value
    real*8, intent(out) :: df           ! function derivative df/dx
    type(aero_particle_t), intent(inout) :: aero_particle ! particle

    real*8, save :: c0, c1, c3, c4, dc0, dc2, dc3
    real*8, save :: A, B
    real*8, save ::  pv
    real*8, save :: M_water, M_solute, rho_water, rho_solute
    real*8, save :: eps, nu, g_water, g_solute

    if (init) then
       ! Start of new Newton loop, compute all constants

       M_water = aero_particle_water_M_w(aero_data)                  ! (kg/mole)
       M_solute = aero_particle_solute_M_w(aero_data, aero_particle) ! (kg/mole)
       nu = aero_particle_solute_nu(aero_data, aero_particle)        ! (1)
       eps = aero_particle_solute_eps(aero_data, aero_particle)      ! (1)
       rho_water = aero_particle_water_rho(aero_data)                ! (kg/m^3)
       rho_solute = aero_particle_solute_rho(aero_data, aero_particle)! (kg/m^3)
       g_water = aero_particle_water_mass(aero_data, aero_particle)  ! (kg)
       g_solute = aero_particle_solute_mass(aero_data, aero_particle) ! (kg)

       pv = aero_particle_volume(aero_particle)

       A = 4d0 * M_water * const%sig / (const%R * env%T * rho_water)
       
       B = nu * eps * M_water * rho_solute &
            * vol2rad(pv)**3.d0 / (M_solute * rho_water)
       
       c4 = log(env%RH) / 8d0
       c3 = A / 8d0
       
       dc3 = log(env%RH) / 2d0
       dc2 = 3d0 * A / 8d0
       
       c1 = B - log(env%RH) * vol2rad(pv)**3d0
       c0 = A * vol2rad(pv)**3d0
       dc0 = c1
    end if
    
    f = c4 * dw**4d0 - c3 * dw**3d0 + c1 * dw + c0
    df = dc3 * dw**3d0 -dc2 * dw**2d0 + dc0

  end subroutine equilibriate_func

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine aero_state_equilibriate(bin_grid, env, aero_data, aero_state)
    
    ! call equilibriate_particle() on each particle in the aerosol

    use mod_bin_grid
    use mod_environ
    use mod_aero_data
    use mod_aero_state
    use mod_aero_particle

    type(bin_grid_t), intent(in) :: bin_grid ! bin grid
    type(environ), intent(inout) :: env ! environment state
    type(aero_data_t), intent(in) :: aero_data ! aerosol data
    type(aero_state_t), intent(inout) :: aero_state ! aerosol state

    integer :: i_bin, i
    
    do i_bin = 1,bin_grid%n_bin
       do i = 1,aero_state%bins(i_bin)%n_part
          call equilibriate_particle(env, aero_data, &
               aero_state%bins(i_bin)%particles(i))
       end do
    end do

  end subroutine aero_state_equilibriate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
end module mod_condensation
