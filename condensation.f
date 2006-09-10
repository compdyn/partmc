! -*- mode: f90;-*-
! Condensation
!

module mod_condensation
contains
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine condense_particles(n_bin, TDV, n_spec, MH, VH, &
       del_t, bin_v, bin_r, bin_g, bin_gs, bin_n, dlnr, env, mat)

    use mod_array
    use mod_array_hybrid
    use mod_bin
    use mod_environ
    use mod_material

    integer, intent(in) :: n_bin ! number of bins
    integer, intent(in) :: TDV   ! second dimension of VH
    integer, intent(in) :: n_spec       ! number of species
    integer, intent(inout) :: MH(n_bin) ! number of particles per bin
    real*8, intent(inout) :: VH(n_bin,TDV,n_spec) ! particle volumes (m^3)
    real*8, intent(in) :: del_t         ! total time to integrate
    real*8, intent(in) :: bin_v(n_bin) ! volume of particles in bins (m^3)
    real*8, intent(in) ::  bin_r(n_bin) ! radius of particles in bins (m)
    real*8, intent(inout) :: bin_g(n_bin) ! mass in bins  
    real*8, intent(inout) :: bin_gs(n_bin,n_spec) ! species mass in bins
    integer, intent(inout) :: bin_n(n_bin)      ! number in bins
    real*8, intent(in) :: dlnr                  ! bin scale factor
    type(environ), intent(inout) :: env  ! environment state
    type(material), intent(in) :: mat    ! material properties
    
    ! local variables
    integer bin, j, new_bin, k
    real*8 pv, pre_water_vol, post_water_vol

    pre_water_vol = sum(bin_gs(:,mat%i_water))

    do bin = 1,n_bin
       do j = 1,MH(bin)
          call condense_particle(n_spec, VH(bin,j,:), del_t, env, mat)
       end do
    end do

    ! We resort the particles in the bins after all particles are
    ! advanced, otherwise we will lose track of which ones have been
    ! advanced and which have not.
    call resort_array_hybrid(n_bin, TDV, n_spec, MH, VH, bin_v, &
         bin_r, dlnr)

    ! update the bin arrays
    call moments_hybrid(n_bin, TDV, n_spec, MH, VH, bin_v, &
         bin_r, bin_g, bin_gs, bin_n, dlnr)

    ! update the environment due to condensation of water
    post_water_vol = sum(bin_gs(:,mat%i_water))
    call change_water_volume(env, mat, post_water_vol - pre_water_vol)

  end subroutine condense_particles

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine condense_particle(n_spec, V, del_t, env, mat)

    ! integrates the condensation growth or decay ODE for total time
    ! del_t

    use mod_util
    use mod_environ
    use mod_material

    integer, intent(in) :: n_spec ! number of species
    real*8, intent(inout) :: V(n_spec) ! particle volumes (m^3)
    real*8, intent(in) :: del_t ! total time to integrate
    type(environ), intent(in) :: env     ! environment state
    type(material), intent(in) :: mat    ! material properties

    real*8 time_step, time
    logical done

    integer i
    real*8 dvdt

    time = 0d0
    done = .false.
    do while (.not. done)
       call condense_step_euler(n_spec, V, del_t - time, &
            time_step, done, env, mat)
       time = time + time_step
    end do
    
  end subroutine condense_particle

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine condense_step_euler(n_spec, V, max_dt, dt, &
       done, env, mat)

    ! Does one timestep (determined by this subroutine) of the
    ! condensation ODE. The timestep will not exceed max_dt, but might
    ! be less. If we in fact step all the way to max_dt then done will
    ! be true.
    
    use mod_environ
    use mod_material

    integer, intent(in) :: n_spec ! number of species
    real*8, intent(inout) :: V(n_spec) ! particle volumes (m^3)
    real*8, intent(in) :: max_dt ! maximum timestep to integrate
    real*8, intent(out) :: dt ! actual timestep used
    logical, intent(out) :: done ! whether we reached the maximum timestep
    type(environ), intent(in) :: env     ! environment state
    type(material), intent(in) :: mat    ! material properties

    real*8 dvdt

    done = .false.
    call find_condense_timestep_variable(n_spec, V, dt, env, mat)
    if (dt .ge. max_dt) then
       dt = max_dt
       done = .true.
    end if

    call cond_growth_rate(n_spec, V, dvdt, env, mat)
    V(mat%i_water) = V(mat%i_water) + dt * dvdt
    V(mat%i_water) = max(0d0, V(mat%i_water))
   
  end subroutine condense_step_euler

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine condense_step_rk_fixed(n_spec, V, max_dt, &
       dt, done, env, mat)

    ! Does one timestep (determined by this subroutine) of the
    ! condensation ODE. The timestep will not exceed max_dt, but might
    ! be less. If we in fact step all the way to max_dt then done will
    ! be true.
    
    use mod_environ
    use mod_material

    integer, intent(in) :: n_spec ! number of species
    real*8, intent(inout) :: V(n_spec) ! particle volumes (m^3)
    real*8, intent(in) :: max_dt ! maximum timestep to integrate
    real*8, intent(out) :: dt ! actual timestep used
    logical, intent(out) :: done ! whether we reached the maximum timestep
    type(environ), intent(in) :: env     ! environment state
    type(material), intent(in) :: mat    ! material properties

    done = .false.
    call find_condense_timestep_variable(n_spec, V, dt, env, mat)
    if (dt .ge. max_dt) then
       dt = max_dt
       done = .true.
    end if

    call condense_step_rk(n_spec, V, dt, env, mat)
   
  end subroutine condense_step_rk_fixed

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine condense_step_rk(n_spec, V, dt, env, mat)

    ! Does one fixed timestep of RK4.

    use mod_environ
    use mod_material

    integer, intent(in) :: n_spec ! number of species
    real*8, intent(inout) :: V(n_spec) ! particle volumes (m^3)
    real*8, intent(out) :: dt ! timestep
    type(environ), intent(in) :: env     ! environment state
    type(material), intent(in) :: mat    ! material properties

    ! local variables
    real*8 k1, k2, k3, k4
    real*8 V_tmp(n_spec)

    V_tmp = V

    ! step 1
    call cond_growth_rate(n_spec, V, k1, env, mat)

    ! step 2
    V_tmp(mat%i_water) = V(mat%i_water) + dt * k1 / 2d0
    V_tmp(mat%i_water) = max(0d0, V_tmp(mat%i_water))
    call cond_growth_rate(n_spec, V_tmp, k2, env, mat)

    ! step 3
    V_tmp(mat%i_water) = V(mat%i_water) + dt * k2 / 2d0
    V_tmp(mat%i_water) = max(0d0, V_tmp(mat%i_water))
    call cond_growth_rate(n_spec, V_tmp, k3, env, mat)

    ! step 4
    V_tmp(mat%i_water) = V(mat%i_water) + dt * k3
    V_tmp(mat%i_water) = max(0d0, V_tmp(mat%i_water))
    call cond_growth_rate(n_spec, V_tmp, k4, env, mat)

    V(mat%i_water) = V(mat%i_water) &
         + dt * (k1 / 6d0 + k2 / 3d0 + k3 / 3d0 + k4 / 6d0)

    V(mat%i_water) = max(0d0, V(mat%i_water))
   
  end subroutine condense_step_rk

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine find_condense_timestep_constant(n_spec, V, dt, env, mat)

    ! constant timestep

    use mod_array
    use mod_environ
    use mod_material

    integer, intent(in) :: n_spec ! number of species
    real*8, intent(in) :: V(n_spec) ! particle volumes (m^3)
    real*8, intent(out) :: dt ! timestep to use
    type(environ), intent(in) :: env     ! environment state
    type(material), intent(in) :: mat    ! material properties

    dt = 5d-3

  end subroutine find_condense_timestep_constant

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine find_condense_timestep_variable(n_spec, V, dt, env, mat)

    ! timestep is proportional to V / (dV/dt)

    use mod_array
    use mod_environ
    use mod_material

    integer, intent(in) :: n_spec ! number of species
    real*8, intent(in) :: V(n_spec) ! particle volumes (m^3)
    real*8, intent(out) :: dt ! timestep to use
    type(environ), intent(in) :: env     ! environment state
    type(material), intent(in) :: mat    ! material properties

    ! parameters
    real*8 scale
    parameter (scale = 0.1d0) ! scale factor for timestep

    real*8 pv, dvdt

    call particle_vol_base(n_spec, V, pv)
    call cond_growth_rate(n_spec, V, dvdt, env, mat)
    dt = abs(scale * pv / dvdt)

  end subroutine find_condense_timestep_variable
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
  subroutine cond_growth_rate(n_spec, V, dvdt, env, mat)
    
    ! Find the water volume growth rate due to condensation.

    use mod_environ
    use mod_material

    integer, intent(in) :: n_spec     ! number of species
    real*8, intent(in) :: V(n_spec)   ! particle volumes (m^3)
    real*8, intent(out) :: dvdt       ! dv/dt (m^3 s^{-1})
    type(environ), intent(in) :: env  ! environment state
    type(material), intent(in) :: mat ! material properties

    real*8 :: dmdt_rel_tol = 1d-8        ! relative dm/dt convergence tolerance
    real*8, parameter :: f_tol = 1d-15   ! function convergence tolerance
    integer, parameter :: iter_max = 100 ! maximum number of iterations

    real*8 dmdt, pm, dmdt_tol

    pm = particle_mass(V, mat)
    dmdt_tol = pm * dmdt_rel_tol

    dmdt = 0d0
    call cond_newt(n_spec, V, dmdt, env, mat, cond_growth_rate_func, &
         dmdt_tol, f_tol, iter_max)

    dvdt = dmdt / mat%rho(mat%i_water)

  end subroutine cond_growth_rate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
  subroutine cond_newt(n_spec, V, x, env, mat, func, x_tol, f_tol, iter_max)
    
    ! Scalar Newton's method for condensation functions.

    use mod_environ
    use mod_material

    integer, intent(in) :: n_spec     ! number of species
    real*8, intent(in) :: V(n_spec)   ! particle volumes (m^3)
    real*8, intent(inout) :: x        ! variable (set to inital value on call)
    type(environ), intent(in) :: env  ! environment state
    type(material), intent(in) :: mat ! material properties
    real*8, intent(in) :: x_tol       ! x convergence tolerance
    real*8, intent(in) :: f_tol       ! f convergence tolerance
    integer, intent(in) :: iter_max   ! maximum number of iterations

    interface
       subroutine func(n_spec, V, env, mat, init, x, f, df)
         integer, intent(in) :: n_spec     ! number of species
         real*8, intent(in) :: V(n_spec)   ! particle volumes (m^3)
         type(environ), intent(in) :: env  ! environment state
         type(material), intent(in) :: mat ! material properties
         logical, intent(in) :: init       ! true if first Newton loop
         real*8, intent(in) :: x           ! independent variable to solve for
         real*8, intent(out) :: f          ! function to solve
         real*8, intent(out) :: df         ! derivative df/dx
       end subroutine func
    end interface
    
    integer iter, k
    real*8 delta_f, delta_x, f, old_f, df

    call func(n_spec, V, env, mat, .true., x, f, df)
    old_f = f

    iter = 0
    do
       iter = iter + 1

       delta_x = f / df
       x = x - delta_x
       call func(n_spec, V, env, mat, .false., x, f, df)
       delta_f = f - old_f
       old_f = f
       
       if (iter .ge. iter_max) then
          write(0,*) 'ERROR: Newton iteration failed to terminate'
          write(0,*) 'iter_max = ', iter_max
          write(0,*) 'x = ', x
          write(0,*) 'V = ', V
          call exit(1)
       end if
       
       if ((abs(delta_x) .lt. x_tol) .and. (abs(delta_f) .lt. f_tol)) exit
    enddo

  end subroutine cond_newt

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine cond_growth_rate_func(n_spec, V, env, mat, init, dmdt, f, df)

    ! Return the error function value and its derivative.

    use mod_array
    use mod_util
    use mod_environ
    use mod_material
    use mod_constants

    integer, intent(in) :: n_spec ! number of species
    real*8, intent(in) :: V(n_spec) ! particle volumes (m^3)
    type(environ), intent(in) :: env     ! environment state
    type(material), intent(in) :: mat    ! material properties
    logical, intent(in) :: init       ! true if first Newton loop
    real*8, intent(in) :: dmdt   ! mass growth rate dm/dt (kg s^{-1})
    real*8, intent(out) :: f  ! error
    real*8, intent(out) :: df ! derivative of error with respect to x
    
    ! local variables
    real*8, save :: k_a, k_ap, k_ap_div, D_v, D_v_div, D_vp, d_p, pv
    real*8, save :: rat, fact1, fact2, c1, c2, c3, c4, c5
    real*8, save :: M_water, M_solute, rho_water, rho_solute
    real*8, save :: eps, nu, g_water, g_solute

    real*8 T_a ! droplet temperature (K), determined as part of solve

! DEBUG
    real*8 t1, t2, t3
! DEBUG

    if (init) then
       ! Start of new Newton loop, compute all constants

!       write(*,*) 'V = ', V

       M_water = average_water_quantity(V, mat, mat%M_w)     ! (kg mole^{-1})
       M_solute = average_solute_quantity(V, mat, mat%M_w)   ! (kg mole^{-1})
       nu = average_solute_quantity(V, mat, dble(mat%nu))    ! (1)
       eps = average_solute_quantity(V, mat, mat%eps)        ! (1)
       rho_water = average_water_quantity(V, mat, mat%rho)   ! (kg m^{-3})
       rho_solute = average_solute_quantity(V, mat, mat%rho) ! (kg m^{-3})
       g_water = total_water_quantity(V, mat, mat%rho)       ! (kg)
       g_solute = total_solute_quantity(V, mat, mat%rho)     ! (kg)

!       write(*,*) 'M_water = ', M_water
!       write(*,*) 'M_solute = ', M_solute
!       write(*,*) 'nu = ', nu
!       write(*,*) 'eps = ', eps
!       write(*,*) 'rho_water = ', rho_water
!       write(*,*) 'rho_solute = ', rho_solute
!       write(*,*) 'g_water = ', g_water
!       write(*,*) 'g_solute = ', g_solute
       
       call particle_vol_base(n_spec, V, pv) ! m^3
       d_p = vol2diam(pv) ! m
       
       ! molecular diffusion coefficient uncorrected
       D_v = 0.211d-4 / (env%p / const%atm) &
            * (env%T / 273d0)**1.94d0 ! m^2 s^{-1}
       ! molecular diffusion coefficient corrected for non-continuum effects
       D_v_div = 1d0 + (2d0 * D_v * 1d-4 / (const%alpha * d_p)) &
            * (2d0 * const%pi * M_water / (const%R * env%T))**0.5d0
       D_vp = D_v / D_v_div
       
!       write(*,*) 'D_v = ', D_v
!       write(*,*) 'D_vp = ', D_vp
       
       ! TEST: use the basic expression for D_vp
       ! D_vp = D_v                ! m^2 s^{-1}
       ! FIXME: check whether we can reinstate the correction
       
       ! thermal conductivity uncorrected
       k_a = 1d-3 * (4.39d0 + 0.071d0 * env%T) ! J m^{-1} s^{-1} K^{-1}
       k_ap_div = 1d0 + 2d0 &
            * k_a / (const%alpha * d_p * const%rho_a * const%cp) &
            * (2d0 * const%pi * const%M_a / (const%R * env%T))**0.5d0 ! dim-less
       ! thermal conductivity corrected
       k_ap = k_a / k_ap_div     ! J m^{-1} s^{-1} K^{-1}

!       write(*,*) 'k_a = ', k_a
!       write(*,*) 'k_ap = ', k_ap
       
       rat = sat_vapor_pressure(env) / (const%R * env%T)
       fact1 = const%L_v * M_water / (const%R * env%T)
       fact2 = const%L_v / (2d0 * const%pi * d_p * k_ap * env%T)

!       write(*,*) 'sat_vapor_pressure = ', sat_vapor_pressure(env)
!       write(*,*) 'sig = ', const%sig
!       write(*,*) 'L_v = ', const%L_v
!       write(*,*) 'R = ', const%R

       t1 = rho_water * const%R * env%T &
            / (sat_vapor_pressure(env) * D_vp * M_water)
       t2 = const%L_v * rho_water / (k_ap * env%T)
       t3 = t2 * (fact1 - 1d0)
 
!       write(*,*) 't1 = ', t1
!       write(*,*) 't2 = ', t2
!       write(*,*) 't3 = ', t3
      
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
    
!    write(*,*) '******************************************'
!    write(*,*) 'dmdt = ', dmdt
!    write(*,*) 'T_a = ', T_a

  end subroutine cond_growth_rate_func
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine equilibriate_particle(n_spec, V, env, mat)

    ! add water to the particle until it is in equilibrium

    use mod_util
    use mod_array
    use mod_environ
    use mod_material
    use mod_constants

    integer, intent(in) :: n_spec      ! number of species
    real*8, intent(inout) :: V(n_spec) ! particle volumes (m^3)
    type(environ), intent(in) :: env   ! environment state
    type(material), intent(in) :: mat  ! material properties

    ! parameters
    integer, parameter :: it_max = 400     ! maximum iterations
    real*8, parameter :: pv_rel_tol = 1d-6 ! dw relative convergence tolerance
    real*8, parameter :: f_tol = 1d-15     ! function convergence tolerance
    integer, parameter :: iter_max = 100   ! maximum number of iterations
    real*8, parameter :: dw_init = 1d0     ! initial value
    
    real*8 dw ! wet diameter of particle
    real*8 dw_tol, pv

    pv = sum(V)
    dw = dw_init
    ! convert volume relative tolerance to diameter absolute tolerance
    dw_tol = vol2diam(pv * (1d0 + pv_rel_tol)) - vol2diam(pv)
    call cond_newt(n_spec, V, dw, env, mat, equilibriate_func, &
         dw_tol, f_tol, iter_max)

    V(mat%i_water) = diam2vol(dw) - pv

  end subroutine equilibriate_particle

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine equilibriate_func(n_spec, V, env, mat, init, dw, f, df)

    use mod_util
    use mod_array
    use mod_environ
    use mod_material
    use mod_constants

    integer, intent(in) :: n_spec     ! number of species
    real*8, intent(in) :: V(n_spec)   ! particle volumes (m^3)
    type(environ), intent(in) :: env  ! environment state
    type(material), intent(in) :: mat ! material properties
    logical, intent(in) :: init       ! true if first Newton loop
    real*8, intent(in) :: dw          ! wet diameter (m)
    real*8, intent(out) :: f          ! function value
    real*8, intent(out) :: df         ! function derivative df/dx

    real*8, save :: c0, c1, c3, c4, dc0, dc2, dc3
    real*8, save :: A, B
    real*8, save ::  pv
    real*8, save :: M_water, M_solute, rho_water, rho_solute
    real*8, save :: eps, nu, g_water, g_solute

    if (init) then
       ! Start of new Newton loop, compute all constants

       M_water = average_water_quantity(V, mat, mat%M_w)     ! (kg mole^{-1})
       M_solute = average_solute_quantity(V, mat, mat%M_w)   ! (kg mole^{-1})
       nu = average_solute_quantity(V, mat, dble(mat%nu))    ! (1)
       eps = average_solute_quantity(V, mat, mat%eps)        ! (1)
       rho_water = average_water_quantity(V, mat, mat%rho)   ! (kg m^{-3})
       rho_solute = average_solute_quantity(V, mat, mat%rho) ! (kg m^{-3})
       g_water = total_water_quantity(V, mat, mat%rho)       ! (kg)
       g_solute = total_solute_quantity(V, mat, mat%rho)     ! (kg)

       call particle_vol_base(n_spec, V, pv)

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
  
end module mod_condensation
