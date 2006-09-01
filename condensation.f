! Condensation
!

module mod_condensation
contains
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine condense_particles(n_bin, TDV, n_spec, MH, VH, rho, i_water, &
       del_t, bin_v, bin_r, bin_g, bin_gs, bin_n, dlnr)

    use mod_array
    use mod_array_hybrid
    use mod_bin

    integer, intent(in) :: n_bin ! number of bins
    integer, intent(in) :: TDV   ! second dimension of VH
    integer, intent(in) :: n_spec       ! number of species
    integer, intent(inout) :: MH(n_bin) ! number of particles per bin
    real*8, intent(inout) :: VH(n_bin,TDV,n_spec) ! particle volumes (m^3)
    real*8, intent(in) :: rho(n_spec) ! density of species (kg m^{-3})
    integer, intent(in) :: i_water ! water species number
    real*8, intent(in) :: del_t         ! total time to integrate
    real*8, intent(in) :: bin_v(n_bin) ! volume of particles in bins (m^3)
    real*8, intent(in) ::  bin_r(n_bin) ! radius of particles in bins (m)
    real*8, intent(inout) :: bin_g(n_bin) ! mass in bins  
    real*8, intent(inout) :: bin_gs(n_bin,n_spec) ! species mass in bins
    integer, intent(inout) :: bin_n(n_bin)      ! number in bins
    real*8, intent(in) :: dlnr                  ! bin scale factor
    
    ! local variables
    integer bin, j, new_bin, k
    real*8 pv

    do bin = 1,n_bin
       do j = 1,MH(bin)
          call condense_particle(n_spec, VH(bin,j,:), rho, i_water, del_t)
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

  end subroutine condense_particles

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine condense_particle(n_spec, V, rho, i_water, del_t)

    ! integrates the condensation growth or decay ODE for total time
    ! del_t

    integer, intent(in) :: n_spec ! number of species
    real*8, intent(inout) :: V(n_spec) ! particle volumes (m^3)
    real*8, intent(in) :: rho(n_spec) ! density of species (kg m^{-3})
    integer, intent(in) :: i_water ! water species number
    real*8, intent(in) :: del_t ! total time to integrate

    real*8 time_step, time
    logical done
    
    time = 0d0
    done = .false.
    do while (.not. done)
       call condense_step_rk_fixed(n_spec, V, rho, i_water, del_t - time, &
            time_step, done)
       time = time + time_step
    end do
    
  end subroutine condense_particle

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine condense_step_euler(n_spec, V, rho, i_water, max_dt, dt, done)

    ! Does one timestep (determined by this subroutine) of the
    ! condensation ODE. The timestep will not exceed max_dt, but might
    ! be less. If we in fact step all the way to max_dt then done will
    ! be true.
    
    integer, intent(in) :: n_spec ! number of species
    real*8, intent(inout) :: V(n_spec) ! particle volumes (m^3)
    real*8, intent(in) :: rho(n_spec) ! density of species (kg m^{-3})  
    integer, intent(in) :: i_water ! water species number
    real*8, intent(in) :: max_dt ! maximum timestep to integrate
    real*8, intent(out) :: dt ! actual timestep used
    logical, intent(out) :: done ! whether we reached the maximum timestep

    real*8 dvdt

    done = .false.
    call find_condense_timestep_variable(n_spec, V, rho, i_water, dt)
    if (dt .ge. max_dt) then
       dt = max_dt
       done = .true.
    end if

    call cond_newt(n_spec, V, rho, i_water, dvdt)
    V(i_water) = V(i_water) + dt * dvdt
    V(i_water) = max(0d0, V(i_water))
   
  end subroutine condense_step_euler

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine condense_step_rk_fixed(n_spec, V, rho, i_water, max_dt, dt, done)

    ! Does one timestep (determined by this subroutine) of the
    ! condensation ODE. The timestep will not exceed max_dt, but might
    ! be less. If we in fact step all the way to max_dt then done will
    ! be true.
    
    integer, intent(in) :: n_spec ! number of species
    real*8, intent(inout) :: V(n_spec) ! particle volumes (m^3)
    real*8, intent(in) :: rho(n_spec) ! density of species (kg m^{-3})  
    integer, intent(in) :: i_water ! water species number
    real*8, intent(in) :: max_dt ! maximum timestep to integrate
    real*8, intent(out) :: dt ! actual timestep used
    logical, intent(out) :: done ! whether we reached the maximum timestep

    done = .false.
    call find_condense_timestep_variable(n_spec, V, rho, i_water, dt)
    if (dt .ge. max_dt) then
       dt = max_dt
       done = .true.
    end if

    call condense_step_rk(n_spec, V, rho, i_water, dt)
   
  end subroutine condense_step_rk_fixed

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine condense_step_rk(n_spec, V, rho, i_water, dt)

    ! Does one fixed timestep of RK4.
    
    integer, intent(in) :: n_spec ! number of species
    real*8, intent(inout) :: V(n_spec) ! particle volumes (m^3)
    real*8, intent(in) :: rho(n_spec) ! density of species (kg m^{-3})
    integer, intent(in) :: i_water ! water species number
    real*8, intent(out) :: dt ! timestep

    ! local variables
    real*8 dmdt, dvdt, k1, k2, k3, k4
    real*8 V_tmp(n_spec)

    V_tmp = V

    ! step 1
    call cond_newt(n_spec, V, rho, i_water, k1)

    ! step 2
    V_tmp(i_water) = V(i_water) + dt * k1 / 2d0
    call cond_newt(n_spec, V_tmp, rho, i_water, k2)

    ! step 3
    V_tmp(i_water) = V(i_water) + dt * k2 / 2d0
    call cond_newt(n_spec, V_tmp, rho, i_water, k3)

    ! step 4
    V_tmp(i_water) = V(i_water) + dt * k3
    call cond_newt(n_spec, V_tmp, rho, i_water, k4)

    V(i_water) = V(i_water) + dt * (k1 / 6d0 + k2 / 3d0 + k3 / 3d0 + k4 / 6d0)

    V(i_water) = max(0d0, V(i_water))
   
  end subroutine condense_step_rk

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine find_condense_timestep_constant(n_spec, V, rho, i_water, dt)

    ! constant timestep

    use mod_array

    integer, intent(in) :: n_spec ! number of species
    real*8, intent(in) :: V(n_spec) ! particle volumes (m^3)
    real*8, intent(in) :: rho(n_spec) ! density of species (kg m^{-3})  
    integer, intent(in) :: i_water ! water species number
    real*8, intent(out) :: dt ! timestep to use

    dt = 5d-3

  end subroutine find_condense_timestep_constant

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine find_condense_timestep_variable(n_spec, V, rho, i_water, dt)

    ! timestep is proportional to V / (dV/dt)

    use mod_array

    integer, intent(in) :: n_spec ! number of species
    real*8, intent(in) :: V(n_spec) ! particle volumes (m^3)
    real*8, intent(in) :: rho(n_spec) ! density of species (kg m^{-3})  
    integer, intent(in) :: i_water ! water species number
    real*8, intent(out) :: dt ! timestep to use

    ! parameters
    real*8 scale
    parameter (scale = 1d0) ! scale factor for timestep

    real*8 pv, dvdt

    call particle_vol_base(n_spec, V, pv)
    call cond_newt(n_spec, V, rho, i_water, dvdt)
    dt = scale * pv / dvdt

  end subroutine find_condense_timestep_variable
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
  subroutine cond_newt(n_spec, V, rho, i_water, dvdt)
    
    ! Newton's method to solve the error equation, determining the
    ! growth rate dm/dt. The extra variable T_a is the local
    ! temperature, which is also implicitly determined, but not
    ! returned at this point.

    use mod_util
    use mod_array

    integer, intent(in) :: n_spec ! number of species
    real*8, intent(in) :: V(n_spec) ! particle volumes (m^3)
    real*8, intent(in) :: rho(n_spec) ! density of species (kg m^{-3})  
    integer, intent(in) :: i_water ! water species number
    real*8, intent(out) :: dvdt  ! dv/dt (m^3 s^{-1})

    ! parameters
    integer iter_max
    real*8 T, RH, p00, T0, p
    real*8 dmdt_min, dmdt_max, dmdt_tol, f_tol

    parameter (T = 298d0)     ! temperature of gas medium (K)
    parameter (RH = 1.01d0)   ! relative humidity (1)
    parameter (p00 = 611d0)   ! equilibrium water vapor pressure at 273 K (Pa)
    parameter (T0 = 273.15d0) ! freezing point of water (K)
    parameter (p = 1d5)       ! ambient pressure (Pa)
    parameter (dmdt_min = 0d0)      ! minimum value of dm/dt (kg s^{-1})
    parameter (dmdt_max = 1d-3)     ! maximum value of dm/dt (kg s^{-1})
    parameter (dmdt_tol = 1d-15) ! dm/dt tolerance for convergence
    parameter (f_tol = 1d-15) ! function tolerance for convergence
    parameter (iter_max = 400)   ! maximum number of iterations

    ! local variables
    integer iter, k
    real*8 g_water, g_solute, pv, p0T
    real*8 dmdt, T_a, delta_f, delta_dmdt, f, old_f, df, d

    ! vapor pressure at temperature T
    p0T = p00 * 10d0**(7.45d0 * (T - T0) / (T - 38d0)) ! Pa

    g_water = V(i_water) * rho(i_water)
    g_solute = 0d0
    do k = 1,n_spec
       if (k .ne. i_water) then
          g_solute = g_solute + V(k) * rho(k)
       end if
    end do

    call particle_vol_base(n_spec, V, pv)
    d = vol2diam(pv)

    dmdt = (dmdt_min + dmdt_max) / 2d0
    call cond_func(dmdt, d, g_water, g_solute, p0T, RH, T, p, f, df, T_a)
    old_f = f

    iter = 0
    do
       iter = iter + 1

       delta_dmdt = f / df
       dmdt = dmdt - delta_dmdt
       call cond_func(dmdt, d, g_water, g_solute, p0T, RH, T, p, f, df, T_a)
       delta_f = f - old_f
       old_f = f
       
       if ((dmdt .lt. dmdt_min) .or. (dmdt .gt. dmdt_max)) then
          write(0,*) 'ERROR: Newton iteration exceeded bounds'
          write(0,'(a15,a15,a15)') 'dmdt', 'lower bound', 'upper bound'
          write(0,'(g15.10,g15.10,g15.10)') dmdt, dmdt_min, dmdt_max
          call exit(1)
       endif

       if (iter .ge. iter_max) then
          write(0,*) 'ERROR: Newton iteration had too many iterations'
          write(0,'(a15,a15)') 'dmdt', 'iter_max'
          write(0,'(g15.10,i15)') dmdt, iter_max
          call exit(2)
       end if
       
       if ((abs(delta_dmdt) .lt. dmdt_tol) &
            .and. (abs(delta_f) .lt. f_tol)) exit
    enddo

    dvdt = dmdt / rho(i_water)

  end subroutine cond_newt

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine cond_func(x, d_p, g1, g2, p0T, RH, T, p, f, df, T_a)

    ! Return the error function value and its derivative.

    real*8, intent(in) :: x   ! mass growth rate dm/dt (kg s^{-1})
    real*8, intent(in) :: d_p ! diameter (m)
    real*8, intent(in) :: g1  ! water mass (kg)
    real*8, intent(in) :: g2  ! solute mass (kg)
    real*8, intent(in) :: p0T ! vapor pressure at temperature T (Pa)
    real*8, intent(in) :: RH  ! relative humidity (1)
    real*8, intent(in) :: T   ! ambient temperature (K)
    real*8, intent(in) :: p   ! ambient pressure (Pa)
    real*8, intent(out) :: f  ! error
    real*8, intent(out) :: df ! derivative of error with respect to x
    real*8, intent(out) :: T_a ! droplet temperature (K)
    
    ! parameters
    real*8 rho, rho_a, rho_n, M_w, M_a, M_s, sig, R, L_v, alpha
    real*8 p00, T0, cp, eps, atm
    integer nu
    parameter (rho = 1000d0)  ! water density (kg m^{-3})
    parameter (rho_a = 1.25d0) ! air density (kg m^{-3})
    parameter (rho_n = 1800d0) ! solute density (kg m^{-3})
    parameter (M_w = 18d-3)   ! molecular weight of water (kg mole^{-1})
    parameter (M_a = 28d-3)   ! molecular weight of air (kg mole^{-1})
    parameter (M_s = 132d-3)  ! molecular weight of solute (kg mole^{-1})
    parameter (sig = 0.073d0) ! surface energy (J m^{-2})
    parameter (R = 8.314d0)   ! universal gas constant (J mole^{-1} K^{-1})
    parameter (L_v = 2.5d6)   ! latent heat (J kg^{-1})
    parameter (alpha = 1d0)   ! accomodation coefficient (the value 0.045 is also used sometimes)
    parameter (p00 = 6.11d2)  ! water saturation vapor pressure at temp T = 273 K (Pa)
    parameter (T0 = 273.15d0) ! freezing point of water (K)
    parameter (cp = 1005d0)   ! specific heat of water (J kg^{-1} K^{-1})
    parameter (nu = 3)        ! number of ions in the solute
    parameter (eps = 0.25d0)  ! solubility of aerosol material (1)
    parameter (atm = 101325d0) ! atmospheric pressure (Pa)

    ! FIXME: nu and eps should be passed as arguments

    ! FIXME: nu, eps, M_s should really be arrays of length n_spec
      
    real*8 pi
    parameter (pi = 3.14159265358979323846d0)

    ! local variables
    real*8 k_a, k_ap, k_ap_div, D_v, D_vp
    real*8 rat, fact1, fact2, c1, c2, c3, c4, c5

    ! molecular diffusion coefficient uncorrected
    D_v = 0.211d-4 / (p / atm) * (T / 273d0)**1.94d0 ! m^2 s^{-1}

    ! molecular diffusion coefficient corrected for non-continuum effects
    ! D_v_div = 1d0 + (2d0 * D_v * 1d-4 / (alpha * d_p)) &
    !      * (2 * pi * M_w / (R * T))**0.5d0
    ! D_vp = D_v / D_v_div

    ! TEST: use the basic expression for D_vp
    D_vp = D_v                ! m^2 s^{-1}
    ! FIXME: check whether we can reinstate the correction

    ! thermal conductivity uncorrected
    k_a = 1d-3 * (4.39d0 + 0.071d0 * T) ! J m^{-1} s^{-1} K^{-1}
    k_ap_div = 1d0 + 2d0 * k_a / (alpha * d_p * rho_a * cp) &
         * (2d0 * pi * M_a / (R * T))**0.5d0 ! dimensionless
    ! thermal conductivity corrected
    k_ap = k_a / k_ap_div     ! J m^{-1} s^{-1} K^{-1}
      
    rat = p0T / (R * T)
    fact1 = L_v * M_w / (R * T)
    fact2 = L_v / (2d0 * pi * d_p * k_ap * T)
    
    c1 = 2d0 * pi * d_p * D_vp * M_w * rat
    c2 = 4d0 * M_w * sig / (R * rho * d_p)
    c3 = c1 * fact1 * fact2
    c4 = L_v / (2d0 * pi * d_p * k_ap)
    ! incorrect expression from Majeed and Wexler:
    ! c5 = nu * eps * M_w * rho_n * r_n**3d0 &
    !     / (M_s * rho * ((d_p / 2)**3d0 - r_n**3))
    ! corrected according to Jim's note:
    c5 = dble(nu) * eps * M_w / M_s * g2 / &
         (g1 + (rho / rho_n) * eps * g2)
    
    T_a = T + c4 * x ! K
    
    f = x - c1 * (RH - exp(c2 / T_a - c5)) &
         / (1d0 + c3 * exp(c2 / T_a - c5))
    
    df = 1d0 + c1 * RH * (1d0 + c3 * exp(c2 / T_a -c5))**(-2d0) * c3 * &
         exp(c2 / T_a - c5) * (-1d0) * c2 * c4 / T_a**2d0 + c1 * &
         (exp(c2 / T_a - c5) * (-1d0) * c2 * c4 / T_a**2d0 * (1d0 + c3 &
         * exp(c2 / T_a -c5))**(-1d0) + exp(c2 / T_a - c5) * (-1d0) * &
         (1d0 + c3 * exp(c2 / T_a -c5))**(-2d0) * c3 * exp(c2 / T_a - &
         c5) * (-1d0) * c2 * c4 / T_a**2d0)
    
  end subroutine cond_func
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine equilibriate_particle(n_spec, V, rho, i_water, &
       nu, eps, M_s)

    ! add water to the particle until it is in equilibrium

    integer, intent(in) :: n_spec ! number of species
    real*8, intent(inout) :: V(n_spec) ! particle volumes (m^3)
    real*8, intent(in) :: rho(n_spec) ! density of species (kg m^{-3})  
    integer, intent(in) :: i_water ! water species number
    integer, intent(in) :: nu      ! number of ions in the solute
    real*8, intent(in) :: eps      ! solubility of aerosol material (1)
    real*8, intent(in) :: M_s      ! molecular weight of solute (kg mole^{-1})

    ! FIXME: nu, eps, M_s should really be arrays of length n_spec

    ! paramters
    real*8 x_min, x_max, x_tol

    real*8 pv

    call particle_vol_base(n_spec, V, pv)

    T = T0
    A = 4d0 * M_w * sig_w / (RR * T * rho_w)
    
    B = nu * eps * M_w * rho_n * vol2rad(pv)**3.d0 / (M_s * rho_w)
    
    c4 = log(RH) / 8.d0
    c3 = A / 8.d0
    
    dc3 = log(RH) / 2.d0
    dc2 = 3.d0 * A / 8.d0
    
    x1 = 0.d0
    x2 = 10.d0
    xacc = 1.d-15
    
    c1 = B - log(RH) * vol2rad(pv)**3d0
    c0 = A * vol2rad(pv)**3d0
    dc0 = c1
    
    call equilibriate_newt(x1, x2, xacc, x, c4, c3, c1, c0, dc3, dc2, dc0)
    
    rw = x / 2d0

  end subroutine equilibriate_particle

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine equilibriate_newt(x1, x2, xacc, x, c4, c3, c1, c0, dc3, &
       dc2, dc0)

    real*8, intent(out) :: x
    real*8, intent(in) :: x1
    real*8, intent(in) :: x2
    real*8, intent(in) :: xacc
    real*8, intent(in) :: c4
    real*8, intent(in) :: c3
    real*8, intent(in) :: c1
    real*8, intent(in) :: c0
    real*8, intent(in) :: dc3
    real*8, intent(in) :: dc2
    real*8, intent(in) :: dc0

    integer jmax
    parameter (jmax=400)
    
    integer j
    real*8 df, dx, f, d 

    x = 0.5d0 * (x1 + x2)
    
    do j = 1,jmax
       call equilibriate_func(x,f,df,d,c4,c3,c1,c0,dc3,dc2,dc0)
       dx = f / df
       x = x - dx
       if((x .lt. x1) .or. (x .gt. x2)) then
          write(6,*)'x1,x2,x ',x1,x2,x
          write(*,*) 'rtnewt jumped out of brackets'
          exit(2)
       endif
       if(abs(dx) .lt. xacc) then
          return
       endif
    enddo
    
    write(*,*) 'rtnewt exceeded maximum iteration '
    exit(2)

  end subroutine equilibriate_newt

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine equlibriate_func(x, f, df, d_p, c4, c3, c1, c0, dc3, dc2, dc0)

    real*8 x, f, df, d_p
    real*8 c4, c3, c1, c0, dc3, dc2, dc0

    f = c4 * x**4d0 - c3 * x**3d0 + c1 * x + c0
    df = dc3 * x**3d0 -dc2 * x**2d0 + dc0

  end subroutine equlibriate_func

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
end module mod_condensation
