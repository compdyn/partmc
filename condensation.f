! Condensation
!

module mod_condensation
contains
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine condense_particles(n_bin, TDV, n_spec, MH, VH, rho, &
       del_t,bin_v, bin_r, bin_g, bin_gs, bin_n, dlnr)

    use mod_array

    integer, intent(in) :: n_bin ! number of bins
    integer, intent(in) :: TDV   ! second dimension of VH
    integer, intent(in) :: n_spec       ! number of species
    integer, intent(inout) :: MH(n_bin) ! number of particles per bin
    real*8, intent(inout) :: VH(n_bin,TDV,n_spec) ! particle volumes (m^3)
    real*8, intent(in) :: rho(n_spec) ! density of species (kg m^{-3})
    real*8, intent(in) :: del_t         ! total time to integrate
    real*8, intent(in) :: bin_v(n_bin) ! volume of particles in bins (m^3)
    real*8, intent(in) ::  bin_r(n_bin) ! radius of particles in bins (m)
    real*8, intent(out) :: bin_g(n_bin) ! mass in bins  
    real*8, intent(out) :: bin_gs(n_bin,n_spec) ! species mass in bins
    integer, intent(inout) :: bin_n(n_bin)      ! number in bins
    real*8, intent(in) :: dlnr                  ! bin scale factor
    
    ! local variables
    integer bin, j, new_bin, k
    real*8 pv
    
    do bin = 1,n_bin
       write(*,*) 'condensation in bin ', bin
       do j = 1,MH(bin)
          ! remove the particle volume from the bin data
          call particle_vol_base(n_spec, VH(bin,j,:), pv)
          bin_g(bin) = bin_g(bin) - pv
          do k = 1,n_spec
             bin_gs(bin,k) = bin_gs(bin,k) - VH(bin,j,k)
          end do
          
          ! do the condensation
          call condense_particle(n_spec, VH(bin,j,:), rho, del_t)
       end do
    end do

    ! re-sort the particles into bins. This has to be done after all
    ! particles are advanced, otherwise we will lose track of which
    ! ones have been advanced and which have not.
    do bin = 1,n_bin
       j = 1
       do while (j .le. MH(bin))
          ! find the new volume
          call particle_vol_base(n_spec, VH(bin,j,:), pv)
          call particle_in_bin(pv, n_bin, bin_v, new_bin)

          ! update bin volumes
          bin_g(new_bin) = bin_g(new_bin) + pv
          do k = 1,n_spec
             bin_gs(new_bin,k) = bin_gs(new_bin,k) + VH(bin,j,k)
          end do
          
          ! if the bin number has changed, move the particle
          if (bin .ne. new_bin) then
!             write(*,*) 'n_bin = ', n_bin, ' TDV = ', TDV, ' n_spec = ', n_spec

!             write(*,*) 'bin = ', bin, ' j = ', j, ' MH(bin) = ', MH(bin)
             ! update the bin structures
             bin_n(bin) = bin_n(bin) - 1
             if (bin_n(bin) .lt. 0) then
                write(*,*) 'ERROR: invalid bin_n'
!                write(*,'(a15,a15,a15)') 'bin','bin_n(bin)', 'new_bin'
!                write(*,'(i15,i15,i15)') bin, bin_n(bin), new_bin
                call exit(2)
             end if
             bin_n(new_bin) = bin_n(new_bin) + 1
!             write(*,*) 'bin = ', bin, ' bin_n(bin) = ', bin_n(bin)
!             write(*,*) 'new_bin = ', new_bin, ' bin_n(new_bin) = ', bin_n(new_bin)
         
             ! move the particle to the new bin, leaving a hole
             MH(new_bin) = MH(new_bin) + 1
             if (MH(new_bin) .gt. TDV) then
                write(*,*) 'ERROR: too many particles in bin ', bin
                call exit(2)
             end if
!             write(*,*) 'new_bin = ', new_bin, ' MH(new_bin) = ', MH(new_bin)
             do k = 1,n_spec
!                write(*,*) 'bin = ', bin, ' j = ', j, ' k = ', k
!                write(*,*) 'new_bin = ', new_bin, ' MH(new_bin) = ', MH(new_bin), ' k = ', k
                VH(new_bin,MH(new_bin),k) = VH(bin,j,k)
             end do
             
             ! copy the last particle in the current bin into the hole
             do k = 1,n_spec
!                write(*,*) 'bin = ', bin, ' MH(bin) = ', MH(bin), ' k = ', k
!                write(*,*) 'bin = ', bin, ' j = ', j, ' k = ', k
                VH(bin,j,k) = VH(bin,MH(bin),k)
             end do
             MH(bin) = MH(bin) - 1
!             write(*,*) 'bin = ', bin, ' MH(bin) = ', MH(bin)
             if (MH(bin) .lt. 0) then
                write(*,*) 'ERROR: invalid MH'
                call exit(2)
             end if

             ! in this case, don't advance j, so that we will still
             ! process the particle we just moved into the hole
          else
             ! if we didn't move the particle, advance j to process
             ! the next particle
             j = j + 1
          end if
       end do
    end do

    ! FIXME: the approach above is inefficient because we might
    ! reprocess particles. For example, if we are doing bin 1 and we
    ! shift a particle up to bin 2, when we do bin 2 we will reprocess
    ! it. It seems to be more trouble than it's worth to worry about
    ! this yet, however.

  end subroutine condense_particles

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine condense_particle(n_spec, V, rho, del_t)

    integer, intent(in) :: n_spec ! number of species
    real*8, intent(inout) :: V(n_spec) ! particle volumes (m^3)
    real*8, intent(in) :: rho(n_spec) ! density of species (kg m^{-3})
    real*8, intent(in) :: del_t ! total time to integrate

    real*8 time_step, time
    logical done
    
    time = 0d0
    done = .false.
    do while (.not. done)
       call condense_step(n_spec, V, rho, del_t - time, time_step, done)
       time = time + time_step
    end do
    
  end subroutine condense_particle

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine condense_step(n_spec, V, rho, max_dt, dt, done)

    integer, intent(in) :: n_spec ! number of species
    real*8, intent(inout) :: V(n_spec) ! particle volumes (m^3)
    real*8, intent(in) :: rho(n_spec) ! density of species (kg m^{-3})  
    real*8, intent(in) :: max_dt ! maximum timestep to integrate
    real*8, intent(out) :: dt ! actual timestep used
    logical, intent(out) :: done ! whether we reached the maximum timestep

    ! parameters
    real*8 T, RH, p00, T0, p
    integer i_water
      
    parameter (T = 298d0)     ! Temperature of gas medium (K)
    parameter (RH = 1.01d0)   ! Relative humidity (1)
    parameter (p00 = 611d0)   ! equilibrium water vapor pressure at 273 K (Pa)
    parameter (T0 = 273.15d0) ! Freezing point of water (K)
    parameter (p = 1d5)       ! ambient pressure (Pa)
    parameter (i_water = 3)   ! water species number

    ! local variables
    real*8 dmdt, p0T, dvdt

    done = .false.
    call find_condense_timestep(n_spec, V, rho, dt)
    if (dt .ge. max_dt) then
       dt = max_dt
       done = .true.
    end if

    ! vapor pressure at temperature T
    p0T = p00 * 10d0**(7.45d0 * (T - T0) / (T - 38d0)) ! Pa

    call cond_newt(n_spec, V, rho, p0T, RH, T, p, i_water, dmdt)

    ! forward-Euler for now
    dvdt = dmdt / rho(i_water) ! m^3 s^{-1}
    V(i_water) = V(i_water) + dt * dvdt
    V(i_water) = max(0d0, V(i_water))
   
  end subroutine condense_step

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine find_condense_timestep(n_spec, V, rho, dt)

    use mod_array

    integer, intent(in) :: n_spec ! number of species
    real*8, intent(in) :: V(n_spec) ! particle volumes (m^3)
    real*8, intent(in) :: rho(n_spec) ! density of species (kg m^{-3})  
    real*8, intent(out) :: dt ! timestep to use

    real*8 pv

    call particle_vol_base(n_spec, V, pv)
    dt = 5d-3

  end subroutine find_condense_timestep

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
  subroutine cond_newt(n_spec, V, rho, p0T, RH, T, p, i_water, dmdt)
    
    ! Newton's method to solve the error equation, determining the
    ! growth rate dm/dt. The extra variable T_a is the local
    ! temperature, which is also implicitly determined, but not
    ! returned at this point.

    use mod_util
    use mod_array

    integer, intent(in) :: n_spec ! number of species
    real*8, intent(in) :: V(n_spec) ! particle volumes (m^3)
    real*8, intent(in) :: rho(n_spec) ! density of species (kg m^{-3})  
    real*8, intent(in) :: p0T ! vapor pressure at temperature T (Pa)
    real*8, intent(in) :: RH  ! relative humidity (???)
    real*8, intent(in) :: T   ! ambient temperature (K)
    real*8, intent(in) :: p   ! ambient pressure (Pa)
    integer, intent(in) :: i_water ! water species number
    real*8, intent(out) :: dmdt  ! dm/dt (kg s^{-1})

    ! parameters
    integer iter_max
    real*8 dmdt_min, dmdt_max, dmdt_tol, f_tol

    parameter (dmdt_min = 0d0)      ! minimum value of dm/dt (kg s^{-1})
    parameter (dmdt_max = 1d-3)     ! maximum value of dm/dt (kg s^{-1})
    parameter (dmdt_tol = 1d-15) ! dm/dt tolerance for convergence
    parameter (f_tol = 1d-15) ! function tolerance for convergence
    parameter (iter_max = 400)   ! maximum number of iterations

    ! local variables
    integer iter, k
    real*8 g_water, g_solute, pv
    real*8 T_a, delta_f, delta_dmdt, f, old_f, df, d

    g_water = V(i_water) * rho(i_water)
    g_solute = 0d0
    do k = 1,n_spec
       if (k .ne. i_water) then
          g_solute = g_solute + V(k) * rho(k)
       end if
    end do

    call particle_vol_base(n_spec, V, pv)
    call vol2diam(pv, d)

    dmdt = (dmdt_min + dmdt_max) / 2d0
    call cond_func(dmdt, d, g_water, g_solute, p0T, RH, T, p, f, df, T_a)
    old_f = f
    
    do iter = 1,iter_max
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
       
       if ((abs(delta_dmdt) .lt. dmdt_tol) &
            .and. (abs(delta_f) .lt. f_tol)) then
          return ! successful termination of Newton's method
       endif
    enddo
    
    write(0,*) 'ERROR: Newton iteration had too many iterations'
    write(0,'(a15,a17)') 'dmdt', 'max iterations'
    write(0,'(g15.10,i17)') dmdt, iter_max
    call exit(2)
    
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
  
end module mod_condensation
