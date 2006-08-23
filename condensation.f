C Condensation
C

      module condensation
      contains

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine cond(n_bin, TDV, n_spec, MH, VH, rho_p)

      integer, intent(in) :: n_bin ! number of bins
      integer, intent(in) :: TDV ! second dimension of VH
      integer, intent(in) :: n_spec ! number of species
      integer, intent(in) :: MH(n_bin) ! number of particles per bin
      real*8,  intent(in) :: rho_p(n_spec) ! density of species (kg m^{-3})
      real*8, intent(inout) :: VH(n_bin,TDV,n_spec) ! particle volumes (m^3)

C     parameters
      real*8 T, RH, p00, T0, pres
      integer i_water
      
      parameter (T = 298d0)     ! Temperature of gas medium (K)
      parameter (RH = 1.01d0)   ! Relative humidity (1)
      parameter (p00 = 611d0)   ! equilibrium water vapor pressure at 273 K (Pa)
      parameter (T0 = 273.15d0) ! Freezing point of water (K)
      parameter (pres = 1d5)    ! ambient pressure (Pa)
      parameter (i_water = 3)   ! water species number
      
      real*8 pi
      parameter (pi = 3.14159265358979323846d0)

C     local variables
      real*8 dmdt(n_bin,TDV)    ! growth rates (kg s^{-1})
      real*8 histot(n_bin,TDV)  ! normalized growth rate: histot=dmdt/m (s^{-1})
      real*8 rho(n_spec)        ! densities (kg m^{-3})
      real*8 p0T

C     vapor pressure at temperature T
      p0T = p00 * 10d0**(7.45d0 * (T - T0) / (T - 38d0)) ! Pa

C      write(6,*) 'p0T ', p0T, p00, RH
      
      call cond_rate(n_bin, TDV, n_spec, MH, T, RH, pres, i_water, rho,
     $     p0T, dmdt, histot, VH)

C     dmdt(i) and histot(i) are output
C     dmdt is growth rate of one droplet (kg s^{-1})
C     histot =  dmdt * 1d6 / e(i) (s^{-1})

      write(6,*)'ende condensation'
      stop

      end subroutine

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 

      subroutine cond_rate(n_bin, TDV, n_spec, MH, T, RH, p, i_water,
     $     rho, p0T, dmdt, histot, VH)

C     Calculation of the term dm/dt according to Majeed and Wexler,
C     Atmos. Env. (2001). Since Eq. (7) in this paper is an implicit
C     equation (T_a depends on dm/dt), a Newton solver is applied.

      use array_hybrid

      integer, intent(in) :: n_bin ! number of bins
      integer, intent(in) :: TDV ! second dimension of VH
      integer, intent(in) :: n_spec ! number of species
      integer, intent(in) :: MH(n_bin) ! number of particles per bin
      real*8, intent(in) :: T   ! ambient temperature (K)
      real*8, intent(in) :: RH  ! relative humidity (1)
      real*8, intent(in) :: p   ! ambient pressure (Pa)
      integer, intent(in) :: i_water ! water species number
      real*8, intent(in) :: rho(n_spec) ! densities of each species (kg m^{-3})
      real*8, intent(in) :: p0T ! vapor pressure at temperature T (Pa)
      real*8, intent(out) :: dmdt(n_bin, TDV) ! growth rates (kg s^{-1}) / (m^3 s^{-1})
      real*8, intent(out) :: histot(n_bin, TDV) ! growth constants (s^{-1})
      real*8, intent(in) :: VH(n_bin, TDV, n_spec) ! particle volumes (m^3)

C     parameters
      real*8 RR, M_w, sig_w, M_s, nu, x1, x2, x_tol, f_tol

      parameter (RR = 8.314d0)  ! universal gas constant (J kg^{-1} K{-1})
      parameter (M_w = 18d-3)   ! molecular weight of water (kg mol^{-1})
      parameter (sig_w = 0.073d0) ! surface tension of water  (J m^{-2})
      parameter (M_s = 132d-3)  ! molecular weight of soluble substance (kg mol^{-1})
      parameter (nu = 3d0)      ! number of ions resulting from the dissociation 
                                ! of one solute molecule (1)
      parameter (x1 = 0d0)      ! minimum value of dm/dt (kg s^{-1})
      parameter (x2 = 1d-3)     ! maximum value of dm/dt (kg s^{-1})
      parameter (x_tol = 1d-15) ! dm/dt tolerance for convergence
      parameter (f_tol = 1d-15) ! function tolerance for convergence

      real*8 pi
      parameter (pi = 3.14159265358979323846d0)

C     local variables
      real*8 g1, g2, x, d, pv
      integer i, j, k

      do i = 1,n_bin
         do j = 1,TDV
            dmdt(i,j) = 0d0
         enddo
      enddo
      
      do i = 1,n_bin
         write(*,*) 'i', i
         do j = 1,MH(i)
            write(6,*) 'VH ', i, j, VH(i,j,1), VH(i,j,2), VH(i,j,3)
            
            g1 = VH(i, j, i_water) * rho(i_water)
            g2 = 0d0
            do k = 1,n_spec
               if (k .ne. i_water) then
                  g2 = g2 + VH(i,j,k) * rho(k)
               endif
            enddo
            
            call particle_vol_hybrid(n_bin, TDV, n_spec, VH, i, j, pv)
            
            d = (6d0 / pi * pv)**(1d0/3d0) ! rn in m, d in m 
            write(6,*)'g1,g2,pv ',g1,g2,pv,d
            
            if ((g1 .ne. 0d0) .or. (g2 .ne. 0d0)) then
               call cond_newt(x1, x2, x_tol, f_tol, d, g1, g2, p0T, RH,
     $              T, p, x)
               dmdt(i,j) = x
               histot(i,j) = dmdt(i,j) / (g1 + g2) ! histot is normalized growth in s^{-1}
               dmdt(i,j) = dmdt(i,j) / rho(i_water) ! dmdt is now the volume growth in m^3 s-1
            endif
         enddo
      enddo

!     FIXME: is dm/dt supposed to be a mass or volume growth rate?
      
      do i = 1,n_bin
         do j = 1,MH(i)
            write(66,*) 'dmdt ', i, j, dmdt(i,j), histot(i,j),
     &           1d0/histot(i,j)
         enddo
      enddo
      
      stop

      end subroutine

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine cond_newt(x1, x2, x_tol, f_tol, d, g1, g2, p0T,
     &     RH, T, p, x)

C     Newton's method to solve the error equation, determining the
C     growth rate dm/dt. The extra variable T_a is the local
C     temperature, which is also implicitly determined, but not returned
C     at this point.

      real*8, intent(in) :: x1  ! lower bound on dm/dt during solve
      real*8, intent(in) :: x2  ! upper bound on dm/dt during solve
      real*8, intent(in) :: x_tol ! tolerance for delta x
      real*8, intent(in) :: f_tol ! tolerance for delta f
      real*8, intent(in) :: d   ! diameter (m)
      real*8, intent(in) :: g1  ! water mass (kg)
      real*8, intent(in) :: g2  ! solute mass (kg)
      real*8, intent(in) :: p0T ! vapor pressure at temperature T (Pa)
      real*8, intent(in) :: RH  ! relative humidity (???)
      real*8, intent(in) :: T   ! ambient temperature (K)
      real*8, intent(in) :: p   ! ambient pressure (Pa)
      real*8, intent(out) :: x  ! dm/dt (kg s^{-1})

      integer j_max
      parameter (j_max = 400)   ! maximum number of iterations

      external funcd            ! error function to solve

      integer j
      real*8 T_a, delta_f, delta_x, f, old_f, df

      x = (x1 + x2) / 2d0
      call cond_func(x, d, g1, g2, p0T, RH, T, p, f, df, T_a)
      old_f = f

C      write(6,*)'rtnewt1 ',x1,x2,x_tol,x
C      write(6,*)'rtnewt2 ',d,g1,g2
C      write(6,*)'rtnewt3 ',p0T,RH,T,p
       
      do j = 1,j_max
         delta_x = f / df
         x = x - delta_x
         call cond_func(x, d, g1, g2, p0T, RH, T, p, f, df, T_a)
         delta_f = f - old_f
         old_f = f

         if ((x .lt. x1) .or. (x .gt. x2)) then
            write(0,*) 'ERROR: Newton iteration exceeded bounds'
            write(0,'(a15,a15,a15)') 'x', 'lower bound', 'upper bound'
            write(0,'(g15.10,g15.10,g15.10)') x, x1, x2
            call exit(1)
         endif

         if ((abs(delta_x) .lt. x_tol)
     &        .and. (abs(delta_f) .lt. f_tol)) then
            return              ! successful termination of Newton's method
         endif
      enddo

      write(0,*) 'ERROR: Newton iteration had too many iterations'
      write(0,'(a15,a17)') 'x', 'max iterations'
      write(0,'(g15.10,i17)') x, j_max
      call exit(1)

      end subroutine
 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine cond_func(x, d_p, g1, g2, p0T, RH, T, p, f, df,
     $     T_a)

C     Return the error function value and its derivative.

      real*8, intent(in) :: x   ! mass growth rate dm/dt (kg s^{-1})
      real*8, intent(in) :: d_p ! diameter (m)
      real*8, intent(in) :: g1  ! water mass (kg)
      real*8, intent(in) :: g2  ! solute mass (kg)
      real*8, intent(in) :: p0T ! vapor pressure at temperature T (Pa)
      real*8, intent(in) :: RH  ! relative humidity (1)
      real*8, intent(in) :: T   ! ambient temperature (K)
      real*8, intent(in) :: p   ! ambient pressure (Pa)
      real*8, intent(out) :: f  ! error
      real*8, intent(out) :: df ! derivative of error
      real*8, intent(out) :: T_a ! droplet temperature (K)
      
C     parameters
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

      real*8 pi
      parameter (pi = 3.14159265358979323846d0)

C     local variables
      real*8 k_a, k_ap, k_ap1, D_v, D_vp
      real*8 rat, fact1, fact2, c1, c2, c3, c4, c5

C     molecular diffusion coefficient (basic)
      D_v = 0.211d0 / (p / atm) * (T / 273d0)**1.94d0 ! cm^2 s^{-1}

C     diffusion coefficient (corrected for non-continuum effects)
C     D_v_de = 1d0 + (2d0 * D_v * 1d-4 / (alpha * d_p))
C     &     * (2 * pi * M_w / (R * T))**0.5d0
C     D_vp = D_v / D_v_de

C     TEST: use the basic expression for D_vp
      D_vp = D_v * 1d-4         ! m^2 s^{-1}

C     thermal conductivity uncorrected and corrected
      k_a = 1d-3 * (4.39d0 + 0.071d0 * T) ! J m^{-1} s^{-1} K^{-1}
      k_ap1 = 1d0 + 2d0 * k_a / (alpha * d_p * rho_a * cp)
     &     * (2d0 * pi * M_a / (R * T))**0.5d0 ! J m^{-1} s^{-1} K^{-1}
      k_ap = k_a / k_ap1        ! J m^{-1} s^{-1} K^{-1}
      
      rat = p0T / (R * T)
      fact1 = L_v * M_w / (R * T)
      fact2 = L_v / (2d0 * pi * d_p * k_ap * T)
      
      c1 = 2d0 * pi * d_p * D_vp * M_w * rat
      c2 = 4d0 * M_w * sig / (R * rho * d_p)
      c3 = c1 * fact1 * fact2
      c4 = L_v / (2d0 * pi * d_p * k_ap)
C     c5 = nu * eps * M_w * rho_n * r_n**3d0
C     &     / (M_s * rho * ((d_p / 2)**3d0 - r_n**3))
C     c5 = nu * eps * M_w / M_s * g2 / g1
C     corrected according to Jim's note:
      c5 = nu * eps * M_w / M_s * g2 / (g1 + (rho / rho_n) * eps * g2)

      f = x - c1 * (RH - exp(c2 / (T + c4 * x) -c5))
     &     / (1d0 + c3 * exp(c2 / (T + c4 * x) - c5))
      
      df = 1d0 + c1 * RH
     &     * (1d0 + c3 * exp(c2 / (T + c4 * x) -c5))**(-2d0)
     &     * c3 * exp(c2 / (T + c4 * x) -c5)
     &     * (-1d0) * c2 * c4 / (T + c4 * x)**2d0
     &     + c1 * (exp(c2 / (T + c4 * x) - c5)
     &     * (-1d0) * c2 * c4 / (T + c4 * x)**2d0
     &     * (1d0 + c3 * exp(c2 / (T + c4 * x) -c5))**(-1d0)
     &     + exp(c2 / (T + c4 * x) -c5)
     &     * (-1d0) * (1d0 + c3 * exp(c2 / (T + c4 * x) -c5))**(-2d0)
     &     * c3 * exp(c2 / (T + c4 * x) - c5)
     &     * (-1d0) * c2 * c4 / (T + c4 * x)**2d0)
      
      T_a = T + c4 * x
      
      end subroutine

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      end module
