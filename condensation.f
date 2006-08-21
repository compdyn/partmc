C Condensation
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine condensation(n_bin, TDV, n_spec, MH, VH)

      integer n_spec, n_bin, TDV
      integer i_water
      real*8 VH(n_bin,TDV,n_spec)  ! INPUT/OUTPUT: volume of particle before/after condensation
      integer MH(n_bin)

      real*8 T, rho(n_spec), RH, pres, pmv, p0T
      real*8 dmdt(n_bin,TDV)       ! growth rate [kg s-1]
      real*8 histot(n_bin,TDV)        ! normalized growth rate: histot=dmdt/m [s-1]
      real*8 p00, T0
      
      parameter (T = 298d0)      ! Temperature of gas medium in K
      parameter (RH= 1.01d0)       ! Relative humidity
      parameter (p00 = 611d0) ! equilibrium water vapor pressure at 273 K in Pa
      parameter (T0  = 273.15d0) ! in K
      parameter (pres = 1000d0) ! ambient pressure in hPa
      parameter (i_water = 3)

      real*8 pi
      parameter (pi = 3.14159265358979323846d0)

      rho(1) = 1800d0
      rho(2) = 1800d0
      rho(3) = 1000d0

      p0T = p00 * 10d0**(7.45d0 * (T - T0) / (T - 38d0))
      pmv = p0T * RH
      write(6,*)'p0T,pmv ',p0T,p00,pmv,RH
  
! 
      call kond(n_bin,TDV,n_spec,MH,T, RH, pres
     $     ,i_water,rho, pmv, p0T,dmdt, histot, 1,n_bin,VH)

! dmdt(i) and histot(i) are output
! dmdt is growth rate of one droplet in kg s-1
! histot =  dmdt * 1.e6 / e(i) in s-1

      write(6,*)'ende condensation'
      stop
      return
      end

cn ****************************************************************
 
      subroutine kond(n_bin,TDV,n_spec,MH,T,RH,p
     $     ,i_water,rho, pmv,p0T,dmdt,histot,ia,ie,VH)

cn *** Calculation of the term dm/dt according to Majeed and Wexler, Atmos. Env. (2001)
cn *** Since Eq. (7) in this paper is an implicit equation (T_a depends on dm/dt), a Newton
cn *** solver is applied.

      integer   i,n_bin,TDV,n_spec,j,k
      integer   i_water
      integer   ia,ie
      real*8    xacc,x1,x2

      real*8    VH(n_bin,TDV,n_spec)
      integer    MH(n_bin)
      real*8    x,d
      real*8    dmdt(n_bin,TDV)                       ! growth rate in kg s-1
      real*8    e(n_bin,TDV)
      real*8    pv

      real*8    pmv
      real*8    histot(n_bin,TDV)                     ! growth constant in s-1
      real*8    p0T
      real*8    T                                ! ambient temperature in K
      real*8    p                                ! ambient pressure  in hPa
      real*8    RH
      real*8    RR,M_w,sig_w,M_s,rho(n_spec),pi,nu
      real*8    g1,g2

      real*8    gmin

      parameter (RR = 8.314d0, M_w = 18d-3)
      parameter (sig_w = 0.073d0, M_s = 132d-3)
      parameter (nu = 3d0)
      parameter (gmin = 0d0)

      parameter (pi = 3.14159265358979323846d0)

      do i=ia,ie
         do j=1,TDV
            dmdt(i,j) = 0d0
         enddo
      enddo

      x1 = 0d0
      x2 = 1000d-7
      xacc = 1d-15

      g1 = 0d0
      g2 = 0d0

      write(6,*)'ia, ie ',ia,ie
      do i=1,160
         write(*,*)'i', i
         do j=1,MH(i)
            write(6,*)'VH ',i,j,VH(i,j,1),VH(i,j,2),VH(i,j,3)
         
            g1 = VH(i,j,i_water)*rho(i_water)
            g2 = 0d0
            do k=1,n_spec-1
               g2 = g2+VH(i,j,k)*rho(k)
            enddo

            call particle_vol_hybrid(n_bin,TDV,n_spec,VH,i,j,pv)
               
            d= (6d0/pi*pv)**(1d0/3d0)                ! rn in m, d in m 
            write(6,*)'g1,g2,pv ',g1,g2,pv,d

            if (g1 .ne. 0. .or. g2.ne. 0.) then

               call rtnewt(x1,x2,xacc,x,d,g1,g2,pmv,p0T,RH,T,p)
               dmdt(i,j) = x
               histot(i,j) = dmdt(i,j)/(g1+g2)     ! histot is normalized growth in s-1 
               dmdt(i,j) = dmdt(i,j)/rho(i_water)  ! dmdt is now the volume growth in m^3 s-1

            endif
         enddo
      enddo

      do i=ia,ie
         do j=1,MH(i)
            write(66,*)'dmdt ', i,j,dmdt(i,j),histot(i,j),1/histot(i,j)
         enddo
      enddo

      STOP

      end

cn ********************************************************************************

      subroutine rtnewt(x1,x2,xacc,x,d,g1,g2,pmv,p0T,RH,T,p)
      integer jmax
      real*8 x,x1,x2,xacc 
      real*8 T,T_a,pmv,p0T,RH,p
      external funcd
      parameter (jmax = 400)
      integer j
      real*8 g1,g2
      real*8 df, dx, f, d

      x = 0.5d0 * (x1 + x2)
      write(6,*)'rtnewt1 ',x1,x2,xacc,x
      write(6,*)'rtnewt2 ',d,g1,g2,pmv
      write(6,*)'rtnewt3 ',p0T,RH,T,p

      do j=1,jmax
         call funcd(x,f,df,d,g1,g2,T_a,pmv,p0T,RH,T,p)
         dx=f/df
         x=x-dx
         if((x1 - x) * (x - x2) .lt. 0d0) then
            write(6,*)'ende bei ',g1,g2,RH,T
            pause 'rtnewt jumped out of brackets'
         endif
         if(abs(dx) .lt. xacc) then

         return
         endif
      enddo
      write(6,*)'ende bei ',g1,g2,RH,T
      pause 'rtnewt exceeded maximum iteration '

      end
 
cn **************************************************************************************
      subroutine funcd(x,f,df,d_p,g1,g2,T_a,pmv,p0T,RH,T,pin)

      real*8 x      ! INPUT: dm/dt
      real*8 f      ! OUTPUT: error
      real*8 df     ! OUTPUT: derivative of error
      real*8 d_p    ! diameter (m)
      real*8 g1
      real*8 g2
      real*8 T_a    ! OUTPUT: droplet temperature (K)
      real*8 pmv    ! vapor pressure of the ambient atmosphere (Pa)
      real*8 p0T    ! vapor pressure at temperature T  (Pa)
      real*8 RH     ! relative humidity
      real*8 T      ! ambient temperature (K)
      real*8 pin    ! ambient pressure (hPa)

C     parameters
      real*8 rho, rho_a, rho_n, M_w, M_a, M_s, sig, R, L_v, alpha
      real*8 p00, T0, cp, nu, eps, atm

      parameter (rho = 1000d0)     ! water density
      parameter (rho_a = 1.25d0)
      parameter (rho_n = 1800d0)
      parameter (M_w = 18d-3)      ! molecular weight of water
      parameter (M_a = 28d-3)
      parameter (M_s = 132d-3)     ! molecular weight of solute
      parameter (sig = 0.073d0)    ! surface energy
      parameter (R = 8.314d0)      ! universal gas constant
      parameter (L_v = 2.5d6)      ! latent heat
      parameter (alpha = 1d0)      ! the value 0.045 is also used sometimes
      parameter (p00 = 6.11d2)     ! vapor pressure at temp T = 273 K (Pa)
      parameter (T0 = 273.15d0)
      parameter (cp = 1005d0)
      parameter (nu = 3d0)         ! number of ions in the solute
      parameter (eps = 0.25d0)     ! solubility of aerosol material
      parameter (atm = 1013.25d0)  ! atmospheric pressure (hPa)

      real*8 pi
      parameter (pi = 3.14159265358979323846d0)

C     local variables
      real*8 D_v             !molecular diffusivity (m2/s)
      real*8 D_vp             !molecular diffusivity (m2/s)
      real*8 k_a, k_ap,k_ap1       ! thermal conductivity uncorrected and corrected (J m^{-1} s^{-1} K^{-1})
      real*8 rat,fact1,fact2,c1,c2,c3,c4,c5

cn ***  diffusion coefficient (basic)
      D_v = 0.211d0 / (pin / atm) * (T / 273d0)**1.94d0 ! in cm^2 s^{-1}

cn ***  diffusion coefficient (corrected for non-continuum effects)
cn      D_v_de = 1+(2*D_v*1.e-04/(alpha*d_p))*(2*pi*M_w/(R*T))**0.5
cn      D_vp = D_v / D_v_de

cn *** TEST: use the basic expression for D_vp
      D_vp = D_v * 1d-4         ! in m^2 s^{-1}
      
      k_a = 1d-3 * (4.39d0 + 0.071d0 * T)
      k_ap1 = 1d0 + 2d0 * k_a / (alpha * d_p * rho_a * cp)
     &     * (2d0 * pi * M_a / (R * T))**0.5d0
      k_ap = k_a / k_ap1
      
      rat = p0T / (R * T)
      fact1 = L_v * M_w / (R * T)
      fact2 = L_v / (2d0 * pi * d_p * k_ap * T)
      
      c1 = 2d0 * pi * d_p * D_vp * M_w * rat
      c2 = 4d0 * M_w * sig / (R * rho * d_p)
      c3 = c1 * fact1 * fact2
      c4 = L_v / (2d0 * pi * d_p * k_ap)
c      c5 = nu*eps*M_w*rho_n*r_n**3/(M_s*rho*((d_p/2)**3-r_n**3))
c      c5 = nu*eps*M_w/M_s * g2/g1
c ** corrected according to Jim's note:
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
      
      end
