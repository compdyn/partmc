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
      
      parameter (T = 298. )      ! Temperature of gas medium in K
      parameter (RH= 1.01 )       ! Relative humidity
      parameter (p00 = 611.    ) ! equilibrium water vapor pressure at 273 K in Pa
      parameter (T0  = 273.15  ) ! in K
      parameter (pres = 1000.  ) ! ambient pressure in hPa
      parameter (i_water = 3)

      real*8 pi
      parameter (pi = 3.14159265358979323846d0)

      write(6,*)'in condensation ',n_bin,TDV,n_spec
      stop
      rho(1) = 1800.
      rho(2) = 1800.
      rho(3) = 1000.

      p0T = p00 *10**(7.45*(T-T0)/(T-38.))
      pmv = p0T * RH

  
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

      parameter (RR=8.314, M_w = 18.e-03)
      parameter (sig_w = 0.073,M_s = 132.e-03)
      parameter (nu=3)
      parameter (gmin = 0.)

      parameter (pi = 3.14159265358979323846d0)

      write(6,*)'in kond'

      stop

      do i=ia,ie
         do j=1,TDV
            dmdt(i,j) = 0.
         enddo
      enddo

      x1 = -1000.e-7
      x2 = 1000.e-7
      xacc = 1.e-15

      g1 = 0.
      g2 = 0.

      do i=ia,ie
         do j=1,MH(i)
          
            g1 = VH(i,j,i_water)*rho(i_water)

            do k=1,n_spec-1
               g2 = g2+VH(i,j,k)*rho(k)
            enddo

            call particle_vol_hybrid(n_bin,TDV,n_spec,VH,j,i,pv)
               
            d= (6/pi*pv)**(1./3.)                ! rn in m, d in m 

            if (g1 .ne. 0. .or. g2.ne. 0.) then

               call rtnewt(x1,x2,xacc,x,d,g1,g2,pmv,p0T,RH,T,p)
               dmdt(i,j) = x

            endif
         enddo
      enddo

cn ** Calculation of histot
      do i=ia,ie
         do j=1,MH(i)
c            histot(i,j) = dmdt(i,j)*1.e06/e(i,j)
            write(6,*)'dmdt ', i,j,dmdt(i,j)
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
      parameter (jmax=400)
      integer j
      real*8 g1,g2
      real*8 df, dx, f, d

      x = .5*(x1+x2)

      do j=1,jmax
         call funcd(x,f,df,d,g1,g2,T_a,pmv,p0T,RH,T,p)
         dx=f/df
         x=x-dx
         if((x1-x)*(x-x2).lt.0) then
            write(6,*)'ende bei ',g1,g2,RH,T
            pause 'rtnewt jumped out of brackets'
         endif
         if(abs(dx).lt.xacc) then

         return
         endif
      enddo
      write(6,*)'ende bei ',g1,g2,RH,T
      pause 'rtnewt exceeded maximum iteration '

      end
 
cn **************************************************************************************
      subroutine funcd(x,f,df,d_p,g1,g2,T_a,pmv,p0T,RH,T,pin)

      real*8 d_p                   !diameter (m)
      real*8 D_v, D_vp             !molecular diffusivity (m2/s)
      real*8 M_w, R                !Molecular weight (water), universal gas constant
      real*8 M_a
      real*8 sig, M_s              !surface energy, Molecular weight (solute)
      real*8 T0
      real*8 T                     !ambient temperature (K)
      real*8 T_a                   !droplet temperature (K)
      real*8 k_a, k_ap,k_ap1       !thermal conductivity uncorrected and corrected (J/(msK)
      real*8 rho, RH               !water density, relative humidity
      real*8 p0T                   !vapor pressure at temperature T  (Pa)
      real*8 p00                   !vapor pressure at T=273 K (Pa)
      real*8 p                     !ambient pressure (atm)
      real*8 pin                   !ambient pressure (hPa)
      real*8 pmv                   !vapor pressure of the ambient atmosphere (Pa)
      real*8 L_v                   !Latent heat
      real*8 pi
      real*8 alpha
      real*8 rho_a, rho_n, g1, g2
      real*8 cp
      real*8 eps, nu
      real*8 f,df,x

      real*8 rat,fact1,fact2,c1,c2,c3,c4,c5
      parameter (rho = 1000., rho_a = 1.25, rho_n=1800.)
      parameter (M_w = 18.*1.e-03, M_a = 28.*1.e-3, M_s=132.*1.e-03)
      parameter (sig = 0.073) 
      parameter (R = 8.314)
      parameter (L_v = 2.5e+6)
      parameter (alpha = 1.)                   ! the value 0.045 is also used sometimes ...
      parameter (p00 = 6.11*1.e2, T0=273.15)
      parameter (cp = 1005)
      parameter (nu = 3)                       ! number of ions in the solute
      parameter (eps = 0.25)                   ! solubility of aerosol material

      parameter (pi = 3.14159265358979323846d0)
cn      parameter (eps = 0.9)

cn ***  conversion hPa in atm

      p = pin/1013.25

cn ***  diffusion coefficient (basic)
      D_v = 0.211/p * (T/273.)**1.94 ! in cm**2/s**-1

cn ***  diffusion coefficient (corrected for non-continuum effects)
cn      D_v_de = 1+(2*D_v*1.e-04/(alpha*d_p))*(2*pi*M_w/(R*T))**0.5
cn      D_vp = D_v / D_v_de

cn *** TEST: use the basic expression for D_vp
      D_vp = D_v*1.e-04               ! in m**2/s**-1
 
      k_a = 1.e-03*(4.39+0.071*T)
      k_ap1=1+2*k_a/(alpha*d_p*rho_a*cp)*(2*pi*M_a/(R*T))**0.5
      k_ap=k_a/k_ap1

      rat = p0T/(R*T)
      fact1 = L_v*M_w/(R*T)
      fact2 = L_v/(2*pi*d_p*k_ap*T)
     
      c1 = 2*pi*d_p*D_vp*M_w*rat * 1. 
      c2 = 4*M_w*sig/(R*rho*d_p)
      c3 = c1 * fact1 * fact2
      c4 = L_v/(2*pi*d_p*k_ap)
c      c5 = nu*eps*M_w*rho_n*r_n**3/(M_s*rho*((d_p/2)**3-r_n**3))

      c5 = nu*eps*M_w/M_s * g2/g1


      f = x- c1*(RH-exp(c2/(T+c4*x)-c5))/(1+c3*exp(c2/(T+c4*x)-c5))

      df = 1+c1*RH*(1+c3*exp(c2/(T+c4*x)-c5))**(-2)
     &       *c3*exp(c2/(T+c4*x)-c5)*(-1)*c2*c4/(T+c4*x)**2
     &       +c1*(exp(c2/(T+c4*x)-c5)*(-1)*c2*c4/(T+c4*x)**2
     &                      *(1+c3*exp(c2/(T+c4*x)-c5))**(-1)
     &       +exp(c2/(T+c4*x)-c5)*(-1)*(1+c3*exp(c2/(T+c4*x)-c5))**(-2)
     &        *c3*exp(c2/(T+c4*x)-c5)*(-1)*c2*c4/(T+c4*x)**2)

      T_a = T + c4 * x
          
         end
