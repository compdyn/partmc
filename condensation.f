C Condensation
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine condensation()

      real*8 a  ! INPUT: volume of particle before condensation
      real*8 k  ! OUTPUT: volume of particle after condensation

      real*8 T, rho_p, RH, pres, pmv, p0T
      real*8 dmdt, histot
      real*8 p00, T0

      parameter ( T = 298. )     ! Temperature of gas medium in K
      parameter (rho_p = 1800. ) ! particle density in kg m-3
      parameter (p00 = 611.    ) ! equilibrium water vapor pressure at 273 K in Pa
      parameter (T0  = 273.15  ) ! in K
      parameter (pres = 1000.  ) ! ambient pressure in hPa

      real*8 pi
      parameter (pi = 3.14159265358979323846d0)

! r_h(i)   : radius of droplet (water + aerosol )in bin i (m)
! gg_h(i,2): mass of water (gg_h(i,1)) and aerosol (gg_h(i,2)) in bin i (kg m-3 = mg cm-3)
! e_h(i)   : mass of dry aerosol (mg)
! RH       : relative humidity

      pmv = p0T * RH
      p0T = p00 *10**(7.45*(T-T0)/(T-38.))

      call kond(r_h, gg_h, e_h, T, RH, pres,
     &     pmv, p0T,
     &     dmdt, histot, ij, i, i, k)

! dmdt(i) and histot(i) are output
! dmdt is growth rate of one droplet in kg s-1
! histot =  dmdt * 1.e6 / e(i) in s-1

      return
      end

cn ****************************************************************
 
      SUBROUTINE kond(r,gg,e,T,RH,p,
     &                pmv,p0T,dmdt,histot,ij,ia,ie,k)

cn *** Calculation of the term dm/dt according to Majeed and Wexler, Atmos. Env. (2001)
cn *** Since Eq. (7) in this paper is an implicit equation (T_a depends on dm/dt), a Newton
cn *** solver is applied.

      integer i,imax
      parameter(imax = 100)

      real*8    dlnr,dt
      real*8    r(imax)                          ! total radius of droplet/aerosol m
      real*8    rn(imax)                         ! radius of aerosol core          m
      real*8    xacc,x1,x2,hh(imax)
      real*8    gg(imax,2),g(imax),e(imax)       ! mass size distribution in kg m-3 bzw. mg cm-3
      real*8    numln(imax)                      ! number size distribution in cm-3 (ln-based)
      real*8    x,d,f,df
      real*8    dmdt(imax)                       ! growth rate in kg s-1
      real*8    dp(imax)                         ! diameter in m
      real*8    rdry                             ! radius of aerosol core in m
      real*8    pmv
      real*8    histot(imax)                     ! growth constant in s-1
      real*8    mv,ml                            ! water vapor concentration in kg m-3 
      real*8    p0T,hm,rmax,rh
      real*8    T                                ! ambient temperature in K
      real*8    p                                ! ambient pressure  in hPa
      real      RR,M_w,sig_w,M_s,rho_w,pi,nu
      real      t1,t2
      real*8    taudiff(imax)
      real*8    g1,g2
      parameter (RR=8.314, M_w = 18.e-03)
      parameter (sig_w = 0.073,M_s = 132.e-03, rho_w = 1000.)
      parameter (nu=3)
      parameter (gmin = 0.)

      parameter (pi = 3.14159265358979323846d0)

      do i=ia,ie
         dmdt(i) = 0.
      enddo

      x1 = -1000.e-7
      x2 = 1000.e-7
      xacc = 1.e-15

      do i=ia,ie

         g1 = gg(i,1)
         g2 = gg(i,2)
         d=2*r(i)   ! rn in m, d in m 

         if (g1 .ne. 0. .or. g2.ne. 0.) then
         call rtnewt(x1,x2,xacc,x,d,g1,g2,pmv,p0T,rh,T,p,ij,i,k)

         dmdt(i) = x

         endif

       enddo

cn ** Calculation of histot
       do i=ia,ie
          histot(i) = dmdt(i)*1.e06/e(i)
       enddo

      end

cn ********************************************************************************

      subroutine rtnewt(x1,x2,xacc,x,d,g1,g2,pmv,p0T,rh,T,p,ij,i,k)
      integer jmax
      real*8 x,x1,x2,xacc 
      real*8 T,T_a,pmv,p0T,rh,p
      external funcd
      parameter (jmax=400)
      integer j
      real*8 g1,g2
      real*8 df, dx, f, d, rdry

      x = .5*(x1+x2)

      do j=1,jmax
         call funcd(x,f,df,d,g1,g2,T_a,pmv,p0T,rh,T,p)
         dx=f/df
         x=x-dx
         if((x1-x)*(x-x2).lt.0) then
            write(6,*)'ende bei ',k,ij,i,d/2.,g1,g2,rh,T
            pause 'rtnewt jumped out of brackets'
         endif
         if(abs(dx).lt.xacc) then

         return
         endif
      enddo
      write(6,*)'ende bei ',k,ij,i,d/2.,g1,g2,rh,T
      pause 'rtnewt exceeded maximum iteration '

      end
 
cn **************************************************************************************
      subroutine funcd(x,f,df,d_p,g1,g2,T_a,pmv,p0T,rh,T,pin)

      real*8 d_p                   !diameter (m)
      real*8 D_v, D_vp             !molecular diffusivity (m2/s)
      real*8 M_w, R                !Molecular weight (water), universal gas constant
      real*8 sig, M_s              !surface energy, Molecular weight (solute)
      real*8 T                     !ambient temperature (K)
      real*8 T_a                   !droplet temperature (K)
      real*8 k_a, k_ap             !thermal conductivity uncorrected and corrected (J/(msK)
      real*8 rho, RH               !water density, relative humidity
      real*8 p0T_a                 !vapor pressure at temperature T_a (flat surface) (Pa)
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
      real*8 eps, nu,r_n
      real*8 f,df,x

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
