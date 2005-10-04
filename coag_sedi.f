C coag_sedi.f
C
C sedimentation coagulation kernel

C &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

      subroutine coag_kernel(a, b, k)

      real*8 a  ! INPUT: size of first particle
      real*8 b  ! INPUT: size of second particle
      real*8 k  ! OUTPUT: kernel k(a,b)

      real*8 pi,const,onethird
      real*8 r1,r2,winf1,winf2,ec
      parameter (pi = 3.14159265358979323846)

      const = 3./(4.*pi)
      onethird  = 1./3.
      r1 = (const*a)**onethird
      r2 = (const*b)**onethird
      call fallg(r1,winf1)
      call fallg(r2,winf2)
      call effic(r1,r2,ec)
      k = ec *pi* (r1+r2)*(r1+r2)*abs(winf1-winf2)
      return
      end

C &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

      subroutine fallg(r,winf)

c terminal velocity of falling drops
      real*8 r, winf
      real*8 b(7),c(6)
      data b /-0.318657e1,0.992696,-0.153193e-2,-0.987059e-3,
     &        -0.578878e-3,0.855176e-4,-0.327815e-5/
      data c /-0.500015e1,0.523778e1,-0.204914e1,0.475294,-0.542819e-1,
     &         0.238449e-2/
      real*8 eta, xlamb, rhow, rhoa, grav, cunh, t0, sigma
      real*8 stok, stb, phy, py, rr, x, y, xrey, bond
      integer i

      eta=1.818e-4
      xlamb=6.62e-6
      rhow=1.
      rhoa=1.225e-3
      grav=980.665
      cunh=1.257*xlamb
      t0=273.15
      sigma=76.1-0.155*(293.15-t0)
      stok=2.*grav*(rhow-rhoa)/(9.*eta)
      stb=32.*rhoa*(rhow-rhoa)*grav/(3.*eta*eta)
      phy=sigma*sigma*sigma*rhoa*rhoa/((eta**4)*grav*(rhow-rhoa))
      py=phy**(1./6.)

c rr: radius in cm-units
      rr=r*1.e+02

      if (rr.le.1.e-3) then
         winf=stok*(rr*rr+cunh*rr)
      elseif (rr.gt.1.e-3.and.rr.le.5.35e-2) then
         x=log(stb*rr*rr*rr)
         y=0.
         do i=1,7
            y=y+b(i)*(x**(i-1))
         enddo
         xrey=(1.+cunh/rr)*exp(y)
         winf=xrey*eta/(2.*rhoa*rr)
      elseif (rr.gt.5.35e-2) then
         bond=grav*(rhow-rhoa)*rr*rr/sigma
         if (rr.gt.0.35) bond=grav*(rhow-rhoa)*0.35*0.35/sigma
         x=log(16.*bond*py/3.)
         y=0.
         do i=1,6
            y=y+c(i)*(x**(i-1))
         enddo
         xrey=py*exp(y)
         winf=xrey*eta/(2.*rhoa*rr)
         if (rr.gt.0.35)  winf=xrey*eta/(2.*rhoa*0.35)
      endif
      winf = winf/100.
      
      return
      end

C &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

      subroutine effic(r1,r2,ec)
      real*8 r1  ! INPUT: radius of first particle
      real*8 r2  ! INPUT: radius of second particle
      real*8 ec  ! OUTPUT: collision efficiency

C     collision efficiencies of hall kernel
      real*8 rat(21),r0(15),ecoll(15,21)
      data r0 /6.,8.,10.,15.,20.,25.,30.,40.,50.,
     &     60.,70.,100.,150.,200.,300./
      data rat /0.,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,
     &     0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0/
      data ecoll /
     &     0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001
     &     ,0.001,0.001,0.001,0.001,0.001,0.003,0.003,0.003,0.004,0.005
     &     ,0.005,0.005,0.010,0.100,0.050,0.200,0.500,0.770,0.870,0.970
     &     ,0.007,0.007,0.007,0.008,0.009,0.010,0.010,0.070,0.400,0.430
     &     ,0.580,0.790,0.930,0.960,1.000,0.009,0.009,0.009,0.012,0.015
     &     ,0.010,0.020,0.280,0.600,0.640,0.750,0.910,0.970,0.980,1.000
     &     ,0.014,0.014,0.014,0.015,0.016,0.030,0.060,0.500,0.700,0.770
     &     ,0.840,0.950,0.970,1.000,1.000,0.017,0.017,0.017,0.020,0.022
     &     ,0.060,0.100,0.620,0.780,0.840,0.880,0.950,1.000,1.000,1.000
     &     ,0.030,0.030,0.024,0.022,0.032,0.062,0.200,0.680,0.830,0.870
     &     ,0.900,0.950,1.000,1.000,1.000,0.025,0.025,0.025,0.036,0.043
     &     ,0.130,0.270,0.740,0.860,0.890,0.920,1.000,1.000,1.000,1.000
     &     ,0.027,0.027,0.027,0.040,0.052,0.200,0.400,0.780,0.880,0.900
     &     ,0.940,1.000,1.000,1.000,1.000,0.030,0.030,0.030,0.047,0.064
     &     ,0.250,0.500,0.800,0.900,0.910,0.950,1.000,1.000,1.000,1.000
     &     ,0.040,0.040,0.033,0.037,0.068,0.240,0.550,0.800,0.900,0.910
     &     ,0.950,1.000,1.000,1.000,1.000,0.035,0.035,0.035,0.055,0.079
     &     ,0.290,0.580,0.800,0.900,0.910,0.950,1.000,1.000,1.000,1.000
     &     ,0.037,0.037,0.037,0.062,0.082,0.290,0.590,0.780,0.900,0.910
     &     ,0.950,1.000,1.000,1.000,1.000,0.037,0.037,0.037,0.060,0.080
     &     ,0.290,0.580,0.770,0.890,0.910,0.950,1.000,1.000,1.000,1.000
     &     ,0.037,0.037,0.037,0.041,0.075,0.250,0.540,0.760,0.880,0.920
     &     ,0.950,1.000,1.000,1.000,1.000,0.037,0.037,0.037,0.052,0.067
     &     ,0.250,0.510,0.770,0.880,0.930,0.970,1.000,1.000,1.000,1.000
     &     ,0.037,0.037,0.037,0.047,0.057,0.250,0.490,0.770,0.890,0.950
     &     ,1.000,1.000,1.000,1.000,1.000,0.036,0.036,0.036,0.042,0.048
     &     ,0.230,0.470,0.780,0.920,1.000,1.020,1.020,1.020,1.020,1.020
     &     ,0.040,0.040,0.035,0.033,0.040,0.112,0.450,0.790,1.010,1.030
     &     ,1.040,1.040,1.040,1.040,1.040,0.033,0.033,0.033,0.033,0.033
     &     ,0.119,0.470,0.950,1.300,1.700,2.300,2.300,2.300,2.300,2.300
     &     ,0.027,0.027,0.027,0.027,0.027,0.125,0.520,1.400,2.300,3.000
     &     ,4.000,4.000,4.000,4.000,4.000/
c     two-dimensional linear interpolation of the collision efficiency
      real*8 rq, r_help, p, q
      integer k, ir, kk, iq, ek

      rq = r1/r2
      if (rq .gt. 1.0) then
         r_help = r1
         r1 = r2
         r2 = r_help
         rq = 1 / rq
      endif

      do k = 2, 15
         if ((r2 .le. r0(k)) .and. (r2 .ge. r0(k-1))) then
            ir = k
         elseif (r2 .gt. r0(15)) then
            ir = 16
         elseif (r2 .lt. r0(1)) then
            ir = 1
         endif
      enddo
      
      do kk=2,21
         if (rq.le.rat(kk).and.rq.gt.rat(kk-1)) iq=kk
      enddo
      
      if (ir.lt.16) then
         if (ir.ge.2) then
            p=(r2-r0(ir-1))/(r0(ir)-r0(ir-1))
            q=(rq-rat(iq-1))/(rat(iq)-rat(iq-1))
            ec=(1.-p)*(1.-q)*ecoll(ir-1,iq-1)+
     &           p*(1.-q)*ecoll(ir,iq-1)+
     &           q*(1.-p)*ecoll(ir-1,iq)+
     &           p*q*ecoll(ir,iq)
         else
            q=(rq-rat(iq-1))/(rat(iq)-rat(iq-1))
            ec=(1.-q)*ecoll(1,iq-1)+q*ecoll(1,iq)
         endif
      else
         q=(rq-rat(iq-1))/(rat(iq)-rat(iq-1))
         ek=(1.-q)*ecoll(15,iq-1)+q*ecoll(15,iq)
         ec=min(ek,1.d0)
      endif
      if (ec.lt.1.e-20) stop 99
      return
      end

