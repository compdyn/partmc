C Sedimentation coagulation kernel.

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine kernel_sedi(v1, v2, k)

      real*8 v1 ! INPUT: volume of first particle (m^3)
      real*8 v2 ! INPUT: volume of second particle (m^3)
      real*8 k  ! OUTPUT: kernel k(a,b) (m^3/s)

      real*8 const, onethird
      real*8 r1, r2, winf1, winf2, ec

      real*8 pi
      parameter (pi = 3.14159265358979323846d0)

      const = 3d0 / (4d0 * pi)
      onethird  = 1d0/3d0
      r1 = (const*v1)**onethird ! m
      r2 = (const*v2)**onethird ! m
      call fall_g(r1, winf1) ! winf1 in m/s
      call fall_g(r2, winf2) ! winf2 in m/s
      call effic(r1 * 1d6, r2 * 1d6, ec) ! ec is dimensionless
      k = ec * pi * (r1 + r2)**2 * abs(winf1 - winf2) 
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine fall_g(r, w_inf)

      real*8 r       ! INPUT: particle radius (m)
      real*8 w_inf   ! OUTPUT: terminal velocity (m/s)

c terminal velocity of falling drops
      real*8 eta, xlamb, rhow, rhoa, grav, cunh, t0, sigma
      real*8 stok, stb, phy, py, rr, x, y, xrey, bond
      integer i
      real*8 b(7),c(6)
      data b /-0.318657d1,0.992696d0,-0.153193d-2,-0.987059d-3,
     &        -0.578878d-3,0.855176d-4,-0.327815d-5/
      data c /-0.500015d1,0.523778d1,-0.204914d1,0.475294d0,
     &         -0.542819d-1,0.238449d-2/

      eta = 1.818d-4
      xlamb = 6.62d-6
      rhow = 1d0
      rhoa = 1.225d-3
      grav = 980.665d0
      cunh = 1.257d0 * xlamb
      t0 = 273.15d0
      sigma = 76.1d0 - 0.155d0 * (293.15d0 - t0)
      stok = 2d0 * grav * (rhow - rhoa) / (9d0 * eta)
      stb = 32d0 * rhoa * (rhow - rhoa) * grav / (3d0 * eta * eta)
      phy = sigma * sigma * sigma * rhoa * rhoa 
     &     / (eta**4 * grav * (rhow - rhoa))
      py = phy**(1d0/6d0)

c rr: radius in cm-units
      rr = r * 1d2

      if (rr .le. 1d-3) then
         w_inf = stok * (rr * rr + cunh * rr)
      elseif (rr .gt. 1d-3 .and. rr .le. 5.35d-2) then
         x = log(stb * rr * rr * rr)
         y = 0d0
         do i = 1,7
            y = y + b(i) * (x**(i - 1))
         enddo
         xrey = (1d0 + cunh/rr) * exp(y)
         w_inf = xrey * eta / (2d0 * rhoa * rr)
      elseif (rr .gt. 5.35d-2) then
         bond = grav * (rhow - rhoa) * rr**2 / sigma
         if (rr .gt. 0.35d0) then
            bond = grav * (rhow - rhoa) * 0.35d0**2 / sigma
         endif
         x = log(16d0 * bond * py / 3d0)
         y = 0d0
         do i = 1,6
            y = y + c(i) * (x**(i - 1))
         enddo
         xrey = py * exp(y)
         w_inf = xrey * eta / (2d0 * rhoa * rr)
         if (rr .gt. 0.35d0) then
            w_inf = xrey * eta / (2d0 * rhoa * 0.35d0)
         endif
      endif
      w_inf = w_inf / 100d0
      
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine effic(r1, r2, ec)
      real*8 r1  ! INPUT: radius of first particle (um)
      real*8 r2  ! INPUT: radius of second particle (um)
      real*8 ec  ! OUTPUT: collision efficiency (dimensionless)

      real*8 r_small, r_big, rq, p, q, ek
      integer k, ir, kk, iq
C     collision efficiencies of hall kernel
      real*8 rat(21),r0(15),ecoll(15,21)
      data r0 /6.,8.,10.,15.,20.,25.,30.,40.,50.,
     &     60.,70.,100.,150.,200.,300./
      data rat /0.,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,
     &     0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0/
c     two-dimensional linear interpolation of the collision efficiency
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

      r_small = min(r1, r2)
      r_big = max(r1, r2)
      rq = r_small / r_big

      ir = 1
      do k = 1, 15
         if (r_big .gt. r0(k)) then
            ir = k + 1
         endif
      enddo

      iq = 1
      do kk = 1,21
         if (rq .gt. rat(kk)) then
            iq = kk + 1
         endif
      enddo
      
      if (ir .lt. 16) then
         if (ir .ge. 2) then
            p = (r_big - r0(ir - 1)) / (r0(ir) - r0(ir - 1))
            q = (rq - rat(iq - 1)) / (rat(iq) - rat(iq - 1))
            ec = (1d0 - p) * (1d0 - q) * ecoll(ir - 1, iq - 1)
     &           + p * (1d0 - q) * ecoll(ir, iq - 1)
     &           + q * (1d0 - p) * ecoll(ir - 1, iq)
     &           + p * q * ecoll(ir, iq)
         else
            q = (rq - rat(iq - 1)) / (rat(iq) - rat(iq - 1))
            ec = (1d0 - q) * ecoll(1, iq - 1) + q * ecoll(1, iq)
         endif
      else
         q = (rq - rat(iq - 1)) / (rat(iq) - rat(iq - 1))
         ek = (1d0 - q) * ecoll(15, iq - 1) + q * ecoll(15, iq)
         ec = min(ek, 1d0)
      endif

      if (ec .lt. 1d-20) stop 99

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
