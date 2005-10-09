C Golovin coagulation kernel.

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine kernel_golovin(a, b, k)

      real*8 a  ! INPUT: volume of first particle
      real*8 b  ! INPUT: volume of second particle
      real*8 k  ! OUTPUT: coagulation kernel

      real*8 beta_1
      parameter (beta_1 = 1000.)

      k = beta_1 * (a + b)

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine soln_golovin_exp(n_bin, vv, rr, dp, dlnr, time,
     &     N_0, V_0, rho_p, g, n_ln)

      integer n_bin       ! INPUT: number of bins
      real*8 vv(n_bin)    ! INPUT: volume of particles in bin
      real*8 rr(n_bin)    ! INPUT; radius of particles in bins
      real*8 dp(n_bin)    ! INPUT: diameter of particles in bins
      real*8 dlnr         ! INPUT: scale factor
      real*8 time         ! INPUT: current time
      real*8 N_0          ! INPUT: particle number concentration (#/m^3)
      real*8 V_0          ! INPUT:
      real*8 rho_p        ! INPUT: particle density (kg/m^3)
      real*8 g(n_bin)     ! OUTPUT: mass of particles in each bin
      real*8 n_ln(n_bin)  ! OUTPUT: number of particles in each bin

      real*8 d_0, beta_1, tau, T, rat_v, nn, b, x
      integer k

      real*8 pi
      parameter (pi = 3.14159265358979323846)

c      write(6,'(3a14)'), 'x', 'bessi1', 'besj0'
c      do k = 0,100
c         d_0 = real(k) / 10d0
c         call bessi1(d_0, b)
c         T = besy0(d_0)
c         write(6,'(3e14.6)'), d_0, b, T
c      enddo
c      call exit(2)

      d_0 = (6d0 * V_0 / pi)**(1d0/3d0)
      call kernel_golovin(1d0, 0d0, beta_1)

      if (time .eq. 0d0) then
         do k = 1,n_bin
            n_ln(k) = pi/2d0 * dp(k)**3d0 * 1d0/V_0 * exp(-(vv(k)/V_0))
         enddo
      else
         tau = N_0 * V_0 * beta_1 * time
         T = 1 - exp(-tau)
         do k = 1,n_bin
            rat_v = vv(k) / V_0
            x = 2d0 * rat_v * sqrt(T)
            if (x .lt. 500d0) then
               call bessi1(x, b)
            else
               b = 0d0
            endif
            nn = N_0/vv(k) * (1d0 - T) / sqrt(T)
     &           * exp(-((1d0 + T) * rat_v)) * b
            n_ln(k) = pi/2d0 * dp(k)**3d0 * nn
c            write(6,'(6a11)'), 'time', 'T', 'b', 'nn', 'n_ln', 'g'
c            write(6,'(6e11.3)'), time, T, b, nn, n_ln(k), g(k)
         enddo
      endif

      do k = 1,n_bin
         g(k) = pi/6d0 * rho_p * dp(k)**3d0 * n_ln(k)
      enddo

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine bessi1(x, r)
      ! bessel function

      real*8 x   ! INPUT: function argument
      real*8 r   ! OUTPUT: function value

      real*8 ax
      real*8 p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,q8,q9,y
      data p1,p2,p3,p4,p5,p6,p7/0.5d0,0.87890594d0,0.51498869d0,
     &     0.15084934d0,0.2658733d-1,0.301532d-2,0.32411d-3/
      data q1,q2,q3,q4,q5,q6,q7,q8,q9/0.39894228d0,-0.3988024d-1,
     &     -0.362018d-2,0.163801d-2,-0.1031555d-1,0.2282967d-1,
     &     -0.2895312d-1,0.1787654d-1,-0.420059d-2/
      
      if(abs(x).lt.3.75) then
         y=(x/3.75)**2
         r = x*(p1+y*(p2+y*(p3+y*(p4+y*(p5+y*(p6+y*p7))))))
      else
         ax = abs(x)
         y=3.75/ax
         r = (exp(ax)/sqrt(ax))*(q1+y*(q2+y*(q3+y*(q4+
     *        y*(q5+y*(q6+y*(q7+y*(q8+y*q9))))))))
         if (x.lt.0.) r = -r
      endif

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
