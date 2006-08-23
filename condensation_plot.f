C Condensation
C

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      program condensation_plot

      integer n_spec, n_bin
      integer i_water
      parameter (n_bin = 120, n_spec = 3)
      real*8 T, rho(n_spec), RH, pres, pmv, p0T
      real*8 dmdt(n_bin)        ! growth rate (kg s^{-1})
      real*8 histot(n_bin)      ! normalized growth rate: histot=dmdt/m (s^{-1})
      real*8 rn(n_bin),rw(n_bin),g(n_bin,2)
      real*8 p00, T0
      
      parameter (T = 290d0)     ! temperature of gas medium (K)
      parameter (RH = 1.03d0)   ! relative humidity (???)
      parameter (p00 = 611d0)   ! equilibrium water vapor pressure at 273 K (Pa)
      parameter (T0  = 273.15d0) ! (K)
      parameter (pres = 1d5)    ! ambient pressure (Pa)
      parameter (i_water = 3)

      integer i

      real*8 pi
      parameter (pi = 3.14159265358979323846d0)

      open (99,file = 'init_eq.d')
      open (66,file = 'dmdt10003c.d')
      do i=1,n_bin
         read(99,'(5e20.11)')rn(i),rw(i),g(i,1),g(i,2)
      enddo

      rho(1) = 1800d0
      rho(2) = 1800d0
      rho(3) = 1000d0

      p0T = p00 * 10d0**(7.45d0 * (T - T0) / (T - 38d0))
      write(6,*)'p0T ',p0T,p00,RH

! 
      call kond_plot(n_bin,n_spec,T, RH, pres,rn,rw,g
     $     ,i_water,rho, p0T,dmdt, histot)

C     dmdt(i) and histot(i) are output
C     dmdt is growth rate of one droplet in kg s^{-1}
C     histot = dmdt * 1.e6 / e(i) in s^{-1}

      write(6,*)'ende condensation'

      contains

cn ****************************************************************
 
      subroutine kond_plot(n_bin,n_spec,T,RH,p,rn,rw,g
     $     ,i_water,rho, p0T,dmdt,histot)

cn *** Calculation of the term dm/dt according to Majeed and Wexler, Atmos. Env. (2001)
cn *** Since Eq. (7) in this paper is an implicit equation (T_a depends on dm/dt), a Newton
cn *** solver is applied.

      use condensation

      integer n_bin         ! number of bins
      integer n_spec        ! number of species
      real*8 T              ! ambient temperature (K)
      real*8 RH             ! relative humidity (???)
      real*8 p              ! ambient pressure  (hPa)
      real*8 rn(n_bin)
      real*8 rw(n_bin)
      real*8 g(n_bin,2)
      integer i_water
      real*8 rho(n_spec)    ! density of each species (???)
      real*8 p0T
      real*8 dmdt(n_bin)    ! OUTPUT: growth rate (kg s^{-1}) 
      real*8 histot(n_bin)  ! OUTPUT: growth constant (s^{-1})

      real*8    x,d

      real*8    RR,M_w,sig_w,M_s
      real*8 pi,nu
      real*8    g1,g2
      real*8    gmin

      parameter (RR = 8.314d0, M_w = 18d-3)
      parameter (sig_w = 0.073d0,M_s = 132d-3)
      parameter (nu = 3d0)
      parameter (gmin = 0d0)

      parameter (pi = 3.14159265358979323846d0)

      integer i
      real*8 x_tol, f_tol, x1, x2

      do i=1,n_bin
            dmdt(i) = 0d0
      enddo

      x1 = 0d0
      x2 = 1000d-7
      x_tol = 1d-15
      f_tol = 1d-15

      g1 = 0d0
      g2 = 0d0

      do i=1,n_bin
            d = 2d0*rw(i)  
            g1 = g(i,1)
            g2 = g(i,2)

            call cond_newt(x1, x2, x_tol, f_tol, d, g1, g2, p0T,
     $           RH, T, p, x)
            dmdt(i) = x
            histot(i) = dmdt(i)/(g1+g2) ! histot is normalized growth in s-1 
            dmdt(i) = dmdt(i)/rho(i_water) ! dmdt is now the volume growth in m^3 s-1
      enddo

      do i=1,n_bin
         write(66,*)i,rn(i)*1d6,rw(i)*1d6,g(i,1)+g(i,2),
     &        dmdt(i),histot(i),1d0/histot(i)
      enddo

      STOP

      end subroutine

cn ********************************************************************************

      end program
