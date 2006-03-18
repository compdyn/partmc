      program coad1d

      integer n, scal, isw
      real*8 dt, rq0b, xmwb, eps, u, rho, emin, tmax
      parameter (n = 220)          ! number of bins
      parameter (scal = 4)         ! bin mesh scale factor
      parameter (isw = 1)          ! kernel (0 = long, 1 = hall, 2 = golovin)
      parameter (dt = 1d0)         ! timestep (s)
      parameter (rq0b = 10d0)      ! mode radius of init dist (um)
      parameter (xmwb = 4.1886d0)  ! total water content (g/m^3)
      parameter (eps = 100d0)      ! epsilon (cm^2/s^3)
      parameter (u = 2.5d0)        ! u' (m/s)
      parameter (rho = 1000d0)     ! particle density (FIXME: units?)
      parameter (emin = 1d-15)     ! 
      parameter (tmax = 600d0)     ! total simulation time (s)

      real*8 dlnr, ax
      real*8 c(n,n), ima(n,n)
      real*8 g(n), r(n), e(n)
      real*8 ck(n,n), ec(n,n)
      real*8 gout(3600,n), hout(3600,n)
      real*8 gin(n), hin(n), bin(n), bout(3600,n)
      real*8 taug(n), taup(n), taul(n), tauu(n)
      real*8 prod(n), ploss(n)

      integer i, j, nt, lmin, ij
      real*8 xn0, xn1, v0, x0, x1
      real*8 tlmin, t, tau, rq0, xmw

      real*8 pi
      parameter (pi = 3.141592654)

C     g   : spectral mass distribution (mg/cm**3)
C     e   : droplet mass grid (mg)
C     r   : droplet radius grid (um)
C     dlnr: constant grid distance of logarithmic grid 
C     xn0 : mean initial droplet mass (mg)
C     xn1 : total initial droplet number concentration (1/cm^3)
C     ax  : growth factor for consecutive masses

C     mass and radius grid
      dlnr = log(2d0) / (3d0 * scal)
      ax = 2d0**(1d0 / scal)
      v0 = 4d0 / 3d0 * pi * rq0**3
      do i = 1,n
         e(i) = emin * 0.5d0 * (ax + 1d0) * ax**(i - 1)
         r(i) = 1000d0 * exp(log(3d0 * e(i) / (4d0 * pi)) / 3d0)
      enddo

C     initial mass distribution
      rq0 = rq0b * 1d-4
      xmw = xmwb * 1d-3
      xn0 = 4d0/3d0 * pi * 1000d0 * exp(log(rq0)*3d0)
      xn1 = xmw / xn0
      x0 = xn1 / xn0
      do i = 1,n
         x1 = e(i)
         g(i) = 3d0 * x1**2 * x0 * exp(-(x1 / xn0))
         if (g(i) .lt. 1d-30) then
            gin(i) = 0d0
         else
         gin(i) = g(i)
         endif
         write(6,*)'init ',i,r(i),g(i)
      enddo
      
      call courant(n, rq0, dlnr, scal, ax, c, ima, g, r, e)
      call trkern (isw, eps, u, n, g, r, e, ck, ec)

C     nt: number of iterations
      nt = int(tmax / dt) + 1
C     multiply kernel with constant timestep and logarithmic grid distance
      do i = 1,n
         do j = 1,n
            ck(i,j) = ck(i,j) * dt * dlnr
         enddo
      enddo
C     time integration
      tlmin = 1d-6
      t = 1d-6
      lmin = 0
      do ij = 2,nt
         t = t + dt
         tlmin = tlmin + dt
         write(6,*)'jetzt kommt coad ',t
C     collision
         call coad(n, dt, taug, taup, taul, tauu, prod, ploss,
     &        c, ima, g, r, e, ck, ec)

C     output for plotting
         if (tlmin .ge. 60d0) then
            tlmin = tlmin - 60d0
            lmin = lmin + 1
            print *,'time in minutes:',lmin
            tau = xn1 * 1500 * v0*t
            write(6,*)'time ', t, tau
            do i = 1,n
               hin(i) = 3d15 * gin(i) / (4d0 * pi * r(i)**3)     ! m^{-3}
               hout(lmin,i) = 3d15 * g(i) / (4d0 * pi * r(i)**3) ! m^{-3}
               gout(lmin,i) = g(i)
               bin(i) = 3d12 * gin(i) / (4d0 * 1000d0 * pi * r(i)**2)
               bout(lmin,i) = 3d12 * gout(lmin,i)
     &              / (4d0 * 1000 * pi * r(i)**2)
            enddo

            call do_mass_balance(n, g, r, dlnr)
         endif
      enddo

      open(15, file = 'sectionm')
      open(25, file = 'sectionn')
      do i = 1,n
         write(15,'(62e14.5)') r(i), gin(i), (gout(j,i), j = 1, lmin)
         write(25,'(62e14.5)') r(i), hin(i), (hout(j,i), j = 1, lmin)
      enddo

      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine coad(n, dt, taug, taup, taul, tauu, prod, ploss,
     &     c, ima, g, r, e, ck, ec)

C     collision subroutine, exponential approach

      integer n
      real*8 dt
      real*8 taug(n)
      real*8 taup(n)
      real*8 taul(n)
      real*8 tauu(n)
      real*8 prod(n)
      real*8 ploss(n)
      real*8 c(n,n)
      real*8 ima(n,n)
      real*8 g(n)
      real*8 r(n)
      real*8 e(n)
      real*8 ck(n,n)
      real*8 ec(n,n)

      real*8 gmin
      parameter (gmin = 1d-60)

      integer i, i0, i1, j, k, kp
      real*8 x0, gsi, gsj, gsk, gk, x1, flux

      do i = 1,n
         prod(i) = 0d0
         ploss(i) = 0d0
      enddo
      
C     lower and upper integration limit i0,i1
      do i0 = 1,(n - 1)
         if (g(i0) .gt. gmin) goto 2000
      enddo
 2000 continue
      do i1 = (n - 1),1,-1
         if (g(i1) .gt. gmin) goto 2010
      enddo
 2010 continue
      
      do i = i0,i1
         do j = i,i1
            k = ima(i,j)
            kp = k + 1
            
            x0 = ck(i,j) * g(i) * g(j)
            x0 = min(x0, g(i) * e(j))
            
            if (j .ne. k) x0 = min(x0, g(j) * e(i))
            gsi = x0 / e(j)
            gsj = x0 / e(i)
            gsk = gsi + gsj
            
C     loss from positions i, j
            ploss(i) = ploss(i) + gsi
            ploss(j) = ploss(j) + gsj
            g(i) = g(i) - gsi
            g(j) = g(j) - gsj
            gk = g(k) + gsk
            
            if (gk .gt. gmin) then
               x1 = log(g(kp) / gk + 1d-60)
               flux = gsk / x1 * (exp(0.5d0 * x1)
     &              - exp(x1 * (0.5d0 - c(i,j))))
               flux = min(flux, gk)
               g(k) = gk - flux
               g(kp) = g(kp) + flux
C     gain for positions i, j
               prod(k) =  prod(k) + gsk - flux           
               prod(kp) = prod(kp) + flux
            endif
         enddo
      enddo
      
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine courant(n, rq0, dlnr, scal, ax, c, ima, g, r, e)

      integer n
      real*8 rq0
      real*8 dlnr
      integer scal
      real*8 ax
      real*8 c(n,n)
      real*8 ima(n,n)
      real*8 g(n)
      real*8 r(n)
      real*8 e(n)

      integer i, j, k, kk
      real*8 x0

      do i = 1,n
         do j = i,n
            x0 = e(i) + e(j)
            do k = j,n
               if ((e(k) .ge. x0) .and. (e(k-1) .lt. x0)) then
                  if (c(i,j) .lt. 1d0 - 1d-08) then
                     kk = k - 1
                     c(i,j) = log(x0 / e(k-1)) / (3d0 * dlnr)
                  else
                     c(i,j) = 0d0
                     kk = k
                  endif
                  ima(i,j) = min(n - 1, kk)
                  goto 2000
               endif
            enddo
 2000       continue
            c(j,i) = c(i,j)
            ima(j,i) = ima(i,j)
         enddo
      enddo
      
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine trkern(isw, eps, u, n, g, r, e, ck, ec)

      integer isw
      real*8 eps
      real*8 u
      integer n
      real*8 g(n)
      real*8 r(n)
      real*8 e(n)
      real*8 ck(n,n)
      real*8 ec(n,n)

      real*8 winf(n)
      real*8 rr(n)
      real*8 cck(n,n)

      integer i, j, im, jm, ip, jp
      real*8 effi, r_tmp, w_tmp, r1, r2, ec_tmp
      real*8 winf1, winf2

      real*8 pi
      parameter (pi = 3.141592654)

      if (isw .eq. 0) then
C long kernel
         do j = 1,n
            r_tmp = r(j) / 1d6
            call fall_g(r_tmp, w_tmp)
            winf(j) = w_tmp * 100d0
            rr(j) = r(j) * 1d-4
         enddo

         do j = 1,n
            do i = 1,j
               if(r(j) .le. 50d0) then
                  effi = 4.5d0 - 4d0 * r(j) * r(j)
     &                 * (1d0 - 3d0 / (max(3d0, dble(r(i))) + 1d-2))
               else
                  effi = 1d0
               endif
               cck(j,i) = pi * (rr(j) + rr(i)) * (rr(j) + rr(i)) * effi
     &              * abs(winf(j) - winf(i))
               cck(i,j) = cck(j,i)
            enddo
         enddo
      elseif (isw .eq. 1) then
C     hall kernel
         do i = 1,n
            do j = 1,n
               r1 = r(i) ! (um)
               r2 = r(j) ! (um)
               call fall_g(r1 / 1d6, winf1) ! winf1 in m/s
               call fall_g(r2 / 1d6, winf2) ! winf2 in m/s
               call effic(r1, r2, ec_tmp) ! ec_tmp dimensionless
               cck(i,j) = 1d-6 * pi * (r1 + r2)**2 * ec_tmp
     &              * abs(winf1 - winf2)
               ! cck in cm^3/s = 1d6 * m^3/s
            enddo
         enddo
      elseif (isw .eq. 2) then
C     golovin kernel
         do j = 1,n
            r_tmp = r(j) / 1d6
            call fall_g(r_tmp, w_tmp)
            winf(j) = w_tmp * 100d0
            rr(j) = r(j) * 1d-4
         enddo

         do j = 1,n
            do i = 1,j
               cck(j,i) = 1d0 * (e(j) + e(i))
               cck(i,j) = cck(j,i)
            enddo
         enddo
      endif

C     two-dimensional linear interpolation of kernel
      do i = 1,n
         do j = 1,n
            jm = max0(j - 1, 1)
            im = max0(i - 1, 1)
            jp = min0(j + 1, n)
            ip = min0(i + 1, n)
            ck(i,j) = 0.125d0 * (cck(i,jm) + cck(im,j)
     &           + cck(ip,j) + cck(i,jp))
     &           + 0.5d0 * cck(i,j)
            if (i .eq. j) ck(i,j) =0.5d0 * ck(i,j)
         enddo
      enddo

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine do_mass_balance(n, g, r, dlnr)

      integer n
      real*8 g(n)
      real*8 r(n)
      real*8 dlnr

      real*8 x0, x1, x_num, hh
      integer i, imax

      real*8 pi
      parameter (pi = 3.141592654)

      x0 = 0d0
      x1 = 0d0
      x_num = 0d0
      do i = 1,n
         x0 = x0 + g(i) * dlnr
         x1 = max(x1, g(i))
         if (abs(x1 - g(i)) .lt. 1d-9) then
            imax = i
         endif
         hh = 3d+15 * g(i) / (4d0 * pi * r(i)**3) * dlnr
         x_num = x_num + hh
      enddo
      
      write(6,*)'mass ', x0, 'max ', x1, 'imax ', imax
      write(6,*)'x_num ', x_num

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
