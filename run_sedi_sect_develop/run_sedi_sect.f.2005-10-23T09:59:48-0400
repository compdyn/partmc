      program coad1d

      integer n, scal, isw
      real*8 dt, rq0b, xmwb, eps, u, rho, emin, tmax
      parameter (n = 160)          ! number of bins
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

c g   : spectral mass distribution (mg/cm**3)
c e   : droplet mass grid (mg)
c r   : droplet radius grid (um)
c dlnr: constant grid distance of logarithmic grid 
c xn0 : mean initial droplet mass (mg)
c xn1 : total initial droplet number concentration (1/cm^3)
c ax  : growth factor for consecutive masses

c mass and radius grid
      dlnr=dlog(2.d0)/(3.*scal)
      ax=2.d0**(1.0/scal)
      v0 = 4.*pi/3. * rq0**3.
      e(1)=emin*0.5*(ax+1.)
      r(1)=1000.*dexp(dlog(3.*e(1)/(4.*pi))/3.)
      do i=2,n
         e(i)=ax*e(i-1)
         r(i)=1000.*dexp(dlog(3.*e(i)/(4.*pi))/3.)
      enddo

c initial mass distribution
      rq0=rq0b*1.e-04
      xmw=xmwb*1.d-3
      xn0=4./3.*pi*1000.*exp(log(rq0)*3.)
      xn1=xmw/xn0
      x0=xn1/xn0
      do i=1,n
         x1=e(i)
         g(i)=3.*x1*x1*x0*dexp(-(x1/xn0))
         if (g(i).lt. 1.e-30) then
            gin(i)=0.
         else
         gin(i) = g(i)
         endif
         write(6,*)'init ',i,r(i),g(i)
      enddo
      
      call courant(n, rq0, dlnr, scal, ax, c, ima, g, r, e)
      call trkern (isw, eps, u, n, g, r, e, ck, ec)

c nt: number of iterations
      nt = int(tmax / dt) + 1
c multiply kernel with constant timestep and logarithmic grid distance
      do i=1,n
         do j=1,n
            ck(i,j)=ck(i,j)*dt*dlnr
         enddo
      enddo
c time integration
      tlmin=1.d-6
      t=1.d-6
      lmin=0
      do ij=2,nt
         t=t+dt
         tlmin=tlmin+dt
         write(6,*)'jetzt kommt coad ',t
c collision
         call coad(n, dt, taug, taup, taul, tauu, prod, ploss,
     &        c, ima, g, r, e, ck, ec)

c output for plotting
         if (tlmin.ge.60.) then
            tlmin=tlmin-60.
            lmin=lmin+1
            print *,'time in minutes:',lmin
            tau = xn1 * 1500 * v0*t
            write(6,*)'time ',t,tau
            do i=1,n
cn ** aufpassen, ob Anzahldichte logarithmisch oder linear def. ist.
cn ** h(i) in cm**(-3)
cn ** lin
c            hin(i) = gin(i)/e(i) * rq0*1.e04/r(i) *1./xn1
c            hout(lmin,i) = g(i)/e(i) * rq0*1.e04/r(i) * 1./xn1
cn ** log
            hin(i) = 3*1.e+15*gin(i)/(4.*pi*r(i)**3.)             ! in m-3
            hout(lmin,i) = 3*1.e+15*g(i)/(4.*pi*r(i)**3.)         ! in m-3
            gout(lmin,i) = g(i)
            bin(i) = 3*1.e+12*gin(i)/(4*1000*pi*r(i)**2)
            bout(lmin,i) = 3*1.e+12*gout(lmin,i)/(4*1000*pi*r(i)**2)
            enddo

            call do_mass_balance(n, g, r, dlnr)

         endif
c     n ** ende output
        
      enddo                     ! ende zeitintegration
      close (16)

      open(15,file='sectionm')
      open(25,file='sectionn')
      do i=1,n
         write(15,'(62e14.5)')r(i),gin(i),(gout(j,i),j=1,lmin)
         write(25,'(62e14.5)')r(i),hin(i),(hout(j,i),j=1,lmin)
         write(35,'(62e14.5)')r(i),bin(i),(bout(j,i),j=1,lmin)
      enddo

      end

cn ****************************************************************************

      subroutine coad(n, dt, taug, taup, taul, tauu, prod, ploss,
     &     c, ima, g, r, e, ck, ec)

c collision subroutine, exponential approach

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

      do i=1,n
         prod(i)=0d0
         ploss(i)=0d0
      enddo
      
      
c     lower and upper integration limit i0,i1
      do i0=1,n-1
         if (g(i0).gt.gmin) go to 2000
      enddo
 2000 continue
      do i1=n-1,1,-1
         if (g(i1).gt.gmin) go to 2010
      enddo
 2010 continue
      
      
      do i=i0,i1
         do j=i,i1
            k=ima(i,j)
            kp=k+1
            
            x0=ck(i,j)*g(i)*g(j)
            x0=min(x0,g(i)*e(j))
            
            if (j.ne.k) x0=min(x0,g(j)*e(i))
            gsi=x0/e(j)
            gsj=x0/e(i)
            gsk=gsi+gsj
            
c     n *** Verlust fuer position i und j merken
            ploss(i) = ploss(i)+gsi
            ploss(j) = ploss(j)+gsj
c     write(6,*)'ploss ',i,j,ploss(i),ploss(j)
c     n ***
            g(i)=g(i)-gsi
            g(j)=g(j)-gsj
            gk=g(k)+gsk
            
            if (gk.gt.gmin) then
               x1=dlog(g(kp)/gk+1.d-60)
               flux=gsk/x1*(dexp(0.5*x1)-dexp(x1*(0.5-c(i,j))))
               flux=min(flux,gk)
               g(k)=gk-flux
               g(kp)=g(kp)+flux
c     n *** Gewinn fuer position k und kp merken
               prod(k) =  prod(k) + gsk - flux           
               prod(kp) = prod(kp)+ flux
               
c     n ***
            endif
         enddo
      enddo
      
      return
      end

cn ****************************************************************************

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

      do i=1,n
         do j=i,n
            x0=e(i)+e(j)
            do k=j,n
               if (e(k).ge.x0.and.e(k-1).lt.x0) then
                  if (c(i,j).lt.1.-1.d-08) then
                     kk=k-1
                     c(i,j)=dlog(x0/e(k-1))/(3.d0*dlnr)
                  else
                     c(i,j)=0.
                     kk=k
                  endif
                  ima(i,j)=min(n-1,kk)
                  go to 2000
               endif
            enddo
 2000       continue
            c(j,i)=c(i,j)
            ima(j,i)=ima(i,j)
         enddo
      enddo
      
      return
      end

cn ********************************************************************************

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

      if (isw.eq.0) then
c long kernel
         do j = 1,n
            r_tmp = r(j) / 1d6
            call fall_g(r_tmp, w_tmp)
            winf(j) = w_tmp * 100d0
            rr(j) = r(j) * 1d-4
         enddo

         do j=1,n
            do i=1,j
               if(r(j).le.50.) then
                  effi=4.5d-4*r(j)*r(j)*
     &                 (1.d0-3.d0/(max(3.d0,dble(r(i)))+1.d-2))
               else
                  effi=1.d0
               endif
               cck(j,i)=pi*(rr(j)+rr(i))*(rr(j)+rr(i))*effi*
     &              abs(winf(j)-winf(i))
               cck(i,j)=cck(j,i)
               write(7,*)'long ',i,j,r(i),r(j),cck(i,j)
            enddo
         enddo
      elseif (isw.eq.1) then
c     hall kernel
         do j = 1,n
         enddo

         do i = 1,n
            do j = 1,n
               r1 = r(i) ! should be in um
               r2 = r(j) ! should be in um
               call fall_g(r1 / 1d6, winf1) ! winf1 in m/s
               call fall_g(r2 / 1d6, winf2) ! winf2 in m/s
               call effic(r1, r2, ec_tmp) ! ec_tmp dimensionless
               cck(i,j) = 1d-6 * pi * (r1 + r2)**2 * ec_tmp
     &              * abs(winf1 - winf2)
               ! cck in cm^3/s = 1d6 * m^3/s
            enddo
         enddo
         
         do i=1,n
            do j=1,n
               write(7,'(5e14.5)')r(i),r(j),cck(i,j)
            enddo
         enddo
      elseif (isw.eq.2) then
c     golovin kernel
         do j = 1,n
            r_tmp = r(j) / 1d6
            call fall_g(r_tmp, w_tmp)
            winf(j) = w_tmp * 100d0
            rr(j) = r(j) * 1d-4
         enddo

         do j=1,n
            do i=1,j
               cck(j,i)=1.0*(e(j)+e(i))
               cck(i,j)=cck(j,i)
            enddo
         enddo
         
      else
      endif

c     two-dimensional linear interpolation of kernel
      do i=1,n
         do j=1,n
            jm=max0(j-1,1)
            im=max0(i-1,1)
            jp=min0(j+1,n)
            ip=min0(i+1,n)
            ck(i,j)=0.125*(cck(i,jm)+cck(im,j)+cck(ip,j)+cck(i,jp))
     &           +.5*cck(i,j)
            if (i.eq.j) ck(i,j)=0.5*ck(i,j)
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

      x0=0.
      x1=0.
      x_num=0.
      do i=1,n
         x0=x0+g(i)*dlnr
         x1=max(x1,g(i))
         if (dabs(x1-g(i)).lt.1.d-9) imax=i
c         gg(i) = g(i)*dlnr
         hh = 3*1.e+15*g(i)/(4.*pi*r(i)**3.)*dlnr
c         mm(i) = 1.e-06*gg(i)/hh(i)
         x_num=x_num+hh
      enddo
      
      print *,'mass ',x0,'max ',x1,'imax ',imax
      write(6,*)'x_num ',x_num

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
