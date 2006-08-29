      
      program init
      
      integer   iz                    ! number of layers
      integer   nmod                  ! number of modes
      real*8    dz                    ! vertical grid size
      parameter (n = 120, scal = 2.,nmod = 3)
      parameter (iz = 3, dz = 5.)
      real*8 dp(n)                    ! wet diameter
      real*8 rn(n),rw(n)              ! dry radius, wet radius
      real*8 x1,x2,xacc,x
      real*8 g(n,iz,2)                ! 1: water mass distribution; 2: aerosol mass distribution
      real*8 numln(n)                 ! number distribution (logarithmic version)
      real*8 s(n)                     ! surface distribution
      real*8 v(n)                     ! volume distribution
      real*8 A, B(n)                  ! needed for equilibrium calculation
      real*8 M_w, M_s                 ! molecular weight of water, aerosol
      real*8 num(n),dg(3),n_0(3),en(n)
      real*8 num_0(3,n),logsg(3)
      real*8 rho_w, rho_n             ! water density, aerosol density
      real*8 c4,c3,c1,c0,dc3,dc2,dc0
      real*8 RR                       ! gas constant
      real*8 RH(iz), T, T0(iz)        ! relative humidity, temperature
      real*8 sig_w, eps, nu           ! surface tension, solubility, number of ions
      real*8 mv(iz)                   ! specific humidity
      real*8 p0T(iz)                  ! partial pressure of water vapor
      real*8 p00, T00, pi, emin
      character *100 outfile

      parameter(M_w=18.d-03,M_s= 132.d-03  )
      parameter(rho_w=1000.d0, rho_n=1800.d0)
      parameter(RR = 8.314d0)
      parameter(p00=611.d0,T00=273.15d0)

      parameter(sig_w =0.073d0 , nu = 3, eps = .25d0)
      parameter(pi = 3.14159265358979323846d0)
      parameter(emin = 1.d-15)


      outfile = 'init_eq.d'

      T0(1) = 290.d0

cn *** Temperature profile
c
      do k=2,iz
         T0(k) = T0(1)
      enddo

      RH(1) = 0.99d0


cn *** Profile of relative humidity
c
      do k=2,iz
         rh(k) = rh(1) 
      enddo
      
      const = 1.d0/sqrt(2.d0*pi)
      const2 = 1.d-12*(pi/6) / 2.303d0
      
      do i=1,n
         rn(i) = 0.d0
         en(i) = 0.d0
         num(i) = 0.d0
         do k=1,iz
            g(i,k,1) = 0.d0
            g(i,k,2) = 0.d0
         enddo
         numln(i) = 0.d0        
      enddo

      do k=1,iz
         p0T(k) = p00 *10**(7.45d0*(T0(k)-T00)/(T0(k)-38.d0))
         mv(k)  = RH(k) * p0T(k)*M_w/(RR*T0(k))
      enddo

      
      dlnr=dlog(2.d0)/(3.d0*dble(scal))     
      ax=2.d0**(1.d0/dble(scal))
        
cn *** mass and radius grid
c
      en(1)=emin*0.5d0*(ax+1.d0)
      rn(1)=1000.d0*dexp(dlog(3.d0*en(1)/(4.d0*pi))/3.d0)
      do i=2,n
         en(i)=ax*en(i-1)
         rn(i)=1000.d0*dexp(dlog(3.d0*en(i)/(4.d0*pi))/3.d0)
      enddo

cn *** convert rn in m

      do i=1,n
         rn(i) = rn(i)*1.d-06
      enddo

cn *** Equilibrium: dm/dt = 0
cn *** Set numerator in eq. (7), Majeed and Wexler, equal to zero and solve for diameter.

      do k=1,iz

         write(6,*)'T ',k,T0(k),RH(k)
         T= T0(k)
         A = 4.d0*M_w*sig_w/(RR*T*rho_w)

         do i=1,n
            B(i) = nu*eps*M_w*rho_n*rn(i)**3.d0/(M_s*rho_w)
         enddo

         c4 = log(RH(k))/8.d0
         c3 = A/8.d0

         dc3 = log(RH(k))/2.d0
         dc2 = 3.d0*A/8.d0

         x1 = 0.d0
         x2 = 10.d0
         xacc = 1.d-15
c *** 
         do i=1,n
            i0=i
            if (rn(i).gt.0.d0) go to 2000
         enddo
 2000    continue
         do i=n,1,-1
            i1=i
            if (rn(i).gt.0.d0) go to 2010
         enddo
 2010    continue

         write(6,*)'i0,i1 ',i0,i1

         do i=i0,i1             ! Find equilibrium 
            c1 = B(i)-log(RH(k))*rn(i)**3.
            c0 = A*rn(i)**3
            dc0 = c1

            call rtnewt(x1,x2,xacc,x,c4,c3,c1,c0,dc3,dc2,dc0)

            dp(i) = x
            rw(i) = dp(i)/2.d0
         enddo

         do i=i0,i1
            g(i,k,1) =
     &           rho_w*4.d0*pi/3.d0*((dp(i)/2.d0)**3 - rn(i)**3)
            g(i,k,2) = rho_n*4.d0*pi/3.d0*rn(i)**3
            v(i) = 1.d+06*
     &           numln(i)*4.d0*pi/3.d0*(dp(i)/2.d0)**3 
         enddo
       enddo !!! K-Loop

cn *** write to output file
c
       open(8,file=outfile)
       do k=1,1
          write(8,*)T0(k),mv(k)
cn *** adjust output, so that it matches input in parcel model!!!
cn *** ( in this case, start with bin number 20) 
          do i=1,n
             write(8,'(5e20.11)')rn(i),rw(i),g(i,k,1),g(i,k,2)
          enddo 
       enddo

       end


cn ************************************************************************

      subroutine rtnewt(x1,x2,xacc,x,c4,c3,c1,c0,dc3,dc2,dc0)
      integer jmax
      real*8 x,x1,x2,xacc 
      real*8 c4,c3,c1,c0,dc3,dc2,dc0
      external funcd
      parameter (jmax=400)
      integer j
      real*8 df, dx, f, d 

      x = .5d0*(x1+x2)

      do j=1,jmax
         call funcd(x,f,df,d,c4,c3,c1,c0,dc3,dc2,dc0)
         dx=f/df
         x=x-dx
         if((x1-x)*(x-x2).lt.0.d0) then
            write(6,*)'x1,x2,x ',x1,x2,x
            pause 'rtnewt jumped out of brackets'
         endif
         if(abs(dx).lt.xacc) then

         return
         endif
      enddo

      pause 'rtnewt exceeded maximum iteration '
      end

C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine funcd(x,f,df,d_p,c4,c3,c1,c0,dc3,dc2,dc0)

      real*8 x,f,df,d_p
      real*8 c4,c3,c1,c0,dc3,dc2,dc0


      f = c4 * x**4.d0 - c3*x**3.d0+c1*x + c0
      df = dc3 * x**3.d0 -dc2*x**2.d0 + dc0


      end
