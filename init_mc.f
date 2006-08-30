      
      program init
      
      integer n, scal, n_spec
      integer MM, M
      integer i_start, i_end
      parameter (n = 160, scal = 3, MM = 100000, n_spec = 2)
      real*8 dp(n)                    ! wet diameter
      real*8 rn(n),rw(n)              ! dry radius, wet radius
      real*8 x1,x2,xacc,x
      real*8 v(n,2)                   ! 1: water mass distribution; 2: aerosol mass distribution
      real*8 numln(n)                 ! number distribution (logarithmic version)
      integer n_ini(n)
      real*8 const
      real*8 A, B(n)                  ! needed for equilibrium calculation
      real*8 M_w, M_s                 ! molecular weight of water, aerosol
      real*8 num(n),bin_v(n)
      real*8 rho_w, rho_n             ! water density, aerosol density
      real*8 c4,c3,c1,c0,dc3,dc2,dc0
      real*8 RR                       ! gas constant
      real*8 RH, T, T0        ! relative humidity, temperature
      real*8 sig_w, eps, nu           ! surface tension, solubility, number of ions
      real*8 mv                   ! specific humidity
      real*8 p0T                  ! partial pressure of water vapor
      real*8 p00, T00, pi, v_min, dlnr

      real*8 V_indiv(MM,n_spec), vol_frac(n_spec, n)
      character *100 outfile

      parameter(M_w=18.d-03,M_s= 132.d-03  )
      parameter(rho_w=1000.d0, rho_n=1800.d0)
      parameter(RR = 8.314d0)
      parameter(p00=611.d0,T00=273.15d0)
      parameter (v_min = 1.d-24)       ! minimum volume (m^3) for making grid
      parameter(sig_w =0.073d0 , nu = 3, eps = .25d0)
      parameter(pi = 3.14159265358979323846d0)

      real*8 dg,logsg
      parameter (dg=0.5e-6,logsg=0.4)

      outfile = 'init_eq.d'

      T0 = 298.d0
      RH = 1.0d0
      const = 1/sqrt(2*pi)

      do i=1,n
         rn(i) = 0.d0
         bin_v(i) = 0.d0
         num(i) = 0.d0
         v(i,1) = 0.d0
         v(i,2) = 0.d0
         numln(i) = 0.d0        
      enddo

       p0T = p00 *10d0**(7.45d0*(T0-T00)/(T0-38.d0))
       mv  = RH * p0T*M_w/(RR*T0)

      call make_grid(n, scal, v_min, bin_v, rn, dlnr)

      do j=1,n 
          num(j) = const * MM /logsg *
     &         dexp(-(dlog10(rn(j)) - dlog10(dg/2))**2./
     &        (2*logsg**2.))
      enddo
         
      do j=1,n     
         numln(j) = num(j)/2.303
      enddo

cn *** Equilibrium: dm/dt = 0
cn *** Set numerator in eq. (7), Majeed and Wexler, equal to zero and solve for diameter.


         T= T0
         A = 4.d0*M_w*sig_w/(RR*T*rho_w)

         do i=1,n
            B(i) = nu*eps*M_w*rho_n*rn(i)**3.d0/(M_s*rho_w)
         enddo

         c4 = log(RH)/8.d0
         c3 = A/8.d0

         dc3 = log(RH)/2.d0
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
            c1 = B(i)-log(RH)*rn(i)**3.
            c0 = A*rn(i)**3
            dc0 = c1

            call rtnewt(x1,x2,xacc,x,c4,c3,c1,c0,dc3,dc2,dc0)

            dp(i) = x
            rw(i) = dp(i)/2.d0
         enddo

         do i=1,n
            n_ini(i) = int(numln(i) * dlnr)
            v(i,1) = numln(i)*
     &           4.d0*pi/3.d0*((dp(i)/2.d0)**3 - rn(i)**3)
            v(i,2) = numln(i)*4.d0*pi/3.d0*rn(i)**3
            vol_frac(1,i) = v(i,1) / (v(i,1) + v(i,2))
            vol_frac(2,i) = v(i,2) / (v(i,1) + v(i,2)) 
         enddo

        i_start = 1
        i_end = MM
        call compute_volumes(n, n_spec, vol_frac, MM, i_start, i_end, 
     &            n_ini, rw, dlnr, V_indiv, M)

cn *** write to output file
c
       open(8,file=outfile)
       write(8,*)T0,mv
       do i=1,n
             write(8,'(10e20.11)')rn(i),rw(i),
     &         numln(i)*dlnr,v(i,1),v(i,2),
     &         v(i,1)/(v(i,1)+v(i,2)), v(i,2)/(v(i,1)+v(i,2)),
     &         (v(i,1)+v(i,2))/numln(i)
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

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine make_grid(n_bin, scal, v_min, bin_v, bin_r,
     $     dlnr)

      integer n_bin        ! INPUT: number of bins
      integer scal         ! INPUT: scale factor
      real*8 bin_v(n_bin)  ! OUTPUT: volume of particles in bins (m^3)
      real*8 bin_r(n_bin)  ! OUTPUT: radius of particles in bins (m)
      real*8 dlnr          ! OUTPUT: scale factor
      real*8 v_min

      integer i
      real*8 ax

      real*8 pi
      parameter (pi = 3.14159265358979323846d0)

      dlnr = dlog(2d0) / (3d0 * dble(scal))
      ax = 2d0**(1d0 / dble(scal)) ! ratio bin_v(i)/bin_v(i-1)


      do i = 1,n_bin
         ! volume (m^3)
         bin_v(i) = v_min * 0.5d0 * (ax + 1d0) * ax**(i - 1)
         ! radius (m)
         bin_r(i) = dexp(dlog(3d0 * bin_v(i) /
     &             (4d0 * pi)) / 3d0)
      enddo

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine compute_volumes(n_bin, n_spec, vol_frac,
     *     MM, i_start, i_end,
     *     n_ini, bin_r, dlnr, V, M)

      integer n_bin        ! INPUT: number of bins
      integer n_spec       ! INPUT: number of species
      real*8 vol_frac(n_spec,n_bin) ! INPUT: composition of particles
      integer MM           ! INPUT: physical size of V
      integer i_start      ! INPUT:
      integer i_end        ! INPUT:
      integer n_ini(n_bin) ! INPUT: initial number distribution
      real*8 bin_r(n_bin)  ! INPUT: wet radius of particles in bin (m)
      real*8 dlnr          ! INPUT: scale factor
      real*8 V(MM,n_spec)  ! OUTPUT: particle volumes  (m^3)
      integer M            ! OUTPUT: logical dimension of V
      
      real*8 total_vol_frac(n_bin)
      integer k, i, j, sum_e, sum_a, delta_n, i_spec
         
      real*8 pi
      parameter (pi = 3.14159265358979323846d0)
      
      sum_e = i_start - 1

      do i=1,n
            total_vol_frac(i) = 0.d0
      enddo

      do i=1,n_bin
         do j=1,n_spec
             total_vol_frac(i) = total_vol_frac(i) + vol_frac(j,i)
         enddo
      enddo
 
      
      do k = 1,n_bin
         delta_n = n_ini(k)
         sum_a = sum_e + 1 
         sum_e = sum_e + delta_n
         do i = sum_a,sum_e
            do i_spec = 1,n_spec
               V(i,i_spec) = vol_frac(i_spec,k)/total_vol_frac(k) *
     *              4d0/3d0 * pi * bin_r(k)**3
            enddo
            write(6,*)'V array ',k,n_ini(k),i,bin_r(k),
     *             V(i,1),V(i,2),V(i,1)+V(i,2)
         enddo
      enddo

      M = sum_e-i_start+1
      write(6,*)'M ', M
      end


