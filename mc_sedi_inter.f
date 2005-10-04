      program MonteCarlo

 
      integer  M, M_initial,N_opt,M_local 
      parameter (MM=10000000)       !WORKING TOTAL NUMBER OF PARTICLES (MC particles)
      integer  i,l,TOPUP
      integer  i_count
      integer  Time_count,lmin

      real*8   TIME, 
     &         del_T,tot_free_max, V_comp, V_bar 
      real*8   TIME_MAX
      real*8   N_tot
      real*8   pi,rho_p
      real*8   V(MM)
C     *** For SPD calculation **//
      integer  n_bin,k
      parameter (n_bin = 160)
      real*8   N_inf,tmc_coll
      real*8   n_ini(n_bin)
      real*8   vv(n_bin),e(n_bin),r(n_bin),dp(n_bin)
      real*8   help(n_bin),n_norm(n_bin),d_norm(n_bin)
C     *** For initialization
      real*8   rr(n_bin),n_ln(n_bin),delta_n(n_bin)
      real*8   dlnr
      integer  sum_a, sum_e
      parameter(TIME_MAX = 300.)

      open(30,file='mc.d')

      call srand(10)

      DO i_loop =1,1

      write(6,*)'START ',i_loop
      write(30,*)'i_loop=',i_loop

      pi=3.1416            
      M=10000000
      M_initial=M               ! ACTUAL INITIAL NUMBER OF PARTICLES
      M_comp=M
      scal = 3

      dlnr=dlog(2.d0)/(3.*scal)
      ax=2.d0**(1.0/scal)
c      V_0 = 1.e-12
      V_0 = 4.1886e-15
      d_0 = (6*V_0/pi)**(1./3.)
c      emin = 1.e-12 
      emin = 1.e-15
      rho_p = 1000.             ! kg m-3, Constant Density of Particle  

   
      N_tot= 1.e+9             !#/m^3, Total particle number concentration in real physical 
                                !volume//
      N_inf =N_tot
      V_comp=M/N_tot
      
c mass and radius grid
      e(1)=emin*0.5*(ax+1.)
      r(1)=1000.*dexp(dlog(3*e(1)/(4.*pi))/3.)
      vv(1)=1.e-06*e(1)/rho_p
      dp(1) = 1.e-06*r(1)
      n_ini(1) = pi*dp(1)**3./2.*M/V_0 * exp(-vv(1)/v_0)
      
      do i=2,n_bin
         e(i)=ax*e(i-1)
         r(i)=1000.*dexp(dlog(3.*e(i)/(4.*pi))/3.)
         vv(i)=1.e-06*e(i)/rho_p
         dp(i)=1.e-06*2.*r(i)
         n_ini(i)=pi/2.* dp(i)**3.*M/v_0*exp(-vv(i)/v_0)
      enddo
      
      do i=1,n_bin
         rr(i) = dp(i)/2.
      enddo

      do i=1,n_bin
         n_norm(i) = n_ini(i)/dp(i) * d_0/M
         d_norm(i) = dp(i)/d_0
      enddo

      sum=0.
      delta_sum = 0.
      sum_mass = 0.
      do i=1,n_bin
         delta_n(i) = int(n_ini(i)*dlnr)
         sum=sum+n_ini(i)*dlnr
         delta_sum = delta_sum + delta_n(i)
         sum_mass = sum_mass + rho_p * pi/6. * dp(i)**3.*n_ini(i)*dlnr
      enddo

      do i=1,n_bin
         help(i) = delta_n(i)*d_0/(dlnr*dp(i)*M)
      enddo
      
      sum_e = 0.
      do k=1,n_bin
         sum_a = sum_e+1
          sum_e = sum_e+delta_n(k)
         do i=sum_a,sum_e
            V(i) = pi/6. * dp(k)**3.
         enddo
      enddo

      sum=0.

      do i=1,MM
        sum=sum+v(i)
      enddo

      TIME=0.                    ! Initialization of Overall Time Frame for Collision Process */
      Time_count = 0
      tlmin = 0.
      lmin = 0
      M_comp = M
      i_count = 0

C *** CRITERIA SET FOR TOPPING UP & REPLICATING THE SUB-SYSTEM ***

      TOPUP=50                   ! TOTAL TOPUPING USED */
      N_opt=M/2                  ! Optimum No.of Particles to be retained in the sub-system for 
                                 ! replicating it.*//

      call coagmax(rr,tot_free_max,n_ini,dlnr)
       
      call moments(V,M,N_tot,M_comp,V_comp,
     &                   TIME,tlmin,del_T,
     &                   V_bar,Time_count,rr,tot_free_max,
     &                   vv,dlnr,dp)
  
      do i_top = 1,TOPUP        ! topping up cycle */

         do l=1,N_opt            ! time step cycle  
            sum=0.
            do i=1,MM
               sum=sum+v(i)
            enddo

C           *** NEXT calculate the collsion probability ***   
            call sub_random(V,M,M_comp,V_comp,tot_free_max,
     &                      del_T,tmc_coll,TIME,tlmin,Time_count,
     &                      i_count)
          
            M_local=M             !CURRENT NUMBER OF PARTICLES IN THE SYSTEM

            sum=0.
            do i=1,MM
               sum=sum+v(i)
            enddo

C           *** CALCULATING MOMENTS *** 

            if (tlmin .gt. 1. ) then
               tlmin = tlmin -1.
               lmin = lmin + 1

               call moments(V,M,N_tot,M_comp,V_comp,
     &                   TIME,tlmin,del_T,
     &                   V_bar,Time_count,rr,tot_free_max,
     &                   vv,dlnr,dp)
             endif

        if (TIME .ge. TIME_MAX) GOTO 2000
        enddo    ! end of l-Loop
      
C *** REPLICATE THE SUB-SYSTEM FOR THE NEXT TOPPING-UP SYSTEM
       sum=0.
       do i=1,MM
          sum=sum+V(i)
       enddo
       call compress(V)

        sum=0.
        do i=1,MM
           sum=sum+V(i)
        enddo

        do i=1,M_local
          M=M+1
          V(M) = V(i)
        enddo
      
        V_comp=2*V_comp         ! UPDATING THE COMPUTATIONAL VOLUME
        M_comp=M     
      enddo                     ! end of topping up loop

 2000 continue
     
      ENDDO

      end

C &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

      subroutine sub_random(V,M,M_comp,V_comp,tot_free_max,
     *                      del_T,tmc_coll,TIME,tlmin,Time_count,
     *                      i_count)

      integer MM,s1,s2,Time_count,n_samp
      parameter (MM=10000000)
      parameter (n_samp = 100)

      real*8 V(MM+1),a,random,V_comp
      real*8 tot_free_max, K_avg
      real*8 del_T, tmc_coll,TIME, pi
      real*8 v_alt(MM+1)
  
      parameter (pi=3.1415)


      i_versuch = 0
 1000 continue
      i_versuch = i_versuch+1
      
C *** NEXT calculate the collsion probability ***   
c *** double(rand())/RAND_MAX: this is the random number in [0,1] ***
          
      random=rand()
      s1=int(random*M_comp)  !CHOOSE A RANDOM  NUMBER s1 in [0,M]        
      random=rand()
      s2=int(random*M_comp)  !CHOOSE A RANDOM  NUMBER s2 in [0,M]

      if ((v(s1) .eq.0) .or. (v(s2) .eq.0)) then
           a = -100.
      else

         if ( s1.ne.0) then             
            if (s2.ne.0) then
                a=coag_kernel(V(s1),V(s2))/tot_free_max; !//collision probability   
                if(abs(s2-s1) .lt. 0.1) a=0.
            else 
              a=-100.
            endif
          else
              a=-100.
          endif
      endif

      if (rand() .gt. a ) goto 1000 
    
C *** END of calculation for Collision Probability and collision analysis 
      
      sum1=0.
      do i=1,M_comp
         sum1=sum1+v(i)
         v_alt(i) = v(i)
      enddo

      i_sum = 0.
      summ = 0.

      do while (i_sum .lt. (n_samp*n_samp))
         random = rand()
         s3 = int(random*M_comp)
         random = rand()
         s4 = int(random*M_comp)
         if (s1.ne. 0.) then
            if(s2.ne.0.) then
                if (V(s3) .ne. 0.) then
                   if (V(s4) .ne. 0.) then
                      summ = summ+coag_kernel(V(s3),V(s4))
                      i_sum = i_sum + 1
                    endif
                 endif
              endif
           endif
        enddo

        K_avg  = summ/(n_samp*(n_samp-1))
      
C *** Calculation for "del_T" for the time increment for each collision 

      del_T=(V_comp*2.)/(K_avg*M*(M-1)) !//time increment

      V(s1) = V(s1)+V(s2)          

C *** Characteristic collision time-MC//

      tmc_coll = (V_comp*2.)/(coag_kernel(V(s1),V(s2))*M)

C *** Removal of the particle s2 by setting its volume to 0.
      
      V(s2) = 0.d+0
      M=M-1
      i_count=i_count+1


      sum2=0.
      do i=1,MM
         sum2=sum2+v(i)
      enddo

C *** If too many zeros in V-array, compress it
 
      if (real(i_count)/M .gt. 0.5) then
         M_comp = M
         call compress(V)
         i_count = 0
      endif
       
      TIME=TIME+del_T            ! Time Updating //
      Time_count = Time_count + 1
      tlmin = tlmin +del_T
      return
      end
      
C &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

      real*8 function coag_kernel(a,b)
      real*8 a,b
      real*8 pi,const,onethird
      real*8 r1,r2,winf1,winf2,ec
      parameter (pi = 3.1415)

      const = 3./(4.*pi)
      onethird  = 1./3.
      r1 = (const*a)**onethird
      r2 = (const*b)**onethird
      call fallg(r1,winf1)
      call fallg(r2,winf2)
      call effic(r1,r2,ec)
          coag_kernel = ec * pi* (r1+r2)*(r1+r2)*abs(winf1-winf2)
      return
      end

C &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

      subroutine fallg(r,winf)

c terminal velocity of falling drops
      implicit double precision (a-h,o-z)
      dimension b(7),c(6),rat(20),r0(15),ecoll(15,20)
      data b /-0.318657e1,0.992696,-0.153193e-2,-0.987059e-3,
     &        -0.578878e-3,0.855176e-4,-0.327815e-5/
      data c /-0.500015e1,0.523778e1,-0.204914e1,0.475294,-0.542819e-1,
     &         0.238449e-2/
      data pi /3.141592654/
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
c collision efficiencies of hall kernel
      implicit double precision (a-h,o-z)
      dimension rat(21),r0(15),ecoll(15,21)
      data r0 /6.,8.,10.,15.,20.,25.,30.,40.,50.,
     &         60.,70.,100.,150.,200.,300./
      data rat /0.,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,
     &          0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0/
      data ecoll /
     &  0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001
     & ,0.001,0.001,0.001,0.001,0.001,0.003,0.003,0.003,0.004,0.005
     & ,0.005,0.005,0.010,0.100,0.050,0.200,0.500,0.770,0.870,0.970
     & ,0.007,0.007,0.007,0.008,0.009,0.010,0.010,0.070,0.400,0.430
     & ,0.580,0.790,0.930,0.960,1.000,0.009,0.009,0.009,0.012,0.015
     & ,0.010,0.020,0.280,0.600,0.640,0.750,0.910,0.970,0.980,1.000
     & ,0.014,0.014,0.014,0.015,0.016,0.030,0.060,0.500,0.700,0.770
     & ,0.840,0.950,0.970,1.000,1.000,0.017,0.017,0.017,0.020,0.022
     & ,0.060,0.100,0.620,0.780,0.840,0.880,0.950,1.000,1.000,1.000
     & ,0.030,0.030,0.024,0.022,0.032,0.062,0.200,0.680,0.830,0.870
     & ,0.900,0.950,1.000,1.000,1.000,0.025,0.025,0.025,0.036,0.043
     & ,0.130,0.270,0.740,0.860,0.890,0.920,1.000,1.000,1.000,1.000
     & ,0.027,0.027,0.027,0.040,0.052,0.200,0.400,0.780,0.880,0.900
     & ,0.940,1.000,1.000,1.000,1.000,0.030,0.030,0.030,0.047,0.064
     & ,0.250,0.500,0.800,0.900,0.910,0.950,1.000,1.000,1.000,1.000
     & ,0.040,0.040,0.033,0.037,0.068,0.240,0.550,0.800,0.900,0.910
     & ,0.950,1.000,1.000,1.000,1.000,0.035,0.035,0.035,0.055,0.079
     & ,0.290,0.580,0.800,0.900,0.910,0.950,1.000,1.000,1.000,1.000
     & ,0.037,0.037,0.037,0.062,0.082,0.290,0.590,0.780,0.900,0.910
     & ,0.950,1.000,1.000,1.000,1.000,0.037,0.037,0.037,0.060,0.080
     & ,0.290,0.580,0.770,0.890,0.910,0.950,1.000,1.000,1.000,1.000
     & ,0.037,0.037,0.037,0.041,0.075,0.250,0.540,0.760,0.880,0.920
     & ,0.950,1.000,1.000,1.000,1.000,0.037,0.037,0.037,0.052,0.067
     & ,0.250,0.510,0.770,0.880,0.930,0.970,1.000,1.000,1.000,1.000
     & ,0.037,0.037,0.037,0.047,0.057,0.250,0.490,0.770,0.890,0.950
     & ,1.000,1.000,1.000,1.000,1.000,0.036,0.036,0.036,0.042,0.048
     & ,0.230,0.470,0.780,0.920,1.000,1.020,1.020,1.020,1.020,1.020
     & ,0.040,0.040,0.035,0.033,0.040,0.112,0.450,0.790,1.010,1.030
     & ,1.040,1.040,1.040,1.040,1.040,0.033,0.033,0.033,0.033,0.033
     & ,0.119,0.470,0.950,1.300,1.700,2.300,2.300,2.300,2.300,2.300
     & ,0.027,0.027,0.027,0.027,0.027,0.125,0.520,1.400,2.300,3.000
     & ,4.000,4.000,4.000,4.000,4.000/
c two-dimensional linear interpolation of the collision efficiency
         rq = r1/r2
         if (rq.gt.1.) then
            r_help = r1
            r1 = r2
            r2 = r_help
            rq=1/rq
         endif
         do k=2,15
            if (r2.le.r0(k).and.r2.ge.r0(k-1)) then
               ir=k
            elseif (r2.gt.r0(15)) then
               ir=16
            elseif (r2.lt.r0(1)) then
               ir=1
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
     &                 p*(1.-q)*ecoll(ir,iq-1)+
     &                 q*(1.-p)*ecoll(ir-1,iq)+
     &                 p*q*ecoll(ir,iq)
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

C &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&


      subroutine moments(V,M,N_tot,M_comp,V_comp,
     &                       TIME,tlmin,del_T,V_bar,
     &                       Time_count,rr,
     &                       tot_free_max,vv,dlnr,dp)

      integer M,MM,n_bin,Time_count,NN_cnt
      real*8 nv_conc,dlnr
      parameter (MM=10000000,n_bin=160)
      real*8 V(MM),N_tot,vv(n_bin),dp(n_bin)
      real*8 V_comp,del_T,TIME,tot_free_max
      real*8 V_bar
      real*8 eta(n_bin),rr(n_bin)
      real*8 psi(n_bin),g(n_bin)
      real*8 n_ln(n_bin)
      real*8 pi,rho_w

      parameter (pi=3.1415)
      parameter (rho_w = 1000.)

      
      V_0 = 1.e-12
      d_0 = (6*V_0/pi)**(1./3.)

      do k=1,n_bin

          NN_cnt = 0
          vv_cnt = 0.
          do i=1,M_comp
              if ((V(i).ge. vv(k-1)) 
     &                .and. (V(i) .lt. vv(k))) then
                    NN_cnt = NN_cnt +1
                    vv_cnt = vv_cnt + v(i)
               endif
           enddo

              eta(k)  = dp(k)/d_0
              nv_conc = NN_cnt/V_comp
              vv_conc = vv_cnt/V_comp
              n_ln(k) = nv_conc/dlnr
              psi(k) =  nv_conc*d_0/(dlnr*dp(k)*MM)
              g(k)   =  vv_conc/dlnr
  
       enddo
          
       call coagmax(rr,tot_free_max,n_ln,dlnr)
c       write(6,*)'tot_free_max nach coagmax ',tot_free_max
                        
c       write(6,*)'gesamtanzahl ',time,M

       write(30,*)'Time = ',TIME
       sum_masse = 0.
       do k=1,n_bin
            write(30,'(i4,6e14.5)')k,eta(k), psi(k),
     &          dp(k)/2.,n_ln(k),g(k)
                sum_masse = sum_masse + g(k)*dlnr
       enddo
c       write(6,*)'masse ',time,sum_masse

       return
       end

C &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

       subroutine coagmax(rr,tot_free_max,n_ln,dlnr)

       implicit double precision (a-h,o-z)
       integer n_bin
       parameter (n_bin=160)
       real*8 v_bin(n_bin), cck, rr(n_bin),tot_free_max 
       real*8 pi,n_ln(n_bin),dlnr
       parameter (pi=3.1415)

       do k=1,n_bin
          v_bin(k)=4./3.*pi*rr(k)**3.
       enddo

       tot_free_max =0.
       do k=1,n_bin
          if(n_ln(k)*dlnr .ge. 1.) then
             do ll=1,k
    
                 cck= coag_kernel(v_bin(k),v_bin(ll))

                 if (cck.gt.tot_free_max) then
                     tot_free_max=cck
                     imax = k
                     jmax= ll
                 endif
              enddo
           endif
       enddo

c       write(6,*)'in coagmax ',imax,jmax,rr(imax),rr(jmax),tot_free_max

       return
       end

C &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

      subroutine compress(V)

      integer MM,i_w,i_v
      parameter (MM = 10000000)
      real*8 V(MM),w(MM)

      sum = 0.
      do i=1,MM
         sum=sum+V(i)
      enddo
c      write(6,*)'totvol in compress1 ',sum

      i_w = 1
      do i_v=1,MM
           if(v(i_v) .ne.0.) then
                 w(i_w)= v(i_v)
                 i_w = i_w+1
           endif
      enddo
      i_count = 0
c     write(6,*)'compressed now ', i_count
      do i=1,MM
           V(i) = W(i)
           W(i) = 0.
       enddo

      sum = 0.
      do i=1,MM
         sum=sum+V(i)
      enddo
c      write(6,*)'totvol in compress2 ',sum
      return
      end
