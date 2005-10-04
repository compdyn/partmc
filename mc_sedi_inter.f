      program MonteCarlo

 
      integer  MM, M, M_initial,N_opt,M_local
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
      real*8   rr(n_bin),delta_n(n_bin)
      real*8   dlnr, tlmin, sum_mass, sum, delta_sum
      real*8   d_0, emin, ax, V_0
      integer  sum_a, sum_e, i_top, M_comp, scal, i_loop
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
      n_ini(1) = pi*dp(1)**3./2.*M/V_0 * exp(-(vv(1)/v_0))
      
      do i=2,n_bin
         e(i)=ax*e(i-1)
         r(i)=1000.*dexp(dlog(3.*e(i)/(4.*pi))/3.)
         vv(i)=1.e-06*e(i)/rho_p
         dp(i)=1.e-06*2.*r(i)
         n_ini(i)=pi/2.* dp(i)**3.*M/v_0*exp(-(vv(i)/v_0))
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
         
         do l=1,N_opt           ! time step cycle  
            sum=0.
            do i=1,MM
               sum=sum+v(i)
            enddo
            
C           *** NEXT calculate the collsion probability ***   
            call sub_random(V,M,M_comp,V_comp,tot_free_max,
     &           del_T,tmc_coll,TIME,tlmin,Time_count,
     &           i_count)
            
            sum=0.
            do i=1,MM
               sum=sum+v(i)
            enddo
            
C     *** CALCULATING MOMENTS *** 
            
            if (tlmin .gt. 1. ) then
               tlmin = tlmin -1.
               lmin = lmin + 1
               
               call moments(V,M,N_tot,M_comp,V_comp,
     &              TIME,tlmin,del_T,
     &              V_bar,Time_count,rr,tot_free_max,
     &              vv,dlnr,dp)
            endif
            
            if (TIME .ge. TIME_MAX) GOTO 2000
         enddo                  ! end of l-Loop
      
C     *** REPLICATE THE SUB-SYSTEM FOR THE NEXT TOPPING-UP SYSTEM
         sum=0.
         do i=1,MM
            sum=sum+V(i)
         enddo
         call compress(V)
         
         sum=0.
         do i=1,MM
            sum=sum+V(i)
         enddo
         
         M_local=M              !CURRENT NUMBER OF PARTICLES IN THE SYSTEM
         do i=1,M_local
            M=M+1
            V(M) = V(i)
         enddo
         
         V_comp=2*V_comp        ! UPDATING THE COMPUTATIONAL VOLUME
         M_comp=M     
      enddo                     ! end of topping up loop

 2000 continue
     
      ENDDO

      end

C &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

      subroutine sub_random(V,M,M_comp,V_comp,tot_free_max,
     *                      del_T,tmc_coll,TIME,tlmin,Time_count,
     *                      i_count)

      integer MM,M,M_comp,s1,s2,Time_count,n_samp, i_count
      parameter (MM=10000000)
      parameter (n_samp = 100)

      real*8 V(MM+1),a,random,V_comp,tlmin
      real*8 tot_free_max, K_avg
      real*8 del_T, tmc_coll,TIME, pi
      real*8 v_alt(MM+1), k, sum1, sum2, summ
      integer s3, s4, i, i_sum, i_versuch
  
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

      if ((v(s1) .eq. 0) .or. (v(s2) .eq. 0)) then
         a = -100.
      else
         if ( s1.ne.0) then             
            if (s2.ne.0) then
               call coag_kernel(V(s1),V(s2),k)
               a = k / tot_free_max; ! collision probability   
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
                      call coag_kernel(V(s3),V(s4),k)
                      summ = summ+k
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

      call coag_kernel(V(s1),V(s2),k)
      tmc_coll = (V_comp*2.)/(k*M)

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
      
