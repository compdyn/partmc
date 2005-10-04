      program MonteCarlo

 
      integer  MM, M, M_comp, N_opt
      parameter (MM=10000)       !WORKING TOTAL NUMBER OF PARTICLES (MC particles)
      integer  i, i_loop
      integer  Time_count,lmin

      real*8   TIME, 
     &         del_T,tot_free_max, V_comp, V_bar 
      real*8   TIME_MAX
      real*8   N_tot
      real*8   pi, rho_p
      real*8   V(MM)
        
      integer  n_bin,k
      parameter (n_bin = 160)
      real*8   N_inf
      real*8   n_ini(n_bin)
      real*8   vv(n_bin),e(n_bin),r(n_bin),dp(n_bin)
      real*8   n_norm(n_bin),d_norm(n_bin)
C     *** For initialization

      real*8   rr(n_bin),delta_n(n_bin)
      real*8   dlnr
      integer  sum_a, sum_e
      parameter(TIME_MAX = 600., del_T = 1.)
      real*8   t1, ax, V_0, d_0, emin, sum
      real*8   delta_sum, sum_mass, tlmin
      integer  scal, n_samp, nt, i_top

      open(30,file='mc.d')

      call srand(10)

      do i_loop = 1,1
      call cpu_time(t1)
      write(6,*)'START ',i_loop, t1
      write(30,*)'i_loop=',i_loop,t1

      pi=3.1416            
      M=MM
      M_comp=M
      scal = 3

      dlnr=dlog(2.d0)/(3.*scal)
      ax=2.d0**(1.0/scal)
      V_0 = 4.1886e-15
c      V_0 = 3.351032e-14
      d_0 = (6*V_0/pi)**(1./3.)
      write(6,*)'d0 ',d_0
      emin = 1.e-15
      rho_p = 1000.               ! kg m-3, Constant Density of Particle  

   
      N_tot= 1.e+9                !#/m^3, Total particle number concentration in real physical 
                                  !volume//
      N_inf =N_tot
      V_comp=M/N_tot
      
c mass and radius grid

      e(1)=emin*0.5*(ax+1.)
      r(1)=1000.*dexp(dlog(3*e(1)/(4.*pi))/3.)
      vv(1)=1.e-06*e(1)/rho_p
      dp(1) = 1.e-06*r(1)
      n_ini(1) = 0. 
      
      do i=2,n_bin
         e(i)=ax*e(i-1)
         r(i)=1000.*dexp(dlog(3.*e(i)/(4.*pi))/3.)
         vv(i)=1.e-06*e(i)/rho_p
         dp(i)=1.e-06*2.*r(i)
         n_ini(i) = 0.
      enddo

c define bidisperse distribution

c      n_ini(97) =  0.99*M/dlnr
c      n_ini(126)  = 0.01*M/dlnr 
      n_ini(97) = (M-1)/dlnr
      n_ini(126) = 1/dlnr

      do i=1,n_bin
         rr(i) = dp(i)/2.
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
         write(6,*)'init ',i,r(i),n_ini(i),delta_n(i)
      enddo

c      write(6,*)'summe ',sum,delta_sum,sum_mass
      
      sum_e = 0.
      do k=1,n_bin
         sum_a = sum_e+1
          sum_e = sum_e+delta_n(k)
         do i=sum_a,sum_e
            V(i) = pi/6. * dp(k)**3.
         enddo
      enddo

      TIME=0.                    ! Initialization of Overall Time Frame for Collision Process */
      Time_count = 0
      tlmin = 0.
      lmin = 0
      M_comp = M

C *** CRITERIA SET FOR TOPPING UP & REPLICATING THE SUB-SYSTEM ***

      N_opt=M/2                  ! Optimum No.of Particles to be retained in the sub-system for 
                                 ! replicating it.*//

      call moments(V,M,N_tot,M_comp,V_comp,
     &                   TIME,tlmin,del_T,
     &                   V_bar,Time_count,rr,tot_free_max,
     &                   vv,dlnr,dp,n_samp)
      tlmin = 0. 
      nt = TIME_MAX/del_T
      do i_top = 1,nt             ! time-step loop
         TIME = real(i_top) / real(nt) * TIME_MAX

C        *** NEXT calculate the collsion probability ***   

         call sub_random(V,M,M_comp,V_comp,N_opt,tot_free_max,
     &                      del_T,TIME,tlmin,Time_count,
     &                      n_samp)
          

C        *** CALCULATING MOMENTS *** 

            call moments(V,M,N_tot,M_comp,V_comp,
     &                   TIME,tlmin,del_T,
     &                   V_bar,Time_count,rr,tot_free_max,
     &                   vv,dlnr,dp,n_samp)

      enddo                     ! end of topping up loop

      ENDDO                     ! end of i-loop

      end

C &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

      subroutine sub_random(V,M,M_comp,V_comp,N_opt,tot_free_max,
     *                      del_T,TIME,tlmin,Time_count,
     *                      n_samp)

      integer MM,M,M_comp,M_local,Time_count,n_samp,N_opt
      parameter (MM=10000)

      real*8 V(MM),V_comp
      real*8 tot_free_max
      real*8 del_T, TIME, pi, tlmin
  
      parameter (pi=3.1415)

      integer i_samp, i

      Time_count = Time_count + 1
      tlmin = tlmin +del_T

      do i_samp = 1,n_samp

         call coag_pair(V, MM, M, M_comp, V_comp, del_T, n_samp)
 
         M_local=M             !CURRENT NUMBER OF PARTICLES IN THE SYSTEM

C ***    REPLICATE THE SUB-SYSTEM FOR THE NEXT TOPPING-UP SYSTEM
         if (M .le. N_opt) then   ! topup
            call compress(V, MM)
            do i=1,M_local
               M=M+1
               V(M) = V(i)
            enddo
            V_comp=2*V_comp         ! UPDATING THE COMPUTATIONAL VOLUME
            M_comp=M     
         endif
      enddo       ! sampling loop
      
C *** If too many zeros in V-array, compress it
 
      if (real(M_comp - M)/M .gt. 0.5) then
         M_comp = M
         call compress(V, MM)
      endif
       
      return
      end

C &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
     
      subroutine coag_pair(V, MM, M, M_comp, V_comp, del_T, n_samp)
      
      integer n_samp, M, MM, M_comp ! INPUT
      real*8 V(MM), V_comp, del_T   ! INPUT

      integer s1, s2
      real*8 expo, p, k

      call find_rand_pair(V, MM, M, M_comp, s1, s2) ! test particles s1, s2

      call coag_kernel(V(s1), V(s2), k)
      expo = k * 1.0/V_comp * del_T * M*(M-1)/n_samp
      p = 1 - exp(-expo) ! probability of coagulation

      if (rand() .lt. p ) then ! coagulate particles s1 and s2
         V(s1) = V(s1) + V(s2)          
         V(s2) = 0.d+0
         M = M - 1
      endif

      return
      end
