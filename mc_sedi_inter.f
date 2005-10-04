      program MonteCarlo

 
      integer  MM, M, M_initial,N_opt
      parameter (MM=10000000)       !WORKING TOTAL NUMBER OF PARTICLES (MC particles)
      integer  i,l,TOPUP
      integer  Time_count,lmin

      real*8   TIME, 
     &         del_T,k_max, V_comp
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
      real*8 g(n_bin), n_ln(n_bin)
C     *** For initialization
      real*8   rr(n_bin),delta_n(n_bin)
      real*8   dlnr, tlmin, sum_mass, delta_sum
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
      n_ini(1) = pi*dp(1)**3./2.*M/V_0 * exp(-(vv(1)/V_0))
      
      do i=2,n_bin
         e(i)=ax*e(i-1)
         r(i)=1000.*dexp(dlog(3.*e(i)/(4.*pi))/3.)
         vv(i)=1.e-06*e(i)/rho_p
         dp(i)=1.e-06*2.*r(i)
         n_ini(i)=pi/2.* dp(i)**3.*M/V_0*exp(-(vv(i)/V_0))
      enddo
      
      do i=1,n_bin
         rr(i) = dp(i)/2.
      enddo

      do i=1,n_bin
         n_norm(i) = n_ini(i)/dp(i) * d_0/M
         d_norm(i) = dp(i)/d_0
      enddo

      delta_sum = 0.
      sum_mass = 0.
      do i=1,n_bin
         delta_n(i) = int(n_ini(i)*dlnr)
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

      TIME=0.                    ! Initialization of Overall Time Frame for Collision Process */
      Time_count = 0
      tlmin = 0.
      lmin = 0
      M_comp = M

C *** CRITERIA SET FOR TOPPING UP & REPLICATING THE SUB-SYSTEM ***

      TOPUP=50                   ! TOTAL TOPUPING USED */
      N_opt=M/2                  ! Optimum No.of Particles to be retained in the sub-system for 
                                 ! replicating it.*//

      call moments(MM, V, n_bin, M_comp, V_comp, vv, dlnr, g, n_ln)
      call coagmax(n_bin, rr, n_ln, dlnr, k_max)
      call print_info(n_bin, TIME, tlmin, dp, g, n_ln)
  
      do i_top = 1,TOPUP        ! topping up cycle */
         do l=1,N_opt           ! time step cycle  
            
C           *** NEXT calculate the collsion probability ***   
            call sub_random(V,M,M_comp,V_comp,k_max,
     &           del_T,tmc_coll,TIME,tlmin,Time_count)
            TIME=TIME+del_T     ! Time Updating //
            Time_count = Time_count + 1
            tlmin = tlmin +del_T
            
C     *** CALCULATING MOMENTS *** 
            
            if (tlmin .gt. 1. ) then
               tlmin = tlmin -1.
               lmin = lmin + 1
               
               call moments(MM, V, n_bin, M_comp, V_comp, vv,
     &              dlnr, g, n_ln)
               call coagmax(n_bin, rr, n_ln, dlnr, k_max)
               call print_info(n_bin, TIME, tlmin, dp, g, n_ln)
            endif
            if (TIME .ge. TIME_MAX) GOTO 2000
         enddo                  ! end of l-Loop
      
C     *** REPLICATE THE SUB-SYSTEM FOR THE NEXT TOPPING-UP SYSTEM
         call double(MM, M_comp, V, V_comp)
      enddo                     ! end of topping up loop

 2000 continue
     
      ENDDO

      end

C &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

      subroutine sub_random(V,M,M_comp,V_comp,k_max,
     *                      del_T,tmc_coll,TIME,tlmin,Time_count)

      integer MM,M,M_comp,s1,s2,Time_count,n_samp
      parameter (MM=10000000)
      parameter (n_samp = 100)

      real*8 V(MM+1),V_comp,tlmin
      real*8 k_max, k_avg
      real*8 del_T, tmc_coll,TIME, pi
      real*8 k
  
      parameter (pi=3.1415)

      external kernel_sedi

      call find_rand_pair_acc_rej(MM, V, M_comp, k_max,
     &     kernel_sedi, s1, s2)
      call kernel_avg(MM, M_comp, V, kernel_sedi, n_samp, k_avg)
      
C *** Calculation for "del_T" for the time increment for each collision 

      del_T=(V_comp*2.)/(k_avg*M*(M-1)) !//time increment

      V(s1) = V(s1)+V(s2)          

C *** Characteristic collision time-MC//

      call kernel_sedi(V(s1),V(s2),k)
      tmc_coll = (V_comp*2.)/(k*M)

C *** Removal of the particle s2 by setting its volume to 0.
      
      V(s2) = 0.d+0
      M=M-1

C *** If too many zeros in V-array, compress it
 
      if (real(M_comp - M)/M .gt. 0.5) then
         call compress(MM, M_comp, V)
      endif
       
      return
      end
      
