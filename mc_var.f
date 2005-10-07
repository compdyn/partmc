C mc_var.f
C
C Monte Carlo with variable timestep.

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine sub_random(MM,V,M,M_comp,V_comp,k_max,
     *                      del_T,tmc_coll,TIME)

      integer MM,M,M_comp,s1,s2,n_samp
      parameter (n_samp = 100)

      real*8 V(MM),V_comp
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
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
