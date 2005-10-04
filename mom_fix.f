C &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      
      subroutine moments(V,M,N_tot,M_comp,V_comp,
     &     TIME,tlmin,del_T,V_bar,
     &     Time_count,rr,
     &     tot_free_max,vv,dlnr,dp,n_samp)
      
      integer M,MM,n_bin,Time_count,NN_cnt, M_comp, n_samp
      real*8 nv_conc,dlnr, tlmin
      parameter (MM=10000,n_bin=160)
      real*8 V(MM),N_tot,vv(n_bin),dp(n_bin)
      real*8 V_comp,del_T,TIME,tot_free_max
      real*8 V_bar
      real*8 rr(n_bin)
      real*8 g(n_bin+1)
      real*8 n_ln(n_bin)
      real*8 pi,rho_w,p_max, V_0, d_0, vv_cnt
      real*8 vv_conc, r_samp, sum_masse
      integer k, i, lmin
      
      parameter (pi=3.1415)
      parameter (rho_w = 1000.)
      parameter (p_max = 0.01)
      
      V_0 = 1.e-12
      d_0 = (6*V_0/pi)**(1./3.)
      
      do k=1,n_bin
         NN_cnt = 0
         vv_cnt = 0.
         do i=1,M_comp
            if ((V(i).ge. vv(k-1)) 
     &           .and. (V(i) .lt. vv(k))) then
               NN_cnt = NN_cnt +1
               vv_cnt = vv_cnt + V(i)
            endif
         enddo
         nv_conc = NN_cnt/V_comp
         vv_conc = vv_cnt/V_comp
         n_ln(k) = nv_conc/dlnr
         g(k)   =  vv_conc/dlnr
      enddo
      
      call coagmax(rr,n_ln,dlnr,tot_free_max)
      r_samp = - (tot_free_max * 1/V_comp *del_T/log(1-p_max))
      n_samp = r_samp * M*(M-1)
      
      if (tlmin.eq.0. .or.tlmin .ge. 60. ) then
         tlmin = tlmin -60.
         lmin = lmin + 1
         write(30,*)'Time = ',TIME
         sum_masse = 0.
         do k=1,n_bin
            write(30,'(i4,6e14.5)')k,
     &           dp(k)/2.,n_ln(k),g(k)
            sum_masse = sum_masse + g(k)*dlnr
         enddo
      endif
      
      return
      end
