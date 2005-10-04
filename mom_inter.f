C &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&


      subroutine moments(V,M,N_tot,M_comp,V_comp,
     &                       TIME,tlmin,del_T,V_bar,
     &                       Time_count,rr,
     &                       tot_free_max,vv,dlnr,dp)

      integer M,MM,n_bin,Time_count,NN_cnt,M_comp
      real*8 nv_conc,dlnr, tlmin
      parameter (MM=10000000,n_bin=160)
      real*8 V(MM),N_tot,vv(n_bin),dp(n_bin)
      real*8 V_comp,del_T,TIME,tot_free_max
      real*8 V_bar
      real*8 eta(n_bin),rr(n_bin)
      real*8 psi(n_bin),g(n_bin)
      real*8 n_ln(n_bin),V_0,d_0
      real*8 pi,rho_w, sum_masse, vv_cnt, vv_conc
      integer k, i

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
       integer imax, jmax, ll, k

       do k=1,n_bin
          v_bin(k)=4./3.*pi*rr(k)**3.
       enddo

       tot_free_max =0.
       do k=1,n_bin
          if(n_ln(k)*dlnr .ge. 1.) then
             do ll=1,k
                 call coag_kernel(v_bin(k), v_bin(ll), cck)
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
