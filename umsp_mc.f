       program umsp_mc
       integer n_bin,i_loop,tmax
       parameter (n_bin=160,i_loop=1,tmax=11)

       real*8 eta(n_bin),g(i_loop,tmax,n_bin),n(i_loop,tmax,n_bin)
       real*8 g_bar(tmax,n_bin),n_bar(tmax,n_bin)
       integer k, l, i, mm
       character*20 dum

       open(20,file='mc_fix.d')
       open(30,file='mc_n.dtest2')
       open(40,file='mc_m.dtest2')
       open(50,file='mcr_n.dtest2')
       open(60,file='mcr_m.dtest2')
      

       do k=1,tmax
          do l=1,n_bin
             n_bar(k,l) = 0.
             g_bar(k,l) = 0.
          enddo
       enddo

       do i=1,i_loop

          read(20,'(a)')dum
          write(6,'(a)')dum
          do k=1,tmax
             read(20,'(a)')dum
             write(6,'(a)')dum
                 do l=1,n_bin
                    read(20,'(i4,5e14.5)')
     *                      mm,eta(l),n(i,k,l),g(i,k,l)
                    write(6,*)mm,eta(l),n(i,k,l),g(i,k,l)
                  enddo
           enddo 
        enddo

       do k=1,tmax
          do l=1,n_bin
             do i=1,i_loop
                    g_bar(k,l) = g_bar(k,l)+g(i,k,l)
                    n_bar(k,l) = n_bar(k,l)+n(i,k,l)
             enddo
          
             g_bar(k,l) = g_bar(k,l)/i_loop
             n_bar(k,l) = n_bar(k,l)/i_loop
          enddo
        enddo

        do k=1,tmax
           write(50,*)'time= ',k-1
           write(60,*)'time= ',k-1
       
           do l=1,n_bin
              write(50,'(i4,100e14.5)')l,eta(l),(n(i,k,l),i=1,i_loop),
     &                 n_bar(k,l)
              write(60,'(i4,100e14.5)')l,eta(l),(g(i,k,l),i=1,i_loop),
     &                 g_bar(k,l)
           enddo
        enddo
 
       do l=1,n_bin
            write(30,'(i4,100e14.5)')l,eta(l)*1.e+06
     &           ,(n_bar(k,l),k=1,tmax)
            write(40,'(i4,100e14.5)')l,eta(l)*1.e+06
     &           ,(g_bar(k,l)*1000.,k=1,tmax)
       enddo


        end
       
     
