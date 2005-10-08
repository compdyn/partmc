      program analy

      integer n,  scal, tmax
      parameter (n=160, scal = 4)
      real*8 nn(3600,n),e(n),r(n),v(n)
      real*8 dlnr, ax
      real   rho, beta_1, T
      real*8 bessi1
      real*8 N_0, v_0, tau
      real*8 n_ini(n), n_i(n), dp(n)
      real*8 n_ln(3600,n),n_ln0(n)
      parameter (rho = 1000., beta_1 = 1000.)
      data emin,pi,tmax/1.e-12,3.141592654,600./
      
      dlnr=dlog(2.d0)/(3.*scal)     
      ax=2.d0**(1.0/scal)
c      const = 4./3.*pi
      const = 1.
c      N_0 = 1.e+09/const
c      v_0 = const*1.e-15
      N_0 = 1.e+09/const
      v_0 = 4.1886e-15 
      d_0 = (6*v_0/pi)**(1./3.) 


c mass and radius grid
      e(1)=emin*0.5*(ax+1.)
      r(1)=1000.*dexp(dlog(3.*e(1)/(4.*pi))/3.)
      v(1)=1.e-06*e(1)/rho
      n_ini(1) = N_0/v_0 * exp(-v(1)/v_0)
      dp(1) = 1.e-06*2*r(1)

      write(6,*)'N_0 in cm-3, v_0 ',N_0*1.e-06, v_0, N_0/v_0

      do i=2,n
         e(i)=ax*e(i-1)
         r(i)=1000.*dexp(dlog(3.*e(i)/(4.*pi))/3.)
         v(i)=1.e-06*e(i)/rho
         n_ini(i)= N_0/v_0 * exp(-v(i)/v_0)
         dp(i) = 1.e-06*2*r(i)
      enddo

      sum=0.
      do i=1,n-1
         n_i(i) = d_0/N_0 * n_ini(i)*pi/2.*dp(i)**2.
         n_ln0(i)=pi/2. * dp(i)**3. * n_ini(i)
         sum = sum+n_ln0(i)*dlnr
      enddo
      write(6,*)'n_tot ',sum


      do ij=1,tmax
         
         t = ij
         
         tau = N_0*v_0*beta_1*t
         T = 1-exp(-tau)
         write(6,*) 
         write(6,*)'tau ',tau,T

         do i=1,n
            rat_v = v(i)/v_0
            nn(ij,i)  = N_0/v(i) * (1-T)/sqrt(T) * exp(-(1+T)*rat_v)*
     *           bessi1(2*rat_v*sqrt(T))
            n_ln(ij,i) = pi/2. * dp(i)**3. * nn(ij,i)

            write(6,*)'nn ', N_0/v(i),bessi1(2*rat_v*sqrt(T))
        enddo

      enddo

      do ij=1,tmax
         do i=1,n-1
            nn(ij,i) = d_0/N_0*nn(ij,i)*pi/2.*dp(i)**2.
         enddo
      enddo

      do i=1,n-1
c     write(99,'(200e14.5)')dp(i)/d_0,n_i(i),(nn(ij,i),ij=1,tmax,5)
         write(98,'(602e14.5)')dp(i)/2.*1.e+06,
     *            n_ln0(i),(n_ln(ij,i),ij=1,tmax)
      enddo

      end

c ***************************************************************      

      function bessi1(x)
      real bessi1, x
      real ax
      double precision p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,
     *                 q7,
     *                 q8,q9,y
      save p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,
     *                 q8,q9
      data p1,p2,p3,p4,p5,p6,p7/0.5d0,0.87890594d0,0.51498869d0,
     *     0.15084934d0,0.2658733d-1,0.301532d-2,0.32411d-3/
      data q1,q2,q3,q4,q5,q6,q7,q8,q9/0.39894228d0,-0.3988024d-1,
     *     -0.362018d-2,0.163801d-2,-0.1031555d-1,0.2282967d-1,
     *     -0.2895312d-1,0.1787654d-1,-0.420059d-2/

      if(abs(x).lt.3.75) then
         y=(x/3.75)**2
         bessi1 = x*(p1+y*(p2+y*(p3+y*(p4+y*(p5+y*(p6+y*p7))))))
      else
         ax = abs(x)
         y=3.75/ax
         bessi1 = (exp(ax)/sqrt(ax))*(q1+y*(q2+y*(q3+y*(q4+
     *        y*(q5+y*(q6+y*(q7+y*(q8+y*q9))))))))
         if (x.lt.0.) bessi1 = -bessi1
      endif
      return
      end
