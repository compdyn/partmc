C Simulation with Golovin kernel and fixed timestepping.

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      program MonteCarlo
 
      integer MM, n_bin, scal
      real*8 t_max, del_t, rho_p, N_tot, p_max, V_0
      parameter (MM = 10000)       ! number of particles
      parameter (n_bin = 160)      ! number of bins
      parameter (scal = 3)         ! scale factor for bins
      parameter (t_max = 600.)     ! total simulation time (seconds)
      parameter (del_t = 1.)       ! timestep (seconds)
      parameter (rho_p = 1000.)    ! particle density (kg m^{-3})
      parameter (N_tot = 1.e+9)    ! particle number concentration (#/m^3)
      parameter (p_max = 0.01)     ! maximum coagulation probability
      parameter (V_0 = 4.1886e-15) !

      integer M, M_comp, i_loop, k
      real*8 V(MM), V_comp, dlnr, t1
      real*8 n_ini(n_bin), vv(n_bin), dp(n_bin), rr(n_bin)

      real*8 pi
      parameter (pi = 3.14159265358979323846)

      external kernel_golovin

      open(30,file='out_golovin_fix.d')
      call srand(10)

      do i_loop = 1,1
         call cpu_time(t1)
         write(6,*)'START ',i_loop, t1
         write(30,*)'i_loop=',i_loop,t1

         M = MM
         M_comp = M
         V_comp = M / N_tot
         
         call make_grid(n_bin, scal, rho_p, vv, dp, rr, dlnr)
         
c     define initial exponential distribution
         do k = 1,n_bin
            n_ini(k) = pi/2.0 * dp(k)**3.0 * M/V_0 * exp(-(vv(k)/V_0))
         enddo

         call compute_volumes(n_bin, MM, n_ini, dp, dlnr, V)

         call mc_fix(MM, M, M_comp, V, V_comp, kernel_golovin, n_bin,
     &        vv, rr, dp, dlnr, t_max, del_t, p_max)

      enddo

      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
