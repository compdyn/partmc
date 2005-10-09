C Simulation with Golovin kernel and variable timestepping.

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      program MonteCarlo
 
      integer MM, n_bin, n_loop, scal, k_avg_samp
      real*8 t_max, rho_p, N_0, t_print, t_k_max, t_k_avg, V_0
      parameter (MM = 10000)        ! number of particles
      parameter (n_bin = 160)       ! number of bins
      parameter (n_loop = 1)         ! number of loops
      parameter (scal = 3)          ! scale factor for bins
      parameter (t_max = 600.)      ! total simulation time (seconds)
      parameter (rho_p = 1000.)     ! particle density (kg/m^3)
      parameter (N_0 = 1d9)         ! particle number concentration (#/m^3)
      parameter (k_avg_samp = 1000) ! number of samples to estimate k_avg
      parameter (t_print = 60)      ! interval between printing (s)
      parameter (t_k_max = 60)      ! interval between estimating k_max (s)
      parameter (t_k_avg = 0.2)     ! interval between estimating k_avg (s)
      parameter (V_0 = 4.1886e-15)   !

      integer M, M_comp, i_loop, k
      real*8 V(MM), V_comp, dlnr, t1
      real*8 n_ini(n_bin), vv(n_bin), dp(n_bin), rr(n_bin)
      real*8 g(n_bin), n_ln(n_bin)

      real*8 pi
      parameter (pi = 3.14159265358979323846)

      external kernel_golovin

      open(30,file='out_golovin_var.d')
      call print_header(n_loop, n_bin, nint(t_max / t_print + 1.))
      call srand(10)

      do i_loop = 1,1
         call cpu_time(t1)
         write(30,*)'i_loop=',i_loop,t1

         M = MM
         M_comp = M
         V_comp = M / N_0
         
         call make_grid(n_bin, scal, rho_p, vv, dp, rr, dlnr)
         
         ! define initial exponential distribution
         do k = 1,n_bin
            n_ini(k) = pi/2.0 * dp(k)**3.0 * M/V_0 * exp(-(vv(k)/V_0))
         enddo

         call compute_volumes(n_bin, MM, n_ini, dp, dlnr, V)

         call mc_var(MM, M, M_comp, V, V_comp, kernel_golovin, n_bin,
     &        vv, rr, dp, g, n_ln, dlnr, t_max, t_print, t_k_max,
     &        t_k_avg, k_avg_samp)

      enddo

      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
