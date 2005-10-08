C Simulation with sedimentation kernel and adaptive timestepping.

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      program MonteCarlo
 
      integer MM, n_bin, scal
      real*8 t_max, rho_p, N_tot, t_print, p_max
      real*8 r_samp_max, del_t_max
      parameter (MM = 10000)         ! number of particles
      parameter (n_bin = 160)        ! number of bins
      parameter (scal = 3)           ! scale factor for bins
      parameter (t_max = 600.)       ! total simulation time (seconds)
      parameter (rho_p = 1000.)      ! particle density (kg/m^3)
      parameter (N_tot = 1.e+9)      ! particle number concentration (#/m^3)
      parameter (t_print = 60)       ! interval between printing (s)
      parameter (p_max = 0.01)       ! maximum coagulation probability
      parameter (r_samp_max = 0.005) ! maximum sampling ratio per timestep
      parameter (del_t_max = 1.0)    ! maximum timestep

      integer M, M_comp, i_loop, k
      real*8 V(MM), V_comp, dlnr, t1
      real*8 n_ini(n_bin), vv(n_bin), dp(n_bin), rr(n_bin)
      real*8 g(n_bin), n_ln(n_bin)

      external kernel_sedi

      open(30,file='out_sedi_adapt.d')
      call srand(10)

      do i_loop = 1,1
         call cpu_time(t1)
         write(6,*)'START ',i_loop, t1
         write(30,*)'i_loop=',i_loop,t1

         M = MM
         M_comp = M
         V_comp = M / N_tot
         
         call make_grid(n_bin, scal, rho_p, vv, dp, rr, dlnr)
         
         ! define bidisperse distribution
         do k = 1,n_bin
            n_ini(k) = 0.
         enddo
         n_ini(97) = (M-1)/dlnr
         n_ini(126) = 1/dlnr

         call compute_volumes(n_bin, MM, n_ini, dp, dlnr, V)

         call mc_adapt(MM, M, M_comp, V, V_comp, kernel_sedi, n_bin,
     &        vv, rr, dp, g, n_ln, dlnr, t_max, t_print,
     &        p_max, r_samp_max, del_t_max)

      enddo

      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
