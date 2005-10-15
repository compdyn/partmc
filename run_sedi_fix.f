C Simulation with sedimentation kernel and fixed timestepping.

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      program MonteCarlo
 
      integer MM, n_bin, n_loop, scal
      real*8 t_max, del_t, rho_p, N_0, p_max, t_print
      parameter (MM = 10000)       ! number of particles
      parameter (n_bin = 160)      ! number of bins
      parameter (n_loop = 1)         ! number of loops
      parameter (scal = 3)         ! scale factor for bins
      parameter (t_max = 600.)     ! total simulation time (seconds)
      parameter (del_t = 1.)       ! timestep (seconds)
      parameter (rho_p = 1000.)    ! particle density (kg m^{-3})
      parameter (N_0 = 1d9)        ! particle number concentration (#/m^3)
      parameter (p_max = 0.01)     ! maximum coagulation probability
      parameter (t_print = 60)       ! interval between printing (s)

      integer M, M_comp, i_loop, k
      real*8 V(MM), V_comp, dlnr
      real*8 n_ini(n_bin), vv(n_bin), dp(n_bin), rr(n_bin)
      real*8 g(n_bin), n_ln(n_bin)

      external kernel_sedi

      open(30,file='out_sedi_fix.d')
      call print_header(n_loop, n_bin, nint(t_max / t_print + 1.))
      call srand(10)

      do i_loop = 1,1

         call make_grid(n_bin, scal, rho_p, vv, dp, rr, dlnr)
         
         ! define bidisperse distribution
         do k = 1,n_bin
            n_ini(k) = 0.
         enddo
         n_ini(97) = (M-1)/dlnr
         n_ini(126) = 1/dlnr

         call compute_volumes(n_bin, MM, n_ini, dp, dlnr, V, M_comp)
         M = M_comp
         V_comp = M / N_0
         
         call mc_fix(MM, M, M_comp, V, V_comp, kernel_sedi, n_bin, vv,
     &        rr, dp, g, n_ln, dlnr, t_max, del_t, p_max, t_print)

      enddo

      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
