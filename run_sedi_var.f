C Simulation with sedimentation kernel and variable timestepping.

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      program MonteCarlo
 
      integer MM, n_bin, n_loop, scal, k_avg_samp
      real*8 t_max, rho_p, N_0, t_print, t_k_max, t_k_avg
      parameter (MM = 10000)        ! number of particles
      parameter (n_bin = 160)       ! number of bins
      parameter (n_loop = 1)        ! number of loops
      parameter (scal = 3)          ! scale factor for bins
      parameter (t_max = 600d0)     ! total simulation time (seconds)
      parameter (rho_p = 1000d0)    ! particle density (kg/m^3)
      parameter (N_0 = 1d9)         ! particle number concentration (#/m^3)
      parameter (k_avg_samp = 1000) ! number of samples to estimate k_avg
      parameter (t_print = 60)      ! interval between printing (s)
      parameter (t_k_max = 60)      ! interval between estimating k_max (s)
      parameter (t_k_avg = 0.2d0)   ! interval between estimating k_avg (s)

      integer M, M_comp, i_loop, k
      real*8 V(MM), V_comp, dlnr
      real*8 n_ini(n_bin), vv(n_bin), rr(n_bin)
      real*8 g(n_bin), n_ln(n_bin)

      external kernel_sedi

      open(30,file='out_sedi_var.d')
      call print_header(n_loop, n_bin, nint(t_max / t_print) + 1)
      call srand(10)

      do i_loop = 1,1

         call make_grid(n_bin, scal, rho_p, vv, rr, dlnr)
         
         ! define bidisperse distribution
         do k = 1,n_bin
            n_ini(k) = 0d0
         enddo
         n_ini(97) = (MM - 1) / dlnr
         n_ini(126) = 1 / dlnr

         call compute_volumes(n_bin, MM, n_ini, rr, dlnr, V, M_comp)
         M = M_comp
         V_comp = M / N_0
         
         call mc_var(MM, M, M_comp, V, V_comp, kernel_sedi, n_bin, vv,
     &        rr, g, n_ln, dlnr, t_max, t_print, t_k_max, t_k_avg,
     &        k_avg_samp)

      enddo

      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
