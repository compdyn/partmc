C Simulation with Golovin kernel and adaptive timestepping.

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      program MonteCarlo
 
      integer MM, n_bin, n_loop, scal
      real*8 t_max, rho_p, N_0, t_print, V_0
      real*8 r_samp_max, del_t_max
      parameter (MM = 10000)           ! number of particles
      parameter (n_bin = 160)          ! number of bins
      parameter (n_loop = 10)          ! number of loops
      parameter (scal = 3)             ! scale factor for bins
      parameter (t_max = 600d0)        ! total simulation time (seconds)
      parameter (rho_p = 1000d0)       ! particle density (kg/m^3)
      parameter (N_0 = 1d9)            ! particle number concentration (#/m^3)
      parameter (t_print = 60d0)       ! interval between printing (s)
      parameter (r_samp_max = 0.005d0) ! maximum sampling ratio per timestep
      parameter (del_t_max = 1d0)      ! maximum timestep
      parameter (V_0 = 4.1886d-15)     ! mean volume of initial distribution

      integer M, i_loop
      real*8 V(MM), V_comp, dlnr
      real*8 bin_v(n_bin),  bin_r(n_bin)
      real*8 bin_g(n_bin)
      real*8 n_den(n_bin)
      integer n_ini(n_bin), bin_n(n_bin)

      external kernel_golovin

      open(30,file='out_golovin_adapt.d')
      call print_header(n_loop, n_bin, nint(t_max / t_print) + 1)
      call srand(10)

      do i_loop = 1,n_loop
         
         call make_grid(n_bin, scal, rho_p, bin_v, bin_r, dlnr)
         call init_exp(V_0, n_bin, bin_v, bin_r, n_den)
         call dist_to_n(MM, dlnr, n_bin, bin_v, bin_r, n_den, bin_n)
         call compute_volumes(n_bin, MM, n_ini, bin_r, dlnr, V, M)
         V_comp = M / N_0

         call mc_adapt(MM, M, V, V_comp,
     &        n_bin, bin_v, bin_r, bin_g, bin_n, dlnr,
     &        kernel_golovin, t_max, t_print, r_samp_max,
     &        del_t_max, i_loop)

      enddo

      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
