C Simulation with Golovin kernel and fixed timestepping.

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      program MonteCarlo
 
      integer MM, MM_1, n_bin, n_spec, n_loop, scal
      real*8 t_max, del_t, rho_p, N_0, t_print
      real*8 V_01, V_02, v_min
      parameter (MM = 10000)       ! total number of particles
      parameter (MM_1 = 2000)      ! number of #1-particles 
      parameter (n_bin = 160)      ! number of bins
      parameter (n_spec = 3)       ! number of species
      parameter (n_loop = 10)      ! number of loops
      parameter (scal = 3)         ! scale factor for bins
      parameter (t_max = 600d0)    ! total simulation time (seconds)
      parameter (del_t = 1d0)      ! timestep (seconds)
      parameter (rho_p = 1000d0)   ! particle density (kg m^{-3})
      parameter (v_min = 1.e-24)   ! minimum volume (m^3) for making grid
      parameter (N_0 = 1d9)        ! particle number concentration (#/m^3)
      parameter (V_01 = 4.1886d-15)! mean volume of #1-initial distribution
      parameter (V_02 = 8*V_01)    ! mean volume of #2-initial distribution
      parameter (t_print = 60)     ! interval between printing (s)

      integer M, M1, M2, i_loop
      real*8 V(MM,n_spec), V_comp, dlnr
      real*8 bin_v(n_bin), bin_r(n_bin)
      real*8 bin_g(n_bin), bin_gs(n_bin,n_spec),vol_frac(n_spec)
      integer n_ini(n_bin), bin_n(n_bin)

      external kernel_golovin

      open(30,file='out_golovin_fix.d')
      call print_header(n_loop, n_bin, n_spec, nint(t_max / t_print) + 1
     $     )
      call srand(10)

      do i_loop = 1,n_loop

         call make_grid(n_bin, scal, v_min, bin_v, bin_r, dlnr)
         call zero_v(MM,n_spec,V)
cn *** initialise first distribution
         vol_frac(1) = 1
         vol_frac(2) = 0
         vol_frac(3) = 0
         call init_exp(MM_1, V_01, dlnr, n_bin, bin_v, bin_r, n_ini)
         call compute_volumes(n_bin, n_spec, vol_frac, MM, 1,MM_1,
     $        n_ini,bin_r, dlnr, V, M1)

cn *** initialise second distribution
         call init_exp(MM-MM_1, V_02, dlnr, n_bin, bin_v, bin_r, n_ini)
         vol_frac(1) = 0
         vol_frac(2) = 0.5
         vol_frac(3) = 0.5
         call compute_volumes(n_bin, n_spec, vol_frac, MM, M1+1,
     $        MM, n_ini, bin_r, dlnr, V, M2)

         M=M1+M2
         V_comp = M / N_0

         call mc_fix(MM, M, V, V_comp,n_spec,
     &        n_bin, bin_v, bin_r, bin_g, bin_gs, bin_n, dlnr,
     &        kernel_golovin, t_max, del_t, t_print, i_loop)

      enddo

      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
