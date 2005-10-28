C Simulation with sedimentation kernel, fixed timestepping and split array.

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      program MonteCarlo

      ! MM = 1e7, t = 480, r_s = 2.2e-5: 32 minutes
      integer MM, n_bin, n_loop, scal
      real*8 t_max, rho_p, N_0, t_print, t_progress
      real*8 del_t, V_0, r_split
      parameter (MM = 100000000)       ! number of particles
      parameter (n_bin = 160)          ! number of bins
      parameter (n_loop = 1)           ! number of loops
      parameter (scal = 3)             ! scale factor for bins
      parameter (t_max = 480d0)        ! total simulation time (seconds)
      parameter (rho_p = 1000d0)       ! particle density (kg/m^3)
      parameter (N_0 = 1d9)            ! particle number concentration (#/m^3)
      parameter (t_print = 60d0)       ! interval between printing (s)
      parameter (t_progress = 1d0)     ! interval between progress (s)
      parameter (del_t = 1d0)          ! timestep (s)
      parameter (V_0 = 4.1886d-15)     ! mean volume of initial distribution (m^3)
      parameter (r_split = 2.2e-5)     ! radius to split at (m)
      !parameter (r_split = 2e-6)       ! radius to split at (m)
      !parameter (r_split = 2e-4)       ! radius to split at (m)

      integer M, i_loop, s_bin, MS
      real*8 V(MM), V_comp, dlnr
      real*8 bin_v(n_bin), bin_r(n_bin)
      real*8 bin_g(n_bin)
      integer n_ini(n_bin), bin_n(n_bin)

      external kernel_sedi

      open(30,file='out_sedi_fix_split.d')
      call print_header(n_loop, n_bin, nint(t_max / t_print) + 1)
      call srand(17)

      do i_loop = 1,n_loop

         call make_grid(n_bin, scal, rho_p, bin_v, bin_r, dlnr)

         call find_split_bin(n_bin, bin_v, r_split, s_bin)
         call init_exp(MM, V_0, dlnr, n_bin, bin_v, bin_r, n_ini)
         !call init_bidisperse(MM, n_bin, n_ini)
         call compute_volumes(n_bin, MM, n_ini, bin_r, dlnr, V, M)
         V_comp = M / N_0
         call find_split_particle(MM, M, V, n_bin, s_bin, bin_v, MS)

         call mc_fix_split(MM, M, MS, V, V_comp, n_bin, s_bin, bin_v,
     $        bin_r, bin_g, bin_n, dlnr, kernel_sedi, t_max, t_print,
     $        t_progress , del_t, i_loop)

      enddo

      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
