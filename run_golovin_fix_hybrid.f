C Simulation with Golovin kernel and fixed timestepping.

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      program run_golovin_fix_hybrid
 
      use mod_bin
      use mod_array
      use mod_init_dist
      use mod_mc_fix_hybrid
      use mod_kernel_golovin
      use mod_condensation
      use mod_environ
      use mod_material
      use mod_constants
      use mod_util

      integer MM, n_bin, n_spec, n_loop, scal, TDV
      real*8 t_max, del_t, rho_p, N_0, t_print, t_progress
      real*8 V_0, v_min
      parameter (MM = 10000)       ! total number of particles
      parameter (TDV =  10000)     ! trailing dimension of VH
      parameter (n_bin = 160)      ! number of bins
      parameter (n_spec = 1)       ! number of species
      parameter (n_loop = 10)      ! number of loops
      parameter (scal = 3)         ! scale factor for bins

      parameter (t_max = 1300d0)   ! total simulation time (seconds)
      parameter (t_print = 100d0)  ! interval between printing (s)
      parameter (t_progress = 1d0) ! interval between progress (s)
      parameter (del_t = 1d0)      ! timestep (s)

      parameter (v_min = 1d-24)    ! minimum volume (m^3) for making grid
      parameter (N_0 = 1d9)        ! particle number concentration (#/m^3)
      parameter (V_0 = 4.1886d-15) ! mean volume of #1-initial distribution

      integer M, i_loop
      real*8 V(MM,n_spec), dlnr, VH(n_bin,TDV,n_spec)
      real*8 bin_v(n_bin), bin_r(n_bin), n_den(n_bin)
      real*8 bin_g(n_bin), bin_gs(n_bin,n_spec), vol_frac(n_spec)
      integer n_ini(n_bin), bin_n(n_bin), MH(n_bin)
      type(environ) :: env
      type(material) :: mat

      call allocate_material(mat, n_spec)
      mat%i_water = 1
      mat%rho = (/ 1000d0 /)
      mat%nu = (/ 0 /)
      mat%eps = (/ 0d0 /)
      mat%M_w = (/ 18d-3 /)

      env%T = 288d0        ! (K)
      env%RH = 0.999d0     ! (1)
      env%p = 1d5          ! (Pa)
      env%dTdt = -0.01d0   ! (K s^{-1})

      open(30,file='out_golovin_fix.d')
      call print_header(n_loop, n_bin, n_spec, nint(t_max / t_print) + 1
     $     )
      call srand(10)

      do i_loop = 1,n_loop

         call make_bin_grid(n_bin, scal, v_min, bin_v, bin_r, dlnr)
         call zero_v(MM, n_spec, V)
         vol_frac(1) = 1d0
         call init_exp(V_0, n_bin, bin_v, bin_r, n_den)
         call dist_to_n(MM, dlnr, n_bin, bin_v, bin_r, n_den, bin_n)
         call compute_volumes(n_bin, n_spec, vol_frac, MM, 1, MM,
     $        n_ini,bin_v, dlnr, V, M)

         env%V_comp = dble(M) / N_0

         call mc_fix_hybrid(MM, M, n_spec, V, n_bin, TDV, MH, VH,
     $        bin_v, bin_r, bin_g, bin_gs, bin_n, dlnr,
     $        kernel_golovin, t_max, t_print, t_progress, del_t, i_loop,
     $        env, mat)

      enddo

      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
