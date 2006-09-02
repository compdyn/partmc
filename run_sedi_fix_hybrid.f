C Simulation with sedimentation kernel, fixed timestepping and hybrid array.

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      program MonteCarlo

      use mod_array
      use mod_init_dist
      use mod_mc_fix_hybrid
      use mod_kernel_sedi
      use mod_condensation
      integer MM, MM_1, TDV, n_bin, n_spec, n_loop, scal, i_water
      real*8 t_max, N_0, t_print, t_progress
      real*8 del_t, del_t_cond, V_01, V_02, v_min
      real*8 d_mean1, d_mean2, log_sigma1, log_sigma2

      parameter (MM =  100000)  ! number of particles
      parameter (TDV =  50000)  ! trailing dimension of VH
      parameter (MM_1 = MM/2)   ! number of #1-particles
      parameter (n_bin = 160)   ! number of bins
      parameter (n_spec = 3)    ! number of species
      parameter (n_loop = 1)    ! number of loops
      parameter (scal = 3)      ! scale factor for bins
      parameter (t_max = 2d0)   ! total simulation time (seconds)
      parameter (v_min = 1.d-24) ! minimum volume (m^3) for making grid
      parameter (N_0 = 1d9)     ! particle number concentration (#/m^3)
      parameter (t_print = 1d0) ! interval between printing (s)
      parameter (t_progress = 1d0) ! interval between progress (s)
      parameter (del_t = 1d0)   ! timestep (s)
      parameter (V_01 = 8.d-2*4.1886d-15) ! mean volume of initial distribution (m^3)
      parameter (V_02 = V_01/8.d0) ! mean volume of #2-initial distribution (m^3)
      parameter (d_mean1 = 0.2d-6) ! mean diameter of #1- initial distribution (m)
      parameter (d_mean2 = 0.2d-6)  ! mean diameter of #2- initial distribution (m)
      parameter (log_sigma1 = 0.25d0) ! log(sigma) of #1- initial distribution
      parameter (log_sigma2 = 0.25d0) ! log(sigma) of #2- initial distribution
      parameter (i_water = 3)   ! water species number

      integer M, M1, M2, i_loop, i
      real*8 V(MM,n_spec), V_comp, dlnr, VH(n_bin,TDV,n_spec)
      real*8 bin_v(n_bin), bin_r(n_bin)
      real*8 bin_g(n_bin), bin_gs(n_bin,n_spec),vol_frac(n_spec)
      real*8 rho_p(n_spec)
      real*8 eps(n_spec), M_w(n_spec)
      real*8 RH_eq
      integer nu(n_spec)
      integer n_ini(n_bin), bin_n(n_bin), MH(n_bin)

      parameter (RH_eq = 0.99d0)                  ! INPUT: equilibrium RH for initial distribution`
      data rho_p / 1800.d0, 1800.d0, 1000.d0 /  ! INPUT: density of species (kg m^{-3})
      data nu / 3, 3, 0 /                          ! INPUT: number of ions in the solute (1)
      data eps / 0.25d0, 0.25d0, 0d0 /               ! INPUT: solubility of solutes (1)
      data M_w / 132d-3, 132d-3, 18d-3 /        ! INPUT: molecular weight of species (kg mole^{-1})

      open(30,file='out_sedi_fix_hybrid.d')
      call print_header(n_loop, n_bin, n_spec, 
     %     nint(t_max / t_print) + 1)
      call srand(17)

      do i_loop = 1,n_loop

         call make_grid(n_bin, scal, v_min, bin_v, bin_r, dlnr)
         call zero_v(MM,n_spec,V)
cn *** initialize first distribution
c         vol_frac(1) = 0.5d0
c         vol_frac(2) = 0.5d0
c         vol_frac(3) = 0.d0
c         call init_exp(MM_1, V_01, dlnr, n_bin, bin_v, bin_r, n_ini)
c         call compute_volumes(n_bin, n_spec, vol_frac, MM, 1,MM_1,
c     &        n_ini, bin_v, dlnr, V, M1)

cn *** initialise second distribution
c         call init_exp(MM-MM_1, V_02, dlnr, n_bin, bin_v, bin_r, n_ini)
c         vol_frac(1) = 0.d0
c         vol_frac(2) = 0.5d0
c         vol_frac(3) = 0.5d0
c         call compute_volumes(n_bin, n_spec, vol_frac, MM, M1+1,
c     $        MM, n_ini, bin_v, dlnr, V, M2)


cn *** initialize first distribution
         call init_log_normal(MM_1, d_mean1, log_sigma1, dlnr, n_bin,
     &     bin_v, bin_r, n_ini)
         vol_frac(1) = 1d0
         vol_frac(2) = 0d0
         vol_frac(3) = 0d0
         call compute_volumes(n_bin, n_spec, vol_frac, MM, 1,MM_1,
     &        n_ini, bin_v, dlnr, V, M1)

cn *** initialise second distribution
         call init_log_normal(MM-MM_1, d_mean2, log_sigma2, dlnr, n_bin,
     &     bin_v, bin_r, n_ini)
         vol_frac(1) = 0d0
         vol_frac(2) = 1d0
         vol_frac(3) = 0d0
         call compute_volumes(n_bin, n_spec, vol_frac, MM, M1+1,
     $        MM, n_ini, bin_v, dlnr, V, M2)

         M=M1+M2
         V_comp = dble(M) / N_0

!     call equlibriate_particle for each particle in V
         do i = 1,M
            call equilibriate_particle(n_spec, V(i,:), rho_p, 
     &           i_water, nu,
     &           eps, M_w, RH_eq)
         enddo

         call mc_fix_hybrid(MM, M, V, n_spec, n_bin, 
     &        TDV, MH, VH, V_comp, bin_v, rho_p, i_water,
     $        bin_r, bin_g, bin_gs, bin_n, dlnr, 
     &        kernel_sedi, t_max, t_print,
     $        t_progress ,del_t, i_loop)

      enddo

      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
