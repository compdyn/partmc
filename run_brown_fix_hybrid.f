! -*- mode: f90; -*-
! Copyright (C) 2005,2006 Nicole Riemer and Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.
!
! Simulation with brownian kernel, fixed timestepping and hybrid array.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program run_brown_fix_hybrid
  
  use mod_bin
  use mod_array
  use mod_array_hybrid
  use mod_init_dist
  use mod_mc_fix_hybrid
  use mod_kernel_brown
  use mod_condensation
  use mod_environ
  use mod_material
  use mod_constants
  use mod_util
  
  ! species #1 is salt, #2 is dust, and #3 is water
  
  integer MM, MM_1, n_bin, n_spec, n_loop, scal, i_water
  real*8 t_max, N_0, t_print, t_state, t_progress
  real*8 del_t, del_t_cond, v_min
  real*8 d_mean1, d_mean2, log_sigma1, log_sigma2
  
  parameter (MM =  100000)  ! number of particles
  parameter (MM_1 = MM/2)   ! number of #1-particles
  parameter (n_bin = 160)   ! number of bins
  parameter (n_spec = 3)    ! number of species
  parameter (n_loop = 1)    ! number of loops
  parameter (scal = 3)      ! scale factor for bins
  parameter (v_min = 1d-24) ! minimum volume (m^3) for making grid
  parameter (N_0 = 1d11)     ! particle number concentration (#/m^3)
  
  parameter (t_max = 3d0*3600d0)  ! total simulation time (seconds)
  parameter (t_print = 600d0) ! interval between output (s)
  parameter (t_state = 0d0)   ! interval between state output (s)
  parameter (t_progress = 1d0) ! interval between printing progress (s)
  parameter (del_t = 1d0)   ! timestep (s)
  parameter (d_mean1 = 0.1d-6) ! mean diameter of #1- initial distribution (m)
  parameter (d_mean2 = 0.05d-6)  ! mean diameter of #2- initial distribution (m)
  parameter (log_sigma1 = 0.21d0) ! log(sigma) of #1- initial distribution
  parameter (log_sigma2 = 0.2d0) ! log(sigma) of #2- initial distribution
  
  integer M, M1, M2, i_loop, i
  real*8 V(MM,n_spec), dlnr
  type(bin_p) VH(n_bin)
  real*8 bin_v(n_bin), n_den(n_bin)
  real*8 bin_g(n_bin), bin_gs(n_bin,n_spec), vol_frac(n_spec)
  integer n_ini(n_bin), bin_n(n_bin), MH(n_bin)
  type(environ) :: env
  type(material) :: mat

  call init_hybrid(n_spec, MH, VH)

  call allocate_material(mat, n_spec)
  mat%i_water = 3
  mat%rho = (/ 2165d0, 2650d0, 1000d0 /)
  mat%nu = (/ 2, 2, 0 /)
  mat%eps = (/ 1d0, 0.5d0, 0d0 /)
  mat%M_w = (/ 58.44d-3, 60.08d-3, 18d-3 /)
  
  env%T = 288d0        ! (K)
  env%RH = 0.70d0      ! (1)
  env%p = 1d5          ! (Pa)
  env%dTdt = -0.00d0   ! (K s^{-1})
  open(30,file='out_brown_fix_hybrid.d')

  call print_header(n_loop, n_bin, n_spec,  &
       nint(t_max / t_print) + 1)
  call srand(17)
  !call srand(time())
  do i_loop = 1,n_loop
     call make_bin_grid(n_bin, scal, v_min, bin_v, dlnr)
     call zero_v(MM,n_spec,V)
     
     ! n *** initialize first distribution
     call init_log_normal(d_mean1, log_sigma1, n_bin, &
          bin_v, n_den)
     call dist_to_n(MM_1, dlnr, n_bin, bin_v, n_den, n_ini)
     vol_frac(1) = 1d0
     vol_frac(2) = 0d0
     vol_frac(3) = 0d0
     call compute_volumes(n_bin, n_spec, vol_frac, MM, 1,MM_1, &
          n_ini, bin_v, dlnr, V, M1)
     ! n *** initialise second distribution
     call init_log_normal(d_mean2, log_sigma2, n_bin, &
          bin_v, n_den)
     call dist_to_n(MM-MM_1, dlnr, n_bin, bin_v, n_den, n_ini)
     vol_frac(1) = 0d0
     vol_frac(2) = 1d0
     vol_frac(3) = 0d0
     call compute_volumes(n_bin, n_spec, vol_frac, MM, M1+1, &
          MM, n_ini, bin_v, dlnr, V, M2)
     
     M = M1 + M2
     env%V_comp = dble(M) / N_0
     ! call equlibriate_particle for each particle in V
     do i = 1,M
        call equilibriate_particle(n_spec, V(i,:), env, mat)
     enddo
     call mc_fix_hybrid(MM, M, n_spec, V, n_bin, MH, VH, &
          bin_v, bin_g, bin_gs, bin_n, dlnr , &
          kernel_brown, t_max, t_print, t_state, t_progress ,del_t, &
          i_loop, env, mat)
     
  enddo
  
end program run_brown_fix_hybrid

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
