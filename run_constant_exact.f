! -*- mode: f90; -*-
! Copyright (C) 2005,2006 Nicole Riemer and Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.
!
! Exact solution with Golovin kernel.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program run_constant_exact
  
  use mod_bin
  use mod_mc_exact
  use mod_kernel_constant
  use mod_array
  use mod_environ
  use mod_material
  
  integer n_bin, n_loop, scal, n_spec
  real*8 t_max, rho_p, N_0, t_print, V_0, V_comp
  real*8 v_min
  parameter (n_spec = 1)
  parameter (n_bin = 160)        ! number of bins
  parameter (n_loop = 1)         ! number of loops
  parameter (scal = 3)           ! scale factor for bins
  parameter (t_max = 600d0)      ! total simulation time (seconds)
  parameter (rho_p = 1000d0)     ! particle density (kg/m^3)
  parameter (v_min = 1.d-24)   ! minimum volume (m^3) for making grid
  parameter (N_0 = 2d8)          ! particle number concentration (#/m^3)
  parameter (t_print = 60d0)     ! interval between printing (s)
  parameter (V_0 = 4.1886d-15)   ! mean volume of initial distribution
  parameter (V_comp = 1d0)       ! computational volume (dummy value)
  
  integer i_loop
  real*8 dlnr
  real*8 bin_v(n_bin), bin_r(n_bin), bin_g(n_bin)
  integer bin_n(n_bin)
  type(environ) :: env
  type(material) :: mat
  
  call allocate_material(mat, n_spec)
  ! FIXME: set rho and other material parameters
  
  ! FIXME: set environment parameters
  
  open(30,file='out_constant_exact.d')
  call print_header(n_loop, n_bin, n_spec,  &
       nint(t_max / t_print) + 1)
  
  do i_loop = 1,n_loop
     
     call make_bin_grid(n_bin, scal, v_min, bin_v, bin_r, dlnr)
     
     call mc_exact(n_bin, bin_v, bin_r, bin_g, bin_n, dlnr, &
          N_0, V_0, rho_p, soln_constant_exp_cond, t_max, t_print, &
          i_loop, V_comp, n_spec, env, mat)
     
  enddo
  
end program run_constant_exact

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
