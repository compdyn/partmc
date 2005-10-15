C Exact solution with Golovin kernel.

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      program MonteCarlo
 
      integer n_bin, n_loop, scal
      real*8 t_max, rho_p, N_0, t_print, V_0
      parameter (n_bin = 160)        ! number of bins
      parameter (n_loop = 1)         ! number of loops
      parameter (scal = 3)           ! scale factor for bins
      parameter (t_max = 600d0)      ! total simulation time (seconds)
      parameter (rho_p = 1000d0)     ! particle density (kg/m^3)
      parameter (N_0 = 1d9)          ! particle number concentration (#/m^3)
      parameter (t_print = 60d0)     ! interval between printing (s)
      parameter (V_0 = 4.1886d-15)   ! mean volume of initial distribution

      integer i_loop
      real*8 dlnr
      real*8 vv(n_bin), rr(n_bin), g(n_bin), n_ln(n_bin)

      external soln_golovin_exp

      open(30,file='out_golovin_exact.d')
      call print_header(n_loop, n_bin, nint(t_max / t_print) + 1)

      do i_loop = 1,n_loop
         
         call make_grid(n_bin, scal, rho_p, vv, rr, dlnr)
         
         call mc_exact(n_bin, vv, rr, g, n_ln, N_0, V_0,
     &        rho_p, soln_golovin_exp, t_max, t_print, i_loop)

      enddo

      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
