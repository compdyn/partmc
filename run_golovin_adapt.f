C Simulation with Golovin kernel and adaptive timestepping.

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      program MonteCarlo
 
      integer MM, n_bin, n_loop, scal
      real*8 t_max, rho_p, N_0, t_print, p_max, V_0
      real*8 r_samp_max, del_t_max
      parameter (MM = 10000)           ! number of particles
      parameter (n_bin = 160)          ! number of bins
      parameter (n_loop = 10)           ! number of loops
      parameter (scal = 3)             ! scale factor for bins
      parameter (t_max = 600d0)        ! total simulation time (seconds)
      parameter (rho_p = 1000d0)       ! particle density (kg/m^3)
      parameter (N_0 = 1d9)            ! particle number concentration (#/m^3)
      parameter (t_print = 60d0)       ! interval between printing (s)
      parameter (p_max = 0.01d0)       ! maximum coagulation probability
      parameter (r_samp_max = 0.005d0) ! maximum sampling ratio per timestep
      parameter (del_t_max = 1d0)      ! maximum timestep
      parameter (V_0 = 4.1886d-15)     ! mean volume of initial distribution

      integer M, M_comp, i_loop, k
      real*8 V(MM), V_comp, dlnr
      real*8 n_ini(n_bin), vv(n_bin),  rr(n_bin)
      real*8 g(n_bin), n_ln(n_bin)

      real*8 pi
      parameter (pi = 3.14159265358979323846d0)

      external kernel_golovin

      open(30,file='out_golovin_adapt.d')
      call print_header(n_loop, n_bin, nint(t_max / t_print) + 1)
      call srand(10)

      do i_loop = 1,n_loop
         
         call make_grid(n_bin, scal, rho_p, vv, rr, dlnr)
         
         ! define initial exponential distribution
         do k = 1,n_bin
            n_ini(k) = pi/2d0 * (2d0*rr(k))**3 * MM/V_0
     &           * exp(-(vv(k) / V_0))
         enddo

         call compute_volumes(n_bin, MM, n_ini, rr, dlnr, V, M_comp)
         M = M_comp
         V_comp = M / N_0

         call mc_adapt(MM, M, M_comp, V, V_comp, kernel_golovin, n_bin,
     &        vv, rr, g, n_ln, dlnr, t_max, t_print,
     &        p_max, r_samp_max, del_t_max)

      enddo

      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
