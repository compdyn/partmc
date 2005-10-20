C Simulation with sedimentation kernel and adaptive timestepping.

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      program MonteCarlo
 
      integer MM, n_bin, n_loop, scal
      real*8 t_max, rho_p, N_0, t_print
      real*8 r_samp_max, del_t_max
      parameter (MM = 10000)           ! number of particles
      parameter (n_bin = 160)          ! number of bins
      parameter (n_loop = 1)           ! number of loops
      parameter (scal = 3)             ! scale factor for bins
      parameter (t_max = 600d0)        ! total simulation time (seconds)
      parameter (rho_p = 1000d0)       ! particle density (kg/m^3)
      parameter (N_0 = 1d9)            ! particle number concentration (#/m^3)
      parameter (t_print = 60)         ! interval between printing (s)
      parameter (r_samp_max = 0.005d0) ! maximum sampling ratio per timestep
      parameter (del_t_max = 1d0)      ! maximum timestep

      integer M, i_loop, k
      real*8 V(MM), V_comp, dlnr
      real*8 n_ini(n_bin), bin_v(n_bin), bin_r(n_bin)
      real*8 bin_g(n_bin)
      integer bin_n(n_bin)

      external kernel_sedi

      open(30,file='out_sedi_adapt.d')
      call print_header(n_loop, n_bin, nint(t_max / t_print) + 1)
      call srand(10)

      do i_loop = 1,n_loop

         call make_grid(n_bin, scal, rho_p, bin_v, bin_r, dlnr)
         
         ! define bidisperse distribution
         do k = 1,n_bin
            n_ini(k) = 0d0
         enddo
         n_ini(97) = (MM - 1) / dlnr
         n_ini(126) = 1 / dlnr

         call compute_volumes(n_bin, MM, n_ini, bin_r, dlnr, V, M)
         V_comp = M / N_0

         do k = 1,MM
            write(*,*)'k,V(k) = ', k, V(k)
         enddo

         call exit(2)
         
         call mc_adapt(MM, M, V, V_comp,
     &        n_bin, bin_v, bin_r, bin_g, bin_n, dlnr,
     &        kernel_sedi, t_max, t_print,
     &        r_samp_max, del_t_max, i_loop)

      enddo

      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
