C Exact solution output.

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine mc_exact(n_bin, bin_v, bin_r, bin_g, bin_n, dlnr,
     &     N_0, V_0, rho_p, soln, t_max, t_print, loop, V_comp, n_spec)
      ! FIXME: N_0 and V_0 are really parameters for the initial value
      ! of the particle distribution. They should be replaced by a n_param,
      ! params() pair.

      integer n_bin        ! INPUT: number of bins
      integer n_spec
      real*8 bin_v(n_bin)  ! INPUT: volume of bins
      real*8 bin_r(n_bin)  ! INPUT: radius of bins
      real*8 bin_g(n_bin)  ! OUTPUT: mass in bins
      integer bin_n(n_bin) ! OUTPUT: number in bins
      real*8 bin_gs(n_bin,n_spec)
      real*8 dlnr          ! INPUT: bin scale factor
      real*8 N_0           ! INPUT: particle number concentration (#/m^3)
      real*8 V_0           ! INPUT:
      real*8 rho_p         ! INPUT: particle density (kg/m^3)
      external soln        ! INPUT: exact solution procedure
      real*8 t_max         ! INPUT: total simulation time
      real*8 t_print       ! INPUT: interval to print info (seconds)
      integer loop         ! INPUT: loop number of run
      real*8 V_comp        ! INPUT: computational volume

      integer i_time, n_time,i
      real*8 time

      n_time = int(t_max / t_print)
      do i_time = 0,n_time
         time = dble(i_time) / dble(n_time) * dble(t_max)
         call soln(n_bin, bin_v, bin_r, bin_g, bin_n, dlnr,
     &        time, N_0, V_0, rho_p, V_comp)

         do i=1,n_bin
            bin_gs(i,1) = bin_g(i)
         enddo

         call print_info(time, V_comp,n_spec,
     &        n_bin, bin_v, bin_r, bin_g, bin_gs,bin_n, dlnr)
      enddo

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
