C Exact solution output.

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine mc_exact(n_bin, vv, rr, dp, g, n_ln, N_0, V_0,
     &     rho_p, soln, t_max, t_print)

      integer n_bin      ! INPUT: number of bins
      real*8 vv(n_bin)   ! INPUT: volume of bins
      real*8 rr(n_bin)   ! INPUT: radius of bins
      real*8 dp(n_bin)   ! INPUT: diameter of bins
      real*8 g(n_bin)    ! OUTPUT: mass in bins
      real*8 n_ln(n_bin) ! OUTPUT: number in bins
      real*8 N_0         ! INPUT: particle number concentration (#/m^3)
      real*8 V_0         ! INPUT:
      real*8 rho_p       ! INPUT: particle density (kg/m^3)
      external soln      ! INPUT: exact solution procedure
      real*8 t_max       ! INPUT: total simulation time
      real*8 t_print     ! INPUT: interval to print info (seconds)

      integer i_time, n_time
      real*8 time

      n_time = int(t_max / t_print)
      do i_time = 0,n_time
         time = dble(i_time) / dble(n_time) * dble(t_max)
         call soln(n_bin, vv, rr, dp, time, N_0, V_0, rho_p,
     &        g, n_ln)
         call print_info(n_bin, time, rr, g, n_ln)
      enddo

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
