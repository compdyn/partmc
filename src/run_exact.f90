! Copyright (C) 2005-2007 Nicole Riemer and Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.
!
! Exact solution output.

module mod_run_exact

  type run_exact_opt_t
    real*8 :: num_conc      ! particle number concentration (#/m^3)
    real*8 :: mean_vol      ! mean init volume (m^3)
    real*8 :: rho_p         ! particle density (kg/m^3)
    real*8 :: t_max         ! total simulation time
    real*8 :: t_output      ! interval to output info (seconds)
    integer :: output_unit  ! unit number to output to
  end type run_exact_opt_t

contains
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine run_exact(bin_grid, env, aero_data, exact_opt, soln)

    ! FIXME: num_conc and mean_vol are really parameters for the
    ! initial value of the particle distribution. They should be
    ! replaced by a n_param, params() pair.

    ! "Run" an exact solution, output data in the same format as
    ! particle-resolved or sectional simulations.
    
    use mod_bin
    use mod_aero_state
    use mod_environ
    use mod_aero_data
    use mod_output_summary
    
    type(bin_grid_t), intent(in) :: bin_grid ! bin grid
    type(environ), intent(inout) :: env ! environment state
    type(aero_data_t), intent(in) :: aero_data ! aerosol data
    type(run_exact_opt_t), intent(in) :: exact_opt ! options
    
    integer :: i_time, n_time, n_spec
    real*8 :: bin_v_den(bin_grid%n_bin), bin_n_den(bin_grid%n_bin)
    real*8 :: bin_vs_den(bin_grid%n_bin,aero_data%n_spec), time
    
    interface
       subroutine soln(bin_grid, bin_v_den, bin_n_den, &
            time, num_conc, mean_vol, rho_p, env)

         use mod_bin
         use mod_environ

         type(bin_grid_t), intent(in) :: bin_grid ! bin grid
         real*8, intent(out) :: bin_v_den(bin_grid%n_bin) ! vol density in bins
         real*8, intent(out) :: bin_n_den(bin_grid%n_bin) ! num density in bins
         real*8, intent(in) :: time              ! current time
         real*8, intent(in) :: num_conc          ! particle number conc (#/m^3)
         real*8, intent(in) :: mean_vol          ! mean init volume (m^3)
         real*8, intent(in) :: rho_p             ! particle density (kg/m^3)
         type(environ), intent(in) :: env        ! environment state
       end subroutine soln
    end interface

    n_spec = 1
    n_time = nint(exact_opt%t_max / exact_opt%t_output)
    call init_environ(env, 0d0)
    bin_vs_den = 0d0
    do i_time = 0,n_time
       time = dble(i_time) / dble(n_time) * exact_opt%t_max
       call update_environ(env, time)
       call soln(bin_grid, bin_v_den, bin_n_den, time, &
            exact_opt%num_conc, exact_opt%mean_vol, exact_opt%rho_p, env)
       bin_vs_den(:,1) = bin_v_den

       call output_summary_density(exact_opt%output_unit, time, &
            bin_grid, aero_data, bin_v_den, bin_vs_den, bin_n_den, env, 1)
    end do
    
  end subroutine run_exact
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
end module mod_run_exact
